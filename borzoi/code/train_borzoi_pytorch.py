#!/usr/bin/env python
"""
Borzoi Fine-Tuning with PyTorch and Hugging Face Accelerate

Based on the official training-borzoi repository:
https://github.com/johahi/training-borzoi
"""

import os
from pathlib import Path

# Suppress TensorFlow warnings (borzoi-pytorch may import TF)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# NCCL settings for better stability
os.environ.setdefault('NCCL_DEBUG', 'WARN')
os.environ.setdefault('NCCL_IB_DISABLE', '0')
os.environ.setdefault('NCCL_P2P_DISABLE', '0')

import torch
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.benchmark = True

import torch.nn as nn

import tqdm
from accelerate import Accelerator
from accelerate.utils import DistributedDataParallelKwargs

# Local imports
from src.config import Config, create_arg_parser
from src.data_loaders import create_data_loaders, create_data_loaders_from_datasets
from src.losses import setup_loss
from src.metrics import get_metrics
from src.utils import set_seed
from src.model import load_borzoi_model, load_borzoi_hires_model, detect_species
from src.training_utils import (
    set_trainable_layers,
    build_optimizer_and_scheduler,
    GracefulKiller,
    save_global_step,
    save_training_stage,
    find_latest_checkpoint,
    detect_checkpoint_stage,
    load_checkpoint_state,
    load_training_stage,
    initialize_wandb,
    training_step,
    handle_oom,
    transition_to_finetune,
    get_metric_direction,
    is_better,
    get_initial_best_value,
    save_metrics_json,
    update_best_model_symlink,
    ThroughputTracker,
)
from src.evaluation import evaluate
from src.batch_manager import AdaptiveBatchManager
from src.inference import run_post_training_inference
from src.logging_utils import log_training_schedule


# =============================================================================
# Setup Functions
# =============================================================================

def setup_accelerator_and_config(args):
    """Initialize accelerator, config, and checkpoint directory.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Tuple of (config, accelerator, ckpt_dir, args).
    """
    config = Config.from_args(args)

    # HiRes models have unused internal encoder heads/layers that cause DDP errors
    # unless we explicitly set find_unused_parameters=True.
    find_unused = config.get('model.hires_resolution') is not None
    ddp_kwargs = DistributedDataParallelKwargs(find_unused_parameters=find_unused)
    
    accelerator = Accelerator(
        log_with="wandb" if config.logging.use_wandb else None,
        step_scheduler_with_optimizer=False,
        mixed_precision="fp16" if config.training.mixed_precision else "no",
        kwargs_handlers=[ddp_kwargs],
    )

    seed = config.get('training.seed', 42)
    set_seed(seed)
    accelerator.print(f"Random seed set to: {seed}")

    accelerator.print("\n" + "=" * 70)
    accelerator.print("Borzoi Fine-Tuning with PyTorch + Accelerate")
    accelerator.print("=" * 70)
    accelerator.print(f"Number of processes: {accelerator.num_processes}")
    accelerator.print(f"Mixed precision: {accelerator.mixed_precision}")

    config.validate()

    ckpt_dir = Path(config.logging.checkpoint_dir)
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    if accelerator.is_main_process:
        config.save(ckpt_dir / "config.yaml")

    return config, accelerator, ckpt_dir


def setup_model_and_data(config, accelerator):
    """Load data and model, configure for species.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance.

    Returns:
        Tuple of (model, train_loader, val_loader, train_dataset, val_dataset,
                  is_human, track_names_for_inference).
    """
    is_mouse, is_human = detect_species(config)

    train_loader, val_loader, train_dataset, val_dataset = create_data_loaders(config, accelerator)

    track_names_for_inference = None
    if hasattr(train_dataset, 'track_names') and train_dataset.track_names:
        track_names_for_inference = list(train_dataset.track_names)

    hires_resolution = config.get('model.hires_resolution', None)
    if hires_resolution is not None:
        model = load_borzoi_hires_model(config, accelerator)
    else:
        model = load_borzoi_model(config, accelerator)

    if hires_resolution is None:
        # Standard Borzoi: validate species-specific head
        if not is_human:
            if not hasattr(model, 'mouse_head'):
                raise RuntimeError("Mouse data detected but model has no mouse_head!")
            accelerator.print(f"Using mouse_head for mouse data")
        else:
            accelerator.print(f"Using human_head for human data")
        accelerator.print(f"Routing model through {'human' if is_human else 'mouse'} head")
    else:
        accelerator.print(f"BorzoiHiRes mode: using dedicated {hires_resolution}bp decoder head")

    model = model.to(accelerator.device)
    accelerator.wait_for_everyone()

    use_sync_bn = config.get('training.use_sync_batchnorm', False)
    if accelerator.num_processes > 1 and use_sync_bn:
        model = nn.SyncBatchNorm.convert_sync_batchnorm(model)
        accelerator.print("Converted BatchNorm to SyncBatchNorm")

    return model, train_loader, val_loader, train_dataset, val_dataset, is_human, track_names_for_inference


def setup_metrics(config, accelerator, track_names=None):
    """Setup loss function and metrics.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance.
        track_names: Optional list of output track names for subset metric selectors.

    Returns:
        Tuple of (loss_fn, metrics, train_metrics).
    """
    loss_fn = setup_loss(config, accelerator)

    # Get PPV config from validation section
    ppv_config = config.get('validation.ppv', {})
    ppv_enabled = ppv_config.get('enabled', False)

    # Get genome-wide Pearson subset config from validation section
    genome_wide_pearson_config = config.get('validation.genome_wide_pearson', {})

    # Get quantile overlap config from validation section
    quantile_overlap_config = config.get('validation.quantile_overlap', {})
    quantile_overlap_enabled = quantile_overlap_config.get('enabled', False)

    # Get AUPRC config from validation section
    auprc_config = config.get('validation.auprc', {})
    auprc_enabled = auprc_config.get('enabled', False)

    # Validation metrics include genome_wide_pearson for cross-resolution comparability
    val_metric_names = ['pearson', 'mse', 'genome_wide_pearson']
    if ppv_enabled:
        val_metric_names.append('ppv')
        accelerator.print(f"\nPPV metric enabled with config: quantile={ppv_config.get('quantile', 0.98)}, "
                         f"smooth_bins={ppv_config.get('smooth_bins', 10)}")

    if quantile_overlap_enabled:
        val_metric_names.append('quantile_overlap')
        accelerator.print(f"\nQuantile overlap metrics enabled with config: "
                         f"quantile={quantile_overlap_config.get('quantile', 0.98)}, "
                         f"smooth_bins={quantile_overlap_config.get('smooth_bins', 1)}, "
                         f"use_log_scale={quantile_overlap_config.get('use_log_scale', False)}")

    if auprc_enabled:
        val_metric_names.append('auprc')
        accelerator.print(f"\nAUPRC metric enabled with config: "
                         f"quantile={auprc_config.get('quantile', 0.98)}, "
                         f"smooth_bins={auprc_config.get('smooth_bins', 1)}, "
                         f"use_log_scale={auprc_config.get('use_log_scale', False)}")

    metrics = get_metrics(
        val_metric_names,
        accelerator=accelerator,
        genome_wide_pearson_config=genome_wide_pearson_config,
        ppv_config=ppv_config,
        quantile_overlap_config=quantile_overlap_config,
        auprc_config=auprc_config,
        track_names=track_names,
    )

    # Training metrics don't need genome_wide_pearson or ppv (computed per-batch anyway)
    train_metrics = get_metrics(['pearson', 'mse'])
    for _, metric in train_metrics:
        metric.to(accelerator.device)

    region_config = config.get('validation.region_metrics', {})
    if region_config.get('enabled', False):
        from src.metrics import RegionR2

        regions_bed = region_config.get('bed_file')
        if regions_bed:
            accelerator.print(f"\nAdding RegionR2 metrics from: {regions_bed}")

            val_chroms = config.get('data.val_chroms', [])
            if val_chroms:
                region_r2_val = RegionR2(
                    regions_bed=regions_bed,
                    bin_size=config.model.bin_size,
                    aggregation_mode=region_config.get('aggregation_mode', 'both'),
                    overlap_mode=region_config.get('overlap_mode', 'strict'),
                    filter_chroms=val_chroms,
                    metric_prefix='region_r2',
                    accelerator=accelerator,
                )
                metrics.append(('region_r2', region_r2_val))
                accelerator.print(f"  Filtering to val chromosomes: {val_chroms}")
                accelerator.print(f"  Aggregation mode: {region_config.get('aggregation_mode', 'both')}")
                accelerator.print(f"  Overlap mode: {region_config.get('overlap_mode', 'strict')}")
        else:
            accelerator.print("Warning: region_metrics enabled but no bed_file specified")

    return loss_fn, metrics, train_metrics


def setup_training_schedule(config, accelerator, train_loader):
    """Calculate training schedule and stage parameters.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance.
        train_loader: Training data loader.

    Returns:
        Dict with training schedule parameters.
    """
    steps_per_epoch = len(train_loader)
    num_epochs = config.training.epochs
    total_steps = steps_per_epoch * num_epochs

    accelerator.print(f"\nTraining configuration:")
    accelerator.print(f"  Epochs: {num_epochs}")
    accelerator.print(f"  Steps per epoch: {steps_per_epoch}")
    accelerator.print(f"  Total steps: {total_steps}")

    lr = config.training.learning_rate
    head_lr = config.get('training.head_learning_rate', lr)
    finetune_lr = config.get('training.finetune_learning_rate', lr)
    warmup_steps = config.training.warmup_steps
    weight_decay = config.training.weight_decay
    weight_decay_transformer = getattr(config.training, 'weight_decay_transformer', 2e-8)
    linear_probe_epochs = config.get('training.linear_probe_epochs', 1)
    if config.get('model.train_from_scratch', False):
        if linear_probe_epochs != 0:
            accelerator.print("train_from_scratch enabled; forcing linear_probe_epochs=0")
        linear_probe_epochs = 0
    allow_unfreeze = not config.model.freeze_trunk
    finetune_batch_size = config.get('training.finetune_batch_size', None)
    finetune_gradient_accumulation_steps = config.get('training.finetune_gradient_accumulation_steps', None)
    num_workers = config.training.num_workers

    if accelerator.num_processes > 1 and num_workers > 0:
        num_workers = max(1, num_workers // accelerator.num_processes)

    linear_probe_steps = min(total_steps, linear_probe_epochs * steps_per_epoch) if allow_unfreeze else total_steps
    if allow_unfreeze and linear_probe_steps >= total_steps:
        linear_probe_steps = max(0, total_steps // 2)
    finetune_steps = max(0, total_steps - linear_probe_steps)

    log_training_schedule(accelerator, linear_probe_steps, finetune_steps, steps_per_epoch, head_lr, finetune_lr,
                          config.training.batch_size, finetune_batch_size, allow_unfreeze)

    return {
        'steps_per_epoch': steps_per_epoch,
        'num_epochs': num_epochs,
        'total_steps': total_steps,
        'lr': lr,
        'head_lr': head_lr,
        'finetune_lr': finetune_lr,
        'warmup_steps': warmup_steps,
        'weight_decay': weight_decay,
        'weight_decay_transformer': weight_decay_transformer,
        'linear_probe_steps': linear_probe_steps,
        'linear_probe_epochs': linear_probe_epochs,  # Store epochs for recalculation after OOM
        'finetune_steps': finetune_steps,
        'allow_unfreeze': allow_unfreeze,
        'finetune_batch_size': finetune_batch_size,
        'finetune_gradient_accumulation_steps': finetune_gradient_accumulation_steps,
        'num_workers': num_workers,
    }


def setup_optimizer_for_stage(model, accelerator, args, config, ckpt_dir, schedule, is_human, resume_info=None):
    """Configure optimizer based on current/resume stage.

    Args:
        model: The model to optimize.
        accelerator: Accelerator instance.
        args: Command-line arguments.
        config: Configuration object.
        ckpt_dir: Checkpoint directory.
        schedule: Training schedule dict from setup_training_schedule.
        is_human: Whether using human data.
        resume_info: Optional precomputed resume info dict.

    Returns:
        Tuple of (optimizer, scheduler, current_stage, resume_info).
    """
    if resume_info is None:
        resume_info = detect_checkpoint_stage(args, config, ckpt_dir, accelerator, schedule['linear_probe_steps'])

    if resume_info is not None and schedule['linear_probe_steps'] == 0 and resume_info['stage'] != 'finetune':
        accelerator.print(
            f"  Warning: checkpoint stage={resume_info['stage']} but linear_probe_steps=0; forcing finetune."
        )
        resume_info['stage'] = 'finetune'

    if resume_info is not None:
        current_stage = resume_info['stage']
        accelerator.print(f"\nResuming in stage: {current_stage}")
        if current_stage == 'finetune':
            accelerator.print("Configuring optimizer for FINETUNE stage (full model)")
            set_trainable_layers(model, accelerator, train_head_only=False, is_human=is_human)
            stage_lr = schedule['finetune_lr']
            remaining_steps = max(0, schedule['total_steps'] - resume_info['global_step'])
            stage_total_steps = max(remaining_steps, 1)
        else:
            accelerator.print("Configuring optimizer for LINEAR PROBE stage (head only)")
            set_trainable_layers(model, accelerator, train_head_only=True, is_human=is_human)
            stage_lr = schedule['head_lr']
            stage_total_steps = schedule['linear_probe_steps']
    else:
        current_stage = "linear_probe" if schedule['linear_probe_steps'] > 0 else "finetune"
        if schedule['linear_probe_steps'] > 0:
            accelerator.print("STAGE 1: Training Head Only (Linear Probing)")
            set_trainable_layers(model, accelerator, train_head_only=True, is_human=is_human)
            stage_lr = schedule['head_lr']
            stage_total_steps = schedule['linear_probe_steps']
        else:
            accelerator.print("No linear probing configured; training all layers")
            set_trainable_layers(model, accelerator, train_head_only=False, is_human=is_human)
            stage_lr = schedule['finetune_lr']
            stage_total_steps = schedule['total_steps']

    optimizer, scheduler = build_optimizer_and_scheduler(
        model, accelerator, lr=stage_lr, weight_decay=schedule['weight_decay'],
        weight_decay_transformer=schedule['weight_decay_transformer'],
        warmup_steps=min(schedule['warmup_steps'], stage_total_steps),
        total_steps=stage_total_steps,
    )

    return optimizer, scheduler, current_stage, resume_info


def setup_batch_manager(config, accelerator):
    """Initialize batch manager for OOM protection.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance.

    Returns:
        AdaptiveBatchManager or None if disabled.
    """
    enable_oom_protection = config.get('training.enable_oom_protection', False)
    if not enable_oom_protection:
        accelerator.print(f"\nOOM protection disabled (training.enable_oom_protection=false)")
        return None

    batch_manager = AdaptiveBatchManager(
        initial_batch_size=config.training.batch_size,
        min_batch_size=config.get('training.min_batch_size', 1),
        reduction_factor=config.get('training.batch_size_reduction_factor', 0.5),
        use_gradient_accumulation=config.get('training.use_gradient_accumulation', True),
        max_accumulation_steps=config.get('training.max_gradient_accumulation_steps', 16),
        recovery_delay=config.get('training.oom_recovery_delay', 2.0),
        max_batch_size=config.get('training.max_batch_size', config.training.batch_size),
        growth_factor=config.get('training.batch_size_growth_factor', 1.25),
        growth_interval=config.get('training.batch_growth_interval', 400),
        enable_growth=config.get('training.enable_auto_batch_growth', True),
    )
    accelerator.print(f"\nOOM protection enabled:")
    accelerator.print(f"  Initial batch size: {batch_manager.current_batch_size}")
    accelerator.print(f"  Min batch size: {batch_manager.min_batch_size}")
    return batch_manager


# =============================================================================
# Training Loop Functions
# =============================================================================

def run_epoch(epoch, train_loader, model, loss_fn, optimizer, scheduler, accelerator,
              is_human, config, schedule, train_metrics, batch_manager, killer,
              train_dataset, val_dataset, val_loader, metrics, ckpt_dir,
              state, throughput_tracker=None):
    """Run a single training epoch.

    Args:
        epoch: Current epoch number.
        train_loader: Training data loader.
        model: The model.
        loss_fn: Loss function.
        optimizer: Optimizer.
        scheduler: Learning rate scheduler.
        accelerator: Accelerator instance.
        is_human: Whether using human data.
        config: Configuration object.
        schedule: Training schedule dict.
        train_metrics: Training metrics.
        batch_manager: Batch manager for OOM protection.
        killer: GracefulKiller instance.
        train_dataset: Training dataset.
        val_dataset: Validation dataset.
        val_loader: Validation data loader.
        metrics: Validation metrics.
        ckpt_dir: Checkpoint directory.
        state: Mutable training state dict.
        throughput_tracker: Optional ThroughputTracker for performance metrics.

    Returns:
        Tuple of (updated train_loader, updated val_loader, optimizer, scheduler, restart_epoch flag).
    """
    clip_grad_norm = config.training.clip_grad_norm
    eval_every_n = config.logging.eval_every_n
    save_every_n = config.logging.save_every_n
    log_every_n = config.logging.log_every_n
    warmup_steps = schedule.get('warmup_steps', 0)

    accelerator.print(f"\nEpoch {epoch + 1}/{schedule['num_epochs']}")
    epoch_loss = 0.0
    epoch_steps = 0
    restart_epoch = False

    # Accumulator for detailed metrics (averaged over log_every_n steps)
    detailed_metrics_accum = {}
    detailed_metrics_count = 0

    progress_bar = tqdm.tqdm(
        train_loader, desc=f"Epoch {epoch + 1}", disable=not accelerator.is_main_process
    )

    for batch_idx, batch_data in enumerate(progress_bar):
        if state['stop_training'] or killer.kill_now:
            break

        # Stage transition
        if (state['current_stage'] == "linear_probe" and schedule['allow_unfreeze'] and
            schedule['finetune_steps'] > 0 and state['global_step'] >= state['stage_transition_step']):

            train_loader, val_loader, optimizer, scheduler, restart_epoch = handle_stage_transition(
                model, accelerator, config, train_dataset, val_dataset, val_loader,
                schedule, ckpt_dir, state, train_loader, batch_manager, is_human, epoch
            )
            if restart_epoch:
                break

        # Training step with optional OOM handling
        # Determine if we should compute detailed metrics this step
        compute_detailed = (state['global_step'] % log_every_n == 0) or (state['global_step'] < 100)

        # Start throughput tracking
        if throughput_tracker is not None:
            throughput_tracker.start_step()

        try:
            loss, _, step_detailed_metrics = training_step(
                model, batch_data, loss_fn, optimizer, scheduler, accelerator,
                is_human, clip_grad_norm, batch_manager, state['global_step'], train_metrics,
                compute_detailed_metrics=compute_detailed,
                warmup_steps=warmup_steps,
            )

            # Accumulate detailed metrics
            if step_detailed_metrics:
                for k, v in step_detailed_metrics.items():
                    if k not in detailed_metrics_accum:
                        detailed_metrics_accum[k] = 0.0
                    detailed_metrics_accum[k] += v
                detailed_metrics_count += 1

        except RuntimeError as e:
            if 'out of memory' in str(e).lower() and batch_manager:
                accelerator.print(f"\n  OOM detected at batch {batch_idx}!")
                train_loader, val_loader, success = handle_oom(
                    batch_manager, optimizer, accelerator, train_dataset, val_dataset, schedule['num_workers']
                )
                if not success:
                    state['stop_training'] = True
                    break

                # Update steps_per_epoch and recalculate stage-related values
                new_steps_per_epoch = len(train_loader)
                state['steps_per_epoch'] = new_steps_per_epoch
                state['current_batch_size'] = batch_manager.current_batch_size

                # Recalculate stage transition step based on epochs (not original steps)
                # This ensures the transition happens at the correct epoch boundary
                linear_probe_epochs = schedule.get('linear_probe_epochs', 1)
                if schedule['allow_unfreeze'] and state['current_stage'] == 'linear_probe':
                    new_linear_probe_steps = min(
                        new_steps_per_epoch * schedule['num_epochs'],
                        linear_probe_epochs * new_steps_per_epoch
                    )
                    state['stage_transition_step'] = new_linear_probe_steps
                    accelerator.print(f"  Recalculated stage_transition_step: {new_linear_probe_steps}")

                # Recalculate total steps
                state['total_steps'] = new_steps_per_epoch * schedule['num_epochs']
                accelerator.print(f"  Recalculated total_steps: {state['total_steps']}")

                restart_epoch = True
                break
            else:
                raise

        if batch_manager:
            batch_manager.register_successful_step()

        # End throughput tracking
        if throughput_tracker is not None:
            batch_size = batch_data[0].shape[0] if hasattr(batch_data[0], 'shape') else 1
            throughput_tracker.end_step(batch_size)

        epoch_loss += loss
        epoch_steps += 1
        state['global_step'] += 1

        current_lr = scheduler.get_last_lr()[0]
        progress_bar.set_postfix({'loss': f'{loss:.4f}', 'lr': f'{current_lr:.2e}'})

        # Logging
        if state['global_step'] % log_every_n == 0:
            log_step = state['global_step'] + state.get('log_step_offset', 0)
            log_dict = {"train/loss": loss, "train/lr": current_lr, "train/epoch": epoch + 1}
            for name, metric in train_metrics:
                log_dict[f"train/{name}"] = metric.compute()
                metric.reset()
            if batch_manager:
                log_dict.update({f"train/{k}": v for k, v in batch_manager.get_status_dict().items()})

            # Add averaged detailed metrics
            if detailed_metrics_count > 0:
                for k, v in detailed_metrics_accum.items():
                    log_dict[f"train/{k}"] = v / detailed_metrics_count
                # Reset accumulators
                detailed_metrics_accum = {}
                detailed_metrics_count = 0

            # Add throughput metrics
            if throughput_tracker is not None:
                throughput_stats = throughput_tracker.get_throughput()
                log_dict.update({f"train/perf/{k}": v for k, v in throughput_stats.items()})

            accelerator.log(log_dict, step=log_step)

        # Periodic evaluation
        if state['global_step'] - state['last_eval_step'] >= eval_every_n:
            # Calculate center-crop bins for cross-resolution comparison
            val_pred_len = config.get('validation.pred_len', config.model.pred_len)
            val_center_crop_bins = None
            if val_pred_len != config.model.pred_len:
                val_center_crop_bins = val_pred_len // config.model.bin_size

            handle_evaluation(accelerator, model, val_loader, loss_fn, metrics, is_human,
                              ckpt_dir, state, config, val_center_crop_bins=val_center_crop_bins)

        # Periodic checkpoint
        if state['global_step'] % save_every_n == 0:
            accelerator.wait_for_everyone()
            if accelerator.is_main_process:
                config_info = f" (config: {config.source_name})" if config.source_name else ""
                accelerator.print(f"\n  Saving checkpoint at step {state['global_step']}{config_info}...")
                step_checkpoint_path = ckpt_dir / f"step_{state['global_step']}"
                accelerator.save_state(output_dir=str(step_checkpoint_path))
                save_global_step(step_checkpoint_path, state['global_step'])
                save_training_stage(step_checkpoint_path, state['current_stage'], state['current_batch_size'])
            accelerator.wait_for_everyone()

        if state['global_step'] >= state['total_steps']:
            state['stop_training'] = True
            break

    if not restart_epoch and epoch_steps > 0:
        avg_epoch_loss = epoch_loss / max(epoch_steps, 1)
        accelerator.print(f"\nEpoch {epoch + 1} completed. Average loss: {avg_epoch_loss:.4f}")

    return train_loader, val_loader, optimizer, scheduler, restart_epoch


def handle_stage_transition(model, accelerator, config, train_dataset, val_dataset, val_loader,
                            schedule, ckpt_dir, state, train_loader, batch_manager, is_human, epoch):
    """Handle transition from linear probe to fine-tuning stage.

    Args:
        Various training components and state.

    Returns:
        Tuple of (new_train_loader, new_val_loader, optimizer, scheduler, restart_epoch).
    """
    # Save checkpoint before transition
    accelerator.wait_for_everyone()
    if accelerator.is_main_process:
        config_info = f" (config: {config.source_name})" if config.source_name else ""
        accelerator.print(f"\n  Saving pre-transition checkpoint at step {state['global_step']}{config_info}...")
        pre_transition_path = ckpt_dir / "pre_finetune_transition"
        accelerator.save_state(output_dir=str(pre_transition_path))
        save_global_step(pre_transition_path, state['global_step'])
        save_training_stage(pre_transition_path, state['current_stage'], state['current_batch_size'])
    accelerator.wait_for_everyone()

    new_train_loader, new_val_loader, steps_per_epoch, total_steps, finetune_steps = transition_to_finetune(
        model, accelerator, config, train_dataset, val_dataset, schedule['num_workers'],
        schedule['finetune_batch_size'], schedule['finetune_lr'], schedule['weight_decay'],
        schedule['weight_decay_transformer'], schedule['warmup_steps'], schedule['finetune_steps'],
        schedule['linear_probe_steps'], is_human, batch_manager, epoch, schedule['num_epochs'],
        current_train_loader=train_loader,
        finetune_gradient_accumulation_steps=schedule['finetune_gradient_accumulation_steps'],
    )
    optimizer, scheduler = build_optimizer_and_scheduler(
        accelerator.unwrap_model(model), accelerator, lr=schedule['finetune_lr'],
        weight_decay=schedule['weight_decay'], weight_decay_transformer=schedule['weight_decay_transformer'],
        warmup_steps=min(schedule['warmup_steps'], finetune_steps), total_steps=max(finetune_steps, 1),
    )
    # Replace the tracked optimizer/scheduler so Accelerate checkpoints only the finetune state.
    accelerator._optimizers = []
    accelerator._schedulers = []
    optimizer, scheduler = accelerator.prepare(optimizer, scheduler)
    state['current_stage'] = "finetune"
    state['steps_per_epoch'] = steps_per_epoch
    state['total_steps'] = total_steps
    model.train()

    # Save checkpoint after successful transition
    accelerator.wait_for_everyone()
    if accelerator.is_main_process:
        config_info = f" (config: {config.source_name})" if config.source_name else ""
        accelerator.print(f"  Saving post-transition checkpoint{config_info}...")
        post_transition_path = ckpt_dir / "finetune_start"
        accelerator.save_state(output_dir=str(post_transition_path))
        save_global_step(post_transition_path, state['global_step'])
        save_training_stage(post_transition_path, state['current_stage'],
                           schedule['finetune_batch_size'] if schedule['finetune_batch_size'] else state['current_batch_size'])
    accelerator.wait_for_everyone()

    if new_train_loader is not None:
        state['current_batch_size'] = schedule['finetune_batch_size'] if schedule['finetune_batch_size'] else state['current_batch_size']
        return new_train_loader, new_val_loader, optimizer, scheduler, True

    return train_loader, val_loader, optimizer, scheduler, False


def handle_evaluation(accelerator, model, val_loader, loss_fn, metrics, is_human,
                      ckpt_dir, state, config, val_center_crop_bins=None):
    """Handle periodic evaluation and checkpointing.

    Args:
        Various training components and state.
        config: Configuration object.
        val_center_crop_bins: If provided, center-crop predictions and targets to this
            many bins for cross-resolution comparison.
    """
    log_step = state['global_step'] + state.get('log_step_offset', 0)
    accelerator.print(f"\n  Evaluating at step {state['global_step']}...")
    val_results = evaluate(accelerator, model, val_loader, loss_fn, metrics, is_human,
                          val_center_crop_bins=val_center_crop_bins)

    accelerator.print(f"  Validation results:")
    for key, value in val_results.items():
        accelerator.print(f"    {key}: {value:.4f}")

    if accelerator.is_main_process:
        accelerator.print(f"  Logging {len(val_results)} metrics to WandB at step {log_step}")
    accelerator.log(val_results, step=log_step)

    # Multi-metric best model tracking
    if accelerator.is_main_process:
        save_metrics = config.get('logging.save_metrics_json', True)
        track_best_metrics = config.get('logging.track_best_metrics', ['val_loss'])
        primary_metric = config.get('training.best_model_metric', 'val_loss')
        primary_mode = config.get('training.best_model_mode', None)

        # Normalize primary metric name
        if not primary_metric.startswith('val_'):
            primary_metric = f'val_{primary_metric}'

        # Auto-detect mode if not specified
        if primary_mode is None:
            primary_mode = get_metric_direction(primary_metric)

        for metric_name in track_best_metrics:
            # Normalize metric name
            normalized_name = metric_name if metric_name.startswith('val_') else f'val_{metric_name}'

            if normalized_name not in val_results:
                continue

            new_value = val_results[normalized_name]
            mode = get_metric_direction(normalized_name)
            old_value = state['best_metrics'].get(normalized_name, get_initial_best_value(mode))

            if is_better(new_value, old_value, mode):
                state['best_metrics'][normalized_name] = new_value

                # Save metric-specific best checkpoint
                checkpoint_name = f"best_model_{normalized_name}"
                checkpoint_path = ckpt_dir / checkpoint_name

                config_info = f" (config: {config.source_name})" if config.source_name else ""
                accelerator.print(f"  New best {normalized_name}: {new_value:.4f} - saving{config_info}")
                accelerator.save_state(output_dir=str(checkpoint_path))
                save_global_step(checkpoint_path, state['global_step'])
                save_training_stage(checkpoint_path, state['current_stage'], state['current_batch_size'])

                if save_metrics:
                    save_metrics_json(
                        checkpoint_path,
                        val_results,
                        state['global_step'],
                        state['current_stage'],
                        state.get('epoch', 0),
                        best_for=[normalized_name],
                    )

                # Update primary best_model symlink if this is the primary metric
                if normalized_name == primary_metric:
                    update_best_model_symlink(ckpt_dir, checkpoint_name)
                    # Keep backward compat with state['best_val_loss']
                    state['best_val_loss'] = new_value

    # Early stopping check
    if state['early_stopping_patience'] is not None:
        metric_name = state.get('early_stopping_metric', 'val_loss')
        metric_mode = state.get('early_stopping_mode', 'min')
        metric_value = val_results.get(metric_name)
        if metric_value is None and not metric_name.startswith('val_'):
            prefixed_name = f"val_{metric_name}"
            metric_value = val_results.get(prefixed_name)
            if metric_value is not None:
                metric_name = prefixed_name
        if metric_value is None:
            accelerator.print(
                f"  Early stopping: metric '{metric_name}' not found in validation results; skipping check"
            )
        else:
            best_metric = state['best_val_metric_for_early_stop']
            if metric_mode == 'max':
                improvement = metric_value - best_metric
            else:
                improvement = best_metric - metric_value

            if improvement > state['min_improvement']:
                state['best_val_metric_for_early_stop'] = metric_value
                state['epochs_without_improvement'] = 0
                accelerator.print(
                    f"  Early stopping: {metric_name} improved by {improvement:.6f}, resetting patience"
                )
            else:
                state['epochs_without_improvement'] += 1
                accelerator.print(
                    f"  Early stopping: no improvement in {metric_name} "
                    f"({state['epochs_without_improvement']}/{state['early_stopping_patience']})"
                )
                if state['epochs_without_improvement'] >= state['early_stopping_patience']:
                    accelerator.print(
                        f"\n  Early stopping triggered after {state['early_stopping_patience']} evaluations "
                        "without improvement"
                    )
                    state['stop_training'] = True

    state['last_eval_step'] = state['global_step']
    model.train()


def finish_training(accelerator, killer, ckpt_dir, best_val_loss, track_names_for_inference, config):
    """Handle training completion and cleanup.

    Args:
        accelerator: Accelerator instance.
        killer: GracefulKiller instance.
        ckpt_dir: Checkpoint directory.
        best_val_loss: Best validation loss achieved.
        track_names_for_inference: Track names for post-training inference.
        config: Configuration object.
    """
    accelerator.print("\n" + "=" * 70)
    if killer.kill_now:
        accelerator.print("Training Interrupted by User!")
    else:
        accelerator.print("Training Completed!")
    accelerator.print(f"Best validation loss: {best_val_loss:.4f}")
    if accelerator.is_main_process:
        accelerator.print(f"Best model saved at: {ckpt_dir / 'best_model'}")
    accelerator.print("=" * 70)

    if not killer.kill_now:
        try:
            run_post_training_inference(config, accelerator, ckpt_dir, track_names_for_inference)
        except Exception as e:
            accelerator.print(f"Post-training inference failed: {e}")

    accelerator.end_training()


# =============================================================================
# Main Training Function
# =============================================================================

def main():
    """Main training entry point."""
    killer = GracefulKiller()

    # Parse arguments and setup
    parser = create_arg_parser()
    args = parser.parse_args()
    config, accelerator, ckpt_dir = setup_accelerator_and_config(args)

    # Load model and data
    model, train_loader, val_loader, train_dataset, val_dataset, is_human, track_names_for_inference = \
        setup_model_and_data(config, accelerator)

    # Setup loss and metrics
    loss_fn, metrics, train_metrics = setup_metrics(config, accelerator, track_names=track_names_for_inference)

    # Detect resume info early to honor saved batch size
    resume_dir = None
    resume_stage_info = None
    if args.resume or args.checkpoint_path:
        if args.checkpoint_path:
            resume_dir = Path(args.checkpoint_path)
        else:
            resume_dir = find_latest_checkpoint(ckpt_dir, accelerator)
        if resume_dir and resume_dir.exists():
            resume_stage_info = load_training_stage(resume_dir)

    desired_batch_size = None
    if config.get('training.force_finetune_on_resume', False) and config.get('training.finetune_batch_size', None):
        desired_batch_size = config.training.finetune_batch_size
    elif resume_stage_info and resume_stage_info.get('batch_size'):
        desired_batch_size = resume_stage_info['batch_size']

    if desired_batch_size and train_loader.batch_size != desired_batch_size:
        accelerator.print(
            f"\nResuming with batch size: {desired_batch_size} "
            f"(current loader batch size: {train_loader.batch_size})"
        )
        config.training.batch_size = desired_batch_size
        train_loader, val_loader = create_data_loaders_from_datasets(
            train_dataset, val_dataset, desired_batch_size, config.training.num_workers, accelerator, verbose=True
        )

    # Prepare dataloaders
    train_loader = accelerator.prepare(train_loader)

    # Calculate training schedule after DDP wrapping for accurate steps/epoch
    schedule = setup_training_schedule(config, accelerator, train_loader)

    # Detect resume info with the finalized schedule
    resume_info = detect_checkpoint_stage(args, config, ckpt_dir, accelerator, schedule['linear_probe_steps'])

    # Setup optimizer for current stage
    optimizer, scheduler, current_stage, resume_info = setup_optimizer_for_stage(
        model, accelerator, args, config, ckpt_dir, schedule, is_human, resume_info=resume_info
    )

    # Prepare model and optimizer
    accelerator.print("\nPreparing model for distributed training...")
    model, optimizer, scheduler = accelerator.prepare(model, optimizer, scheduler)
    accelerator.print("Model prepared successfully!")

    # Load checkpoint state
    starting_epoch, resume_global_step = load_checkpoint_state(resume_info, accelerator, schedule['steps_per_epoch'])
    val_loader = accelerator.prepare(val_loader)

    # Initialize WandB
    wandb_step_offset = initialize_wandb(config, accelerator, ckpt_dir, starting_epoch, resume_global_step)

    # Setup batch manager
    batch_manager = setup_batch_manager(config, accelerator)

    # Initialize training state
    global_step = resume_global_step if resume_global_step > 0 else starting_epoch * schedule['steps_per_epoch']
    early_stopping_metric = config.get('training.early_stopping_metric', 'val_loss')
    early_stopping_mode = config.get('training.early_stopping_mode', 'min')
    if early_stopping_mode not in {'min', 'max'}:
        accelerator.print(
            f"Warning: unknown early_stopping_mode '{early_stopping_mode}', defaulting to 'min'"
        )
        early_stopping_mode = 'min'

    early_stopping_best = float('inf') if early_stopping_mode == 'min' else float('-inf')

    # Initialize multi-metric best tracking
    track_best_metrics = config.get('logging.track_best_metrics', ['val_loss'])
    best_metrics = {}
    for metric_name in track_best_metrics:
        normalized_name = metric_name if metric_name.startswith('val_') else f'val_{metric_name}'
        mode = get_metric_direction(normalized_name)
        best_metrics[normalized_name] = get_initial_best_value(mode)

    state = {
        'global_step': global_step,
        'best_val_loss': float('inf'),
        'best_metrics': best_metrics,
        'current_stage': current_stage,
        'stage_transition_step': schedule['linear_probe_steps'],
        'current_batch_size': resume_info['batch_size'] if resume_info else config.training.batch_size,
        'stop_training': False,
        'last_eval_step': (global_step // config.logging.eval_every_n) * config.logging.eval_every_n,
        'log_step_offset': wandb_step_offset,
        'early_stopping_patience': config.get('training.early_stopping_patience', None),
        'early_stopping_metric': early_stopping_metric,
        'early_stopping_mode': early_stopping_mode,
        'min_improvement': config.get('training.min_improvement', 1e-4),
        'epochs_without_improvement': 0,
        'best_val_metric_for_early_stop': early_stopping_best,
        'steps_per_epoch': schedule['steps_per_epoch'],
        'total_steps': schedule['total_steps'],
    }

    # Training loop
    accelerator.print("\n" + "=" * 70)
    accelerator.print("Starting Training")
    accelerator.print("=" * 70)

    # Initialize throughput tracker for performance metrics
    throughput_tracker = ThroughputTracker(window_size=100)

    model.train()
    epoch = starting_epoch

    while epoch < schedule['num_epochs']:
        if killer.kill_now:
            config_info = f" (config: {config.source_name})" if config.source_name else ""
            accelerator.print(f"\nGraceful shutdown requested - saving checkpoint{config_info}...")
            if accelerator.is_main_process:
                checkpoint_path = ckpt_dir / f"interrupted_epoch_{epoch}"
                accelerator.save_state(output_dir=str(checkpoint_path))
                save_global_step(checkpoint_path, state['global_step'])
                save_training_stage(checkpoint_path, state['current_stage'], state['current_batch_size'])
            state['stop_training'] = True
            break

        train_loader, val_loader, optimizer, scheduler, restart_epoch = run_epoch(
            epoch, train_loader, model, loss_fn, optimizer, scheduler, accelerator,
            is_human, config, schedule, train_metrics, batch_manager, killer,
            train_dataset, val_dataset, val_loader, metrics, ckpt_dir, state,
            throughput_tracker=throughput_tracker,
        )

        if restart_epoch:
            continue

        if state['stop_training'] or killer.kill_now:
            break

        epoch += 1

    # Finish training
    finish_training(accelerator, killer, ckpt_dir, state['best_val_loss'],
                    track_names_for_inference, config)


if __name__ == "__main__":
    main()
