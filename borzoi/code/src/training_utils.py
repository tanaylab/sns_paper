"""Training utilities for optimizer, scheduler, and weight decay configuration."""

import os
import shutil
import signal
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import torch
import torch.nn as nn
import torch.utils.data
from torch.optim.lr_scheduler import SequentialLR, LinearLR, LambdaLR
from accelerate import Accelerator, DistributedType

if TYPE_CHECKING:
    from .batch_manager import AdaptiveBatchManager


# =============================================================================
# Training Metrics Computation
# =============================================================================

def compute_layer_gradient_norms(model: nn.Module) -> Dict[str, float]:
    """Compute gradient norms for different layer groups.

    Args:
        model: PyTorch model with gradients computed.

    Returns:
        Dict with gradient norms for transformer, conv, and head layers.
    """
    grad_norms = {
        'transformer': 0.0,
        'conv': 0.0,
        'head': 0.0,
        'other': 0.0,
    }
    grad_counts = {k: 0 for k in grad_norms}

    for name, param in model.named_parameters():
        if param.grad is None:
            continue

        grad_norm_sq = param.grad.data.norm(2).item() ** 2
        name_lower = name.lower()

        if 'transformer' in name_lower or 'attention' in name_lower:
            grad_norms['transformer'] += grad_norm_sq
            grad_counts['transformer'] += 1
        elif 'conv' in name_lower:
            grad_norms['conv'] += grad_norm_sq
            grad_counts['conv'] += 1
        elif 'head' in name_lower:
            grad_norms['head'] += grad_norm_sq
            grad_counts['head'] += 1
        else:
            grad_norms['other'] += grad_norm_sq
            grad_counts['other'] += 1

    # Take sqrt to get actual norms
    return {k: (v ** 0.5) for k, v in grad_norms.items()}


def compute_parameter_norm(model: nn.Module) -> float:
    """Compute total L2 norm of trainable parameters.

    Args:
        model: PyTorch model.

    Returns:
        L2 norm of all trainable parameters.
    """
    total_norm_sq = 0.0
    for param in model.parameters():
        if param.requires_grad:
            total_norm_sq += param.data.norm(2).item() ** 2
    return total_norm_sq ** 0.5


def compute_output_statistics(outputs: torch.Tensor) -> Dict[str, float]:
    """Compute statistics of model outputs for monitoring.

    Args:
        outputs: Model output tensor of shape (batch, bins, tracks).

    Returns:
        Dict with output statistics.
    """
    with torch.no_grad():
        return {
            'output_mean': outputs.mean().item(),
            'output_std': outputs.std().item(),
            'output_max': outputs.abs().max().item(),
            'output_min': outputs.min().item(),
            'frac_zeros': (outputs.abs() < 1e-6).float().mean().item(),
            'frac_negative': (outputs < 0).float().mean().item(),
        }


def compute_loss_statistics(
    loss_fn: nn.Module,
    outputs: torch.Tensor,
    targets: torch.Tensor,
) -> Dict[str, float]:
    """Compute per-sample loss statistics.

    Args:
        loss_fn: Loss function (must support reduction='none' or we estimate).
        outputs: Model outputs.
        targets: Target values.

    Returns:
        Dict with loss variance and other statistics.
    """
    with torch.no_grad():
        # Compute per-sample MSE as a proxy for loss variance
        per_sample_mse = ((outputs - targets) ** 2).mean(dim=[1, 2])
        return {
            'loss_variance': per_sample_mse.var().item(),
            'loss_max': per_sample_mse.max().item(),
            'loss_min': per_sample_mse.min().item(),
        }


def get_gpu_memory_stats() -> Dict[str, float]:
    """Get GPU memory usage statistics.

    Returns:
        Dict with memory statistics in GB.
    """
    if not torch.cuda.is_available():
        return {}

    return {
        'gpu_memory_allocated_gb': torch.cuda.memory_allocated() / 1e9,
        'gpu_memory_reserved_gb': torch.cuda.memory_reserved() / 1e9,
        'gpu_memory_max_allocated_gb': torch.cuda.max_memory_allocated() / 1e9,
    }


class ThroughputTracker:
    """Track training throughput (samples per second)."""

    def __init__(self, window_size: int = 100):
        """Initialize throughput tracker.

        Args:
            window_size: Number of steps to average over.
        """
        self.window_size = window_size
        self.step_times = []
        self.step_samples = []
        self.last_time = None

    def start_step(self):
        """Mark the start of a training step."""
        self.last_time = time.perf_counter()

    def end_step(self, num_samples: int):
        """Mark the end of a training step.

        Args:
            num_samples: Number of samples processed in this step.
        """
        if self.last_time is None:
            return

        elapsed = time.perf_counter() - self.last_time
        self.step_times.append(elapsed)
        self.step_samples.append(num_samples)

        # Keep only recent history
        if len(self.step_times) > self.window_size:
            self.step_times.pop(0)
            self.step_samples.pop(0)

        self.last_time = None

    def get_throughput(self) -> Dict[str, float]:
        """Get throughput statistics.

        Returns:
            Dict with samples_per_second, steps_per_second, avg_step_time.
        """
        if not self.step_times:
            return {}

        total_time = sum(self.step_times)
        total_samples = sum(self.step_samples)
        num_steps = len(self.step_times)

        return {
            'samples_per_second': total_samples / total_time if total_time > 0 else 0,
            'steps_per_second': num_steps / total_time if total_time > 0 else 0,
            'avg_step_time_ms': (total_time / num_steps) * 1000 if num_steps > 0 else 0,
        }

# =============================================================================
# Metric Direction Mapping
# =============================================================================
# Defines whether each metric should be minimized or maximized for "best" model selection.

METRIC_DIRECTIONS = {
    'val_loss': 'min',
    'val_mse': 'min',
    'val_pearson': 'max',
    'val_genome_wide_pearson': 'max',
    'val_ppv': 'max',
    'val_sensitivity': 'max',
    'val_iou': 'max',
    'val_precision': 'max',
    'val_recall': 'max',
    'val_f1': 'max',
    'val_auprc': 'max',
    'val_r2': 'max',
}


def get_metric_direction(metric_name: str) -> str:
    """Return 'min' or 'max' for the given metric.

    Args:
        metric_name: Name of the metric (with or without 'val_' prefix).

    Returns:
        'min' for metrics that should be minimized (loss, mse),
        'max' for metrics that should be maximized (correlation, ppv).
    """
    normalized = metric_name if metric_name.startswith('val_') else f'val_{metric_name}'
    if normalized in METRIC_DIRECTIONS:
        return METRIC_DIRECTIONS[normalized]

    # Support suffixed metrics like val_iou_atac or val_genome_wide_pearson_subset.
    # Longest keys first so specific base metrics (e.g. genome_wide_pearson) win.
    for base_metric in sorted(METRIC_DIRECTIONS.keys(), key=len, reverse=True):
        if normalized.startswith(f"{base_metric}_"):
            return METRIC_DIRECTIONS[base_metric]

    return 'min'  # Default to min for unknown


def is_better(new_value: float, old_value: float, mode: str) -> bool:
    """Check if new_value is better than old_value given optimization mode.

    Args:
        new_value: The new metric value.
        old_value: The old/best metric value.
        mode: 'min' for minimization, 'max' for maximization.

    Returns:
        True if new_value is better than old_value.
    """
    if mode == 'max':
        return new_value > old_value
    return new_value < old_value


def get_initial_best_value(mode: str) -> float:
    """Get initial best value for comparison.

    Args:
        mode: 'min' for minimization, 'max' for maximization.

    Returns:
        float('inf') for 'min' mode, float('-inf') for 'max' mode.
    """
    return float('-inf') if mode == 'max' else float('inf')


# =============================================================================
# Training Constants
# =============================================================================

# Weight decay for convolutional/other layers.
# Small value to prevent overfitting while allowing effective learning.
DEFAULT_WEIGHT_DECAY = 1e-8

# Weight decay for transformer attention layers.
# Slightly higher than conv layers due to transformer's tendency to overfit.
DEFAULT_WEIGHT_DECAY_TRANSFORMER = 2e-8

# Learning rate warmup start factor.
# Training starts at this fraction of the target LR and linearly increases.
# 0.01 = start at 1% of target LR to avoid large initial updates.
WARMUP_START_FACTOR = 0.01

# Learning rate end factor for linear decay.
# After warmup, LR decays to this fraction of the peak LR.
LR_END_FACTOR = 0.01


def add_weight_decay(
    model: nn.Module,
    accelerator: Accelerator,
    weight_decay: float = DEFAULT_WEIGHT_DECAY,
    weight_decay_transformer: float = DEFAULT_WEIGHT_DECAY_TRANSFORMER,
    skip_list: Tuple[str, ...] = ("bias", "LayerNorm", "layer_norm", "ln_"),
    verbose: bool = False,
) -> List[dict]:
    """Add weight decay to model parameters, following training-borzoi approach.

    Different weight decay for transformer layers vs other layers.
    Skip certain parameters like biases and layer norms.

    Args:
        model: PyTorch model.
        accelerator: Accelerator instance for printing.
        weight_decay: Weight decay for non-transformer layers.
        weight_decay_transformer: Weight decay for transformer layers.
        skip_list: Parameter name patterns to skip weight decay.
        verbose: If True, print each parameter assignment. Default False for cleaner logs.

    Returns:
        List of parameter groups for optimizer.
    """
    decay_transformer = []
    decay_other = []
    no_decay = []
    no_decay_names = []

    for name, param in model.named_parameters():
        if not param.requires_grad:
            continue

        if any(skip_name in name for skip_name in skip_list):
            no_decay.append(param)
            no_decay_names.append(name)
            if verbose:
                accelerator.print(f"  No decay: {name}")
        elif "transformer" in name.lower():
            decay_transformer.append(param)
        else:
            decay_other.append(param)

    # Print summary instead of per-parameter logging
    accelerator.print(f"  Parameter groups:")
    accelerator.print(f"    Transformer layers (wd={weight_decay_transformer}): {len(decay_transformer)} params")
    accelerator.print(f"    Other layers (wd={weight_decay}): {len(decay_other)} params")
    accelerator.print(f"    No decay (bias/norm): {len(no_decay)} params")

    return [
        {"params": decay_transformer, "weight_decay": weight_decay_transformer},
        {"params": decay_other, "weight_decay": weight_decay},
        {"params": no_decay, "weight_decay": 0.0},
    ]


def set_trainable_layers(
    model: nn.Module,
    accelerator: Accelerator,
    train_head_only: bool,
    is_human: bool = True,
) -> None:
    """Toggle requires_grad for linear probing vs full fine-tuning.

    Args:
        model: PyTorch model.
        accelerator: Accelerator instance for printing.
        train_head_only: If True, only train head layers.
        is_human: If True, freeze mouse_head; if False, freeze human_head.
    """
    total = 0
    trainable = 0
    for name, param in model.named_parameters():
        total += 1
        # Always freeze the opposite species head to avoid unused-param in DDP
        if (is_human and "mouse_head" in name) or (not is_human and "human_head" in name):
            param.requires_grad = False
        elif train_head_only and "head" not in name.lower():
            param.requires_grad = False
        else:
            param.requires_grad = True
            trainable += 1
    stage_desc = "head-only (linear probe)" if train_head_only else "full model"
    accelerator.print(f"Trainable parameters ({stage_desc}): {trainable}/{total}")


def build_optimizer_and_scheduler(
    model: nn.Module,
    accelerator: Accelerator,
    lr: float,
    weight_decay: float,
    weight_decay_transformer: float,
    warmup_steps: int,
    total_steps: int,
) -> Tuple[torch.optim.Optimizer, torch.optim.lr_scheduler.LRScheduler]:
    """Construct AdamW optimizer and LR scheduler with warmup + linear decay.

    Args:
        model: PyTorch model.
        accelerator: Accelerator instance.
        lr: Learning rate.
        weight_decay: Weight decay for non-transformer layers.
        weight_decay_transformer: Weight decay for transformer layers.
        warmup_steps: Number of warmup steps.
        total_steps: Total training steps.

    Returns:
        Tuple of (optimizer, scheduler).
    """
    parameters = add_weight_decay(
        model,
        accelerator,
        weight_decay=weight_decay,
        weight_decay_transformer=weight_decay_transformer,
    )

    optimizer = torch.optim.AdamW(parameters, lr=lr)

    warmup_steps = int(min(warmup_steps, total_steps))
    decay_steps = max(total_steps - warmup_steps, 1)

    def warmup_fn(current_step: int) -> float:
        if warmup_steps == 0:
            return 1.0
        # Linear warmup from WARMUP_START_FACTOR to 1.0
        progress = float(current_step) / float(max(1, warmup_steps))
        return WARMUP_START_FACTOR + progress * (1.0 - WARMUP_START_FACTOR)

    warmup_scheduler = LambdaLR(optimizer, lr_lambda=warmup_fn)

    train_scheduler = LinearLR(
        optimizer,
        start_factor=1.0,
        end_factor=LR_END_FACTOR,
        total_iters=decay_steps,
    )

    if warmup_steps > 0:
        scheduler = SequentialLR(
            optimizer,
            schedulers=[warmup_scheduler, train_scheduler],
            milestones=[warmup_steps],
        )
    else:
        scheduler = train_scheduler

    accelerator.print(f"Optimizer: AdamW (lr={lr}, weight_decay={weight_decay})")
    accelerator.print(
        f"Scheduler: Warmup ({warmup_steps} steps) + Linear decay over {decay_steps} steps"
    )

    return optimizer, scheduler


class GracefulKiller:
    """Handle SIGINT and SIGTERM for graceful shutdown."""

    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, signum, frame):
        print(f"\n{'='*70}")
        print(f"Received signal {signum} - will save checkpoint and exit gracefully...")
        print(f"{'='*70}")
        self.kill_now = True


def _atomic_write(target_path: Path, content: str) -> None:
    """Write content to a file atomically.

    Writes to a temporary file first, then atomically renames it to the target.
    This prevents corruption if the process crashes during the write.

    Args:
        target_path: The target file path.
        content: The content to write.
    """
    import tempfile
    # Create temp file in the same directory to ensure atomic rename works
    target_dir = target_path.parent
    fd, temp_path = tempfile.mkstemp(dir=target_dir, prefix='.tmp_', suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        # Atomic rename (on POSIX systems)
        os.replace(temp_path, target_path)
    except Exception:
        # Clean up temp file on failure
        try:
            os.unlink(temp_path)
        except OSError:
            pass
        raise


def save_global_step(checkpoint_dir: Path, global_step: int) -> None:
    """Save global_step to a file in the checkpoint directory.

    This ensures accurate resume even when batch size changes between runs.
    Uses atomic write to prevent corruption on crash.

    Args:
        checkpoint_dir: Directory where checkpoint is saved.
        global_step: Current global step count.
    """
    global_step_file = checkpoint_dir / "global_step.txt"
    _atomic_write(global_step_file, str(global_step))


def save_training_stage(checkpoint_dir: Path, stage: str, current_batch_size: int) -> None:
    """Save training stage information to checkpoint directory.

    Uses atomic write to prevent corruption on crash.

    Args:
        checkpoint_dir: Directory where checkpoint is saved.
        stage: Current training stage ("linear_probe" or "finetune").
        current_batch_size: Current batch size being used.
    """
    import json
    stage_file = checkpoint_dir / "training_stage.json"
    stage_info = {
        "stage": stage,
        "batch_size": current_batch_size,
    }
    _atomic_write(stage_file, json.dumps(stage_info))


def load_training_stage(checkpoint_dir: Path) -> dict:
    """Load training stage information from checkpoint directory.

    Args:
        checkpoint_dir: Directory where checkpoint is saved.

    Returns:
        Dict with 'stage' and 'batch_size', or None if not found.
    """
    import json
    stage_file = checkpoint_dir / "training_stage.json"
    if stage_file.exists():
        try:
            with open(stage_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return None


def load_global_step(checkpoint_dir: Path) -> int:
    """Load global_step from a checkpoint directory.

    Args:
        checkpoint_dir: Directory where checkpoint is saved.

    Returns:
        Global step count, or 0 if not found.
    """
    global_step_file = checkpoint_dir / "global_step.txt"
    if global_step_file.exists():
        try:
            with open(global_step_file, 'r') as f:
                return int(f.read().strip())
        except (ValueError, IOError):
            pass
    return 0


def save_metrics_json(
    checkpoint_dir: Path,
    metrics: dict,
    global_step: int,
    stage: str,
    epoch: int,
    best_for: Optional[List[str]] = None,
) -> None:
    """Save metrics.json to checkpoint directory.

    Args:
        checkpoint_dir: Directory where checkpoint is saved.
        metrics: Dict of metric names to values.
        global_step: Current global step count.
        stage: Current training stage ("linear_probe" or "finetune").
        epoch: Current epoch number.
        best_for: List of metric names this checkpoint is best for.
    """
    import json
    from datetime import datetime

    metrics_data = {
        "global_step": global_step,
        "training_stage": stage,
        "epoch": epoch,
        "timestamp": datetime.now().isoformat(),
        "metrics": metrics,
        "best_for": best_for or [],
    }

    metrics_file = checkpoint_dir / "metrics.json"
    _atomic_write(metrics_file, json.dumps(metrics_data, indent=2))


def load_metrics_json(checkpoint_dir: Path) -> Optional[dict]:
    """Load metrics.json from checkpoint directory.

    Args:
        checkpoint_dir: Directory where checkpoint is saved.

    Returns:
        Dict with metrics data, or None if not found or invalid.
    """
    import json

    metrics_file = checkpoint_dir / "metrics.json"
    if metrics_file.exists():
        try:
            with open(metrics_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return None


def update_best_model_symlink(ckpt_dir: Path, source_name: str) -> None:
    """Update best_model symlink to point to source checkpoint.

    Args:
        ckpt_dir: Parent checkpoint directory.
        source_name: Name of the checkpoint directory to link to.
    """
    best_model_path = ckpt_dir / "best_model"
    source_path = ckpt_dir / source_name

    if not source_path.exists():
        return

    # Remove existing symlink/directory if present
    if best_model_path.is_symlink():
        best_model_path.unlink()
    elif best_model_path.exists():
        shutil.rmtree(best_model_path)

    # Create relative symlink
    best_model_path.symlink_to(source_name)


def is_checkpoint_complete(checkpoint_path: Path) -> bool:
    """Check if a checkpoint directory contains all required files.

    Args:
        checkpoint_path: Path to checkpoint directory.

    Returns:
        True if checkpoint appears complete.
    """
    has_model = (checkpoint_path / "model.safetensors").exists() or \
                (checkpoint_path / "pytorch_model.bin").exists()
    has_optimizer = (checkpoint_path / "optimizer.bin").exists()
    return has_model and has_optimizer


def _parse_epoch_from_checkpoint(checkpoint_name: str, steps_per_epoch: int = 1) -> int:
    """Parse epoch number from checkpoint directory name.

    Handles various checkpoint naming conventions:
    - step_X: Returns X // steps_per_epoch
    - epoch_X: Returns X
    - interrupted_epoch_X: Returns X
    - Other patterns: Returns 0

    Args:
        checkpoint_name: Name of the checkpoint directory.
        steps_per_epoch: Number of steps per epoch (for step-based checkpoints).

    Returns:
        Epoch number extracted from the checkpoint name.
    """
    if 'step_' in checkpoint_name:
        try:
            step_num = int(checkpoint_name.split('step_')[1].split('_')[0])
            return step_num // max(steps_per_epoch, 1)
        except (ValueError, IndexError):
            return 0
    elif 'epoch_' in checkpoint_name:
        try:
            return int(checkpoint_name.split('epoch_')[-1].split('_')[0])
        except (ValueError, IndexError):
            return 0
    return 0


def find_latest_checkpoint(ckpt_dir: Path, accelerator) -> Path:
    """Find the latest valid checkpoint in a directory.

    Args:
        ckpt_dir: Checkpoint directory.
        accelerator: Accelerator instance for printing.

    Returns:
        Path to latest checkpoint, or None if none found.
    """
    all_checkpoints = []

    for pattern in ["epoch_*", "interrupted_epoch_*", "step_*"]:
        for p in ckpt_dir.glob(pattern):
            if is_checkpoint_complete(p):
                try:
                    mtime = p.stat().st_mtime
                    all_checkpoints.append((mtime, p))
                except OSError:
                    pass

    if all_checkpoints:
        all_checkpoints.sort(reverse=True)
        return all_checkpoints[0][1]
    return None


def detect_checkpoint_stage(args, config, ckpt_dir: Path, accelerator, linear_probe_steps: int) -> Optional[dict]:
    """Detect the training stage from checkpoint before loading.

    This must be called BEFORE building the optimizer to ensure the correct
    optimizer configuration is used.

    Args:
        args: Command-line arguments.
        config: Configuration object.
        ckpt_dir: Checkpoint directory.
        accelerator: Accelerator instance.
        linear_probe_steps: Number of linear probe steps configured.

    Returns:
        Dict with 'resume_dir', 'stage', 'batch_size', 'global_step', or None if not resuming.
    """
    if not (args.resume or args.checkpoint_path):
        return None

    if args.checkpoint_path:
        resume_dir = Path(args.checkpoint_path)
    else:
        resume_dir = find_latest_checkpoint(ckpt_dir, accelerator)

    if not resume_dir or not resume_dir.exists():
        return None

    accelerator.print(f"\nDetecting checkpoint stage from: {resume_dir}")

    # Load stage info if available
    stage_info = load_training_stage(resume_dir)
    global_step = load_global_step(resume_dir)

    if stage_info:
        stage = stage_info.get('stage', 'linear_probe')
        batch_size = stage_info.get('batch_size', config.training.batch_size)
        accelerator.print(f"  Found saved stage info: stage={stage}, batch_size={batch_size}")
    else:
        # Infer stage from global_step if stage info not saved
        if global_step >= linear_probe_steps:
            stage = 'finetune'
        else:
            stage = 'linear_probe'
        batch_size = config.training.batch_size
        accelerator.print(f"  Inferred stage from global_step={global_step}: {stage}")

    if config.get('training.force_finetune_on_resume', False):
        if stage != 'finetune':
            accelerator.print("  Forcing resume stage to finetune (override checkpoint stage)")
        stage = 'finetune'

    return {
        'resume_dir': resume_dir,
        'stage': stage,
        'batch_size': batch_size,
        'global_step': global_step,
        'stage_info_found': stage_info is not None,
    }


def _is_optimizer_state_mismatch(error: Exception) -> bool:
    message = str(error).lower()
    return (
        "parameter group" in message and "optimizer" in message
    ) or "size of optimizer's group" in message


def _load_checkpoint_model_only(resume_dir: Path, accelerator: Accelerator) -> None:
    """Load model (and RNG) state while skipping optimizer/scheduler states."""
    if accelerator.distributed_type in {
        DistributedType.DEEPSPEED,
        DistributedType.FSDP,
        DistributedType.MEGATRON_LM,
    }:
        raise RuntimeError(
            "Optimizer-mismatch fallback is not supported for DeepSpeed/FSDP/Megatron-LM."
        )

    from accelerate.checkpointing import load_accelerator_state

    models = list(accelerator._models)
    if not models:
        raise RuntimeError("No models registered with Accelerator to load checkpoint state.")

    map_location = "cpu"
    if accelerator.num_processes > 1 and accelerator.multi_device and accelerator.distributed_type != DistributedType.MULTI_XPU:
        map_location = "on_device"

    override_attributes = load_accelerator_state(
        resume_dir,
        models,
        [],
        [],
        [],
        accelerator.state.process_index,
        accelerator.scaler,
        map_location,
        None,
    )
    if "step" in override_attributes:
        accelerator.step = override_attributes["step"]


def _prepare_alt_optimizer_checkpoint(resume_dir: Path) -> Optional[Path]:
    optimizer_alt = resume_dir / "optimizer_1.bin"
    scheduler_alt = resume_dir / "scheduler_1.bin"
    if not (optimizer_alt.exists() and scheduler_alt.exists()):
        return None

    temp_dir = Path(tempfile.mkdtemp(prefix=".tmp_resume_", dir=resume_dir))
    for item in resume_dir.iterdir():
        if item.name in {"optimizer.bin", "scheduler.bin"}:
            continue
        if item.is_dir():
            continue
        dest = temp_dir / item.name
        try:
            os.symlink(item, dest)
        except OSError:
            shutil.copy2(item, dest)

    for src_name, dest_name in [("optimizer_1.bin", "optimizer.bin"), ("scheduler_1.bin", "scheduler.bin")]:
        src = resume_dir / src_name
        dest = temp_dir / dest_name
        if src.exists():
            try:
                os.symlink(src, dest)
            except OSError:
                shutil.copy2(src, dest)

    return temp_dir


def load_checkpoint_state(resume_info: Optional[dict], accelerator, steps_per_epoch: int) -> Tuple[int, int]:
    """Load checkpoint state after optimizer is properly configured.

    Args:
        resume_info: Dict from detect_checkpoint_stage.
        accelerator: Accelerator instance.
        steps_per_epoch: Number of steps per epoch.

    Returns:
        Tuple of (starting_epoch, resume_global_step).
    """
    if resume_info is None:
        return 0, 0

    resume_dir = resume_info['resume_dir']
    accelerator.print(f"\nLoading checkpoint state from: {resume_dir}")
    try:
        accelerator.load_state(str(resume_dir))
    except (ValueError, RuntimeError) as exc:
        if not _is_optimizer_state_mismatch(exc):
            raise

        if resume_info.get("stage") == "finetune":
            alt_dir = _prepare_alt_optimizer_checkpoint(resume_dir)
            if alt_dir is not None:
                accelerator.print(
                    "  Warning: optimizer state mismatch; retrying with optimizer_1.bin/scheduler_1.bin."
                )
                try:
                    accelerator.load_state(str(alt_dir))
                except (ValueError, RuntimeError) as alt_exc:
                    if _is_optimizer_state_mismatch(alt_exc):
                        accelerator.print(
                            "  Warning: alternate optimizer state mismatch; loading model weights only."
                        )
                        _load_checkpoint_model_only(resume_dir, accelerator)
                    else:
                        raise
                finally:
                    shutil.rmtree(alt_dir, ignore_errors=True)
            else:
                accelerator.print(
                    "  Warning: optimizer state mismatch and no optimizer_1.bin found; "
                    "loading model weights only."
                )
                _load_checkpoint_model_only(resume_dir, accelerator)
        else:
            accelerator.print(
                "  Warning: optimizer state does not match current parameter groups; "
                "loading model weights only."
            )
            _load_checkpoint_model_only(resume_dir, accelerator)

    starting_epoch = _parse_epoch_from_checkpoint(resume_dir.name, steps_per_epoch)

    resume_global_step = resume_info.get('global_step', 0)
    if resume_global_step > 0:
        accelerator.print(f"  Loaded global_step: {resume_global_step}")

    return starting_epoch, resume_global_step


def handle_checkpoint_resume(args, config, ckpt_dir: Path, accelerator, steps_per_epoch: int) -> Tuple[int, int]:
    """Handle checkpoint resumption with direct loading.

    Loads checkpoint state directly without stage detection. For two-stage training
    with automatic stage detection, use detect_checkpoint_stage() + load_checkpoint_state().

    Args:
        args: Command-line arguments.
        config: Configuration object.
        ckpt_dir: Checkpoint directory.
        accelerator: Accelerator instance.
        steps_per_epoch: Number of steps per epoch.

    Returns:
        Tuple of (starting_epoch, resume_global_step).
    """
    if not (args.resume or args.checkpoint_path):
        return 0, 0

    if args.checkpoint_path:
        resume_dir = Path(args.checkpoint_path)
    else:
        resume_dir = find_latest_checkpoint(ckpt_dir, accelerator)

    if not resume_dir or not resume_dir.exists():
        return 0, 0

    accelerator.print(f"\nResuming from checkpoint: {resume_dir}")
    accelerator.load_state(str(resume_dir))

    starting_epoch = _parse_epoch_from_checkpoint(resume_dir.name, steps_per_epoch)
    resume_global_step = load_global_step(resume_dir)

    if resume_global_step > 0:
        accelerator.print(f"  Loaded global_step: {resume_global_step}")

    return starting_epoch, resume_global_step


def initialize_wandb(config, accelerator, ckpt_dir: Path, starting_epoch: int, resume_global_step: int) -> int:
    """Initialize WandB with resume support.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance.
        ckpt_dir: Checkpoint directory.
        starting_epoch: Starting epoch number.
        resume_global_step: Global step from checkpoint (0 if not resuming).

    Returns:
        Step offset to keep W&B logging monotonic when resuming from older checkpoints.
    """
    if not config.logging.use_wandb:
        accelerator.print("WandB logging disabled (use_wandb=False)")
        return 0

    if not accelerator.is_main_process:
        return 0

    wandb_run_id = None
    if starting_epoch > 0 or resume_global_step > 0:
        wandb_id_file = ckpt_dir / "wandb_run_id.txt"
        if wandb_id_file.exists():
            with open(wandb_id_file, 'r') as f:
                wandb_run_id = f.read().strip()
            accelerator.print(f"  Resuming W&B run: {wandb_run_id}")

    wandb_config = {"name": config.logging.wandb_run_name, "config": config.to_dict()}
    if wandb_run_id:
        wandb_config["id"] = wandb_run_id
        wandb_config["resume"] = "allow"

    accelerator.print(f"\nInitializing WandB:")
    accelerator.print(f"  Project: {config.logging.wandb_project}")
    accelerator.print(f"  Run name: {config.logging.wandb_run_name}")

    accelerator.init_trackers(
        project_name=config.logging.wandb_project,
        init_kwargs={"wandb": wandb_config}
    )

    accelerator.print("  WandB initialized successfully!")

    wandb_step_offset = 0
    if wandb_run_id:
        wandb_tracker = accelerator.get_tracker("wandb")
        if wandb_tracker and wandb_tracker.run:
            run_step = getattr(wandb_tracker.run, "step", None)
            if isinstance(run_step, int) and run_step > resume_global_step:
                wandb_step_offset = run_step - resume_global_step
                accelerator.print(
                    f"  Aligning W&B logging steps with offset={wandb_step_offset} (run step={run_step})"
                )

    if wandb_run_id is None:
        wandb_tracker = accelerator.get_tracker("wandb")
        if wandb_tracker and wandb_tracker.run:
            wandb_id_file = ckpt_dir / "wandb_run_id.txt"
            with open(wandb_id_file, 'w') as f:
                f.write(wandb_tracker.run.id)

    return wandb_step_offset


def training_step(
    model: nn.Module,
    batch_data: tuple,
    loss_fn: nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LRScheduler,
    accelerator: Accelerator,
    is_human: bool,
    clip_grad_norm: float,
    batch_manager: Optional["AdaptiveBatchManager"],
    global_step: int,
    train_metrics: List[tuple],
    compute_detailed_metrics: bool = False,
    warmup_steps: int = 0,
) -> Tuple[float, bool, Dict[str, float]]:
    """Execute a single training step.

    Args:
        model: PyTorch model.
        batch_data: Tuple of (inputs, targets) or (inputs, targets, metadata).
        loss_fn: Loss function.
        optimizer: Optimizer.
        scheduler: Learning rate scheduler.
        accelerator: Accelerator instance.
        is_human: Whether using human data.
        clip_grad_norm: Gradient clipping value.
        batch_manager: Batch manager for gradient accumulation.
        global_step: Global step counter (used for accumulation tracking across epochs).
        train_metrics: List of training metrics.
        compute_detailed_metrics: Whether to compute detailed gradient/output metrics.
        warmup_steps: Number of warmup steps (for logging warmup progress).

    Returns:
        Tuple of (loss_value, oom_occurred, metrics_dict).
    """
    if len(batch_data) == 3:
        inputs, targets, metadata = batch_data
    else:
        inputs, targets = batch_data

    accum_steps = batch_manager.gradient_accumulation_steps if batch_manager else 1
    detailed_metrics = {}

    # Only zero gradients on first accumulation step (use global_step for correct tracking across epochs)
    if batch_manager is None or global_step % accum_steps == 0:
        optimizer.zero_grad()

    with accelerator.autocast():
        inputs = inputs.permute(0, 2, 1)  # (batch, seq_len, 4) -> (batch, 4, seq_len)
        outputs = model(inputs, is_human=is_human)

        if isinstance(outputs, tuple):
            outputs = outputs[0]

        outputs = outputs.permute(0, 2, 1)  # (batch, tracks, bins) -> (batch, bins, tracks)
        loss = loss_fn(outputs, targets.float())

        # Update training metrics
        for _, metric in train_metrics:
            metric.update(outputs, targets)

        # Compute detailed output statistics if requested
        if compute_detailed_metrics:
            output_stats = compute_output_statistics(outputs)
            detailed_metrics.update({f'output/{k}': v for k, v in output_stats.items()})

            # Loss statistics
            loss_stats = compute_loss_statistics(loss_fn, outputs, targets)
            detailed_metrics.update({f'loss/{k}': v for k, v in loss_stats.items()})

            # Check for NaN/Inf in loss
            detailed_metrics['loss/is_nan'] = 1.0 if torch.isnan(loss).any() else 0.0
            detailed_metrics['loss/is_inf'] = 1.0 if torch.isinf(loss).any() else 0.0

        if accum_steps > 1:
            loss = loss / accum_steps

    accelerator.backward(loss)

    # Only step optimizer on accumulation boundary (use global_step for correct tracking)
    should_step = batch_manager is None or batch_manager.should_step(global_step)

    grad_norm = None
    if should_step:
        # Compute gradient statistics before clipping (if detailed metrics requested)
        if compute_detailed_metrics:
            # Compute total gradient norm (before clipping)
            total_grad_norm = 0.0
            for p in model.parameters():
                if p.grad is not None:
                    total_grad_norm += p.grad.data.norm(2).item() ** 2
            total_grad_norm = total_grad_norm ** 0.5
            detailed_metrics['grad/norm_before_clip'] = total_grad_norm

            # Layer-wise gradient norms
            layer_grad_norms = compute_layer_gradient_norms(model)
            detailed_metrics['grad/norm_transformer'] = layer_grad_norms['transformer']
            detailed_metrics['grad/norm_conv'] = layer_grad_norms['conv']
            detailed_metrics['grad/norm_head'] = layer_grad_norms['head']
            detailed_metrics['grad/norm_other'] = layer_grad_norms['other']

            # Parameter norm
            detailed_metrics['param/norm'] = compute_parameter_norm(model)

        # Clip gradients and get the actual norm
        if clip_grad_norm > 0:
            grad_norm = accelerator.clip_grad_norm_(model.parameters(), clip_grad_norm)
            if hasattr(grad_norm, 'item'):
                grad_norm = grad_norm.item()

            if compute_detailed_metrics:
                detailed_metrics['grad/norm'] = grad_norm
                detailed_metrics['grad/clipped'] = 1.0 if grad_norm > clip_grad_norm else 0.0
                # Clip ratio: how much the gradients were scaled down
                if grad_norm > clip_grad_norm:
                    detailed_metrics['grad/clip_ratio'] = clip_grad_norm / grad_norm
                else:
                    detailed_metrics['grad/clip_ratio'] = 1.0

        optimizer.step()
        scheduler.step()

        # Warmup progress
        if compute_detailed_metrics and warmup_steps > 0:
            detailed_metrics['lr/warmup_progress'] = min(1.0, global_step / warmup_steps)

    # GPU memory stats (always compute if detailed metrics requested)
    if compute_detailed_metrics:
        gpu_stats = get_gpu_memory_stats()
        detailed_metrics.update({f'memory/{k}': v for k, v in gpu_stats.items()})

    display_loss = loss.item() * accum_steps if accum_steps > 1 else loss.item()
    return display_loss, False, detailed_metrics


def handle_oom(
    batch_manager: "AdaptiveBatchManager",
    optimizer: torch.optim.Optimizer,
    accelerator: Accelerator,
    train_dataset,
    val_dataset,
    num_workers: int,
) -> Tuple[Optional[torch.utils.data.DataLoader], Optional[torch.utils.data.DataLoader], bool]:
    """Handle OOM by reducing batch size and recreating loaders.

    Args:
        batch_manager: Batch manager instance.
        optimizer: Optimizer.
        accelerator: Accelerator instance.
        train_dataset: Training dataset.
        val_dataset: Validation dataset.
        num_workers: Number of workers.

    Returns:
        Tuple of (new_train_loader, new_val_loader, success).
    """
    import time
    from src.data_loaders import create_data_loaders_from_datasets

    optimizer.zero_grad(set_to_none=True)
    torch.cuda.empty_cache()
    time.sleep(batch_manager.recovery_delay)

    accelerator.wait_for_everyone()

    if not batch_manager.reduce_batch_size():
        accelerator.print(f"  ERROR: Cannot reduce batch size below {batch_manager.min_batch_size}")
        return None, None, False

    accelerator.print(f"  Reducing batch size to {batch_manager.current_batch_size}")
    accelerator.print(f"  Gradient accumulation steps: {batch_manager.gradient_accumulation_steps}")
    accelerator.print(f"  Effective batch size: {batch_manager.effective_batch_size}")

    train_loader, val_loader = create_data_loaders_from_datasets(
        train_dataset, val_dataset, batch_manager.current_batch_size, num_workers, accelerator
    )
    train_loader, val_loader = accelerator.prepare(train_loader, val_loader)

    return train_loader, val_loader, True


def transition_to_finetune(
    model: nn.Module,
    accelerator: Accelerator,
    config,
    train_dataset,
    val_dataset,
    num_workers: int,
    finetune_batch_size: Optional[int],
    finetune_lr: float,
    weight_decay: float,
    weight_decay_transformer: float,
    warmup_steps: int,
    finetune_steps: int,
    linear_probe_steps: int,
    is_human: bool,
    batch_manager: Optional["AdaptiveBatchManager"],
    epoch: int,
    num_epochs: int,
    current_train_loader=None,
    finetune_gradient_accumulation_steps: Optional[int] = None,
) -> Tuple[Optional[torch.utils.data.DataLoader], Optional[torch.utils.data.DataLoader], int, int, int]:
    """Transition from linear probing to fine-tuning.

    Args:
        model: PyTorch model.
        accelerator: Accelerator instance.
        config: Configuration object.
        train_dataset: Training dataset.
        val_dataset: Validation dataset.
        num_workers: Number of workers.
        finetune_batch_size: Batch size for fine-tuning.
        finetune_lr: Learning rate for fine-tuning.
        weight_decay: Weight decay for non-transformer layers.
        weight_decay_transformer: Weight decay for transformer layers.
        warmup_steps: Number of warmup steps.
        finetune_steps: Number of fine-tuning steps.
        linear_probe_steps: Number of linear probing steps.
        is_human: Whether using human data.
        batch_manager: Batch manager instance.
        epoch: Current epoch.
        num_epochs: Total number of epochs.
        current_train_loader: Current training loader (used to get actual batch size).
        finetune_gradient_accumulation_steps: Gradient accumulation steps for fine-tuning.
            If provided, maintains effective batch size when physical batch is reduced.

    Returns:
        Tuple of (train_loader, val_loader, steps_per_epoch, total_steps, finetune_steps).
    """
    from src.data_loaders import create_data_loaders_from_datasets

    accelerator.print("\nSTAGE 2: Unfreezing Trunk for Fine-tuning")
    base_model = accelerator.unwrap_model(model)
    set_trainable_layers(base_model, accelerator, train_head_only=False, is_human=is_human)

    # Use current loader length to get accurate steps_per_epoch
    # This accounts for any batch size changes from OOM protection
    if current_train_loader is not None:
        steps_per_epoch = len(current_train_loader)
        current_batch_size = current_train_loader.batch_size if hasattr(current_train_loader, 'batch_size') else config.training.batch_size
    else:
        # Fallback: use batch_manager's current batch size if available
        current_batch_size = batch_manager.current_batch_size if batch_manager else config.training.batch_size
        steps_per_epoch = len(train_dataset) // current_batch_size

    total_steps = linear_probe_steps + finetune_steps

    # Determine the target batch size for finetune stage
    target_batch_size = finetune_batch_size if finetune_batch_size else current_batch_size

    if target_batch_size != current_batch_size:
        accelerator.print(f"  Changing batch size from {current_batch_size} to {target_batch_size}")
        train_loader, val_loader = create_data_loaders_from_datasets(
            train_dataset, val_dataset, target_batch_size, num_workers, accelerator, verbose=True
        )
        train_loader, val_loader = accelerator.prepare(train_loader, val_loader)
        steps_per_epoch = len(train_loader)
        remaining_epochs = num_epochs - epoch
        finetune_steps = steps_per_epoch * remaining_epochs
        total_steps = linear_probe_steps + finetune_steps
        accelerator.print(f"  Updated steps_per_epoch: {steps_per_epoch}")
        accelerator.print(f"  Updated finetune_steps: {finetune_steps}")

        if batch_manager:
            batch_manager.reset_for_stage(
                target_batch_size,
                gradient_accumulation_steps=finetune_gradient_accumulation_steps,
            )
            if finetune_gradient_accumulation_steps and finetune_gradient_accumulation_steps > 1:
                accelerator.print(f"  Gradient accumulation steps: {finetune_gradient_accumulation_steps}")
                accelerator.print(f"  Effective batch size: {batch_manager.effective_batch_size}")
    else:
        train_loader, val_loader = None, None
        # Recalculate finetune_steps based on actual steps_per_epoch
        remaining_epochs = num_epochs - epoch
        finetune_steps = steps_per_epoch * remaining_epochs
        total_steps = linear_probe_steps + finetune_steps
        accelerator.print(f"  Keeping batch size at {current_batch_size}")
        accelerator.print(f"  steps_per_epoch: {steps_per_epoch}")
        accelerator.print(f"  finetune_steps: {finetune_steps}")

        # Still apply gradient accumulation even if batch size unchanged
        if batch_manager and finetune_gradient_accumulation_steps:
            batch_manager.reset_for_stage(
                current_batch_size,
                gradient_accumulation_steps=finetune_gradient_accumulation_steps,
            )
            if finetune_gradient_accumulation_steps > 1:
                accelerator.print(f"  Gradient accumulation steps: {finetune_gradient_accumulation_steps}")
                accelerator.print(f"  Effective batch size: {batch_manager.effective_batch_size}")

    return train_loader, val_loader, steps_per_epoch, total_steps, finetune_steps
