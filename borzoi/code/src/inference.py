"""Post-training inference utilities."""

import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from accelerate import Accelerator

from .config import Config


def run_post_training_inference(
    config: Config,
    accelerator: Accelerator,
    ckpt_dir: Path,
    track_names: Optional[List[str]] = None,
):
    """Run inference after training completes (main process only).

    Supports multiple inference configurations. If the config's inference section
    is a list, each config is run sequentially. If it's a dict (legacy format),
    it's treated as a single inference config.

    Args:
        config: Training configuration.
        accelerator: Accelerator instance.
        ckpt_dir: Checkpoint directory.
        track_names: Track names from training data (auto-detected if None).
    """
    if not accelerator.is_main_process:
        return

    # Get inference configs as list (handles both dict and list formats)
    inference_configs = config.get_inference_configs()

    if not inference_configs:
        return

    # Filter to configs that have run_after_training=True
    configs_to_run = [cfg for cfg in inference_configs if cfg.get('run_after_training', False)]

    if not configs_to_run:
        return

    accelerator.print("\n" + "=" * 70)
    accelerator.print("Running Post-Training Inference")
    accelerator.print(f"  {len(configs_to_run)} inference configuration(s) to run")
    accelerator.print("=" * 70)

    # Get top-level misha_track for fallback
    top_level_misha = config.get('misha_track', {})

    for idx, inference_cfg in enumerate(configs_to_run):
        config_name = inference_cfg.get('name', f'config_{idx + 1}')

        accelerator.print(f"\n--- Inference [{idx + 1}/{len(configs_to_run)}]: {config_name} ---")

        # Determine checkpoint path (supports metric selection)
        checkpoint_option = inference_cfg.get('checkpoint', 'best_model')
        metric_option = inference_cfg.get('metric')
        checkpoint_path = _resolve_checkpoint_path(ckpt_dir, checkpoint_option)

        if not checkpoint_path.exists():
            accelerator.print(f"Warning: Checkpoint not found at {checkpoint_path}, skipping inference config '{config_name}'")
            continue

        # Merge misha_track: per-inference config overrides top-level
        misha_track = _merge_misha_config(top_level_misha, inference_cfg.get('misha_track'))

        # Build inference command
        cmd = _build_inference_command(
            config, inference_cfg, ckpt_dir, checkpoint_path, track_names, accelerator,
            misha_track=misha_track, metric=metric_option
        )

        output_dir = inference_cfg.get('output_dir') or str(ckpt_dir / 'predictions')
        accelerator.print(f"Command: {' '.join(cmd[:8])}...")
        accelerator.print(f"Output directory: {output_dir}")

        _execute_inference(cmd, accelerator)

    accelerator.print("\n" + "=" * 70)
    accelerator.print("All inference configurations completed")
    accelerator.print("=" * 70)


def _merge_misha_config(
    top_level: Dict[str, Any],
    per_inference: Optional[Dict[str, Any]]
) -> Dict[str, Any]:
    """Merge top-level misha_track with per-inference overrides.

    Uses shallow merge: per-inference values completely override top-level
    values at the same key. This is appropriate since misha_track fields
    are flat (no nested structures).

    Args:
        top_level: Top-level misha_track config (fallback).
        per_inference: Per-inference misha_track config (overrides).

    Returns:
        Merged misha_track config.
    """
    if not per_inference:
        return top_level.copy() if top_level else {}

    # Start with top-level, override with per-inference (shallow merge)
    merged = top_level.copy() if top_level else {}
    merged.update(per_inference)
    return merged


def _resolve_checkpoint_path(ckpt_dir: Path, checkpoint_option: str) -> Path:
    """Resolve checkpoint path from option string.

    Args:
        ckpt_dir: Checkpoint directory.
        checkpoint_option: 'best_model', 'last', or explicit path.

    Returns:
        Resolved checkpoint path.
    """
    if checkpoint_option == 'best_model':
        return ckpt_dir / 'best_model'
    elif checkpoint_option == 'last':
        epoch_ckpts = list(ckpt_dir.glob("epoch_*"))
        if epoch_ckpts:
            epoch_nums = []
            for p in epoch_ckpts:
                try:
                    epoch_num = int(p.name.split('_')[1])
                    epoch_nums.append((epoch_num, p))
                except (ValueError, IndexError):
                    pass
            if epoch_nums:
                epoch_nums.sort(key=lambda x: x[0], reverse=True)
                return epoch_nums[0][1]
        return ckpt_dir / 'best_model'
    else:
        return Path(checkpoint_option)


def _build_inference_command(
    config: Config,
    inference_cfg: dict,
    ckpt_dir: Path,
    checkpoint_path: Path,
    track_names: Optional[List[str]],
    accelerator: Accelerator,
    misha_track: Optional[Dict[str, Any]] = None,
    metric: Optional[str] = None,
) -> List[str]:
    """Build the torchrun command for inference.

    Args:
        config: Training configuration.
        inference_cfg: Inference section of config.
        ckpt_dir: Checkpoint directory.
        checkpoint_path: Path to checkpoint.
        track_names: Track names for output.
        accelerator: Accelerator instance.
        misha_track: Merged misha_track config (per-inference + top-level).
        metric: Optional metric name for checkpoint selection.

    Returns:
        Command as list of strings.
    """
    output_dir = inference_cfg.get('output_dir') or str(ckpt_dir / 'predictions')
    batch_size = inference_cfg.get('batch_size') or config.training.batch_size
    num_workers = inference_cfg.get('num_workers') or config.training.num_workers
    mode = inference_cfg.get('mode', 'genome')
    output_format = inference_cfg.get('output_format', ['bigwig'])
    bigwig_jobs = inference_cfg.get('bigwig_jobs', 4)
    prefix = inference_cfg.get('prefix', '')
    genome_fasta = inference_cfg.get('genome_fasta')

    infer_track_names = inference_cfg.get('track_names') or track_names

    # Determine number of GPUs
    cuda_devices = os.environ.get('CUDA_VISIBLE_DEVICES', '')
    if cuda_devices:
        num_gpus = len(cuda_devices.split(','))
    else:
        num_gpus = accelerator.num_processes

    # Build command
    infer_script = Path(__file__).parent.parent / 'infer_borzoi_pytorch.py'
    config_path = ckpt_dir / 'config.yaml'

    # Get inference config name (required to prevent running all configs)
    inference_name = inference_cfg.get('name')

    cmd = [
        'torchrun',
        f'--nproc_per_node={num_gpus}',
        '--master_port=29516',
        str(infer_script),
        '--config', str(config_path),
        '--checkpoint', str(checkpoint_path),
        '--output_dir', output_dir,
        '--mode', mode,
        '--batch_size', str(batch_size),
        '--num_workers', str(num_workers),
        '--bigwig_jobs', str(bigwig_jobs),
        '--no_run_all',  # Only run the specific config we're targeting
    ]

    # Add inference_name to select the specific inference config
    if inference_name:
        cmd.extend(['--inference_name', inference_name])

    # Add metric for checkpoint selection if specified
    if metric:
        cmd.extend(['--metric', metric])

    # Add genome_fasta override if specified
    if genome_fasta:
        cmd.extend(['--genome_fasta', genome_fasta])

    # Add prefix if specified
    if prefix:
        cmd.extend(['--prefix', prefix])

    # Add output format(s)
    if isinstance(output_format, list):
        for fmt in output_format:
            cmd.extend(['--output_format', fmt])
    else:
        cmd.extend(['--output_format', output_format])

    # Add track names
    if infer_track_names:
        cmd.extend(['--track_names'] + list(infer_track_names))

    # Mode-specific arguments
    if mode == 'chromosomes' and inference_cfg.get('chromosomes'):
        cmd.extend(['--chromosomes'] + inference_cfg['chromosomes'])
    elif mode == 'bed' and inference_cfg.get('bed_file'):
        cmd.extend(['--bed_file', inference_cfg['bed_file']])
    elif mode == 'genome' and inference_cfg.get('exclude_chroms'):
        cmd.extend(['--exclude_chroms'] + inference_cfg['exclude_chroms'])

    # Optional flags
    if not inference_cfg.get('use_rc_average', True):
        cmd.append('--no_rc_average')
    if not inference_cfg.get('mixed_precision', True):
        cmd.append('--no_mixed_precision')
    if inference_cfg.get('stride'):
        cmd.extend(['--stride', str(inference_cfg['stride'])])

    # Misha track parameters (from merged config)
    if misha_track and misha_track.get('enabled', False):
        if misha_track.get('groot'):
            cmd.extend(['--misha_groot', str(misha_track['groot'])])
        if misha_track.get('track_prefix'):
            cmd.extend(['--misha_track_prefix', str(misha_track['track_prefix'])])
        suffix = misha_track.get('suffix') or misha_track.get('track_suffix')
        if suffix:
            cmd.extend(['--misha_suffix', str(suffix)])
        if misha_track.get('binsize'):
            cmd.extend(['--misha_binsize', str(misha_track['binsize'])])

    return cmd


def _execute_inference(cmd: List[str], accelerator: Accelerator) -> None:
    """Execute the inference command.

    Args:
        cmd: Command as list of strings.
        accelerator: Accelerator instance for logging.
    """
    try:
        env = os.environ.copy()
        # Remove accelerate-specific env vars that might interfere
        for key in list(env.keys()):
            if key.startswith('ACCELERATE_') or key in ['RANK', 'WORLD_SIZE', 'LOCAL_RANK']:
                del env[key]

        subprocess.run(cmd, check=True, env=env)
        accelerator.print("Inference completed successfully!")
    except subprocess.CalledProcessError as e:
        accelerator.print(f"Inference failed with exit code {e.returncode}")
    except FileNotFoundError:
        accelerator.print("Error: torchrun not found. Please ensure PyTorch is installed correctly.")
    except Exception as e:
        accelerator.print(f"Inference failed: {e}")
