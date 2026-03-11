"""Logging utilities for training progress and configuration."""


def log_training_schedule(
    accelerator,
    linear_probe_steps: int,
    finetune_steps: int,
    steps_per_epoch: int,
    head_lr: float,
    finetune_lr: float,
    batch_size: int,
    finetune_batch_size: int,
    allow_unfreeze: bool
) -> None:
    """Log the training schedule.

    Args:
        accelerator: Accelerator instance for printing.
        linear_probe_steps: Number of linear probing steps.
        finetune_steps: Number of fine-tuning steps.
        steps_per_epoch: Steps per epoch.
        head_lr: Learning rate for head-only training.
        finetune_lr: Learning rate for fine-tuning.
        batch_size: Initial batch size.
        finetune_batch_size: Batch size for fine-tuning.
        allow_unfreeze: Whether unfreezing is allowed.
    """
    accelerator.print("\nStage schedule:")
    if linear_probe_steps > 0:
        linear_probe_epochs = linear_probe_steps / steps_per_epoch
        accelerator.print(f"  Stage 1 (head only): {linear_probe_steps} steps ({linear_probe_epochs:.1f} epochs) at lr={head_lr}")
    if allow_unfreeze and finetune_steps > 0:
        finetune_epochs = finetune_steps / steps_per_epoch
        batch_info = f"batch_size={finetune_batch_size}" if finetune_batch_size else f"batch_size={batch_size}"
        accelerator.print(f"  Stage 2 (full fine-tune): {finetune_steps} steps ({finetune_epochs:.1f} epochs) at lr={finetune_lr}, {batch_info}")
