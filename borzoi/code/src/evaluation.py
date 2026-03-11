"""Evaluation utilities for Borzoi fine-tuning."""

import inspect
import warnings
from typing import Dict, List, Optional

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import tqdm
from accelerate import Accelerator


def _center_crop_bins(tensor: torch.Tensor, target_bins: int) -> torch.Tensor:
    """Center-crop tensor along the bins dimension.

    Args:
        tensor: Tensor of shape (batch, bins, tracks)
        target_bins: Number of bins to keep from center

    Returns:
        Center-cropped tensor of shape (batch, target_bins, tracks)
    """
    current_bins = tensor.shape[1]
    if target_bins >= current_bins:
        return tensor

    start = (current_bins - target_bins) // 2
    end = start + target_bins
    return tensor[:, start:end, :]


def _metric_supports_metadata(metric: nn.Module) -> bool:
    """Check if a metric's update method accepts a metadata parameter.

    Uses instance attribute caching to avoid repeated signature introspection.
    This is safe for long-running jobs since the cache is stored on the metric
    object itself and will be garbage collected with it.

    Args:
        metric: The metric module to check.

    Returns:
        True if the metric's update method accepts 'metadata' parameter.
    """
    # Check for cached result on the metric instance itself
    if hasattr(metric, '_supports_metadata_cached'):
        return metric._supports_metadata_cached

    # Introspect signature
    if hasattr(metric, 'update'):
        sig = inspect.signature(metric.update)
        supports_metadata = 'metadata' in sig.parameters
    else:
        supports_metadata = False

    # Cache result on the metric instance (safe - tied to object lifetime)
    try:
        metric._supports_metadata_cached = supports_metadata
    except (AttributeError, TypeError):
        # Some modules don't allow setting attributes; skip caching
        pass

    return supports_metadata


def evaluate(
    accelerator: Accelerator,
    model: nn.Module,
    val_loader: DataLoader,
    loss_fn: nn.Module,
    metrics: List[tuple],
    is_human: bool,
    val_center_crop_bins: Optional[int] = None,
) -> Dict[str, float]:
    """Run evaluation on validation set.

    Args:
        accelerator: Accelerator instance.
        model: Model to evaluate.
        val_loader: Validation data loader.
        loss_fn: Loss function.
        metrics: List of (name, metric) tuples.
        is_human: Whether to route through human head.
        val_center_crop_bins: If provided, center-crop predictions and targets
            to this many bins before computing loss and metrics. This allows
            comparing models with different native pred_len using the same
            effective pred_len for evaluation.

    Returns:
        Dictionary of metric names to values. Returns NaN values if validation
        loader is empty.
    """
    model.eval()

    # Handle empty validation loader
    if len(val_loader) == 0:
        warnings.warn(
            "Validation loader is empty - skipping evaluation. "
            "Check your validation data configuration.",
            RuntimeWarning
        )
        results = {"val_loss": float('nan')}
        for name, _ in metrics:
            results[f"val_{name}"] = float('nan')
        return results

    # Reset metrics
    for _, metric in metrics:
        metric.reset()

    total_loss = 0.0
    num_batches = 0
    num_samples = 0

    with torch.no_grad():
        for batch_data in tqdm.tqdm(val_loader, desc="Evaluating", disable=not accelerator.is_main_process):
            # Handle both 2-tuple (backward compat) and 3-tuple (with metadata)
            if len(batch_data) == 3:
                inputs, targets, metadata = batch_data
            else:
                inputs, targets = batch_data
                metadata = None

            with accelerator.autocast():
                # Transpose input: (batch, seq_len, 4) -> (batch, 4, seq_len)
                inputs = inputs.permute(0, 2, 1)

                outputs = model(inputs, is_human=is_human)

                # Handle different output formats
                if isinstance(outputs, tuple):
                    outputs = outputs[0]

                # Model outputs (batch, tracks, bins); transpose to (batch, bins, tracks)
                outputs = outputs.permute(0, 2, 1)

                # Center-crop if validation pred_len differs from model pred_len
                if val_center_crop_bins is not None:
                    outputs = _center_crop_bins(outputs, val_center_crop_bins)
                    targets = _center_crop_bins(targets, val_center_crop_bins)

                loss = loss_fn(outputs, targets.float())

            total_loss += loss.item()
            num_batches += 1
            num_samples += targets.shape[0]

            # Update metrics - pass metadata if available (with cached signature check)
            for _, metric in metrics:
                if metadata is not None and _metric_supports_metadata(metric):
                    metric.update(outputs, targets, metadata=metadata)
                else:
                    metric.update(outputs, targets)

    # Reduce loss/metrics across processes
    loss_tensor = torch.tensor(total_loss, device=accelerator.device)
    batches_tensor = torch.tensor(num_batches, device=accelerator.device)
    samples_tensor = torch.tensor(num_samples, device=accelerator.device)
    loss_tensor = accelerator.reduce(loss_tensor, reduction="sum")
    batches_tensor = accelerator.reduce(batches_tensor, reduction="sum")
    samples_tensor = accelerator.reduce(samples_tensor, reduction="sum")

    results = {"val_loss": (loss_tensor / torch.clamp(batches_tensor, min=1)).item(), "val_num_batches": batches_tensor.item(), "val_num_samples": samples_tensor.item()}

    for name, metric in metrics:
        value = metric.compute()
        # Handle both scalar and dict returns
        if isinstance(value, dict):
            for key, val in value.items():
                value_tensor = torch.tensor(val, device=accelerator.device)
                value_tensor = accelerator.reduce(value_tensor, reduction="mean")
                results[f"val_{key}"] = value_tensor.item()
        else:
            value_tensor = torch.tensor(value, device=accelerator.device)
            value_tensor = accelerator.reduce(value_tensor, reduction="mean")
            results[f"val_{name}"] = value_tensor.item()

    model.train()
    return results
