"""
PyTorch loss implementations for Borzoi fine-tuning.
All losses return scalar values for proper distributed training.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import warnings
from typing import Optional

# =============================================================================
# Numerical Constants
# =============================================================================

# Numerical stability constant for division and logarithm operations.
# Standard in deep learning: prevents log(0)=-inf, division by zero,
# while being small enough to not affect gradient magnitude.
EPSILON = 1e-8

# Maximum value for clipping predictions and targets.
# Handles extreme genomic coverage values and prevents numerical overflow.
# Typical bulk ChIP-seq/ATAC-seq coverage rarely exceeds 1e5.
MAX_COVERAGE_VAL = 1e6

# Maximum value for log-space operations (post-log1p transformation).
# log1p(3e6) ≈ 15, typical max for bulk genomic data after squash transform.
MAX_LOG_SPACE_VAL = 15.0

# Clipping range for raw (non-log) MSE operations.
# Wider range than MAX_COVERAGE_VAL since MSE can handle larger values
# without the overflow risk of log/exp operations.
MAX_RAW_MSE_VAL = 1e9

# Default weight for multinomial term in Poisson-multinomial loss.
# Empirically chosen to balance magnitude and shape learning.
# Higher values emphasize profile shape over total magnitude.
POISSON_MULTINOMIAL_ALPHA = 5.0

# Global flag to enable/disable NaN checking in loss functions.
# Set to False in production for ~5-10% performance improvement.
# Can be configured via set_nan_checking_enabled() function.
_NAN_CHECKING_ENABLED = True


def set_nan_checking_enabled(enabled: bool) -> None:
    """Enable or disable NaN checking in loss functions.

    Disabling NaN checking can provide ~5-10% performance improvement
    in production training, but should be kept enabled during debugging.

    Args:
        enabled: If True, check for NaN/Inf in loss values and warn.
                If False, skip checking for performance.
    """
    global _NAN_CHECKING_ENABLED
    _NAN_CHECKING_ENABLED = enabled


def is_nan_checking_enabled() -> bool:
    """Return whether NaN checking is currently enabled."""
    return _NAN_CHECKING_ENABLED


def _check_for_nan(loss: torch.Tensor, loss_name: str) -> torch.Tensor:
    """Check for NaN/Inf in loss and warn if found.

    Warns on non-finite values but lets them propagate to trigger
    PyTorch's anomaly detection rather than silently masking issues.

    This check can be disabled globally for performance via
    set_nan_checking_enabled(False).

    Args:
        loss: Loss tensor to check
        loss_name: Name of the loss for warning message

    Returns:
        Original loss tensor (unchanged)
    """
    if _NAN_CHECKING_ENABLED and not torch.isfinite(loss).all():
        warnings.warn(
            f"{loss_name}: Non-finite values detected in loss! "
            f"This indicates numerical instability (NaN or Inf). "
            f"Loss value: {loss.item() if loss.numel() == 1 else loss}",
            RuntimeWarning,
            stacklevel=3
        )
    return loss


def _center_crop(
    tensor: torch.Tensor,
    crop_bins: Optional[int]
) -> torch.Tensor:
    """Center-crop bins dimension if requested."""
    if crop_bins is None:
        return tensor

    if tensor.shape[1] < crop_bins:
        raise ValueError(
            f"Cannot center-crop to {crop_bins} bins; tensor only has {tensor.shape[1]} bins"
        )

    trim = tensor.shape[1] - crop_bins
    start = trim // 2
    end = start + crop_bins
    return tensor[:, start:end]


class CenterCropLoss(nn.Module):
    """Wrap another loss to crop predictions/targets before computing loss."""

    def __init__(self, base_loss: nn.Module, crop_bins: Optional[int] = None):
        super().__init__()
        self.base_loss = base_loss
        self.crop_bins = crop_bins

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        y_pred = _center_crop(y_pred, self.crop_bins)
        y_true = _center_crop(y_true, self.crop_bins)
        return self.base_loss(y_pred, y_true)


def poisson_multinomial_torch(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    alpha: float = POISSON_MULTINOMIAL_ALPHA,
    epsilon: float = EPSILON,
    max_val: float = MAX_COVERAGE_VAL,
    normalize_multinomial: bool = False
) -> torch.Tensor:
    """
    Combined Poisson + Multinomial profile loss.

    This is the original Borzoi loss function. Combines:
    - Poisson term: Models total signal magnitude per track
    - Multinomial term: Models signal distribution/shape across bins

    Args:
        y_pred: Predictions of shape (batch, bins, tracks)
        y_true: Targets of shape (batch, bins, tracks)
        alpha: Weight for multinomial term (see POISSON_MULTINOMIAL_ALPHA)
        epsilon: Numerical stability constant (see EPSILON)
        max_val: Maximum value for clipping (see MAX_COVERAGE_VAL)
        normalize_multinomial: If True, normalize multinomial term by number of
            bins to prevent high-coverage regions from dominating the loss.
            Default False for backward compatibility.

    Returns:
        Scalar loss value

    Note on asymmetric clipping:
        y_true is clamped to [0, max_val] while y_pred is clamped to [epsilon, max_val].
        This asymmetry is intentional:
        - Targets (y_true) can legitimately be zero (no coverage at a position)
        - Predictions (y_pred) are clamped to a minimum of epsilon to prevent:
          1. log(0) = -inf in the Poisson term: log(y_pred + epsilon)
          2. Division by zero in the multinomial term: y_pred / total_pred
        Setting y_pred minimum to epsilon (1e-8) has negligible effect on loss values
        but ensures numerical stability.
    """
    # Asymmetric clipping: targets can be 0, predictions need minimum epsilon
    # to avoid log(0) and division by zero (see docstring)
    y_true = torch.clamp(y_true, 0.0, max_val)
    y_pred = torch.clamp(y_pred, epsilon, max_val)

    num_bins = y_true.shape[1]

    # Poisson term on total counts per track
    total_true = y_true.sum(dim=1)  # (batch, tracks)
    total_pred = y_pred.sum(dim=1)  # (batch, tracks)

    poisson_term = total_pred - total_true * torch.log(total_pred + epsilon)
    poisson_term = poisson_term.sum(dim=1)  # (batch,)

    # Multinomial term on shape/profile
    # Reuse total_true/total_pred for normalization to avoid redundant sum
    pred_probs = y_pred / (total_pred.unsqueeze(1) + epsilon)
    true_probs = y_true / (total_true.unsqueeze(1) + epsilon)

    multinomial_term = -true_probs * torch.log(pred_probs + epsilon)
    multinomial_term = multinomial_term.sum(dim=[1, 2])  # (batch,)

    # Reuse total_true for total counts (sum over tracks)
    total_counts = total_true.sum(dim=1)  # (batch,)
    multinomial_term = multinomial_term * total_counts

    # Optional normalization to prevent high-coverage regions from dominating
    if normalize_multinomial:
        # Normalize by number of bins to make per-bin contribution equal
        multinomial_term = multinomial_term / num_bins
        # Also normalize Poisson term for consistency
        poisson_term = poisson_term / num_bins

    total_loss = alpha * multinomial_term + poisson_term

    return _check_for_nan(total_loss.mean(), "PoissonMultinomial")


def poisson_loss(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    epsilon: float = EPSILON,
    max_val: float = MAX_COVERAGE_VAL
) -> torch.Tensor:
    """
    Poisson negative log-likelihood loss.

    Args:
        y_pred: Predictions of shape (batch, bins, tracks)
        y_true: Targets of shape (batch, bins, tracks)
        epsilon: Numerical stability constant (see EPSILON)
        max_val: Maximum value for clipping (see MAX_COVERAGE_VAL)

    Returns:
        Scalar loss value
    """
    y_true = torch.clamp(y_true, 0.0, max_val)
    y_pred = torch.clamp(y_pred, epsilon, max_val)

    loss = y_pred - y_true * torch.log(y_pred)

    return _check_for_nan(loss.mean(), "Poisson")


def mse_loss(y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
    """
    Mean Squared Error loss.

    Note: Raw MSE on genomic data often leads to exploding gradients.
    Consider using mse_log_loss or cosine_mse_log_loss instead.
    """
    y_true = torch.clamp(y_true, -MAX_RAW_MSE_VAL, MAX_RAW_MSE_VAL)
    y_pred = torch.clamp(y_pred, -MAX_RAW_MSE_VAL, MAX_RAW_MSE_VAL)

    squared_diff = (y_true - y_pred) ** 2

    return _check_for_nan(squared_diff.mean(), "MSE")


def mse_log_loss(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    epsilon: float = EPSILON
) -> torch.Tensor:
    """
    MSE performed on log-transformed data (log1p).
    Much more stable than standard MSE for genomic tracks.
    """
    y_true = torch.clamp(y_true, 0.0, MAX_RAW_MSE_VAL)
    y_pred = torch.clamp(y_pred, epsilon, MAX_RAW_MSE_VAL)

    target_log = torch.log1p(y_true)
    pred_log = torch.log1p(y_pred)

    squared_diff = (target_log - pred_log) ** 2

    return _check_for_nan(squared_diff.mean(), "MSELog")


def cosine_mse_log_loss(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    mse_weight: float = 1.0,
    cosine_weight: float = 1.0,
    epsilon: float = EPSILON
) -> torch.Tensor:
    """
    Robust Cosine + MSE Loss in Log Space (CRESTED-style).

    Equivalent to CRESTED's CosineMSELogLoss.

    Fixes gradient explosion by:
    1. Compressing inputs with log1p(x)
    2. Combining MSE (magnitude) with Cosine Distance (shape)
    3. Clipping inputs to finite ranges
    """
    # Clip for numerical stability
    y_true = torch.clamp(y_true, 0.0, MAX_RAW_MSE_VAL)
    y_pred = torch.clamp(y_pred, epsilon, MAX_RAW_MSE_VAL)
    
    # Transform to log space
    target_log = torch.log1p(y_true)
    pred_log = torch.log1p(y_pred)
    
    # MSE component (in log space)
    mse = ((target_log - pred_log) ** 2).mean()
    
    # Cosine component - computed over flattened bins and tracks
    # Shape: (batch, bins, tracks)
    
    # Dot product over bins and tracks
    dot_product = (target_log * pred_log).sum(dim=[1, 2])  # (batch,)
    
    # L2 norms
    target_norm = torch.sqrt((target_log ** 2).sum(dim=[1, 2]) + epsilon)  # (batch,)
    pred_norm = torch.sqrt((pred_log ** 2).sum(dim=[1, 2]) + epsilon)  # (batch,)
    
    # Cosine similarity per sample
    cosine_sim = dot_product / (target_norm * pred_norm + epsilon)
    
    # Convert to distance (1 - similarity)
    cosine_dist = 1.0 - cosine_sim
    cosine_loss = cosine_dist.mean()
    
    # Combine
    total_loss = mse_weight * mse + cosine_weight * cosine_loss

    return _check_for_nan(total_loss, "CosineMSELog")


def cosine_mse_raw_loss(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    mse_weight: float = 1.0,
    cosine_weight: float = 1.0,
    epsilon: float = EPSILON,
    max_val: float = MAX_LOG_SPACE_VAL,
) -> torch.Tensor:
    """
    Cosine + MSE Loss for data that is already log-transformed or normalized.

    Args:
        y_pred: Predictions of shape (batch, bins, tracks).
        y_true: Targets of shape (batch, bins, tracks).
        mse_weight: Weight for MSE component.
        cosine_weight: Weight for cosine distance component.
        epsilon: Numerical stability constant (see EPSILON).
        max_val: Maximum value for clipping (see MAX_LOG_SPACE_VAL).

    Returns:
        Scalar loss value.
    """
    # Clip for safety - max_val corresponds to reasonable range for log-transformed data
    y_true = torch.clamp(y_true, 0.0, max_val)
    y_pred = torch.clamp(y_pred, 0.0, max_val)
    
    # MSE component
    mse = ((y_true - y_pred) ** 2).mean()
    
    # Cosine component
    dot_product = (y_true * y_pred).sum(dim=[1, 2])
    true_norm = torch.sqrt((y_true ** 2).sum(dim=[1, 2]) + epsilon)
    pred_norm = torch.sqrt((y_pred ** 2).sum(dim=[1, 2]) + epsilon)
    
    cosine_sim = dot_product / (true_norm * pred_norm + epsilon)
    cosine_dist = 1.0 - cosine_sim
    cosine_loss = cosine_dist.mean()

    total_loss = mse_weight * mse + cosine_weight * cosine_loss

    return _check_for_nan(total_loss, "CosineMSERaw")


def pearson_loss(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    epsilon: float = EPSILON
) -> torch.Tensor:
    """
    Negative Pearson correlation as loss (1 - correlation).
    """
    y_true = torch.clamp(y_true, -MAX_RAW_MSE_VAL, MAX_RAW_MSE_VAL)
    y_pred = torch.clamp(y_pred, -MAX_RAW_MSE_VAL, MAX_RAW_MSE_VAL)
    
    batch_size = y_true.shape[0]
    num_tracks = y_true.shape[2]
    
    # Flatten to (batch * tracks, bins)
    y_true_flat = y_true.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
    y_pred_flat = y_pred.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
    
    # Center the data
    y_true_centered = y_true_flat - y_true_flat.mean(dim=-1, keepdim=True)
    y_pred_centered = y_pred_flat - y_pred_flat.mean(dim=-1, keepdim=True)
    
    # Compute correlation with robust denominator handling
    covariance = (y_true_centered * y_pred_centered).sum(dim=-1)
    true_std = torch.sqrt((y_true_centered ** 2).sum(dim=-1) + epsilon)
    pred_std = torch.sqrt((y_pred_centered ** 2).sum(dim=-1) + epsilon)

    # Robust division: if both stds are near zero, set correlation to 0
    # (constant tracks have undefined correlation - treat as uncorrelated)
    denom = true_std * pred_std
    # Where denom is too small, use 1.0 to avoid division instability
    # This effectively sets correlation to covariance, which will be ~0 for constant tracks
    safe_denom = torch.where(denom > epsilon, denom, torch.ones_like(denom))
    correlation = covariance / safe_denom
    # Mask out invalid correlations (where original denom was too small)
    correlation = torch.where(denom > epsilon, correlation, torch.zeros_like(correlation))
    correlation = torch.clamp(correlation, -1.0, 1.0)

    loss = 1.0 - correlation

    return _check_for_nan(loss.mean(), "Pearson")


class PoissonMultinomialLoss(nn.Module):
    """Module wrapper for poisson_multinomial_torch loss."""

    def __init__(
        self,
        alpha: float = POISSON_MULTINOMIAL_ALPHA,
        epsilon: float = EPSILON,
        normalize_multinomial: bool = False
    ):
        super().__init__()
        self.alpha = alpha
        self.epsilon = epsilon
        self.normalize_multinomial = normalize_multinomial

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        return poisson_multinomial_torch(
            y_pred, y_true, self.alpha, self.epsilon,
            normalize_multinomial=self.normalize_multinomial
        )


class PoissonLoss(nn.Module):
    """Module wrapper for Poisson loss."""

    def __init__(self, epsilon: float = EPSILON):
        super().__init__()
        self.epsilon = epsilon

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        return poisson_loss(y_pred, y_true, self.epsilon)


class MSELogLoss(nn.Module):
    """Module wrapper for MSE log loss."""

    def __init__(self, epsilon: float = EPSILON):
        super().__init__()
        self.epsilon = epsilon

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        return mse_log_loss(y_pred, y_true, self.epsilon)


class CosineMSELogLoss(nn.Module):
    """Module wrapper for Cosine + MSE log loss (CRESTED-style)."""

    def __init__(self, mse_weight: float = 1.0, cosine_weight: float = 1.0, epsilon: float = EPSILON):
        super().__init__()
        self.mse_weight = mse_weight
        self.cosine_weight = cosine_weight
        self.epsilon = epsilon

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        return cosine_mse_log_loss(y_pred, y_true, self.mse_weight, self.cosine_weight, self.epsilon)


class CosineMSERawLoss(nn.Module):
    """Module wrapper for Cosine + MSE raw loss."""

    def __init__(self, mse_weight: float = 1.0, cosine_weight: float = 1.0, epsilon: float = EPSILON):
        super().__init__()
        self.mse_weight = mse_weight
        self.cosine_weight = cosine_weight
        self.epsilon = epsilon

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        return cosine_mse_raw_loss(y_pred, y_true, self.mse_weight, self.cosine_weight, self.epsilon)


class PearsonLoss(nn.Module):
    """Module wrapper for Pearson loss."""

    def __init__(self, epsilon: float = EPSILON):
        super().__init__()
        self.epsilon = epsilon

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        return pearson_loss(y_pred, y_true, self.epsilon)


class MaskedLoss(nn.Module):
    """Wrap a loss to ignore positions where targets are NaN.

    Useful for sparse data (e.g. methylation, scATAC) where many positions
    are unobserved.
    """

    def __init__(self, base_loss: nn.Module):
        super().__init__()
        self.base_loss = base_loss

    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        mask = torch.isfinite(y_true)
        if not mask.any():
            # If everything is masked, return 0 to avoid NaN gradients
            return y_pred.sum() * 0.0

        # Create masked versions
        # We replace NaNs with 0 in the tensors passed to the base loss,
        # BUT we expect the base loss to handle them or we'll need to
        # use a loss that supports weights/masks.
        # For MSE/Poisson/Cosine, we can often just mask the per-element loss.

        # If base_loss is a known standard loss, we can apply the mask
        y_true_filled = torch.where(mask, y_true, torch.zeros_like(y_true))
        y_pred_masked = torch.where(mask, y_pred, torch.zeros_like(y_pred))

        # This simple masking works for MSE/Poisson because (0-0)^2 = 0
        # and Poisson(0,0) = 0. For others like Cosine, it's more complex.
        loss = self.base_loss(y_pred_masked, y_true_filled)

        # Scale by fraction of valid elements to keep loss magnitude consistent
        scale = mask.numel() / (mask.float().sum() + EPSILON)
        return loss * scale


def get_loss(loss_name: str, center_crop_bins: Optional[int] = None, mask_nan: bool = False, **kwargs) -> nn.Module:
    """
    Return a loss module by name.
    
    Args:
        loss_name: 'poisson', 'poisson_multinomial', 'mse', 'mse_log', 
                   'pearson', 'crested', 'cosine_mse'
        center_crop_bins: If set, crop predictions/targets to this many bins
                          along the sequence dimension before computing loss.
        mask_nan: If True, ignore positions where target is NaN.
        **kwargs: Parameters for the loss
    
    Returns:
        nn.Module loss function
    """
    name = loss_name.lower()
    
    if name == 'poisson':
        base_loss = PoissonLoss(**kwargs)
    
    elif name == 'poisson_multinomial':
        base_loss = PoissonMultinomialLoss(**kwargs)
    
    elif name == 'mse':
        base_loss = nn.MSELoss()
    
    elif name == 'mse_log':
        base_loss = MSELogLoss(**kwargs)
    
    elif name == 'pearson':
        base_loss = PearsonLoss(**kwargs)
    
    elif name in ('crested', 'cosine_mse', 'cosine_mse_log'):
        base_loss = CosineMSELogLoss(**kwargs)
    
    elif name == 'cosine_mse_raw':
        base_loss = CosineMSERawLoss(**kwargs)
    
    else:
        raise ValueError(
            f"Unknown loss: {loss_name}. "
            "Options: 'poisson', 'poisson_multinomial', 'mse', 'mse_log', 'pearson', 'crested'"
        )
    
    if mask_nan:
        base_loss = MaskedLoss(base_loss)

    if center_crop_bins is not None:
        return CenterCropLoss(base_loss, crop_bins=center_crop_bins)

    return base_loss


def setup_loss(config, accelerator) -> nn.Module:
    """Setup loss function with optional center cropping and masking.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance for printing.

    Returns:
        Loss function module.
    """
    loss_name = config.training.loss
    accelerator.print(f"\nUsing loss: {loss_name}")

    mask_nan = config.get('training.mask_nan', False)
    if mask_nan:
        accelerator.print("  NaN masking enabled (ignoring missing targets in loss)")

    crop_bp = config.get('training.loss_center_crop_bp', None)
    crop_bins_cfg = config.get('training.loss_center_crop_bins', None)

    if crop_bp is not None and crop_bins_cfg is not None:
        expected_bins = crop_bp // config.model.bin_size
        if expected_bins != crop_bins_cfg:
            accelerator.print(f"Warning: loss_center_crop_bp and loss_center_crop_bins disagree")
        crop_bins = crop_bins_cfg
    elif crop_bins_cfg is not None:
        crop_bins = crop_bins_cfg
    elif crop_bp is not None:
        crop_bins = crop_bp // config.model.bin_size
    else:
        crop_bins = None

    if crop_bins is not None:
        accelerator.print(f"  Center-cropping loss to {crop_bins} bins ({crop_bins * config.model.bin_size} bp)")

    return get_loss(loss_name, center_crop_bins=crop_bins, mask_nan=mask_nan)
