"""
PyTorch metrics for Borzoi fine-tuning.
Includes Pearson and Spearman correlation for evaluating profile predictions.
"""

import numpy as np
import torch
import torch.nn as nn
from typing import Any, Dict, List, Optional, Set, Tuple
import warnings
from bisect import bisect_left
import re

# Numerical stability constant for division operations
EPSILON = 1e-8

# Maximum value for clipping coverage predictions/targets.
# Handles extreme genomic coverage values and prevents numerical overflow.
MAX_COVERAGE_VAL = 1e6

# Minimum variance threshold for R² and correlation calculations.
# Samples with variance below this are considered constant and excluded.
MIN_VARIANCE_THRESHOLD = 1e-6


class PearsonCorrelation(nn.Module):
    """
    Pearson correlation coefficient metric.
    
    Computes correlation between predicted and true profiles
    across bins for each track, then averages.
    """
    
    def __init__(self, epsilon: float = 1e-8):
        super().__init__()
        self.epsilon = epsilon
        self.reset()
        
    def reset(self):
        """Reset accumulated state."""
        self.total_correlation = 0.0
        self.count = 0
        
    def update(self, y_pred: torch.Tensor, y_true: torch.Tensor):
        """
        Update metric state with batch of predictions.

        Args:
            y_pred: Predictions of shape (batch, bins, tracks)
            y_true: Targets of shape (batch, bins, tracks)
        """
        with torch.no_grad():
            # Conditional type conversion to avoid unnecessary overhead
            if y_true.dtype != torch.float32:
                y_true = y_true.float()
            if y_pred.dtype != torch.float32:
                y_pred = y_pred.float()

            # Identify NaN/Inf BEFORE clamping — torch.clamp may convert NaN to
            # min on some CUDA devices, which would defeat the isfinite check.
            finite_mask = torch.isfinite(y_true) & torch.isfinite(y_pred)

            # Clip for numerical stability (genomic coverage can be large)
            y_true = torch.clamp(y_true, 0.0, MAX_COVERAGE_VAL)
            y_pred = torch.clamp(y_pred, 0.0, MAX_COVERAGE_VAL)

            batch_size = y_true.shape[0]
            num_tracks = y_true.shape[2]

            # Reshape to (batch * tracks, bins)
            y_true_flat = y_true.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
            y_pred_flat = y_pred.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)

            # Use pre-clamp mask (reshaped to match)
            valid_mask = finite_mask.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
            
            # Compute means only over valid positions
            # We use a safe denominator to avoid div by zero for fully masked rows
            valid_counts = valid_mask.float().sum(dim=-1, keepdim=True)
            sum_true = torch.where(valid_mask, y_true_flat, torch.zeros_like(y_true_flat)).sum(dim=-1, keepdim=True)
            sum_pred = torch.where(valid_mask, y_pred_flat, torch.zeros_like(y_pred_flat)).sum(dim=-1, keepdim=True)
            
            mean_true = sum_true / (valid_counts + self.epsilon)
            mean_pred = sum_pred / (valid_counts + self.epsilon)
            
            y_true_c = torch.where(valid_mask, y_true_flat - mean_true, torch.zeros_like(y_true_flat))
            y_pred_c = torch.where(valid_mask, y_pred_flat - mean_pred, torch.zeros_like(y_pred_flat))
            
            # Covariance / Std
            covariance = (y_true_c * y_pred_c).sum(dim=-1)
            true_std = torch.sqrt((y_true_c ** 2).sum(dim=-1) + self.epsilon)
            pred_std = torch.sqrt((y_pred_c ** 2).sum(dim=-1) + self.epsilon)
            
            corr = covariance / (true_std * pred_std + self.epsilon)
            
            # Mask rows that were fully invalid or had zero variance
            mask = (valid_counts.squeeze(-1) > 1) & torch.isfinite(corr)
            valid_corr = torch.where(mask, corr, torch.zeros_like(corr))
            num_valid = mask.float().sum()
            
            if num_valid > 0:
                avg_corr = valid_corr.sum() / num_valid
                self.total_correlation += avg_corr.item()
                self.count += 1
                
    def compute(self) -> float:
        """Compute final metric value."""
        if self.count == 0:
            return 0.0
        return self.total_correlation / self.count
    
    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        """Compute Pearson correlation for a single batch (differentiable)."""
        y_true = y_true.float()
        y_pred = y_pred.float()
        
        batch_size = y_true.shape[0]
        num_tracks = y_true.shape[2]
        
        y_true_flat = y_true.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
        y_pred_flat = y_pred.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
        
        mean_true = y_true_flat.mean(dim=-1, keepdim=True)
        mean_pred = y_pred_flat.mean(dim=-1, keepdim=True)
        
        y_true_c = y_true_flat - mean_true
        y_pred_c = y_pred_flat - mean_pred
        
        covariance = (y_true_c * y_pred_c).sum(dim=-1)
        true_std = torch.sqrt((y_true_c ** 2).sum(dim=-1) + self.epsilon)
        pred_std = torch.sqrt((y_pred_c ** 2).sum(dim=-1) + self.epsilon)
        
        corr = covariance / (true_std * pred_std + self.epsilon)
        corr = torch.clamp(corr, -1.0, 1.0)
        
        mask = torch.isfinite(corr)
        valid_corr = torch.where(mask, corr, torch.zeros_like(corr))
        num_valid = mask.float().sum()
        
        return valid_corr.sum() / torch.clamp(num_valid, min=1.0)


class SpearmanCorrelation(nn.Module):
    """
    Spearman rank correlation coefficient metric.

    Uses Pearson on ranks. Note: This is an approximation since
    proper ranking with tie handling is complex in PyTorch.
    """

    def __init__(self, epsilon: float = 1e-8):
        super().__init__()
        self.epsilon = epsilon
        self.reset()

    def reset(self):
        """Reset accumulated state."""
        self.total_correlation = 0.0
        self.count = 0

    def _masked_rank_transform(self, x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        """Rank transformation that ignores masked values.
        
        Masked values are assigned a dummy rank that won't affect correlation
        as they will be masked again during the Pearson computation.
        """
        ranks = torch.zeros_like(x)
        for i in range(x.shape[0]):
            m = mask[i]
            if m.any():
                valid_x = x[i][m]
                # Double argsort for ranks
                r = torch.argsort(torch.argsort(valid_x)).float()
                ranks[i][m] = r
        return ranks

    def update(self, y_pred: torch.Tensor, y_true: torch.Tensor):
        """Update metric state with batch of predictions."""
        with torch.no_grad():
            # Conditional type conversion to avoid unnecessary overhead
            if y_true.dtype != torch.float32:
                y_true = y_true.float()
            if y_pred.dtype != torch.float32:
                y_pred = y_pred.float()
            
            batch_size = y_true.shape[0]
            num_tracks = y_true.shape[2]
            
            y_true_flat = y_true.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
            y_pred_flat = y_pred.permute(0, 2, 1).reshape(batch_size * num_tracks, -1)
            
            # Mask NaNs in targets
            valid_mask = torch.isfinite(y_true_flat)
            
            # Convert to ranks only for valid positions
            y_true_ranks = self._masked_rank_transform(y_true_flat, valid_mask)
            y_pred_ranks = self._masked_rank_transform(y_pred_flat, valid_mask)
            
            # Compute Pearson on ranks only over valid positions
            valid_counts = valid_mask.float().sum(dim=-1, keepdim=True)
            sum_true = torch.where(valid_mask, y_true_ranks, torch.zeros_like(y_true_ranks)).sum(dim=-1, keepdim=True)
            sum_pred = torch.where(valid_mask, y_pred_ranks, torch.zeros_like(y_pred_ranks)).sum(dim=-1, keepdim=True)
            
            mean_true = sum_true / (valid_counts + self.epsilon)
            mean_pred = sum_pred / (valid_counts + self.epsilon)
            
            y_true_c = torch.where(valid_mask, y_true_ranks - mean_true, torch.zeros_like(y_true_ranks))
            y_pred_c = torch.where(valid_mask, y_pred_ranks - mean_pred, torch.zeros_like(y_pred_ranks))
            
            covariance = (y_true_c * y_pred_c).sum(dim=-1)
            true_std = torch.sqrt((y_true_c ** 2).sum(dim=-1) + self.epsilon)
            pred_std = torch.sqrt((y_pred_c ** 2).sum(dim=-1) + self.epsilon)
            
            corr = covariance / (true_std * pred_std + self.epsilon)
            
            mask = (valid_counts.squeeze(-1) > 1) & torch.isfinite(corr)
            valid_corr = torch.where(mask, corr, torch.zeros_like(corr))
            num_valid = mask.float().sum()
            
            if num_valid > 0:
                avg_corr = valid_corr.sum() / num_valid
                self.total_correlation += avg_corr.item()
                self.count += 1
                
    def compute(self) -> float:
        """Compute final metric value."""
        if self.count == 0:
            return 0.0
        return self.total_correlation / self.count
    
    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        """Compute Spearman correlation for a single batch."""
        batch_size = y_true.shape[0]
        num_tracks = y_true.shape[2]
        
        y_true_flat = y_true.permute(0, 2, 1).reshape(batch_size * num_tracks, -1).float()
        y_pred_flat = y_pred.permute(0, 2, 1).reshape(batch_size * num_tracks, -1).float()
        
        y_true_ranks = self._rank_transform(y_true_flat)
        y_pred_ranks = self._rank_transform(y_pred_flat)
        
        mean_true = y_true_ranks.mean(dim=-1, keepdim=True)
        mean_pred = y_pred_ranks.mean(dim=-1, keepdim=True)
        
        y_true_c = y_true_ranks - mean_true
        y_pred_c = y_pred_ranks - mean_pred
        
        covariance = (y_true_c * y_pred_c).sum(dim=-1)
        true_std = torch.sqrt((y_true_c ** 2).sum(dim=-1) + self.epsilon)
        pred_std = torch.sqrt((y_pred_c ** 2).sum(dim=-1) + self.epsilon)
        
        corr = covariance / (true_std * pred_std + self.epsilon)
        
        mask = torch.isfinite(corr)
        valid_corr = torch.where(mask, corr, torch.zeros_like(corr))
        num_valid = mask.float().sum()
        
        return valid_corr.sum() / torch.clamp(num_valid, min=1.0)


class R2Score(nn.Module):
    """
    R-squared (coefficient of determination) metric.
    """
    
    def __init__(self, epsilon: float = 1e-8):
        super().__init__()
        self.epsilon = epsilon
        self.reset()
        
    def reset(self):
        """Reset accumulated state."""
        self.total_r2 = 0.0
        self.count = 0
        
    def update(self, y_pred: torch.Tensor, y_true: torch.Tensor):
        """Update metric state with batch of predictions."""
        with torch.no_grad():
            # Conditional type conversion to avoid unnecessary overhead
            if y_true.dtype != torch.float32:
                y_true = y_true.float()
            if y_pred.dtype != torch.float32:
                y_pred = y_pred.float()

            # Identify NaN/Inf BEFORE clamping — torch.clamp may convert NaN to
            # min on some CUDA devices, which would defeat the isfinite check.
            finite_mask = torch.isfinite(y_true) & torch.isfinite(y_pred)

            # Clip for numerical stability (genomic coverage can be large)
            y_true = torch.clamp(y_true, 0.0, MAX_COVERAGE_VAL)
            y_pred = torch.clamp(y_pred, 0.0, MAX_COVERAGE_VAL)

            # Flatten spatial dims: (batch, bins, tracks) -> (batch, bins * tracks)
            y_true_flat = y_true.reshape(y_true.shape[0], -1)
            y_pred_flat = y_pred.reshape(y_pred.shape[0], -1)

            # Use pre-clamp mask (reshaped to match)
            valid_mask = finite_mask.reshape(y_true.shape[0], -1)
            
            # Sum of Squares Residuals (only over valid positions)
            ss_res = torch.where(valid_mask, (y_true_flat - y_pred_flat) ** 2, torch.zeros_like(y_true_flat)).sum(dim=-1)
            
            # Total Sum of Squares (only over valid positions)
            valid_counts = valid_mask.float().sum(dim=-1, keepdim=True)
            sum_true = torch.where(valid_mask, y_true_flat, torch.zeros_like(y_true_flat)).sum(dim=-1, keepdim=True)
            y_mean = sum_true / (valid_counts + self.epsilon)
            
            ss_tot = torch.where(valid_mask, (y_true_flat - y_mean) ** 2, torch.zeros_like(y_true_flat)).sum(dim=-1)
            
            # Filter low variance
            variance_mask = (ss_tot > MIN_VARIANCE_THRESHOLD) & (valid_counts.squeeze(-1) > 1)
            
            r2_values = 1.0 - (ss_res / (ss_tot + self.epsilon))
            
            valid_r2 = torch.where(variance_mask, r2_values, torch.zeros_like(r2_values))
            num_valid = variance_mask.float().sum()
            
            if num_valid > 0:
                batch_r2 = valid_r2.sum() / num_valid
                self.total_r2 += batch_r2.item()
                self.count += 1
                
    def compute(self) -> float:
        """Compute final metric value."""
        if self.count == 0:
            return 0.0
        return self.total_r2 / self.count
    
    def forward(self, y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        """Compute R2 for a single batch."""
        y_true = y_true.float()
        y_pred = y_pred.float()
        
        y_true_flat = y_true.reshape(y_true.shape[0], -1)
        y_pred_flat = y_pred.reshape(y_pred.shape[0], -1)
        
        ss_res = ((y_true_flat - y_pred_flat) ** 2).sum(dim=-1)
        y_mean = y_true_flat.mean(dim=-1, keepdim=True)
        ss_tot = ((y_true_flat - y_mean) ** 2).sum(dim=-1)
        
        variance_mask = ss_tot > MIN_VARIANCE_THRESHOLD
        r2_values = 1.0 - (ss_res / (ss_tot + self.epsilon))
        
        valid_r2 = torch.where(variance_mask, r2_values, torch.zeros_like(r2_values))
        num_valid = variance_mask.float().sum()
        
        return valid_r2.sum() / torch.clamp(num_valid, min=1.0)


class MSE(nn.Module):
    """Mean Squared Error metric."""

    def __init__(self):
        super().__init__()
        self.reset()

    def reset(self):
        self.total_mse = 0.0
        self.count = 0

    def update(self, y_pred: torch.Tensor, y_true: torch.Tensor):
        with torch.no_grad():
            mask = torch.isfinite(y_true)
            if not mask.any():
                return
                
            squared_diff = (y_true - y_pred) ** 2
            # Mean only over valid positions
            mse = squared_diff[mask].mean().item()
            self.total_mse += mse
            self.count += 1

    def compute(self) -> float:
        if self.count == 0:
            return 0.0
        return self.total_mse / self.count


class GenomeWidePearson(nn.Module):
    """
    Genome-wide Pearson correlation metric.

    Unlike PearsonCorrelation which computes correlation within each window
    then averages, this metric accumulates all bins across the entire validation
    set and computes ONE correlation per track at the end.

    This makes the metric comparable across different resolutions (1k, 2k, 4k, etc.)
    since it measures correlation at the same genomic scale regardless of window size.

    Uses online/streaming formula to avoid storing all predictions:
    r = (n*sum(xy) - sum(x)*sum(y)) / sqrt((n*sum(x^2) - sum(x)^2) * (n*sum(y^2) - sum(y)^2))
    """

    def __init__(
        self,
        epsilon: float = 1e-8,
        accelerator=None,
        track_indices: Optional[List[int]] = None,
    ):
        super().__init__()
        self.epsilon = epsilon
        self.accelerator = accelerator
        self.track_indices = list(track_indices) if track_indices is not None else None
        self.reset()

    def reset(self):
        """Reset accumulated state."""
        # Running sums per track for online Pearson computation
        # Will be initialized on first update when we know num_tracks
        self.sum_x = None      # sum of predictions
        self.sum_y = None      # sum of targets
        self.sum_xy = None     # sum of pred * target
        self.sum_x2 = None     # sum of pred^2
        self.sum_y2 = None     # sum of target^2
        self.count = None      # count of bins per track
        self.num_tracks = None
        self._resolved_track_indices = None
        self._track_selector_initialized = False

    def _resolve_track_subset(self, total_tracks: int) -> Optional[List[int]]:
        """Resolve and validate selected track indices for this run."""
        if self._track_selector_initialized:
            return self._resolved_track_indices

        self._track_selector_initialized = True
        if self.track_indices is None:
            self._resolved_track_indices = list(range(total_tracks))
            return self._resolved_track_indices

        valid_indices = sorted({i for i in self.track_indices if 0 <= i < total_tracks})
        dropped = sorted({i for i in self.track_indices if i < 0 or i >= total_tracks})
        if dropped:
            warnings.warn(
                f"Dropping out-of-range track indices for {self.__class__.__name__}: {dropped}",
                stacklevel=2,
            )

        if not valid_indices:
            warnings.warn(
                f"{self.__class__.__name__}: no valid track indices selected; metric will stay at 0.0",
                stacklevel=2,
            )
            self._resolved_track_indices = None
            return None

        self._resolved_track_indices = valid_indices
        return self._resolved_track_indices

    def update(self, y_pred: torch.Tensor, y_true: torch.Tensor):
        """
        Update metric state with batch of predictions.

        Args:
            y_pred: Predictions of shape (batch, bins, tracks)
            y_true: Targets of shape (batch, bins, tracks)
        """
        with torch.no_grad():
            # Conditional type conversion
            if y_true.dtype != torch.float32:
                y_true = y_true.float()
            if y_pred.dtype != torch.float32:
                y_pred = y_pred.float()

            # Identify NaN/Inf BEFORE clamping — torch.clamp may convert NaN to
            # min on some CUDA devices, which would defeat the isfinite check.
            finite_mask = torch.isfinite(y_true) & torch.isfinite(y_pred)

            # Clip for numerical stability
            y_true = torch.clamp(y_true, 0.0, MAX_COVERAGE_VAL)
            y_pred = torch.clamp(y_pred, 0.0, MAX_COVERAGE_VAL)

            _, _, num_tracks = y_true.shape
            selected_track_indices = self._resolve_track_subset(num_tracks)
            if not selected_track_indices:
                return
            selected_tracks = len(selected_track_indices)

            # Initialize accumulators on first call
            if self.sum_x is None:
                self.num_tracks = selected_tracks
                device = y_pred.device
                self.sum_x = torch.zeros(selected_tracks, device=device, dtype=torch.float64)
                self.sum_y = torch.zeros(selected_tracks, device=device, dtype=torch.float64)
                self.sum_xy = torch.zeros(selected_tracks, device=device, dtype=torch.float64)
                self.sum_x2 = torch.zeros(selected_tracks, device=device, dtype=torch.float64)
                self.sum_y2 = torch.zeros(selected_tracks, device=device, dtype=torch.float64)
                self.count = torch.zeros(selected_tracks, device=device, dtype=torch.float64)

            # Reshape to (batch * bins, tracks) for accumulation
            y_pred_flat = y_pred.reshape(-1, num_tracks).double()
            y_true_flat = y_true.reshape(-1, num_tracks).double()
            valid_mask = finite_mask.reshape(-1, num_tracks)
            if selected_tracks != num_tracks:
                track_index_tensor = torch.tensor(selected_track_indices, device=y_pred.device, dtype=torch.long)
                y_pred_flat = y_pred_flat.index_select(dim=1, index=track_index_tensor)
                y_true_flat = y_true_flat.index_select(dim=1, index=track_index_tensor)
                valid_mask = valid_mask.index_select(dim=1, index=track_index_tensor)

            # Zero out invalid values for safe accumulation
            y_pred_masked = torch.where(valid_mask, y_pred_flat, torch.zeros_like(y_pred_flat))
            y_true_masked = torch.where(valid_mask, y_true_flat, torch.zeros_like(y_true_flat))

            # Accumulate sums per track
            self.sum_x += y_pred_masked.sum(dim=0)
            self.sum_y += y_true_masked.sum(dim=0)
            self.sum_xy += (y_pred_masked * y_true_masked).sum(dim=0)
            self.sum_x2 += (y_pred_masked ** 2).sum(dim=0)
            self.sum_y2 += (y_true_masked ** 2).sum(dim=0)
            self.count += valid_mask.float().sum(dim=0).double()

    def compute(self) -> float:
        """Compute final genome-wide Pearson correlation (averaged across tracks)."""
        if self.sum_x is None or self.count is None:
            return 0.0

        # Handle distributed training: reduce sums across all processes
        if self.accelerator is not None and self.accelerator.num_processes > 1:
            self.sum_x = self.accelerator.reduce(self.sum_x, reduction="sum")
            self.sum_y = self.accelerator.reduce(self.sum_y, reduction="sum")
            self.sum_xy = self.accelerator.reduce(self.sum_xy, reduction="sum")
            self.sum_x2 = self.accelerator.reduce(self.sum_x2, reduction="sum")
            self.sum_y2 = self.accelerator.reduce(self.sum_y2, reduction="sum")
            self.count = self.accelerator.reduce(self.count, reduction="sum")

        # Compute Pearson correlation per track using the streaming formula
        # r = (n*sum(xy) - sum(x)*sum(y)) / sqrt((n*sum(x^2) - sum(x)^2) * (n*sum(y^2) - sum(y)^2))
        n = self.count

        numerator = n * self.sum_xy - self.sum_x * self.sum_y
        denominator_x = n * self.sum_x2 - self.sum_x ** 2
        denominator_y = n * self.sum_y2 - self.sum_y ** 2

        # Handle edge cases
        denominator = torch.sqrt(torch.clamp(denominator_x, min=0) * torch.clamp(denominator_y, min=0))

        # Compute correlation, avoiding division by zero
        valid_mask = (denominator > self.epsilon) & (n > 1)
        corr = torch.zeros_like(numerator)
        corr[valid_mask] = numerator[valid_mask] / (denominator[valid_mask] + self.epsilon)

        # Clamp to valid range
        corr = torch.clamp(corr, -1.0, 1.0)

        # Average across tracks (only valid ones)
        if valid_mask.any():
            avg_corr = corr[valid_mask].mean().item()
        else:
            avg_corr = 0.0

        return avg_corr


class GenomeWideHistogramMetric(nn.Module):
    """
    Base class for genome-wide metrics using histogram-based accumulation.

    Accumulates 1D histograms for predictions and targets, plus a 2D joint histogram,
    across the entire validation set. This enables memory-efficient computation of
    quantile-based metrics genome-wide.

    Subclasses implement compute() to derive specific metrics from the histograms.
    """

    def __init__(
        self,
        quantile: float = 0.98,
        smooth_bins: int = 1,
        num_hist_bins: int = 1000,
        value_max: float = 1000.0,
        use_log_scale: bool = False,
        epsilon: float = 1e-8,
        accelerator=None,
        track_indices: Optional[List[int]] = None,
        metric_suffix: Optional[str] = None,
    ):
        """
        Args:
            quantile: Quantile threshold for defining "high" regions (e.g., 0.98 = top 2%)
            smooth_bins: Number of bins to smooth over (moving average window size).
                        Set to 1 for no smoothing (default).
            num_hist_bins: Number of histogram bins for quantile estimation
            value_max: Maximum value for histogram (values are clipped)
            use_log_scale: Use log-scale bins. Set to False (default) if signal is
                          already log-transformed. Set to True for raw signal values.
            epsilon: Small constant for numerical stability
            accelerator: Accelerator instance for distributed training
            track_indices: Optional subset of track indices to evaluate.
            metric_suffix: Optional suffix appended to output keys (e.g., 'atac').
        """
        super().__init__()
        self.quantile = quantile
        self.smooth_bins = smooth_bins
        self.num_hist_bins = num_hist_bins
        self.value_max = value_max
        self.use_log_scale = use_log_scale
        self.epsilon = epsilon
        self.accelerator = accelerator
        self.track_indices = list(track_indices) if track_indices is not None else None
        self.metric_suffix = metric_suffix.strip() if metric_suffix else None

        # Pre-compute bin edges for histogram
        if use_log_scale:
            # Log-scale bins with a small offset to handle zeros
            log_min = np.log1p(0)  # = 0
            log_max = np.log1p(value_max)
            self.bin_edges = torch.from_numpy(
                np.expm1(np.linspace(log_min, log_max, num_hist_bins + 1))
            ).float()
        else:
            self.bin_edges = torch.linspace(0, value_max, num_hist_bins + 1)

        self.reset()

    def reset(self):
        """Reset accumulated state."""
        self.pred_hist = None       # (num_tracks, num_hist_bins)
        self.target_hist = None     # (num_tracks, num_hist_bins)
        self.joint_hist = None      # (num_tracks, num_hist_bins, num_hist_bins)
        self.num_tracks = None
        self.total_samples = 0
        self._resolved_track_indices = None
        self._track_selector_initialized = False

    def _smooth(self, x: torch.Tensor) -> torch.Tensor:
        """
        Apply moving average smoothing along the bins dimension.

        Args:
            x: Tensor of shape (batch, bins, tracks)

        Returns:
            Smoothed tensor of shape (batch, bins - smooth_bins + 1, tracks)
        """
        if self.smooth_bins <= 1:
            return x

        batch_size, num_bins, num_tracks = x.shape

        if num_bins < self.smooth_bins:
            # If fewer bins than smooth window, just average all
            return x.mean(dim=1, keepdim=True)

        # Transpose to (batch, tracks, bins) for conv1d-style operation
        x_t = x.permute(0, 2, 1)  # (batch, tracks, bins)

        # Use avg_pool1d for efficient smoothing
        x_flat = x_t.reshape(-1, 1, num_bins)
        smoothed = torch.nn.functional.avg_pool1d(
            x_flat,
            kernel_size=self.smooth_bins,
            stride=1,
            padding=0
        )

        # Reshape back to (batch, tracks, new_bins) then to (batch, new_bins, tracks)
        new_bins = smoothed.shape[-1]
        smoothed = smoothed.reshape(batch_size, num_tracks, new_bins)
        smoothed = smoothed.permute(0, 2, 1)  # (batch, new_bins, tracks)

        return smoothed

    def _digitize(self, values: torch.Tensor) -> torch.Tensor:
        """Assign values to histogram bins."""
        # Clip values to valid range
        values = torch.clamp(values, 0, self.value_max - self.epsilon)

        # Move bin_edges to same device
        bin_edges = self.bin_edges.to(values.device)

        # Use bucketize (PyTorch's digitize equivalent)
        bins = torch.bucketize(values, bin_edges[1:-1])
        return bins

    def _format_output_key(self, key: str) -> str:
        """Format output key with optional metric suffix."""
        if self.metric_suffix:
            return f"{key}_{self.metric_suffix}"
        return key

    def _resolve_track_subset(self, total_tracks: int) -> Optional[List[int]]:
        """Resolve and validate selected track indices for this run."""
        if self._track_selector_initialized:
            return self._resolved_track_indices

        self._track_selector_initialized = True
        if self.track_indices is None:
            self._resolved_track_indices = list(range(total_tracks))
            return self._resolved_track_indices

        valid_indices = sorted({i for i in self.track_indices if 0 <= i < total_tracks})
        dropped = sorted({i for i in self.track_indices if i < 0 or i >= total_tracks})
        if dropped:
            warnings.warn(
                f"Dropping out-of-range track indices for {self.__class__.__name__}: {dropped}",
                stacklevel=2,
            )

        if not valid_indices:
            suffix = f" ({self.metric_suffix})" if self.metric_suffix else ""
            warnings.warn(
                f"{self.__class__.__name__}{suffix}: no valid track indices selected; metric will stay at 0.0",
                stacklevel=2,
            )
            self._resolved_track_indices = None
            return None

        self._resolved_track_indices = valid_indices
        return self._resolved_track_indices

    def update(self, y_pred: torch.Tensor, y_true: torch.Tensor):
        """
        Update metric state with batch of predictions.

        Args:
            y_pred: Predictions of shape (batch, bins, tracks)
            y_true: Targets of shape (batch, bins, tracks)
        """
        with torch.no_grad():
            if y_true.dtype != torch.float32:
                y_true = y_true.float()
            if y_pred.dtype != torch.float32:
                y_pred = y_pred.float()

            # Apply smoothing (no-op if smooth_bins <= 1)
            y_pred_smooth = self._smooth(y_pred)
            y_true_smooth = self._smooth(y_true)

            _, _, num_tracks = y_pred_smooth.shape
            device = y_pred.device
            selected_track_indices = self._resolve_track_subset(num_tracks)
            if not selected_track_indices:
                return

            selected_tracks = len(selected_track_indices)

            # Initialize histograms on first call
            if self.pred_hist is None:
                self.num_tracks = selected_tracks
                self.pred_hist = torch.zeros(selected_tracks, self.num_hist_bins, device=device, dtype=torch.long)
                self.target_hist = torch.zeros(selected_tracks, self.num_hist_bins, device=device, dtype=torch.long)
                self.joint_hist = torch.zeros(selected_tracks, self.num_hist_bins, self.num_hist_bins, device=device, dtype=torch.long)

            # Reshape to (batch * bins, tracks)
            y_pred_flat = y_pred_smooth.reshape(-1, num_tracks)
            y_true_flat = y_true_smooth.reshape(-1, num_tracks)
            if selected_tracks != num_tracks:
                track_index_tensor = torch.tensor(selected_track_indices, device=device, dtype=torch.long)
                y_pred_flat = y_pred_flat.index_select(dim=1, index=track_index_tensor)
                y_true_flat = y_true_flat.index_select(dim=1, index=track_index_tensor)

            # Mask out NaN/Inf values to avoid corrupting histograms
            valid_mask = torch.isfinite(y_pred_flat) & torch.isfinite(y_true_flat)  # (N, tracks)

            # Replace non-finite values with 0 so _digitize doesn't produce undefined bins
            y_pred_safe = torch.where(valid_mask, y_pred_flat, torch.zeros_like(y_pred_flat))
            y_true_safe = torch.where(valid_mask, y_true_flat, torch.zeros_like(y_true_flat))

            # Digitize values
            pred_bins = self._digitize(y_pred_safe)  # (N, tracks)
            target_bins = self._digitize(y_true_safe)  # (N, tracks)

            # Update histograms for each track (only valid positions)
            for t in range(selected_tracks):
                track_valid = valid_mask[:, t]
                if not track_valid.any():
                    continue

                valid_pred = pred_bins[:, t][track_valid]
                valid_target = target_bins[:, t][track_valid]

                # 1D histograms
                pred_counts = torch.bincount(valid_pred, minlength=self.num_hist_bins)
                target_counts = torch.bincount(valid_target, minlength=self.num_hist_bins)
                self.pred_hist[t] += pred_counts[:self.num_hist_bins]
                self.target_hist[t] += target_counts[:self.num_hist_bins]

                # 2D joint histogram: joint_hist[pred_bin, target_bin]
                joint_indices = valid_pred * self.num_hist_bins + valid_target
                joint_counts = torch.bincount(joint_indices, minlength=self.num_hist_bins * self.num_hist_bins)
                self.joint_hist[t] += joint_counts[:self.num_hist_bins * self.num_hist_bins].reshape(
                    self.num_hist_bins, self.num_hist_bins
                )

            self.total_samples += y_pred_flat.shape[0]

    def _reduce_histograms(self):
        """Reduce histograms across processes for distributed training."""
        if self.accelerator is not None and self.accelerator.num_processes > 1:
            self.pred_hist = self.accelerator.reduce(self.pred_hist, reduction="sum")
            self.target_hist = self.accelerator.reduce(self.target_hist, reduction="sum")
            self.joint_hist = self.accelerator.reduce(self.joint_hist, reduction="sum")

    def _find_quantile_bin(self, hist: torch.Tensor, quantile: float) -> int:
        """Find the bin index corresponding to a quantile threshold."""
        cumsum = torch.cumsum(hist, dim=0).float()
        total = cumsum[-1]
        if total == 0:
            return self.num_hist_bins - 1
        threshold_count = quantile * total
        bin_idx = torch.searchsorted(cumsum, threshold_count).item()
        return min(bin_idx, self.num_hist_bins - 1)

    def compute(self) -> Dict[str, float]:
        """Compute metrics from accumulated histograms. Override in subclasses."""
        raise NotImplementedError("Subclasses must implement compute()")


class GenomeWidePPV(GenomeWideHistogramMetric):
    """
    Genome-wide Positive Predictive Value (PPV) metric.

    Computes PPV (precision), sensitivity (recall), and F1 based on quantile thresholds:
    - Target positive: value >= quantile threshold of all targets (per-track)
    - Prediction positive: value >= quantile threshold of all predictions (per-track)

    PPV = TP / (TP + FP)
    Sensitivity = TP / (TP + FN)
    F1 = 2 × PPV × Sensitivity / (PPV + Sensitivity)

    Optionally applies smoothing before threshold computation.
    """

    def __init__(
        self,
        quantile: float = 0.98,
        smooth_bins: int = 10,
        num_hist_bins: int = 1000,
        value_max: float = 1000.0,
        use_log_scale: bool = False,
        epsilon: float = 1e-8,
        accelerator=None,
        track_indices: Optional[List[int]] = None,
        metric_suffix: Optional[str] = None,
    ):
        """
        Args:
            quantile: Quantile threshold for defining positives (e.g., 0.98 = top 2%)
            smooth_bins: Number of bins to smooth over (moving average). Default 10.
            num_hist_bins: Number of histogram bins for quantile estimation
            value_max: Maximum value for histogram (values are clipped)
            use_log_scale: Use log-scale bins. Default False (signal often already log-scaled).
            epsilon: Small constant for numerical stability
            accelerator: Accelerator instance for distributed training
            track_indices: Optional subset of track indices to evaluate.
            metric_suffix: Optional suffix appended to output keys.
        """
        super().__init__(
            quantile=quantile,
            smooth_bins=smooth_bins,
            num_hist_bins=num_hist_bins,
            value_max=value_max,
            use_log_scale=use_log_scale,
            epsilon=epsilon,
            accelerator=accelerator,
            track_indices=track_indices,
            metric_suffix=metric_suffix,
        )

    def compute(self) -> Dict[str, float]:
        """
        Compute PPV, sensitivity, and F1 metrics.

        Returns:
            Dictionary with 'ppv', 'sensitivity', 'f1'.
        """
        if self.pred_hist is None or self.total_samples == 0:
            return {
                self._format_output_key('ppv'): 0.0,
                self._format_output_key('sensitivity'): 0.0,
                self._format_output_key('f1'): 0.0,
            }

        self._reduce_histograms()

        ppv_per_track = []
        sensitivity_per_track = []
        f1_per_track = []

        for t in range(self.num_tracks):
            pred_threshold_bin = self._find_quantile_bin(self.pred_hist[t], self.quantile)
            target_threshold_bin = self._find_quantile_bin(self.target_hist[t], self.quantile)

            joint = self.joint_hist[t]

            # TP: pred >= threshold AND target >= threshold
            tp = joint[pred_threshold_bin:, target_threshold_bin:].sum().float()
            # FP: pred >= threshold AND target < threshold
            fp = joint[pred_threshold_bin:, :target_threshold_bin].sum().float()
            # FN: pred < threshold AND target >= threshold
            fn = joint[:pred_threshold_bin, target_threshold_bin:].sum().float()

            ppv = tp / (tp + fp + self.epsilon)
            sensitivity = tp / (tp + fn + self.epsilon)
            f1 = 2 * ppv * sensitivity / (ppv + sensitivity + self.epsilon)

            ppv_per_track.append(ppv.item())
            sensitivity_per_track.append(sensitivity.item())
            f1_per_track.append(f1.item())

        return {
            self._format_output_key('ppv'): np.mean(ppv_per_track),
            self._format_output_key('sensitivity'): np.mean(sensitivity_per_track),
            self._format_output_key('f1'): np.mean(f1_per_track),
        }


class GenomeWideQuantileOverlap(GenomeWideHistogramMetric):
    """
    Genome-wide quantile overlap metrics (IoU, Precision, Recall, F1).

    Computes overlap between high-signal regions in predictions and targets:
    - "High" is defined as bins above the quantile threshold (default: 98th percentile)

    Two sets of metrics are computed:

    1. IoU (Jaccard) - using SEPARATE thresholds for pred and target:
       - pred_high = pred >= pred's quantile threshold
       - target_high = target >= target's quantile threshold
       - IoU = |intersection| / |union|

    2. Precision/Recall/F1 - using TARGET's threshold for both:
       - threshold = target's quantile threshold (applied to both pred and target)
       - Precision = |intersection| / |pred_high| (of predicted peaks, what fraction are true?)
       - Recall = |intersection| / |target_high| (of true peaks, what fraction did we find?)
       - F1 = harmonic mean of precision and recall

    Optionally applies smoothing before threshold computation.
    """

    def __init__(
        self,
        quantile: float = 0.98,
        smooth_bins: int = 1,
        num_hist_bins: int = 1000,
        value_max: float = 1000.0,
        use_log_scale: bool = False,
        epsilon: float = 1e-8,
        accelerator=None,
        track_indices: Optional[List[int]] = None,
        metric_suffix: Optional[str] = None,
    ):
        """
        Args:
            quantile: Quantile threshold for defining "high" regions (e.g., 0.98 = top 2%)
            smooth_bins: Number of bins to smooth over. Default 1 (no smoothing).
            num_hist_bins: Number of histogram bins for quantile estimation
            value_max: Maximum value for histogram (values are clipped)
            use_log_scale: Use log-scale bins. Default False (signal often already log-scaled).
            epsilon: Small constant for numerical stability
            accelerator: Accelerator instance for distributed training
            track_indices: Optional subset of track indices to evaluate.
            metric_suffix: Optional suffix appended to output keys.
        """
        super().__init__(
            quantile=quantile,
            smooth_bins=smooth_bins,
            num_hist_bins=num_hist_bins,
            value_max=value_max,
            use_log_scale=use_log_scale,
            epsilon=epsilon,
            accelerator=accelerator,
            track_indices=track_indices,
            metric_suffix=metric_suffix,
        )

    def compute(self) -> Dict[str, float]:
        """
        Compute IoU, Precision, Recall, and F1 metrics.

        Returns:
            Dictionary with:
            - 'iou': Jaccard index using separate thresholds for pred/target
            - 'precision': using target's threshold for both
            - 'recall': using target's threshold for both
            - 'f1': harmonic mean of precision and recall
        """
        if self.pred_hist is None or self.total_samples == 0:
            return {
                self._format_output_key('iou'): 0.0,
                self._format_output_key('precision'): 0.0,
                self._format_output_key('recall'): 0.0,
                self._format_output_key('f1'): 0.0,
            }

        self._reduce_histograms()

        iou_per_track = []
        precision_per_track = []
        recall_per_track = []
        f1_per_track = []

        for t in range(self.num_tracks):
            # Get thresholds
            pred_threshold_bin = self._find_quantile_bin(self.pred_hist[t], self.quantile)
            target_threshold_bin = self._find_quantile_bin(self.target_hist[t], self.quantile)

            joint = self.joint_hist[t]

            # =================================================================
            # IoU: using SEPARATE thresholds for pred and target
            # =================================================================
            # pred_high = pred >= pred_threshold, target_high = target >= target_threshold
            intersection_iou = joint[pred_threshold_bin:, target_threshold_bin:].sum().float()
            pred_high_count_iou = joint[pred_threshold_bin:, :].sum().float()
            target_high_count_iou = joint[:, target_threshold_bin:].sum().float()
            union_iou = pred_high_count_iou + target_high_count_iou - intersection_iou

            iou = intersection_iou / (union_iou + self.epsilon)

            # =================================================================
            # Precision/Recall/F1: using TARGET's threshold for BOTH
            # =================================================================
            # Both pred_high and target_high use target_threshold_bin
            # intersection: pred >= target_threshold AND target >= target_threshold
            intersection_prf = joint[target_threshold_bin:, target_threshold_bin:].sum().float()

            # pred_high: pred >= target_threshold (count across all target bins)
            pred_high_count_prf = joint[target_threshold_bin:, :].sum().float()

            # target_high: target >= target_threshold (count across all pred bins)
            target_high_count_prf = joint[:, target_threshold_bin:].sum().float()

            precision = intersection_prf / (pred_high_count_prf + self.epsilon)
            recall = intersection_prf / (target_high_count_prf + self.epsilon)
            f1 = 2 * precision * recall / (precision + recall + self.epsilon)

            iou_per_track.append(iou.item())
            precision_per_track.append(precision.item())
            recall_per_track.append(recall.item())
            f1_per_track.append(f1.item())

        return {
            self._format_output_key('iou'): np.mean(iou_per_track),
            self._format_output_key('precision'): np.mean(precision_per_track),
            self._format_output_key('recall'): np.mean(recall_per_track),
            self._format_output_key('f1'): np.mean(f1_per_track),
        }


class GenomeWideAUPRC(GenomeWideHistogramMetric):
    """
    Genome-wide area under precision-recall curve (AUPRC).

    Defines positives as target bins above a per-track quantile threshold, then
    sweeps prediction thresholds across histogram bins to compute a PR curve.
    AUPRC is computed per track and averaged across tracks.
    """

    def __init__(
        self,
        quantile: float = 0.98,
        smooth_bins: int = 1,
        num_hist_bins: int = 1000,
        value_max: float = 1000.0,
        use_log_scale: bool = False,
        epsilon: float = 1e-8,
        accelerator=None,
        track_indices: Optional[List[int]] = None,
        metric_suffix: Optional[str] = None,
    ):
        super().__init__(
            quantile=quantile,
            smooth_bins=smooth_bins,
            num_hist_bins=num_hist_bins,
            value_max=value_max,
            use_log_scale=use_log_scale,
            epsilon=epsilon,
            accelerator=accelerator,
            track_indices=track_indices,
            metric_suffix=metric_suffix,
        )

    def compute(self) -> Dict[str, float]:
        """
        Compute AUPRC using target-quantile positives and a sweep over pred bins.

        Returns:
            Dictionary with key 'auprc' (or suffixed key for subset metrics).
        """
        if self.pred_hist is None or self.total_samples == 0:
            return {self._format_output_key('auprc'): 0.0}

        self._reduce_histograms()

        auprc_per_track = []
        for t in range(self.num_tracks):
            target_threshold_bin = self._find_quantile_bin(self.target_hist[t], self.quantile)
            joint = self.joint_hist[t].float()

            # Positives are bins where target >= target quantile threshold.
            positive_total = joint[:, target_threshold_bin:].sum()
            if positive_total <= 0:
                auprc_per_track.append(0.0)
                continue

            # For each prediction threshold bin i:
            # TP(i): pred >= i and target >= target_threshold
            # FP(i): pred >= i and target < target_threshold
            tp_by_pred_bin = joint[:, target_threshold_bin:].sum(dim=1)
            fp_by_pred_bin = joint[:, :target_threshold_bin].sum(dim=1)

            # Sweep thresholds from high to low prediction bins.
            tp_cum = torch.cumsum(torch.flip(tp_by_pred_bin, dims=[0]), dim=0)
            fp_cum = torch.cumsum(torch.flip(fp_by_pred_bin, dims=[0]), dim=0)

            precision = tp_cum / (tp_cum + fp_cum + self.epsilon)
            recall = tp_cum / (positive_total + self.epsilon)

            # Add the origin point for stable trapezoidal integration.
            recall_curve = torch.cat([torch.zeros(1, device=recall.device), recall], dim=0)
            precision_curve = torch.cat([precision[:1], precision], dim=0)

            auprc = torch.trapz(precision_curve, recall_curve)
            auprc = torch.clamp(auprc, 0.0, 1.0)
            auprc_per_track.append(auprc.item())

        return {self._format_output_key('auprc'): np.mean(auprc_per_track)}


class RegionR2(nn.Module):
    """
    Region-based R² metric for diagnostic purposes.

    Computes R² by summing predictions/targets within genomic regions
    and calculating R² across all regions.

    Args:
        regions_bed: Path to BED file with regions (chrom, start, end)
        bin_size: Bin size in bp (must match data), default 32
        aggregation_mode: 'per_track', 'aggregate', or 'both' (default)
        overlap_mode: 'strict' (full containment only) or 'partial' (weighted by overlap)
        filter_chroms: Optional list of chromosomes to filter regions to
        metric_prefix: Prefix for metric names in output dict (e.g., 'val' or 'train')
        epsilon: Numerical stability constant
        accelerator: Optional Accelerator instance for distributed gather
    """

    def __init__(
        self,
        regions_bed: str,
        bin_size: int = 32,
        aggregation_mode: str = 'both',
        overlap_mode: str = 'strict',
        filter_chroms: Optional[List[str]] = None,
        metric_prefix: str = 'region_r2',
        epsilon: float = 1e-8,
        accelerator=None,
    ):
        super().__init__()
        self.bin_size = bin_size
        self.aggregation_mode = aggregation_mode
        self.overlap_mode = overlap_mode
        self.metric_prefix = metric_prefix
        self.epsilon = epsilon
        self.accelerator = accelerator

        # Load regions from BED file
        self.regions = self._load_regions(regions_bed)

        # Filter regions by chromosome if specified
        if filter_chroms is not None:
            self.active_region_ids = self._filter_regions_by_chroms(filter_chroms)
        else:
            self.active_region_ids = set(r['id'] for r in self.regions)

        # Build bin-to-region mapping (chromosome -> interval tree)
        self.region_index, self.region_starts = self._build_region_index()

        # Reset accumulation state
        self.reset()

    def _load_regions(self, bed_file: str) -> List[Dict]:
        """Load regions from BED file."""
        regions = []
        with open(bed_file, 'r') as f:
            for line_idx, line in enumerate(f):
                if line.startswith('#'):
                    continue
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                regions.append({
                    'id': line_idx,
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2])
                })
        return regions

    def _build_region_index(self) -> Tuple[Dict[str, List[Dict]], Dict[str, List[int]]]:
        """Build chromosome -> sorted regions mapping for fast lookup."""
        from collections import defaultdict
        index = defaultdict(list)
        starts_index = defaultdict(list)

        # Only index active regions
        for region in self.regions:
            if region['id'] in self.active_region_ids:
                index[region['chrom']].append(region)

        # Sort by start position for binary search
        for chrom in index:
            index[chrom].sort(key=lambda r: r['start'])
            starts_index[chrom] = [r['start'] for r in index[chrom]]

        return dict(index), dict(starts_index)

    def _filter_regions_by_chroms(self, chroms: List[str]) -> Set[int]:
        """Return set of region IDs that are on the specified chromosomes."""
        chrom_set = set(chroms)
        return {r['id'] for r in self.regions if r['chrom'] in chrom_set}

    def _find_overlapping_regions(
        self,
        chrom: str,
        bin_start: int,
        bin_end: int
    ) -> List[Tuple[int, float]]:
        """
        Find regions overlapping a bin.

        Args:
            chrom: Chromosome name
            bin_start: Bin start position (bp)
            bin_end: Bin end position (bp)

        Returns:
            List of (region_id, overlap_fraction) tuples
        """
        if chrom not in self.region_index:
            return []

        regions = self.region_index[chrom]
        region_starts = self.region_starts[chrom]
        overlaps = []

        # Use binary search to jump to the first region that could overlap this bin.
        # We walk backwards from there until regions can no longer overlap, which keeps
        # the per-bin work proportional to the local overlap density instead of the
        # total number of regions.
        start_idx = bisect_left(region_starts, bin_end)
        for idx in range(start_idx - 1, -1, -1):
            region = regions[idx]
            if region['end'] <= bin_start:
                break  # no earlier region can overlap

            # Compute overlap
            overlap_start = max(bin_start, region['start'])
            overlap_end = min(bin_end, region['end'])
            overlap_len = overlap_end - overlap_start

            if overlap_len > 0:
                if self.overlap_mode == 'strict':
                    # Only include if bin is fully contained
                    if bin_start >= region['start'] and bin_end <= region['end']:
                        overlaps.append((region['id'], 1.0))
                else:  # partial
                    fraction = overlap_len / self.bin_size
                    overlaps.append((region['id'], fraction))

        return overlaps

    def reset(self):
        """Reset accumulated state."""
        # Store per-region sums: region_id -> tensor of shape (num_tracks,)
        self.region_sums_pred = {}
        self.region_sums_target = {}
        self.num_tracks = None

    def _find_overlapping_regions_batch(
        self,
        chrom: str,
        bin_starts: np.ndarray,
        bin_ends: np.ndarray
    ) -> List[List[Tuple[int, float]]]:
        """
        Find regions overlapping multiple bins efficiently.

        Args:
            chrom: Chromosome name
            bin_starts: Array of bin start positions
            bin_ends: Array of bin end positions

        Returns:
            List of lists of (region_id, overlap_fraction) tuples, one per bin
        """
        if chrom not in self.region_index:
            return [[] for _ in range(len(bin_starts))]

        regions = self.region_index[chrom]
        region_starts = self.region_starts[chrom]
        num_regions = len(regions)

        if num_regions == 0:
            return [[] for _ in range(len(bin_starts))]

        # Pre-extract region data for faster access
        region_ends = np.array([r['end'] for r in regions])
        region_starts_arr = np.array(region_starts)
        region_ids = [r['id'] for r in regions]

        results = []
        for i in range(len(bin_starts)):
            bin_start = bin_starts[i]
            bin_end = bin_ends[i]
            overlaps = []

            # Binary search for first region that could overlap
            start_idx = bisect_left(region_starts, bin_end)
            for idx in range(start_idx - 1, -1, -1):
                if region_ends[idx] <= bin_start:
                    break

                # Compute overlap
                overlap_start = max(bin_start, region_starts_arr[idx])
                overlap_end = min(bin_end, region_ends[idx])
                overlap_len = overlap_end - overlap_start

                if overlap_len > 0:
                    if self.overlap_mode == 'strict':
                        if bin_start >= region_starts_arr[idx] and bin_end <= region_ends[idx]:
                            overlaps.append((region_ids[idx], 1.0))
                    else:
                        fraction = overlap_len / self.bin_size
                        overlaps.append((region_ids[idx], fraction))

            results.append(overlaps)

        return results

    def update(
        self,
        y_pred: torch.Tensor,
        y_true: torch.Tensor,
        metadata: Optional[Dict] = None
    ):
        """
        Update metric with batch predictions.

        Uses vectorized bin coordinate computation for improved performance.

        Args:
            y_pred: Predictions (batch, bins, tracks)
            y_true: Targets (batch, bins, tracks)
            metadata: Dict with 'chrom', 'start', 'end', 'bin_size'
                     Values can be single items or lists/tensors (one per sample in batch)
        """
        if metadata is None:
            warnings.warn("RegionR2: No metadata provided, skipping batch", stacklevel=2)
            return

        with torch.no_grad():
            y_pred = y_pred.float()
            y_true = y_true.float()

            batch_size, num_bins, num_tracks = y_pred.shape

            # Initialize num_tracks on first call
            if self.num_tracks is None:
                self.num_tracks = num_tracks

            # Extract metadata - handle both single values and batched values
            chroms = metadata.get('chrom')
            starts = metadata.get('start')
            bin_size = metadata.get('bin_size', self.bin_size)

            if chroms is None or starts is None:
                return

            # Handle bin_size - could be tensor or scalar
            if isinstance(bin_size, torch.Tensor):
                bin_size = bin_size[0].item() if bin_size.numel() > 0 else self.bin_size
            elif isinstance(bin_size, (list, tuple)):
                bin_size = bin_size[0] if len(bin_size) > 0 else self.bin_size

            # Convert starts to numpy array for vectorized operations
            if isinstance(starts, torch.Tensor):
                starts_arr = starts.cpu().numpy()
            elif isinstance(starts, (list, tuple)):
                starts_arr = np.array(starts)
            else:
                starts_arr = np.array([starts] * batch_size)

            # Convert chroms to list if not already
            if isinstance(chroms, str):
                chroms = [chroms] * batch_size
            elif not isinstance(chroms, (list, tuple)):
                chroms = list(chroms)

            # Pre-compute bin offsets (vectorized)
            bin_offsets = np.arange(num_bins, dtype=np.int64) * int(bin_size)

            # Process each sample in batch
            for sample_idx in range(batch_size):
                sample_chrom = chroms[sample_idx] if sample_idx < len(chroms) else chroms[0]
                sample_start = starts_arr[sample_idx] if sample_idx < len(starts_arr) else starts_arr[0]

                # Skip if chromosome not in our filtered regions
                if sample_chrom not in self.region_index:
                    continue

                # Compute all bin coordinates for this sample (vectorized)
                bin_starts = (sample_start + bin_offsets).astype(np.int64)
                bin_ends = bin_starts + int(bin_size)

                # Find overlapping regions for all bins at once
                all_overlaps = self._find_overlapping_regions_batch(
                    sample_chrom, bin_starts, bin_ends
                )

                # Process bins with overlaps
                for bin_idx, overlaps in enumerate(all_overlaps):
                    if not overlaps:
                        continue

                    # Get predictions and targets for this bin
                    bin_pred = y_pred[sample_idx, bin_idx, :]
                    bin_target = y_true[sample_idx, bin_idx, :]

                    # Skip bins where any track has NaN target
                    if not torch.isfinite(bin_target).all():
                        continue

                    # Add to region sums
                    for region_id, fraction in overlaps:
                        weighted_pred = bin_pred * fraction
                        weighted_target = bin_target * fraction

                        if region_id not in self.region_sums_pred:
                            self.region_sums_pred[region_id] = torch.zeros(
                                num_tracks, device=y_pred.device
                            )
                            self.region_sums_target[region_id] = torch.zeros(
                                num_tracks, device=y_true.device
                            )

                        self.region_sums_pred[region_id] += weighted_pred
                        self.region_sums_target[region_id] += weighted_target

    def compute(self) -> Dict[str, float]:
        """
        Compute final R² values.

        Returns:
            Dictionary with metric_prefix + '_per_track' and/or '_aggregate'
        """
        if not self.region_sums_pred:
            # No data accumulated, return zeros
            results = {}
            if self.aggregation_mode in ('per_track', 'both'):
                results[f'{self.metric_prefix}_per_track'] = 0.0
            if self.aggregation_mode in ('aggregate', 'both'):
                results[f'{self.metric_prefix}_aggregate'] = 0.0
            return results

        # Stack all region sums into tensors
        region_ids = sorted(self.region_sums_pred.keys())

        sum_preds = torch.stack([self.region_sums_pred[rid] for rid in region_ids])  # (num_regions, tracks)
        sum_targets = torch.stack([self.region_sums_target[rid] for rid in region_ids])  # (num_regions, tracks)

        # Handle distributed training: gather region sums from all processes
        if self.accelerator is not None and self.accelerator.num_processes > 1:
            # Gather region IDs and sums from all processes
            # Each process may have different regions, so we need to aggregate
            all_region_ids = self._gather_all_region_ids(region_ids)
            sum_preds, sum_targets = self._gather_and_aggregate_sums(
                all_region_ids, region_ids, sum_preds, sum_targets
            )

        results = self._compute_r2(sum_preds, sum_targets, self.metric_prefix)

        return results

    def _gather_all_region_ids(self, local_region_ids: List[int]) -> List[int]:
        """Gather all unique region IDs from all processes."""
        # Convert to tensor for gathering
        local_ids_tensor = torch.tensor(local_region_ids, dtype=torch.long,
                                         device=self.accelerator.device)
        # Pad to max length across processes
        max_len = torch.tensor([len(local_region_ids)], device=self.accelerator.device)
        all_max_lens = self.accelerator.gather(max_len)
        global_max_len = all_max_lens.max().item()

        # Pad local tensor
        padded_ids = torch.full((global_max_len,), -1, dtype=torch.long,
                                 device=self.accelerator.device)
        padded_ids[:len(local_region_ids)] = local_ids_tensor

        # Gather from all processes
        all_padded_ids = self.accelerator.gather(padded_ids)  # (num_processes * global_max_len,)

        # Get unique valid IDs
        valid_ids = all_padded_ids[all_padded_ids >= 0].unique().cpu().tolist()
        return sorted(valid_ids)

    def _gather_and_aggregate_sums(
        self,
        all_region_ids: List[int],
        local_region_ids: List[int],
        local_preds: torch.Tensor,
        local_targets: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Gather and aggregate sums from all processes for each region."""
        num_tracks = local_preds.shape[1] if local_preds.numel() > 0 else self.num_tracks
        device = self.accelerator.device

        # Create tensors for all regions (fill with zeros for regions not on this process)
        num_regions = len(all_region_ids)
        all_preds = torch.zeros((num_regions, num_tracks), dtype=torch.float32, device=device)
        all_targets = torch.zeros((num_regions, num_tracks), dtype=torch.float32, device=device)

        # Map local region IDs to global indices
        local_id_set = set(local_region_ids)
        for global_idx, region_id in enumerate(all_region_ids):
            if region_id in local_id_set:
                local_idx = local_region_ids.index(region_id)
                all_preds[global_idx] = local_preds[local_idx]
                all_targets[global_idx] = local_targets[local_idx]

        # Sum across processes (each process contributes its partial sums)
        gathered_preds = self.accelerator.reduce(all_preds, reduction="sum")
        gathered_targets = self.accelerator.reduce(all_targets, reduction="sum")

        return gathered_preds, gathered_targets

    def _compute_r2(
        self,
        preds: torch.Tensor,
        targets: torch.Tensor,
        prefix: str
    ) -> Dict[str, float]:
        """
        Helper to compute both per-track and aggregate R².

        Args:
            preds: Predictions (num_regions, num_tracks)
            targets: Targets (num_regions, num_tracks)
            prefix: Prefix for metric names

        Returns:
            Dictionary with R² metrics
        """
        results = {}

        # Per-track R²
        if self.aggregation_mode in ('per_track', 'both'):
            # Compute R² independently for each track
            mean_target = targets.mean(dim=0, keepdim=True)  # (1, tracks)
            ss_res = ((targets - preds) ** 2).sum(dim=0)  # (tracks,)
            ss_tot = ((targets - mean_target) ** 2).sum(dim=0)  # (tracks,)

            # Avoid division by zero and clamp to [0, 1] for interpretability
            variance_mask = ss_tot > MIN_VARIANCE_THRESHOLD
            r2_per_track = torch.zeros_like(ss_res)
            r2_per_track[variance_mask] = 1 - ss_res[variance_mask] / (ss_tot[variance_mask] + self.epsilon)
            r2_per_track = torch.clamp(r2_per_track, 0.0, 1.0)

            # Average across tracks
            results[f'{prefix}_per_track'] = r2_per_track.mean().item()

        # Aggregate R² (all tracks together)
        if self.aggregation_mode in ('aggregate', 'both'):
            # Flatten tracks dimension
            preds_flat = preds.reshape(-1)  # (num_regions * tracks,)
            targets_flat = targets.reshape(-1)

            mean_target = targets_flat.mean()
            ss_res = ((targets_flat - preds_flat) ** 2).sum()
            ss_tot = ((targets_flat - mean_target) ** 2).sum()

            if ss_tot > MIN_VARIANCE_THRESHOLD:
                r2_aggregate = 1 - ss_res / (ss_tot + self.epsilon)
                r2_aggregate = torch.clamp(r2_aggregate, 0.0, 1.0)
            else:
                r2_aggregate = torch.tensor(0.0, device=preds.device)
            results[f'{prefix}_aggregate'] = r2_aggregate.item()

        return results


def get_metrics(
    metric_names: Optional[List[str]] = None,
    accelerator=None,
    genome_wide_pearson_config: Optional[Dict] = None,
    ppv_config: Optional[Dict] = None,
    quantile_overlap_config: Optional[Dict] = None,
    auprc_config: Optional[Dict] = None,
    track_names: Optional[List[str]] = None,
) -> List[nn.Module]:
    """
    Factory function to get a list of metrics by name.

    Args:
        metric_names: List of metric names. Options: 'pearson', 'spearman', 'r2', 'mse',
                     'genome_wide_pearson', 'ppv', 'quantile_overlap', 'auprc'
        accelerator: Optional Accelerator instance for distributed training support
                    (required for genome_wide_pearson and histogram metrics in multi-GPU setups)
        genome_wide_pearson_config: Optional dict with genome-wide Pearson config:
                   - track_subsets: optional list of selectors for additional
                     subset-specific metrics:
                     - name: suffix for output keys (e.g., 'rna')
                     - regex: regex string or list of regex strings
                     - tracks: exact track name string or list of names
        ppv_config: Optional dict with PPV metric configuration:
                   - quantile: float (default 0.98)
                   - smooth_bins: int (default 10)
                   - num_hist_bins: int (default 1000)
                   - value_max: float (default 1000.0)
                   - use_log_scale: bool (default False, set True for raw signal)
        quantile_overlap_config: Optional dict with quantile overlap metric configuration:
                   - quantile: float (default 0.98)
                   - smooth_bins: int (default 1, no smoothing)
                   - num_hist_bins: int (default 1000)
                   - value_max: float (default 1000.0)
                   - use_log_scale: bool (default False, set True for raw signal)
                   - track_subsets: optional list of selectors for additional
                     subset-specific metrics:
                     - name: suffix for output keys (e.g., 'atac')
                     - regex: regex string or list of regex strings
                     - tracks: exact track name string or list of names
        auprc_config: Optional dict with AUPRC metric configuration:
                   - quantile: float (default 0.98, defines positives from targets)
                   - smooth_bins: int (default 1, no smoothing)
                   - num_hist_bins: int (default 1000)
                   - value_max: float (default 1000.0)
                   - use_log_scale: bool (default False, set True for raw signal)
                   - track_subsets: optional list of selectors for additional
                     subset-specific metrics:
                     - name: suffix for output keys (e.g., 'atac')
                     - regex: regex string or list of regex strings
                     - tracks: exact track name string or list of names

    Returns:
        List of metric modules
    """
    if metric_names is None:
        metric_names = ['pearson']

    if genome_wide_pearson_config is None:
        genome_wide_pearson_config = {}

    if ppv_config is None:
        ppv_config = {}

    if quantile_overlap_config is None:
        quantile_overlap_config = {}

    if auprc_config is None:
        auprc_config = {}

    def _as_list(value: Any) -> List[Any]:
        if value is None:
            return []
        if isinstance(value, (list, tuple)):
            return list(value)
        return [value]

    def _sanitize_suffix(raw_name: str) -> str:
        clean = re.sub(r'[^0-9a-zA-Z_]+', '_', str(raw_name)).strip('_').lower()
        return clean or "subset"

    def _resolve_subset_config(
        selector: Dict[str, Any],
        available_track_names: Optional[List[str]],
        metric_label: str,
    ) -> Optional[Tuple[str, List[int]]]:
        if not isinstance(selector, dict):
            warnings.warn(f"{metric_label}: invalid track_subsets entry (expected dict), skipping")
            return None

        if not available_track_names:
            warnings.warn(
                f"{metric_label}: track_subsets configured but track names unavailable; skipping subset metric",
            )
            return None

        subset_name = _sanitize_suffix(selector.get('name', 'subset'))
        exact_tracks = set(_as_list(selector.get('tracks', selector.get('include_tracks'))))
        regex_patterns = _as_list(selector.get('regex', selector.get('include_regex')))

        if not exact_tracks and not regex_patterns:
            warnings.warn(
                f"{metric_label}.{subset_name}: no 'tracks' or 'regex' provided; skipping subset metric",
            )
            return None

        compiled_patterns = []
        for pattern in regex_patterns:
            try:
                compiled_patterns.append(re.compile(str(pattern)))
            except re.error as e:
                warnings.warn(
                    f"{metric_label}.{subset_name}: invalid regex '{pattern}': {e}; skipping this pattern"
                )

        selected_indices: List[int] = []
        for idx, track_name in enumerate(available_track_names):
            in_exact = track_name in exact_tracks
            in_regex = any(rx.search(track_name) for rx in compiled_patterns)
            if in_exact or in_regex:
                selected_indices.append(idx)

        if not selected_indices:
            warnings.warn(
                f"{metric_label}.{subset_name}: selector matched no tracks; skipping subset metric"
            )
            return None

        return subset_name, selected_indices

    metrics = []
    for name in metric_names:
        name_lower = name.lower()
        if name_lower == 'pearson':
            metrics.append(('pearson', PearsonCorrelation()))
        elif name_lower == 'spearman':
            metrics.append(('spearman', SpearmanCorrelation()))
        elif name_lower == 'r2':
            metrics.append(('r2', R2Score()))
        elif name_lower == 'mse':
            metrics.append(('mse', MSE()))
        elif name_lower == 'genome_wide_pearson':
            metrics.append(('genome_wide_pearson', GenomeWidePearson(accelerator=accelerator)))
            for selector in genome_wide_pearson_config.get('track_subsets', []):
                resolved = _resolve_subset_config(
                    selector,
                    track_names,
                    metric_label='genome_wide_pearson',
                )
                if resolved is None:
                    continue
                subset_name, subset_indices = resolved
                metrics.append((
                    f'genome_wide_pearson_{subset_name}',
                    GenomeWidePearson(
                        accelerator=accelerator,
                        track_indices=subset_indices,
                    ),
                ))
        elif name_lower == 'ppv':
            base_metric = GenomeWidePPV(
                quantile=ppv_config.get('quantile', 0.98),
                smooth_bins=ppv_config.get('smooth_bins', 10),
                num_hist_bins=ppv_config.get('num_hist_bins', 1000),
                value_max=ppv_config.get('value_max', 1000.0),
                use_log_scale=ppv_config.get('use_log_scale', False),
                accelerator=accelerator
            )
            metrics.append(('ppv', base_metric))

            for selector in ppv_config.get('track_subsets', []):
                resolved = _resolve_subset_config(selector, track_names, metric_label='ppv')
                if resolved is None:
                    continue
                subset_name, subset_indices = resolved
                metrics.append((
                    f'ppv_{subset_name}',
                    GenomeWidePPV(
                        quantile=ppv_config.get('quantile', 0.98),
                        smooth_bins=ppv_config.get('smooth_bins', 10),
                        num_hist_bins=ppv_config.get('num_hist_bins', 1000),
                        value_max=ppv_config.get('value_max', 1000.0),
                        use_log_scale=ppv_config.get('use_log_scale', False),
                        accelerator=accelerator,
                        track_indices=subset_indices,
                        metric_suffix=subset_name,
                    ),
                ))
        elif name_lower == 'quantile_overlap':
            overlap_metric = GenomeWideQuantileOverlap(
                quantile=quantile_overlap_config.get('quantile', 0.98),
                smooth_bins=quantile_overlap_config.get('smooth_bins', 1),
                num_hist_bins=quantile_overlap_config.get('num_hist_bins', 1000),
                value_max=quantile_overlap_config.get('value_max', 1000.0),
                use_log_scale=quantile_overlap_config.get('use_log_scale', False),
                accelerator=accelerator
            )
            metrics.append(('quantile_overlap', overlap_metric))

            for selector in quantile_overlap_config.get('track_subsets', []):
                resolved = _resolve_subset_config(selector, track_names, metric_label='quantile_overlap')
                if resolved is None:
                    continue
                subset_name, subset_indices = resolved
                metrics.append((
                    f'quantile_overlap_{subset_name}',
                    GenomeWideQuantileOverlap(
                        quantile=quantile_overlap_config.get('quantile', 0.98),
                        smooth_bins=quantile_overlap_config.get('smooth_bins', 1),
                        num_hist_bins=quantile_overlap_config.get('num_hist_bins', 1000),
                        value_max=quantile_overlap_config.get('value_max', 1000.0),
                        use_log_scale=quantile_overlap_config.get('use_log_scale', False),
                        accelerator=accelerator,
                        track_indices=subset_indices,
                        metric_suffix=subset_name,
                    ),
                ))
        elif name_lower == 'auprc':
            auprc_metric = GenomeWideAUPRC(
                quantile=auprc_config.get('quantile', 0.98),
                smooth_bins=auprc_config.get('smooth_bins', 1),
                num_hist_bins=auprc_config.get('num_hist_bins', 1000),
                value_max=auprc_config.get('value_max', 1000.0),
                use_log_scale=auprc_config.get('use_log_scale', False),
                accelerator=accelerator,
            )
            metrics.append(('auprc', auprc_metric))

            for selector in auprc_config.get('track_subsets', []):
                resolved = _resolve_subset_config(selector, track_names, metric_label='auprc')
                if resolved is None:
                    continue
                subset_name, subset_indices = resolved
                metrics.append((
                    f'auprc_{subset_name}',
                    GenomeWideAUPRC(
                        quantile=auprc_config.get('quantile', 0.98),
                        smooth_bins=auprc_config.get('smooth_bins', 1),
                        num_hist_bins=auprc_config.get('num_hist_bins', 1000),
                        value_max=auprc_config.get('value_max', 1000.0),
                        use_log_scale=auprc_config.get('use_log_scale', False),
                        accelerator=accelerator,
                        track_indices=subset_indices,
                        metric_suffix=subset_name,
                    ),
                ))
        else:
            print(f"Warning: Unknown metric '{name}', skipping")

    return metrics
