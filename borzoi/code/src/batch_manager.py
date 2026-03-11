"""Adaptive batch size management with OOM handling."""

import math
from typing import Optional


class AdaptiveBatchManager:
    """Manages adaptive batch size with optional gradient accumulation for OOM recovery.

    This class provides automatic batch size reduction when CUDA OOM errors occur,
    with optional gradient accumulation to maintain effective batch size.
    """

    def __init__(
        self,
        initial_batch_size: int,
        min_batch_size: int = 1,
        reduction_factor: float = 0.5,
        use_gradient_accumulation: bool = True,
        max_accumulation_steps: int = 16,
        recovery_delay: float = 2.0,
        max_batch_size: Optional[int] = None,
        growth_factor: float = 1.2,
        growth_interval: int = 500,
        enable_growth: bool = True,
    ):
        """Initialize the adaptive batch manager.

        Args:
            initial_batch_size: Starting batch size.
            min_batch_size: Minimum allowed batch size.
            reduction_factor: Factor to reduce batch size on OOM (e.g., 0.5 = halve).
                Must be in range (0, 1) exclusive.
            use_gradient_accumulation: Whether to use gradient accumulation to maintain
                effective batch size after reduction.
            max_accumulation_steps: Maximum gradient accumulation steps.
            recovery_delay: Seconds to wait after OOM before retrying.
            max_batch_size: Maximum batch size for growth (defaults to initial).
            growth_factor: Factor to increase batch size after stability. Must be > 1.
            growth_interval: Steps of stability before attempting growth.
            enable_growth: Whether to allow batch size growth after OOM recovery.

        Raises:
            ValueError: If reduction_factor is not in (0, 1) or growth_factor <= 1.
        """
        # Validate reduction_factor
        if not (0 < reduction_factor < 1):
            raise ValueError(
                f"reduction_factor must be in range (0, 1), got {reduction_factor}. "
                f"Use values like 0.5 to halve batch size on OOM."
            )

        # Validate growth_factor
        if growth_factor <= 1:
            raise ValueError(
                f"growth_factor must be > 1, got {growth_factor}. "
                f"Use values like 1.2 to grow batch size by 20%."
            )

        # Validate batch sizes
        if initial_batch_size < 1:
            raise ValueError(f"initial_batch_size must be >= 1, got {initial_batch_size}")
        if min_batch_size < 1:
            raise ValueError(f"min_batch_size must be >= 1, got {min_batch_size}")

        self.original_batch_size = initial_batch_size
        self.current_batch_size = initial_batch_size
        self.min_batch_size = max(1, min_batch_size)
        self.reduction_factor = reduction_factor
        self.use_gradient_accumulation = use_gradient_accumulation
        self.max_accumulation_steps = max_accumulation_steps
        self.recovery_delay = recovery_delay
        self.max_batch_size = max_batch_size or initial_batch_size
        self.growth_factor = growth_factor
        self.growth_interval = growth_interval
        self.enable_growth = enable_growth

        self.gradient_accumulation_steps = 1
        self.oom_count = 0
        self.needs_loader_rebuild = False
        self.steps_since_oom = 0

    @property
    def effective_batch_size(self) -> int:
        """Effective batch size after accounting for gradient accumulation."""
        return self.current_batch_size * self.gradient_accumulation_steps

    def reduce_batch_size(self) -> bool:
        """Reduce batch size on OOM.

        Returns:
            True if reduction was successful, False if already at minimum.
        """
        new_batch_size = max(
            self.min_batch_size,
            int(self.current_batch_size * self.reduction_factor)
        )

        if new_batch_size >= self.current_batch_size:
            return False

        self.current_batch_size = new_batch_size
        self.oom_count += 1
        self.needs_loader_rebuild = True
        self.steps_since_oom = 0

        # Adjust gradient accumulation to maintain effective batch size
        if self.use_gradient_accumulation:
            target_effective = self.original_batch_size
            new_accum = min(
                self.max_accumulation_steps,
                max(1, target_effective // self.current_batch_size)
            )
            self.gradient_accumulation_steps = new_accum

        return True

    def register_successful_step(self):
        """Track stable steps since the last OOM event."""
        self.steps_since_oom += 1

    def maybe_increase_batch_size(self) -> Optional[int]:
        """Cautiously grow batch size after sustained stability.

        Returns:
            New batch size if growth occurred, None otherwise.
        """
        if not self.enable_growth:
            return None

        if self.current_batch_size >= self.max_batch_size:
            return None

        if self.steps_since_oom < self.growth_interval:
            return None

        proposed_batch = min(
            self.max_batch_size,
            max(self.current_batch_size + 1, int(math.ceil(self.current_batch_size * self.growth_factor)))
        )

        if proposed_batch <= self.current_batch_size:
            return None

        self.current_batch_size = proposed_batch
        self.needs_loader_rebuild = True
        self.steps_since_oom = 0

        if self.use_gradient_accumulation:
            target_effective = self.original_batch_size
            new_accum = min(
                self.max_accumulation_steps,
                max(1, target_effective // self.current_batch_size)
            )
            self.gradient_accumulation_steps = new_accum

        return self.current_batch_size

    def should_step(self, step: int) -> bool:
        """Check if optimizer should step at this step.

        Args:
            step: Current global step (0-indexed, will be incremented after this call).

        Returns:
            True if optimizer should step (accumulation boundary reached).
        """
        # step is 0-indexed; we want to step on the last step of each accumulation window
        # e.g., with accum_steps=4: step at 3, 7, 11, ... (i.e., when (step+1) % accum == 0)
        return (step + 1) % self.gradient_accumulation_steps == 0

    def reset_for_stage(
        self,
        new_batch_size: int,
        gradient_accumulation_steps: Optional[int] = None,
        target_effective_batch_size: Optional[int] = None,
    ):
        """Reset manager for a new training stage (e.g., fine-tuning).

        Args:
            new_batch_size: New initial batch size for the stage.
            gradient_accumulation_steps: Explicit accumulation steps to use.
                If provided, overrides automatic calculation.
            target_effective_batch_size: Target effective batch size to achieve
                via gradient accumulation. If provided and gradient_accumulation_steps
                is not set, calculates accumulation steps as:
                target_effective_batch_size // new_batch_size.
                Capped by max_accumulation_steps.
        """
        self.original_batch_size = new_batch_size
        self.current_batch_size = new_batch_size
        self.max_batch_size = max(self.max_batch_size, new_batch_size)
        self.needs_loader_rebuild = False
        self.steps_since_oom = 0

        # Determine gradient accumulation steps
        if gradient_accumulation_steps is not None:
            # Explicit accumulation steps provided
            self.gradient_accumulation_steps = min(
                self.max_accumulation_steps,
                max(1, gradient_accumulation_steps)
            )
        elif target_effective_batch_size is not None:
            # Calculate from target effective batch size
            self.gradient_accumulation_steps = min(
                self.max_accumulation_steps,
                max(1, target_effective_batch_size // new_batch_size)
            )
        else:
            # Default: no accumulation
            self.gradient_accumulation_steps = 1

    def get_status_dict(self) -> dict:
        """Get current status for logging.

        Returns:
            Dict with batch_size, grad_accum_steps, effective_batch_size, oom_count.
        """
        return {
            "batch_size": self.current_batch_size,
            "grad_accum_steps": self.gradient_accumulation_steps,
            "effective_batch_size": self.effective_batch_size,
            "oom_count": self.oom_count,
        }
