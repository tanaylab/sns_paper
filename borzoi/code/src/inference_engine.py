"""
Inference Engine for Borzoi Model Predictions.

This module provides a clean separation between inference and perturbation/data generation.
The InferenceEngine handles all model predictions (single, batch, distributed) while
being agnostic to the source of the input sequences.

Architecture:
    InferenceEngine - handles model predictions
    AdaptiveBatchManager - handles OOM recovery with adaptive batch sizing
    SequencePerturbationEngine - handles sequence modifications (in perturbation_operations.py)
    Analysis classes - orchestrate the two

Usage:
    # Create engine
    engine = InferenceEngine(
        model=model,
        device=device,
        use_rc_average=True,
        mixed_precision=True,
        use_compile=True,  # Optional torch.compile() acceleration
    )

    # Single prediction
    pred = engine.predict(sequence)  # (seq_len, 4) -> (bins, tracks)

    # Batch prediction
    preds = engine.predict_batch(sequences)  # (batch, seq_len, 4) -> (batch, bins, tracks)

    # Batch prediction with adaptive OOM handling
    preds = engine.predict_batch_adaptive(sequences, batch_manager)

    # Distributed batch prediction (multi-GPU)
    preds = engine.predict_batch_distributed(sequences, rank, world_size)
"""

import time
import warnings
from typing import Optional, Union, List
import numpy as np
import torch
import torch.distributed as dist

# Enable optimized backends for faster inference
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.benchmark = True


# =============================================================================
# Adaptive Batch Size Manager (OOM Handling)
# =============================================================================

class AdaptiveBatchManager:
    """
    Manages adaptive batch size for OOM recovery during inference.

    This class tracks batch size adjustments and provides methods to reduce
    batch size when CUDA out-of-memory errors occur.

    Usage:
        manager = AdaptiveBatchManager(initial_batch_size=64)

        while not done:
            try:
                predictions = engine.predict_batch(sequences[:manager.current_batch_size])
                done = True
            except RuntimeError as e:
                if 'out of memory' in str(e).lower():
                    if not manager.reduce_batch_size():
                        raise  # Can't reduce further
                    torch.cuda.empty_cache()
                    time.sleep(manager.recovery_delay)
                else:
                    raise
    """

    def __init__(
        self,
        initial_batch_size: int,
        min_batch_size: int = 1,
        reduction_factor: float = 0.5,
        recovery_delay: float = 2.0,
    ):
        """
        Initialize the adaptive batch manager.

        Args:
            initial_batch_size: Starting batch size
            min_batch_size: Minimum batch size (won't reduce below this)
            reduction_factor: Factor to multiply batch size on OOM (e.g., 0.5 = halve)
            recovery_delay: Seconds to wait after OOM before retry
        """
        self.original_batch_size = initial_batch_size
        self.current_batch_size = initial_batch_size
        self.min_batch_size = max(1, min_batch_size)
        self.reduction_factor = reduction_factor
        self.recovery_delay = recovery_delay
        self.oom_count = 0

    def reduce_batch_size(self) -> bool:
        """
        Reduce batch size on OOM.

        Returns:
            True if batch size was reduced, False if already at minimum.
        """
        new_batch_size = max(
            self.min_batch_size,
            int(self.current_batch_size * self.reduction_factor)
        )

        if new_batch_size >= self.current_batch_size:
            # Can't reduce further
            return False

        self.current_batch_size = new_batch_size
        self.oom_count += 1
        return True

    def reset(self) -> None:
        """Reset to original batch size."""
        self.current_batch_size = self.original_batch_size
        self.oom_count = 0


# =============================================================================
# Inference Engine
# =============================================================================

class InferenceEngine:
    """
    Unified inference engine for Borzoi model predictions.

    Handles:
    - Single sequence prediction
    - Batch prediction with automatic OOM handling
    - Reverse complement averaging
    - Mixed precision inference
    - Multi-GPU distributed inference
    - Optional torch.compile() acceleration with warmup

    The engine is agnostic to the source of input sequences - they can come from
    perturbation operations, raw genomic data, or any other source.

    OOM Handling Strategy:
        The engine uses a single, unified OOM handling strategy in predict_batch().
        When OOM occurs, it automatically retries with progressively smaller batch
        sizes (batch/2, batch/4, ..., 1). The deprecated predict_batch_adaptive()
        method is kept for backwards compatibility but delegates to predict_batch().
    """

    def __init__(
        self,
        model: torch.nn.Module,
        device: torch.device,
        use_rc_average: bool = True,
        mixed_precision: bool = True,
        is_human: bool = True,
        use_compile: bool = False,
        compile_mode: str = "reduce-overhead",
        warmup: bool = True,
        warmup_seq_len: int = 524288,
    ):
        """
        Initialize the inference engine.

        Args:
            model: Loaded Borzoi model (already on device).
            device: Torch device for inference.
            use_rc_average: Whether to average forward and reverse complement predictions.
            mixed_precision: Whether to use mixed precision (fp16) inference.
            is_human: Whether the model is trained on human data (affects some model internals).
            use_compile: Whether to use torch.compile() for acceleration (10-30% speedup).
            compile_mode: torch.compile mode ('default', 'reduce-overhead', 'max-autotune').
            warmup: Whether to run a warmup forward pass (recommended for compiled models).
            warmup_seq_len: Sequence length for warmup pass (default: 524288 for Borzoi).
        """
        self.device = device
        self.use_rc_average = use_rc_average
        self.mixed_precision = mixed_precision
        self.is_human = is_human

        # Ensure model is in eval mode
        model.eval()

        # Optional torch.compile() for acceleration
        if use_compile and hasattr(torch, 'compile'):
            try:
                self.model = torch.compile(model, mode=compile_mode)
                self._compiled = True
            except Exception as e:
                warnings.warn(f"torch.compile() failed, using eager mode: {e}")
                self.model = model
                self._compiled = False
        else:
            self.model = model
            self._compiled = False

        # Pre-compute RC indices on device (avoids tensor creation per batch)
        self._rc_indices = torch.tensor([3, 2, 1, 0], device=device)

        # Track if warmed up
        self._warmed_up = False

        # Run warmup if requested (especially important for compiled models)
        if warmup and self._compiled:
            self._run_warmup(warmup_seq_len)

    def _run_warmup(self, seq_len: int) -> None:
        """
        Run a warmup forward pass to trigger JIT compilation.

        This is especially important for torch.compile() which does lazy
        compilation on the first forward pass. Running warmup here ensures
        the first real inference call is fast.

        Args:
            seq_len: Sequence length for the warmup pass.
        """
        try:
            # Create dummy input (batch=1, channels=4, seq_len)
            dummy_input = torch.zeros(1, 4, seq_len, device=self.device, dtype=torch.float32)

            with torch.inference_mode():
                if self.mixed_precision:
                    with torch.amp.autocast('cuda'):
                        _ = self.model(dummy_input, is_human=self.is_human)
                else:
                    _ = self.model(dummy_input, is_human=self.is_human)

            # Clear GPU cache after warmup
            torch.cuda.empty_cache()
            self._warmed_up = True
        except Exception as e:
            warnings.warn(f"Warmup failed: {e}. First inference call may be slow.")

    def predict(self, sequence: np.ndarray) -> np.ndarray:
        """
        Run inference on a single sequence.

        Args:
            sequence: One-hot encoded sequence of shape (seq_len, 4).

        Returns:
            Predictions of shape (num_bins, num_tracks).
        """
        batch = sequence[np.newaxis, ...]  # (1, seq_len, 4)
        predictions = self.predict_batch(batch)
        return predictions[0]  # (num_bins, num_tracks)

    def predict_batch(
        self,
        sequences: Union[np.ndarray, torch.Tensor],
        handle_oom: bool = True,
    ) -> np.ndarray:
        """
        Run batched inference with optional reverse complement averaging.

        Includes automatic OOM handling with adaptive batch size reduction.
        On OOM, tries progressively smaller batches (batch/2, batch/4, ...)
        before falling back to batch-size-1.

        Args:
            sequences: Batch of sequences, shape (batch, seq_len, 4).
            handle_oom: Whether to automatically handle OOM by reducing batch size.
                       Set to False if you want to handle OOM yourself.

        Returns:
            Predictions of shape (batch, num_bins, num_tracks).
        """
        # Convert NumPy to torch and pin memory for efficient non-blocking transfer
        if isinstance(sequences, np.ndarray):
            # Convert to tensor and immediately pin memory
            # This enables efficient non_blocking=True transfers to GPU
            sequences = torch.from_numpy(sequences).pin_memory()
        elif sequences.device.type == 'cpu':
            # Only pin if not already pinned (calling pin_memory on pinned tensor is a no-op but wasteful)
            if not sequences.is_pinned():
                sequences = sequences.pin_memory()

        # Try batch inference first, with adaptive retry on OOM if enabled
        try:
            return self._predict_batch_inner(sequences)
        except RuntimeError as e:
            if handle_oom and 'out of memory' in str(e).lower():
                return self._predict_batch_with_adaptive_retry(sequences)
            else:
                raise

    def _predict_batch_with_adaptive_retry(
        self,
        sequences: torch.Tensor,
        max_attempts: int = 5,
    ) -> np.ndarray:
        """
        Retry inference with progressively smaller batch sizes on OOM.

        Tries: full_batch → batch/2 → batch/4 → ... → batch_size=1

        Args:
            sequences: Batch of sequences (batch, seq_len, 4)
            max_attempts: Maximum number of reduction attempts

        Returns:
            Predictions of shape (batch, num_bins, num_tracks).
        """
        batch_size = sequences.shape[0]
        torch.cuda.empty_cache()

        # Try progressively smaller batches
        for attempt in range(max_attempts):
            try_batch_size = max(1, batch_size // (2 ** attempt))

            # Stop if we've reached batch_size=1 or batch too small to divide
            if try_batch_size < 4 and attempt > 0:
                # Just use batch_size=1 for final attempt
                try_batch_size = 1

            try:
                # Process in smaller batches
                results = []
                for start_idx in range(0, batch_size, try_batch_size):
                    end_idx = min(start_idx + try_batch_size, batch_size)
                    batch_chunk = sequences[start_idx:end_idx]
                    pred = self._predict_batch_inner(batch_chunk)
                    results.append(pred)

                    # Clear cache between chunks for memory efficiency
                    if try_batch_size > 1:
                        torch.cuda.empty_cache()

                return np.concatenate(results, axis=0)

            except RuntimeError as e:
                if 'out of memory' in str(e).lower():
                    torch.cuda.empty_cache()
                    # If this was already batch_size=1, we're truly out of options
                    if try_batch_size == 1:
                        raise RuntimeError(
                            f"OOM even with batch_size=1. Cannot process sequences of this length."
                        ) from e
                    # Otherwise, continue to next smaller batch size
                    continue
                else:
                    raise

        # Should never reach here due to the batch_size=1 raise above
        raise RuntimeError(f"Inference failed after {max_attempts} adaptive retry attempts")

    @torch.inference_mode()
    def _predict_batch_inner(
        self,
        sequences: torch.Tensor,
    ) -> np.ndarray:
        """
        Inner prediction function without OOM handling.

        Uses torch.inference_mode() which is faster than torch.no_grad()
        because it also disables view tracking and version counting.

        Args:
            sequences: Batch of sequences as torch tensor, shape (batch, seq_len, 4).
                      Should be in pinned memory for efficient non-blocking transfer.

        Returns:
            Predictions of shape (batch, num_bins, num_tracks).
        """
        # Non-blocking transfer for better GPU utilization
        # This is efficient because predict_batch() ensures input is in pinned memory
        sequences = sequences.to(self.device, non_blocking=True)
        sequences_cf = sequences.permute(0, 2, 1)

        if self.mixed_precision:
            with torch.amp.autocast('cuda'):
                predictions = self._forward_with_rc(sequences_cf)
        else:
            predictions = self._forward_with_rc(sequences_cf)

        # Handle tuple returns (some model variants return auxiliary outputs)
        if isinstance(predictions, tuple):
            predictions = predictions[0]

        # Model outputs (batch, tracks, bins); transpose to (batch, bins, tracks)
        predictions = predictions.permute(0, 2, 1)

        return predictions.cpu().numpy()

    def _forward_with_rc(self, sequences_cf: torch.Tensor) -> torch.Tensor:
        """
        Forward pass with optional reverse complement averaging.

        Args:
            sequences_cf: Channel-first sequences (batch, 4, seq_len)

        Returns:
            Predictions tensor
        """
        predictions = self.model(sequences_cf, is_human=self.is_human)
        if isinstance(predictions, tuple):
            predictions = predictions[0]

        if not self.use_rc_average:
            return predictions

        # Reverse complement: flip sequence and swap channels using pre-computed indices
        sequences_rc = sequences_cf.flip(2).index_select(1, self._rc_indices)
        predictions_rc = self.model(sequences_rc, is_human=self.is_human)
        if isinstance(predictions_rc, tuple):
            predictions_rc = predictions_rc[0]
        predictions_rc = predictions_rc.flip(2)

        return (predictions + predictions_rc) * 0.5

    def predict_batch_distributed(
        self,
        sequences: np.ndarray,
        rank: int,
        world_size: int,
    ) -> Optional[np.ndarray]:
        """
        Distributed batch prediction across multiple GPUs.

        Each GPU processes a shard of the sequences, then results are gathered
        to rank 0 using NCCL all_gather.

        Args:
            sequences: Full batch of sequences, shape (total_batch, seq_len, 4).
            rank: Current GPU rank.
            world_size: Total number of GPUs.

        Returns:
            On rank 0: Full predictions array, shape (total_batch, num_bins, num_tracks).
            On other ranks: None.
        """
        total_batch = len(sequences)

        # Distribute sequences across ranks
        seqs_per_rank = (total_batch + world_size - 1) // world_size
        start_idx = rank * seqs_per_rank
        end_idx = min(start_idx + seqs_per_rank, total_batch)

        # Process local shard
        local_sequences = sequences[start_idx:end_idx]
        local_num = len(local_sequences)

        if local_num > 0:
            local_predictions = self.predict_batch(local_sequences)
        else:
            # Edge case: this rank has no sequences
            # Need to determine output shape from a dummy prediction
            dummy_pred = self.predict(sequences[0])
            local_predictions = np.zeros((0,) + dummy_pred.shape, dtype=np.float32)

        # Gather results from all ranks
        if not dist.is_initialized():
            # Single GPU mode - just return local results
            return local_predictions

        # Pad local predictions for uniform tensor size in all_gather
        num_bins = local_predictions.shape[1] if local_num > 0 else 0
        num_tracks = local_predictions.shape[2] if local_num > 0 else 0

        # Broadcast shape info from rank 0
        if local_num > 0:
            shape_tensor = torch.tensor([num_bins, num_tracks], device=self.device)
        else:
            shape_tensor = torch.tensor([0, 0], device=self.device)

        # Gather shape info to handle edge cases
        gathered_shapes = [torch.zeros(2, dtype=torch.long, device=self.device) for _ in range(world_size)]
        dist.all_gather(gathered_shapes, shape_tensor)

        # Find the correct shape from a rank that has data
        for shape in gathered_shapes:
            if shape[0] > 0:
                num_bins = int(shape[0])
                num_tracks = int(shape[1])
                break

        # Pad to uniform size
        local_padded = np.zeros((seqs_per_rank, num_bins, num_tracks), dtype=np.float32)
        local_padded[:local_num] = local_predictions

        # Convert to GPU tensor for NCCL
        local_tensor = torch.from_numpy(local_padded).to(self.device)
        gathered_tensors = [torch.zeros_like(local_tensor) for _ in range(world_size)]

        dist.all_gather(gathered_tensors, local_tensor)

        # Only rank 0 assembles the result
        if rank != 0:
            return None

        # Concatenate and trim to actual size
        all_predictions = []
        for r in range(world_size):
            r_start = r * seqs_per_rank
            r_end = min(r_start + seqs_per_rank, total_batch)
            r_num = r_end - r_start
            if r_num > 0:
                all_predictions.append(gathered_tensors[r][:r_num].cpu().numpy())

        return np.concatenate(all_predictions, axis=0)

    def predict_pairs(
        self,
        reference_sequences: np.ndarray,
        perturbed_sequences: np.ndarray,
    ) -> tuple:
        """
        Predict both reference and perturbed sequences in a single batch call.

        This is more efficient than separate calls because it maximizes GPU utilization.

        Args:
            reference_sequences: Reference sequences, shape (n, seq_len, 4).
            perturbed_sequences: Perturbed sequences, shape (n, seq_len, 4).

        Returns:
            Tuple of (reference_predictions, perturbed_predictions),
            each of shape (n, num_bins, num_tracks).
        """
        n = len(reference_sequences)
        assert len(perturbed_sequences) == n, "Must have same number of reference and perturbed sequences"

        # Stack all sequences for single batch call
        all_sequences = np.concatenate([reference_sequences, perturbed_sequences], axis=0)

        # Single prediction call
        all_predictions = self.predict_batch(all_sequences)

        # Split back
        ref_predictions = all_predictions[:n]
        pert_predictions = all_predictions[n:]

        return ref_predictions, pert_predictions

    def predict_pairs_distributed(
        self,
        reference_sequences: np.ndarray,
        perturbed_sequences: np.ndarray,
        rank: int,
        world_size: int,
    ) -> Optional[tuple]:
        """
        Distributed prediction of reference and perturbed sequence pairs.

        Args:
            reference_sequences: Reference sequences, shape (n, seq_len, 4).
            perturbed_sequences: Perturbed sequences, shape (n, seq_len, 4).
            rank: Current GPU rank.
            world_size: Total number of GPUs.

        Returns:
            On rank 0: Tuple of (reference_predictions, perturbed_predictions).
            On other ranks: None.
        """
        n = len(reference_sequences)

        # Stack for efficient distributed processing
        all_sequences = np.concatenate([reference_sequences, perturbed_sequences], axis=0)

        # Distributed prediction
        all_predictions = self.predict_batch_distributed(all_sequences, rank, world_size)

        if all_predictions is None:
            return None

        # Split back
        ref_predictions = all_predictions[:n]
        pert_predictions = all_predictions[n:]

        return ref_predictions, pert_predictions

    def predict_batch_adaptive(
        self,
        sequences: Union[np.ndarray, torch.Tensor],
        batch_manager: Optional[AdaptiveBatchManager] = None,
        show_progress: bool = False,
        max_retries: int = 10,
    ) -> np.ndarray:
        """
        Run batched inference with adaptive OOM recovery.

        .. deprecated::
            This method is deprecated. Use predict_batch() instead, which now
            includes automatic OOM handling. This method is kept for backwards
            compatibility and delegates to predict_batch().

        Automatically reduces batch size on OOM errors and retries,
        providing robust inference even with limited GPU memory.

        Args:
            sequences: Batch of sequences, shape (batch, seq_len, 4).
            batch_manager: AdaptiveBatchManager instance (optional, for external batch control).
            show_progress: Whether to show progress bar.
            max_retries: Maximum number of OOM recovery attempts (unused, kept for compatibility).

        Returns:
            Predictions of shape (batch, num_bins, num_tracks).

        Raises:
            RuntimeError: If inference fails after max_retries attempts.
        """
        warnings.warn(
            "predict_batch_adaptive() is deprecated. Use predict_batch() instead, "
            "which now includes automatic OOM handling.",
            DeprecationWarning,
            stacklevel=2,
        )

        if isinstance(sequences, np.ndarray):
            sequences = torch.from_numpy(sequences)

        total_seqs = len(sequences)

        # If batch_manager provided, use it for chunking (for progress bar)
        if batch_manager is not None and show_progress:
            from tqdm import tqdm
            all_predictions = []
            batch_size = batch_manager.current_batch_size

            iterator = tqdm(
                range(0, total_seqs, batch_size),
                desc=f"Inference (batch={batch_size})",
                total=(total_seqs + batch_size - 1) // batch_size,
            )

            for start_idx in iterator:
                end_idx = min(start_idx + batch_size, total_seqs)
                batch_seqs = sequences[start_idx:end_idx]
                # Use predict_batch which handles OOM automatically
                batch_preds = self.predict_batch(batch_seqs)
                all_predictions.append(batch_preds)

            return np.concatenate(all_predictions, axis=0)
        else:
            # Just delegate to predict_batch which handles OOM
            return self.predict_batch(sequences)
