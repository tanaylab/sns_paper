"""
Simplified perturbation engine for Borzoi models.

This module provides a clean, table-based API for sequence perturbation analysis:
1. Region table: Defines regions of interest (ROIs)
2. Splice perturbation table: Defines transplantations (sequence replacements)
3. Sequence perturbation table: Defines fine-grained edits applied after transplantation

Workflow:
    regions = load_regions("regions.bed")
    splice_perts = load_splice_pert("splice.pert")
    seq_perts = load_seq_pert("seq.pert")

    analyzer = SimplePerturbationAnalyzer(model, genome_cache, device, ...)
    df = analyzer.compute(regions, splice_perts, seq_perts, track_indices, track_names)

    # Output: Parquet with columns like baseline, splice1, double_mut, observed

Key features:
- All coordinates are genomic (no relative coordinates)
- Aggregated perturbations: Group multiple edits by splice_id
- Wide-format output with {splice_id} columns
- Multi-track support with column suffixes
"""

import warnings
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
import pysam
import torch
import torch.nn as nn
from tqdm import tqdm

from .data_loaders import one_hot_encode
from .inference_engine import InferenceEngine


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class Region:
    """
    Defines a genomic region for analysis.

    Attributes:
        chrom: Chromosome name (e.g., 'chr1')
        start: Start position (0-based, inclusive) - genomic coordinate
        end: End position (0-based, exclusive) - genomic coordinate
        region_id: Unique identifier for this region
        normalize_to: If set, extract window of this size centered on region
                     (e.g., 524288 for Borzoi)
    """
    chrom: str
    start: int
    end: int
    region_id: str
    normalize_to: Optional[int] = None

    @property
    def center(self) -> int:
        """Center position of the region."""
        return (self.start + self.end) // 2

    @property
    def length(self) -> int:
        """Length of the region in bp."""
        return self.end - self.start

    def get_extraction_window(self, seq_len: int) -> Tuple[int, int]:
        """
        Get the genomic window to extract for this region.

        Args:
            seq_len: Model sequence length (e.g., 524288)

        Returns:
            (seq_start, seq_end) tuple of genomic coordinates
        """
        center = self.center

        if self.normalize_to:
            # Use normalized size
            half_len = self.normalize_to // 2
        else:
            # Use max of region length or seq_len
            effective_len = max(self.length, seq_len)
            half_len = effective_len // 2

        return (center - half_len, center + half_len)

    def __repr__(self) -> str:
        norm_str = f", norm={self.normalize_to}" if self.normalize_to else ""
        return f"Region({self.chrom}:{self.start}-{self.end}, id={self.region_id}{norm_str})"


@dataclass
class SplicePert:
    """
    Defines a transplantation operation - replacing a region with donor sequence.

    Attributes:
        chrom: Chromosome of target region (genomic)
        start: Start position of target (genomic, 0-based)
        end: End position of target (genomic, 0-based)
        id: Unique identifier for this splice operation (renamed from splice_id)
        chrom_src: Source chromosome (if using genomic donor)
        start_src: Source start position (if using genomic donor)
        end_src: Source end position (if using genomic donor)
        seq: Literal DNA sequence (alternative to genomic donor)
        region_id: Optional filter - only apply to regions with this ID
    """
    chrom: str
    start: int
    end: int
    id: str
    chrom_src: Optional[str] = None
    start_src: Optional[int] = None
    end_src: Optional[int] = None
    seq: Optional[str] = None
    region_id: Optional[str] = None

    def __post_init__(self):
        """Validate that exactly one donor specification is provided."""
        has_genomic = (
            self.chrom_src is not None
            and self.start_src is not None
            and self.end_src is not None
        )
        has_literal = self.seq is not None

        if not (has_genomic or has_literal):
            raise ValueError(
                f"SplicePert {self.id}: Must specify either genomic donor "
                "(chrom_src, start_src, end_src) or literal sequence (seq)"
            )
        if has_genomic and has_literal:
            raise ValueError(
                f"SplicePert {self.id}: Cannot specify both genomic donor and literal seq"
            )

    def get_donor_length(self) -> int:
        """Get the length of the donor sequence."""
        if self.seq is not None:
            return len(self.seq)
        else:
            return self.end_src - self.start_src

    @property
    def target_length(self) -> int:
        """Get the length of the target region."""
        return self.end - self.start

    def __repr__(self) -> str:
        if self.seq:
            donor = f"seq={self.seq[:20]}..."
        else:
            donor = f"{self.chrom_src}:{self.start_src}-{self.end_src}"
        return f"SplicePert({self.id}: {self.chrom}:{self.start}-{self.end} ← {donor})"


@dataclass
class SeqPert:
    """
    Defines a sequence perturbation applied after transplantation.

    Attributes:
        chrom: Chromosome of target region (genomic)
        start: Start position of target (genomic, 0-based)
        end: End position of target (genomic, 0-based)
        id: Which variant this applies to (renamed from splice_id)
        chrom_src: Source chromosome (if using genomic donor)
        start_src: Source start position (if using genomic donor)
        end_src: Source end position (if using genomic donor)
        seq: Literal DNA sequence (alternative to genomic donor)
        region_id: Optional filter - only apply to regions with this ID
    """
    chrom: str
    start: int
    end: int
    id: str
    chrom_src: Optional[str] = None
    start_src: Optional[int] = None
    end_src: Optional[int] = None
    seq: Optional[str] = None
    region_id: Optional[str] = None

    def __post_init__(self):
        """Validate that exactly one donor specification is provided."""
        has_genomic = (
            self.chrom_src is not None
            and self.start_src is not None
            and self.end_src is not None
        )
        has_literal = self.seq is not None

        if not (has_genomic or has_literal):
            pid = f"{self.id}:{self.start}"
            raise ValueError(
                f"SeqPert {pid}: Must specify either genomic donor "
                "(chrom_src, start_src, end_src) or literal sequence (seq)"
            )
        if has_genomic and has_literal:
            pid = f"{self.id}:{self.start}"
            raise ValueError(
                f"SeqPert {pid}: Cannot specify both genomic donor and literal seq"
            )

    def __repr__(self) -> str:
        if self.seq:
            donor = f"seq={self.seq[:20]}..."
        else:
            donor = f"{self.chrom_src}:{self.start_src}-{self.end_src}"
        return f"SeqPert({self.id}: {self.chrom}:{self.start}-{self.end} ← {donor})"


@dataclass
class VariantSpec:
    """
    Specification for a single variant to compute.

    Attributes:
        variant_id: Full identifier (maps 1:1 with input id, or "baseline")
        region: Region this variant belongs to
        splice_ops: List of splice operations to apply
        seq_ops: List of sequence perturbations to apply (after splices)
    """
    variant_id: str
    region: Region
    splice_ops: List[SplicePert] = field(default_factory=list)
    seq_ops: List[SeqPert] = field(default_factory=list)

    def __repr__(self) -> str:
        return f"VariantSpec({self.variant_id}, {len(self.splice_ops)} splices, {len(self.seq_ops)} perts)"


# =============================================================================
# Table Loaders
# =============================================================================

def load_regions(
    bed_file: str,
    normalize_to: Optional[int] = None,
) -> List[Region]:
    """
    Load region table from BED or CSV file.

    Expected format (tab or comma-separated):
        chrom  start  end  region_id  [normalize_to]

    Args:
        bed_file: Path to BED/CSV file
        normalize_to: Global normalization size (can be overridden per-region in 5th column)

    Returns:
        List of Region objects
    """
    # Detect separator
    sep = '\t' if bed_file.endswith('.bed') else ','

    # Read file - try with header first, fall back to no header
    df = pd.read_csv(bed_file, sep=sep, comment='#')

    # Check if first row looks like a header (has 'chrom' or 'start' as values)
    # If columns are named properly, use them; otherwise treat as headerless
    if 'chrom' in df.columns or 'start' in df.columns:
        # Has header - columns are properly named
        pass
    else:
        # No header - re-read without header
        df = pd.read_csv(bed_file, sep=sep, header=None, comment='#')

    if len(df.columns) < 4:
        raise ValueError(
            f"Region table must have at least 4 columns (chrom, start, end, region_id). "
            f"Found {len(df.columns)} columns."
        )

    regions = []
    for idx, row in df.iterrows():
        # Handle both named columns and positional access
        if 'chrom' in df.columns:
            chrom = str(row['chrom'])
            start = int(row['start'])
            end = int(row['end'])
            region_id = str(row['region_id'])
        else:
            chrom = str(row.iloc[0])
            start = int(row.iloc[1])
            end = int(row.iloc[2])
            region_id = str(row.iloc[3])

        # Optional per-region normalization (5th column if present)
        norm = normalize_to
        if 'chrom' in df.columns:
            # Named columns - check for normalize_to column
            if 'normalize_to' in df.columns and pd.notna(row.get('normalize_to')):
                norm = int(row['normalize_to'])
        else:
            # Positional - check 5th column
            if len(row) >= 5 and pd.notna(row.iloc[4]):
                norm = int(row.iloc[4])

        regions.append(Region(
            chrom=chrom,
            start=start,
            end=end,
            region_id=region_id,
            normalize_to=norm,
        ))

    return regions


def load_splice_pert(tsv_file: str) -> List[SplicePert]:
    """
    Load splice perturbation table from TSV file.

    Expected format (tab-separated with header):
        chrom  start  end  id  [chrom_src  start_src  end_src]  [seq]  [region_id]

    Must have either (chrom_src, start_src, end_src) OR seq, not both.

    Args:
        tsv_file: Path to TSV file

    Returns:
        List of SplicePert objects
    """
    df = pd.read_csv(tsv_file, sep='\t', comment='#')
    df.columns = [c.lower().strip() for c in df.columns]

    # Validate required columns
    required = {'chrom', 'start', 'end', 'id'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Splice table missing required columns: {missing}")

    splice_perts = []
    for idx, row in df.iterrows():
        # Detect source type
        has_genomic = 'chrom_src' in row.index and pd.notna(row.get('chrom_src'))
        has_seq = 'seq' in row.index and pd.notna(row.get('seq'))

        splice_perts.append(SplicePert(
            chrom=str(row['chrom']),
            start=int(row['start']),
            end=int(row['end']),
            id=str(row['id']),
            chrom_src=str(row['chrom_src']) if has_genomic else None,
            start_src=int(row['start_src']) if has_genomic else None,
            end_src=int(row['end_src']) if has_genomic else None,
            seq=str(row['seq']).upper() if has_seq else None,
            region_id=str(row['region_id']) if 'region_id' in row.index and pd.notna(row.get('region_id')) else None,
        ))

    return splice_perts


def load_seq_pert(tsv_file: str) -> List[SeqPert]:
    """
    Load sequence perturbation table from TSV file.

    Expected format (tab-separated with header):
        chrom  start  end  id  [chrom_src  start_src  end_src]  [seq]  [region_id]

    Must have either (chrom_src, start_src, end_src) OR seq, not both.
    id serves as the grouping identifier (variant ID).

    Args:
        tsv_file: Path to TSV file

    Returns:
        List of SeqPert objects
    """
    df = pd.read_csv(tsv_file, sep='\t', comment='#')
    df.columns = [c.lower().strip() for c in df.columns]

    # Validate required columns
    required = {'chrom', 'start', 'end', 'id'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Seq pert table missing required columns: {missing}")

    seq_perts = []
    for idx, row in df.iterrows():
        # Detect source type
        has_genomic = 'chrom_src' in row.index and pd.notna(row.get('chrom_src'))
        has_seq = 'seq' in row.index and pd.notna(row.get('seq'))

        seq_perts.append(SeqPert(
            chrom=str(row['chrom']),
            start=int(row['start']),
            end=int(row['end']),
            id=str(row['id']),
            chrom_src=str(row['chrom_src']) if has_genomic else None,
            start_src=int(row['start_src']) if has_genomic else None,
            end_src=int(row['end_src']) if has_genomic else None,
            seq=str(row['seq']).upper() if has_seq else None,
            region_id=str(row['region_id']) if 'region_id' in row.index and pd.notna(row.get('region_id')) else None,
        ))

    return seq_perts


# =============================================================================
# Genome Cache (Thread-Safe FASTA Access)
# =============================================================================

class GenomeCache:
    """
    Thread-safe shared genome file handle with pooled access.

    Uses multiple file handles for parallel access from worker threads.
    Uses thread-local storage to reduce lock contention.

    Supports context manager protocol for reliable cleanup:
        with GenomeCache(genome_fasta) as cache:
            seq = cache.fetch('chr1', 0, 1000)
    """

    def __init__(
        self,
        genome_fasta: str,
        num_handles: int = 8,
    ):
        """
        Initialize genome cache.

        Args:
            genome_fasta: Path to indexed FASTA file (.fa or .fasta with .fai index)
            num_handles: Number of parallel file handles for thread pool
        """
        self.genome_fasta = genome_fasta
        self.num_handles = num_handles
        self._handles: List[pysam.FastaFile] = []
        self._handle_locks: List[Lock] = []
        self._handle_idx = 0
        self._idx_lock = Lock()
        self._initialized = False
        self._closed = False

        # Thread-local storage for better performance
        import threading
        self._local = threading.local()

    def __enter__(self) -> "GenomeCache":
        """Context manager entry - ensures handles are initialized."""
        self._ensure_initialized()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        """Context manager exit - ensures handles are closed."""
        self.close()
        return False  # Don't suppress exceptions

    def _ensure_initialized(self):
        """Lazy initialization of file handles."""
        if self._closed:
            raise RuntimeError("GenomeCache has been closed and cannot be reused")
        if not self._initialized:
            with self._idx_lock:
                if not self._initialized:
                    for _ in range(self.num_handles):
                        self._handles.append(pysam.FastaFile(self.genome_fasta))
                        self._handle_locks.append(Lock())
                    self._initialized = True

    def _get_handle(self) -> Tuple[pysam.FastaFile, Lock]:
        """Get next available file handle (round-robin with thread affinity)."""
        self._ensure_initialized()

        # Try to use thread-local cached handle first
        if hasattr(self._local, 'handle_idx'):
            idx = self._local.handle_idx
            return self._handles[idx], self._handle_locks[idx]

        # Assign handle to this thread (reduces lock contention)
        with self._idx_lock:
            idx = self._handle_idx
            self._handle_idx = (self._handle_idx + 1) % self.num_handles
            self._local.handle_idx = idx

        return self._handles[idx], self._handle_locks[idx]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """
        Fetch sequence from genome (thread-safe).

        Args:
            chrom: Chromosome name
            start: Start position (0-based, inclusive)
            end: End position (0-based, exclusive)

        Returns:
            DNA sequence string (uppercase)
        """
        handle, lock = self._get_handle()
        with lock:
            return handle.fetch(chrom, max(0, start), end).upper()

    def fetch_batch(
        self,
        regions: List[Tuple[str, int, int]],
    ) -> List[str]:
        """
        Fetch multiple sequences efficiently (thread-safe).

        For regions on the same chromosome, this can coalesce fetches
        to reduce I/O overhead.

        Args:
            regions: List of (chrom, start, end) tuples

        Returns:
            List of DNA sequences (uppercase)
        """
        if not regions:
            return []

        # Group by chromosome for potential optimization
        from collections import defaultdict
        by_chrom: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
        for idx, (chrom, start, end) in enumerate(regions):
            by_chrom[chrom].append((idx, start, end))

        results = [None] * len(regions)

        for chrom, chrom_regions in by_chrom.items():
            # Sort by start position
            chrom_regions.sort(key=lambda x: x[1])

            # Fetch each region (could be optimized to coalesce adjacent regions)
            for idx, start, end in chrom_regions:
                results[idx] = self.fetch(chrom, start, end)

        return results

    def close(self):
        """Close all file handles."""
        if self._closed:
            return
        for handle in self._handles:
            try:
                handle.close()
            except Exception:
                pass  # Ignore errors during cleanup
        self._handles = []
        self._handle_locks = []
        self._initialized = False
        self._closed = True

    def __del__(self):
        """Destructor - closes handles if not already closed."""
        try:
            self.close()
        except Exception:
            pass  # Ignore errors during destruction


# =============================================================================
# Sequence Builder
# =============================================================================

class SequenceBuilder:
    """
    Builds variant sequences by applying splice and sequence perturbations.
    """

    def __init__(self, genome_cache: GenomeCache):
        """
        Initialize sequence builder.

        Args:
            genome_cache: GenomeCache instance for fetching sequences
        """
        self.genome_cache = genome_cache

    def build_variant_sequence(
        self,
        variant_spec: VariantSpec,
        seq_len: int,
        chrom_sizes: Dict[str, int],
    ) -> np.ndarray:
        """
        Build one-hot encoded sequence for a variant.

        Steps:
        1. Extract base sequence (centered on region, length=seq_len)
        2. Apply splice operations (replace regions with donor sequences)
        3. Apply seq operations (replace regions with donor sequences)
        4. One-hot encode

        Args:
            variant_spec: Variant specification
            seq_len: Model sequence length (e.g., 524288)
            chrom_sizes: Chromosome sizes for boundary handling

        Returns:
            One-hot encoded sequence (seq_len, 4)
        """
        region = variant_spec.region

        # 1. Determine extraction window
        seq_start, seq_end = region.get_extraction_window(seq_len)

        # 2. Fetch base sequence with explicit left/right boundary handling
        # Handle negative seq_start (left boundary)
        left_pad_len = 0
        fetch_start = seq_start
        if seq_start < 0:
            left_pad_len = abs(seq_start)
            fetch_start = 0

        # Handle seq_end exceeding chromosome (right boundary)
        chrom_size = chrom_sizes.get(region.chrom, float('inf'))
        fetch_end = min(seq_end, int(chrom_size))

        # Fetch sequence (pysam will handle fetch_end > chrom_size gracefully)
        seq_str = self.genome_cache.fetch(region.chrom, fetch_start, fetch_end)

        # Convert string directly to uint8 array (ASCII codes)
        # This avoids slow Python string/list operations on 500k+ elements
        seq_arr = np.frombuffer(seq_str.upper().encode('ascii'), dtype=np.uint8).copy()

        # Pad left if needed (region starts before chromosome start)
        if left_pad_len > 0:
            # Create array of 'N' characters (ASCII 78)
            left_pad = np.full(left_pad_len, ord('N'), dtype=np.uint8)
            seq_arr = np.concatenate([left_pad, seq_arr])

        # Pad right if needed (region extends past chromosome end)
        expected_len = seq_len
        if len(seq_arr) < expected_len:
            right_pad_len = expected_len - len(seq_arr)
            right_pad = np.full(right_pad_len, ord('N'), dtype=np.uint8)
            seq_arr = np.concatenate([seq_arr, right_pad])
        elif len(seq_arr) > expected_len:
            # Trim if we got too much (shouldn't happen with correct logic above)
            seq_arr = seq_arr[:expected_len]

        # 3. Apply splice operations
        for splice in variant_spec.splice_ops:
            donor_str = self._get_donor_sequence(splice)
            donor_arr = np.frombuffer(donor_str.upper().encode('ascii'), dtype=np.uint8)

            # Convert genomic coords to seq-relative coords
            rel_start = splice.start - seq_start
            rel_end = splice.end - seq_start

            # Validate bounds
            if rel_start < 0 or rel_end > len(seq_arr):
                warnings.warn(
                    f"Splice {splice.id} at {splice.chrom}:{splice.start}-{splice.end} "
                    f"is outside extraction window {region.chrom}:{seq_start}-{seq_end}. Skipping."
                )
                continue

            # NumPy slice assignment (much faster than list slicing)
            target_len = rel_end - rel_start
            seq_arr[rel_start:rel_end] = donor_arr[:target_len]

        # 4. Apply seq operations
        for seq_op in variant_spec.seq_ops:
            donor_str = self._get_donor_sequence(seq_op)
            donor_arr = np.frombuffer(donor_str.upper().encode('ascii'), dtype=np.uint8)

            rel_start = seq_op.start - seq_start
            rel_end = seq_op.end - seq_start

            # Validate bounds
            if rel_start < 0 or rel_end > len(seq_arr):
                warnings.warn(
                    f"SeqPert {seq_op.id} at {seq_op.chrom}:{seq_op.start}-{seq_op.end} "
                    f"is outside extraction window {region.chrom}:{seq_start}-{seq_end}. Skipping."
                )
                continue

            # NumPy slice assignment
            target_len = rel_end - rel_start
            seq_arr[rel_start:rel_end] = donor_arr[:target_len]

        # 5. One-hot encode directly from NumPy array
        # Use the ONE_HOT_LOOKUP table from data_loaders
        from .data_loaders import ONE_HOT_LOOKUP
        one_hot = ONE_HOT_LOOKUP[seq_arr]

        return one_hot

    def _get_donor_sequence(self, pert: Union[SplicePert, SeqPert]) -> str:
        """
        Get donor sequence from genomic coords or literal seq.

        Args:
            pert: SplicePert or SeqPert object

        Returns:
            DNA sequence string (uppercase)
        """
        if pert.seq is not None:
            return pert.seq.upper()
        else:
            return self.genome_cache.fetch(pert.chrom_src, pert.start_src, pert.end_src)


# =============================================================================
# Variant Graph Builder
# =============================================================================

def _check_overlaps(
    operations: List[Union[SplicePert, SeqPert]],
    operation_type: str,
) -> List[Tuple[Union[SplicePert, SeqPert], Union[SplicePert, SeqPert]]]:
    """
    Check for overlapping perturbation operations.

    Args:
        operations: List of SplicePert or SeqPert objects
        operation_type: Description for warning messages ('splice' or 'seq')

    Returns:
        List of overlapping operation pairs
    """
    if len(operations) < 2:
        return []

    # Sort by start position
    sorted_ops = sorted(operations, key=lambda x: (x.chrom, x.start))
    overlaps = []

    for i in range(len(sorted_ops) - 1):
        op1 = sorted_ops[i]
        op2 = sorted_ops[i + 1]

        # Check if on same chromosome and overlapping
        if op1.chrom == op2.chrom and op1.end > op2.start:
            overlaps.append((op1, op2))

    return overlaps


def build_variant_specs(
    region: Region,
    splice_perts: List[SplicePert],
    seq_perts: List[SeqPert],
    seq_len: int = 524288,
    validate_overlaps: bool = True,
) -> List[VariantSpec]:
    """
    Build all variant specifications for a region.

    Creates variants based on aggregated perturbations:
    1. Baseline (no perturbations)
    2. {id} (contains splice op + ALL seq perts with this id)

    Args:
        region: Region to build variants for
        splice_perts: All splice perturbations
        seq_perts: All sequence perturbations
        seq_len: Model sequence length for validation (default: 524288)
        validate_overlaps: If True, warn about overlapping perturbations

    Returns:
        List of VariantSpec objects
    """
    variants = []

    # 1. Baseline variant (no perturbations)
    variants.append(VariantSpec(
        variant_id="baseline",
        region=region,
        splice_ops=[],
        seq_ops=[],
    ))

    # Filter applicable operations for this region
    applicable_splices = [
        s for s in splice_perts
        if s.region_id is None or s.region_id == region.region_id
    ]

    # Map id -> List[SplicePert] (allow multiple splices per id)
    from collections import defaultdict
    splice_map = defaultdict(list)
    for s in applicable_splices:
        splice_map[s.id].append(s)

    applicable_seq_perts = [
        p for p in seq_perts
        if p.region_id is None or p.region_id == region.region_id
    ]

    # Map id -> List[SeqPert]
    seq_map = defaultdict(list)
    for p in applicable_seq_perts:
        seq_map[p.id].append(p)

    # Identify all unique variant IDs
    all_ids = set(splice_map.keys()) | set(seq_map.keys())
    all_ids.discard("baseline")

    # Get extraction window for validation using provided seq_len
    seq_start, seq_end = region.get_extraction_window(seq_len=seq_len)

    for var_id in sorted(list(all_ids)):
        current_splice_ops = splice_map.get(var_id, [])
        current_seq_ops = seq_map.get(var_id, [])

        # Validate splices are within extraction window
        for splice_op in current_splice_ops:
            if not (seq_start <= splice_op.start < splice_op.end <= seq_end):
                warnings.warn(
                    f"Splice {splice_op.id} at {splice_op.chrom}:{splice_op.start}-{splice_op.end} "
                    f"is outside region {region.region_id} extraction window {region.chrom}:{seq_start}-{seq_end}. "
                    f"This may cause issues."
                )

        # Validate overlaps within this variant
        if validate_overlaps:
            # Check for overlapping splice operations
            if len(current_splice_ops) > 1:
                splice_overlaps = _check_overlaps(current_splice_ops, 'splice')
                for op1, op2 in splice_overlaps:
                    warnings.warn(
                        f"Overlapping splice perturbations in variant {var_id}: "
                        f"{op1.chrom}:{op1.start}-{op1.end} overlaps with "
                        f"{op2.chrom}:{op2.start}-{op2.end}. "
                        f"Later operations will overwrite earlier ones."
                    )

            # Check for overlapping seq operations
            if current_seq_ops:
                seq_overlaps = _check_overlaps(current_seq_ops, 'seq')
                for op1, op2 in seq_overlaps:
                    warnings.warn(
                        f"Overlapping seq perturbations in variant {var_id}: "
                        f"{op1.chrom}:{op1.start}-{op1.end} overlaps with "
                        f"{op2.chrom}:{op2.start}-{op2.end}. "
                        f"Later operations will overwrite earlier ones."
                    )

                # Check if seq ops overlap with splice ops
                for splice_op in current_splice_ops:
                    for seq_op in current_seq_ops:
                        if (seq_op.chrom == splice_op.chrom and
                            seq_op.start < splice_op.end and seq_op.end > splice_op.start):
                            # This is expected - seq ops modify within splice region
                            pass

        variants.append(VariantSpec(
            variant_id=var_id,
            region=region,
            splice_ops=current_splice_ops,
            seq_ops=current_seq_ops,
        ))

    return variants


# =============================================================================
# Main Analyzer
# =============================================================================

class SimplePerturbationAnalyzer:
    """
    Simplified perturbation analyzer with table-based API.

    Workflow:
    1. Load regions, splice.pert, seq.pert tables
    2. For each region, build all variant specifications
    3. Build sequences for all variants
    4. Run inference in batches
    5. Format results as wide DataFrame
    """

    def __init__(
        self,
        model: nn.Module,
        genome_cache: GenomeCache,
        chrom_sizes: Dict[str, int],
        device: torch.device,
        seq_len: int = 524288,
        pred_len: int = 196608,
        model_bin_size: int = 32,
        is_human: bool = True,
        use_rc_average: bool = True,
        mixed_precision: bool = True,
        use_compile: bool = False,
        compile_mode: str = "default",
        seed: int = 42,
    ):
        """
        Initialize analyzer.

        Args:
            model: Borzoi model (nn.Module)
            genome_cache: GenomeCache for sequence fetching
            chrom_sizes: Chromosome sizes dict
            device: torch device
            seq_len: Model sequence length (default: 524288)
            pred_len: Prediction window length (default: 196608)
            model_bin_size: Bin size in bp (default: 32)
            is_human: Whether genome is human (for RC averaging)
            use_rc_average: Average predictions from forward and reverse complement
            mixed_precision: Use mixed precision inference
            use_compile: Enable torch.compile()
            compile_mode: Compile mode if use_compile=True
            seed: Random seed
        """
        self.model = model
        self.genome_cache = genome_cache
        self.chrom_sizes = chrom_sizes
        self.device = device
        self.seq_len = seq_len
        self.pred_len = pred_len
        self.model_bin_size = model_bin_size
        self.is_human = is_human
        self.seed = seed

        # Create inference engine
        self.inference_engine = InferenceEngine(
            model=model,
            device=device,
            is_human=is_human,
            use_rc_average=use_rc_average,
            mixed_precision=mixed_precision,
            use_compile=use_compile,
            compile_mode=compile_mode,
        )

        # Create sequence builder
        self.sequence_builder = SequenceBuilder(genome_cache)

    def compute(
        self,
        regions: List[Region],
        splice_perts: List[SplicePert],
        seq_perts: List[SeqPert],
        track_indices: List[int],
        track_names: List[str],
        data_matrix_path: Optional[str] = None,
        inference_batch_size: int = 64,
        num_workers: int = 16,
        show_progress: bool = True,
        profile: bool = False,
    ) -> pd.DataFrame:
        """
        Compute perturbation effects for all regions and variants.

        Args:
            regions: List of regions to analyze
            splice_perts: List of splice perturbations
            seq_perts: List of sequence perturbations
            track_indices: Track indices to extract from model output
            track_names: Track names (for column naming)
            data_matrix_path: Optional path to observed data (parquet)
            inference_batch_size: Batch size for inference
            num_workers: Number of workers for sequence loading
            show_progress: Show progress bar
            profile: Measure and print timing statistics

        Returns:
            DataFrame with columns:
            - region_id, readout_position
            - baseline[_{track}], {splice}.baseline[_{track}], {splice}.{pert}[_{track}]
            - observed_{track} (if data_matrix_path provided)
        """
        # 1. Build variant specifications
        t0 = time.time()
        all_variants = []
        for region in regions:
            variants = build_variant_specs(
                region, splice_perts, seq_perts,
                seq_len=self.seq_len,
                validate_overlaps=True,
            )
            all_variants.extend(variants)
        t_vars = time.time()

        if show_progress:
            print(f"Built {len(all_variants)} variants across {len(regions)} regions")

        # 2. Build sequences for all variants (parallel)
        sequences = self._build_sequences_parallel(
            all_variants,
            num_workers=num_workers,
            show_progress=show_progress,
        )
        t_seqs = time.time()

        # 3. Run inference in batches
        results = self._run_inference(
            all_variants,
            sequences,
            track_indices,
            inference_batch_size,
            show_progress,
        )
        t_infer = time.time()

        # 4. Format results as DataFrame
        df = self._format_results(results, track_names)

        # 5. Add observed data if provided
        if data_matrix_path:
            df = self._add_observed_data(df, data_matrix_path, track_names)
        
        if profile:
            total_seqs = len(all_variants)
            build_time = t_seqs - t_vars
            infer_time = t_infer - t_seqs
            total_compute_time = t_infer - t0
            
            print("\n--- Profiling Report ---")
            print(f"Total variants processed: {total_seqs}")
            print(f"Variant specs build time: {t_vars - t0:.2f}s")
            print(f"Sequence build time:      {build_time:.2f}s ({build_time/total_seqs*1000:.2f} ms/seq)")
            print(f"Inference time:           {infer_time:.2f}s ({infer_time/total_seqs*1000:.2f} ms/seq)")
            print(f"Total compute time:       {total_compute_time:.2f}s ({total_compute_time/total_seqs*1000:.2f} ms/seq)")
            print("------------------------\n")

        return df

    def _build_sequences_parallel(
        self,
        variants: List[VariantSpec],
        num_workers: int,
        show_progress: bool,
    ) -> List[np.ndarray]:
        """Build sequences for all variants in parallel."""
        sequences = [None] * len(variants)

        def build_one(idx):
            variant = variants[idx]
            seq = self.sequence_builder.build_variant_sequence(
                variant, self.seq_len, self.chrom_sizes
            )
            return idx, seq

        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(build_one, i) for i in range(len(variants))]

            iterator = as_completed(futures)
            if show_progress:
                iterator = tqdm(iterator, total=len(futures), desc="Building sequences")

            for future in iterator:
                idx, seq = future.result()
                sequences[idx] = seq

        return sequences

    def _run_inference(
        self,
        variants: List[VariantSpec],
        sequences: List[np.ndarray],
        track_indices: List[int],
        batch_size: int,
        show_progress: bool,
    ) -> List[Dict]:
        """Run inference on all sequences in batches."""
        results = []

        # Prepare batches
        num_batches = (len(sequences) + batch_size - 1) // batch_size

        iterator = range(num_batches)
        if show_progress:
            iterator = tqdm(iterator, desc="Running inference")

        for batch_idx in iterator:
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(sequences))

            # Get batch
            batch_variants = variants[start_idx:end_idx]
            batch_seqs = sequences[start_idx:end_idx]

            # Stack into batch array (NumPy)
            # Let inference_engine handle pinning and GPU transfer for efficiency
            batch_array = np.stack(batch_seqs, axis=0)

            # Run inference (inference_engine will pin memory and transfer to GPU)
            with torch.no_grad():
                preds = self.inference_engine.predict_batch(batch_array)  # (B, bins, tracks)

            # Extract tracks - preds is already numpy array from inference_engine
            preds = preds[:, :, track_indices]  # (B, bins, len(track_indices))

            # Store results
            for i, variant in enumerate(batch_variants):
                region = variant.region

                # Calculate readout positions (start coordinate of each 32bp bin)
                seq_start, seq_end = region.get_extraction_window(self.seq_len)
                pred_start_offset = (self.seq_len - self.pred_len) // 2
                pred_start_genomic = seq_start + pred_start_offset

                num_bins = preds.shape[1]
                readout_positions = pred_start_genomic + np.arange(num_bins) * self.model_bin_size

                results.append({
                    'region_id': region.region_id,
                    'chrom': region.chrom,
                    'variant_id': variant.variant_id,
                    'readout_positions': readout_positions,
                    'predictions': preds[i],  # (bins, tracks)
                })

        return results

    def _format_results(
        self,
        results: List[Dict],
        track_names: List[str],
    ) -> pd.DataFrame:
        """
        Format results as wide DataFrame using vectorized operations.

        Schema:
        - region_id, readout_position
        - baseline (single track) or baseline_{track} (multi-track)
        - {splice}.baseline, {splice}.{pert}, etc.

        Optimized using pandas concat + pivot for 2-5x speedup over nested loops.
        """
        if not results:
            return pd.DataFrame()

        # Vectorized approach: build list of DataFrames and concat
        all_dfs = []
        single_track = len(track_names) == 1

        for result in results:
            region_id = result['region_id']
            chrom = result['chrom']
            variant_id = result['variant_id']
            positions = result['readout_positions']
            predictions = result['predictions']  # (bins, tracks)

            num_bins = len(positions)

            # Create base DataFrame for this result
            df_chunk = pd.DataFrame({
                'region_id': [region_id] * num_bins,
                'chrom': [chrom] * num_bins,
                'readout_position': positions.astype(np.int64),
                'variant_id': variant_id,
            })

            # Add prediction columns for each track
            for t, track_name in enumerate(track_names):
                if single_track:
                    df_chunk['value'] = predictions[:, t]
                else:
                    df_chunk[f'value_{track_name}'] = predictions[:, t]

            all_dfs.append(df_chunk)

        # Concatenate all chunks
        long_df = pd.concat(all_dfs, ignore_index=True)

        if single_track:
            # Pivot: rows are (region_id, chrom, readout_position), columns are variant_id
            wide_df = long_df.pivot_table(
                index=['region_id', 'chrom', 'readout_position'],
                columns='variant_id',
                values='value',
                aggfunc='first',
            ).reset_index()

            # Flatten column names
            wide_df.columns.name = None
        else:
            # For multi-track: need to pivot each track's values
            # First, get unique (region_id, chrom, readout_position) tuples
            index_df = long_df[['region_id', 'chrom', 'readout_position']].drop_duplicates()

            # Get unique variant IDs
            variant_ids = long_df['variant_id'].unique()

            # Build wide format by merging each variant's data
            wide_df = index_df.copy()

            for variant_id in variant_ids:
                variant_data = long_df[long_df['variant_id'] == variant_id].copy()
                value_cols = [c for c in variant_data.columns if c.startswith('value_')]

                # Rename value columns to include variant_id
                rename_map = {}
                for col in value_cols:
                    track_name = col.replace('value_', '')
                    rename_map[col] = f"{variant_id}_{track_name}"

                variant_data = variant_data[['region_id', 'chrom', 'readout_position'] + value_cols]
                variant_data = variant_data.rename(columns=rename_map)

                wide_df = wide_df.merge(
                    variant_data,
                    on=['region_id', 'chrom', 'readout_position'],
                    how='left',
                )

        # Sort by region_id and position
        wide_df = wide_df.sort_values(['region_id', 'readout_position']).reset_index(drop=True)

        return wide_df

    def _add_observed_data(
        self,
        df: pd.DataFrame,
        data_matrix_path: str,
        track_names: List[str],
    ) -> pd.DataFrame:
        """
        Add observed data columns to DataFrame.

        Adds columns: observed_{track_name}

        The observed data is matched by aligning readout_position to the observed
        data's 32bp bin grid. Since readout_position is the bin start coordinate,
        we snap it to the nearest observed bin.
        """
        # Load observed data
        obs_df = pd.read_parquet(data_matrix_path)

        # Expected columns: region_id, readout_position, {track_name}_observed or observed_{track_name}
        # Try both naming conventions
        for track_name in track_names:
            obs_col_name = f"observed_{track_name}"

            # Try finding the column (might be named differently)
            if obs_col_name in obs_df.columns:
                src_col = obs_col_name
            elif f"{track_name}_observed" in obs_df.columns:
                src_col = f"{track_name}_observed"
            elif track_name in obs_df.columns:
                src_col = track_name
            else:
                warnings.warn(f"Could not find observed data for track {track_name}")
                continue

            # Merge strategy depends on observed data format
            if 'region_id' in obs_df.columns and 'readout_position' in obs_df.columns:
                # Direct match on region_id + readout_position
                merge_cols = ['region_id', 'readout_position']
                df = df.merge(
                    obs_df[merge_cols + [src_col]].rename(columns={src_col: obs_col_name}),
                    on=merge_cols,
                    how='left',
                )
            elif 'chrom' in obs_df.columns and 'start' in obs_df.columns:
                # Observed data has chrom/start/end format (genome-wide 32bp bins)
                # Need to extract chrom from results and match by chrom + start

                # Get chromosome info - try to extract from region or use provided
                if 'chrom' not in df.columns:
                    # Try to get chrom from region_id or other source
                    # For now, we need to iterate per-track
                    warnings.warn(
                        f"Cannot merge observed data: results missing 'chrom' column. "
                        f"Observed data has chrom/start format but results only have readout_position."
                    )
                    continue

                # The observed 'start' is on a 32bp grid from chromosome start
                # Our readout_position is the bin start but may not align to the same grid
                # Snap readout_position to nearest 32bp boundary
                bin_size = self.model_bin_size
                df['_obs_start'] = (df['readout_position'] // bin_size) * bin_size

                # Merge on chrom + snapped start
                obs_subset = obs_df[['chrom', 'start', src_col]].copy()
                obs_subset = obs_subset.rename(columns={'start': '_obs_start', src_col: obs_col_name})
                obs_subset['_obs_start'] = obs_subset['_obs_start'].astype(np.int64)

                df = df.merge(
                    obs_subset,
                    on=['chrom', '_obs_start'],
                    how='left',
                )

                # Clean up temp column
                df = df.drop(columns=['_obs_start'])
            else:
                warnings.warn(f"Could not merge observed data: missing region_id/readout_position or chrom/start columns")
                continue

        return df
