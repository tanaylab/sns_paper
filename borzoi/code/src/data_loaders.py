"""PyTorch data loading utilities for Borzoi fine-tuning."""

import json
import os
import shutil
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
import pysam
from tqdm import tqdm
import torch
from torch.utils.data import Dataset, DataLoader, ConcatDataset, get_worker_info

# Numerical stability constant for power operations on near-zero values
SQUASH_EPSILON = 1e-6

# Threshold for squash transformation where sqrt scaling kicks in.
# Borzoi's squash applies y^0.75, then sqrt(excess) for values above threshold.
SQUASH_THRESHOLD = 384.0

ONE_HOT_LOOKUP = np.zeros((256, 4), dtype=np.float32)
for base, vec in (('A', [1, 0, 0, 0]),
                  ('C', [0, 1, 0, 0]),
                  ('G', [0, 0, 1, 0]),
                  ('T', [0, 0, 0, 1]),
                  ('N', [0, 0, 0, 0])):
    ONE_HOT_LOOKUP[ord(base)] = vec
    ONE_HOT_LOOKUP[ord(base.lower())] = vec


def one_hot_encode(sequence: str) -> np.ndarray:
    """Convert a DNA sequence to one-hot encoding.

    Args:
        sequence: DNA string (A/C/G/T/N).

    Returns:
        Array of shape (len(sequence), 4).
    """
    seq_array = np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)
    return ONE_HOT_LOOKUP[seq_array]


def reverse_complement_onehot(onehot: np.ndarray) -> np.ndarray:
    """Reverse complement a one-hot encoded sequence."""
    return onehot[::-1, [3, 2, 1, 0]].copy()


def squash_coverage(y: np.ndarray, threshold: float = SQUASH_THRESHOLD) -> np.ndarray:
    """Apply Borzoi's squash transformation to coverage values."""
    y_exp = np.power(y.astype(np.float32) + SQUASH_EPSILON, 0.75)
    mask = y_exp > threshold
    result = y_exp.copy()
    result[mask] = threshold + np.sqrt(y_exp[mask] - threshold)
    return result


def calculate_padding(seq_len: int, pred_len: int) -> int:
    """Calculate the padding/context needed on each side."""
    assert (seq_len - pred_len) % 2 == 0
    return (seq_len - pred_len) // 2


def calculate_target_bins(pred_len: int, bin_size: int = 32) -> int:
    """Calculate the number of target bins."""
    assert pred_len % bin_size == 0
    return pred_len // bin_size


def load_chrom_sizes(chrom_sizes_file: str) -> Dict[str, int]:
    """Load chromosome sizes from a file."""
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def load_regions(bed_file: str) -> pd.DataFrame:
    """Load genomic regions from a BED file."""
    df = pd.read_csv(bed_file, sep='\t', header=None,
                     usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    return df


def generate_genome_wide_regions(
    chrom_sizes: Dict[str, int],
    tile_size: int = 1000000,
    exclude_chroms: Optional[List[str]] = None,
    pred_len: Optional[int] = None
) -> pd.DataFrame:
    """Create tiled regions across the genome.

    Args:
        chrom_sizes: Map of chromosome to length.
        tile_size: Window size for tiling in bp.
        exclude_chroms: Chromosomes to skip.
        pred_len: Prediction window length. If provided, tile_size will be aligned
                  to be a multiple of pred_len to prevent data loss.

    Returns:
        DataFrame with columns chrom/start/end.
    """
    if exclude_chroms is None:
        exclude_chroms = []

    # Align tile_size to be a multiple of pred_len to prevent data loss
    # Example: tile_size=1,000,000, pred_len=24,576 -> aligned_tile_size=999,424 (40 * 24,576)
    if pred_len is not None and pred_len > 0:
        num_windows = tile_size // pred_len
        if num_windows == 0:
            num_windows = 1
        aligned_tile_size = num_windows * pred_len
        if aligned_tile_size != tile_size:
            print(f"INFO: Aligning tile_size from {tile_size} to {aligned_tile_size} "
                  f"({num_windows} x {pred_len} bp) to prevent data loss")
        tile_size = aligned_tile_size

    regions = []
    for chrom, chrom_len in chrom_sizes.items():
        # Skip excluded chromosomes
        if chrom in exclude_chroms:
            continue

        # Generate tiles across the chromosome
        for start in range(0, chrom_len, tile_size):
            end = min(start + tile_size, chrom_len)
            regions.append({'chrom': chrom, 'start': start, 'end': end})

    return pd.DataFrame(regions)


def build_parquet_cache_streaming(
    input_file: str,
    output_dir: str,
    bin_size: int = 32,
    track_cols: Optional[List[str]] = None,
    batch_size: int = 200000,
) -> Tuple[List[str], int]:
    """Build a memory-mappable NPY-dir cache from Parquet using bounded RAM.

    This avoids materializing the full Parquet table in memory.

    Args:
        input_file: Path to source Parquet file.
        output_dir: Directory to write per-chromosome .npy files + metadata.json.
        bin_size: Expected bin size; adjusted if Parquet 'end-start' indicates otherwise.
        track_cols: Optional subset of track columns to cache.
        batch_size: Number of rows per Arrow batch.

    Returns:
        Tuple of (track_names, detected_bin_size).
    """
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError("pyarrow required for Parquet. pip install pyarrow") from exc

    parquet_path = Path(input_file)
    out_path = Path(output_dir)
    os.makedirs(out_path, exist_ok=True)

    pf = pq.ParquetFile(str(parquet_path))
    schema_names = pf.schema.names

    if "chrom" not in schema_names or "start" not in schema_names:
        raise ValueError(
            f"Parquet must include 'chrom' and 'start' columns. Found columns: {schema_names}"
        )

    exclude_cols = {'chrom', 'start', 'end', 'region', 'chr', 'chromosome', 'intervalID'}
    available_tracks = [col for col in schema_names if col not in exclude_cols]
    if not available_tracks:
        raise ValueError("No track columns found in Parquet data matrix")

    if track_cols is None:
        track_names = available_tracks
        print(f"  Using all {len(track_names)} tracks: {', '.join(track_names[:5])}" +
              (f" ... (+{len(track_names)-5} more)" if len(track_names) > 5 else ""))
    else:
        missing_tracks = set(track_cols) - set(available_tracks)
        if missing_tracks:
            raise ValueError(
                f"Selected tracks not found in data: {missing_tracks}\n"
                f"Available tracks: {available_tracks}"
            )
        track_names = list(track_cols)
        print(f"  Selected {len(track_names)} tracks (from {len(available_tracks)} available):")
        print(f"    {', '.join(track_names)}")

    has_end = "end" in schema_names
    first_pass_cols = ["chrom", "start"] + (["end"] if has_end else [])

    print("  Streaming Parquet metadata pass...")
    chrom_order: List[str] = []
    chrom_seen = set()
    chrom_counts: Dict[str, int] = {}
    chrom_max_start: Dict[str, int] = {}
    chrom_min_start: Dict[str, int] = {}
    detected_bin_size: Optional[int] = None
    delta_counts: Dict[int, int] = {}
    dropped_invalid_rows = 0
    row_count = 0

    for batch in tqdm(
        pf.iter_batches(columns=first_pass_cols, batch_size=batch_size),
        desc="  [pass1] scanning",
        ncols=100,
    ):
        df = batch.to_pandas()
        if df.empty:
            continue

        valid_mask = pd.notna(df["chrom"]) & pd.notna(df["start"])
        if has_end:
            valid_mask &= pd.notna(df["end"])
        invalid_count = int((~valid_mask).sum())
        if invalid_count:
            dropped_invalid_rows += invalid_count
        df = df.loc[valid_mask]
        if df.empty:
            continue

        chroms = df["chrom"].to_numpy(copy=False)
        starts = df["start"].to_numpy(copy=False)
        row_count += len(df)

        if np.any(starts < 0):
            raise ValueError("Found negative start coordinates in data matrix")

        if has_end:
            ends = df["end"].to_numpy(copy=False)
            deltas = ends - starts
            if np.any(deltas <= 0):
                raise ValueError("Found non-positive interval lengths in data matrix")
            delta_int = np.rint(deltas).astype(np.int64)
            if not np.allclose(deltas, delta_int, atol=1e-6):
                raise ValueError("Found non-integer interval lengths in data matrix")
            unique_delta, counts = np.unique(delta_int, return_counts=True)
            for d, c in zip(unique_delta, counts):
                delta_counts[int(d)] = delta_counts.get(int(d), 0) + int(c)

        unique_chroms = pd.unique(chroms)
        for chrom in unique_chroms:
            mask = chroms == chrom
            chrom_starts = starts[mask]
            if chrom_starts.size == 0:
                continue
            chrom_str = str(chrom)
            if chrom_str not in chrom_seen:
                chrom_seen.add(chrom_str)
                chrom_order.append(chrom_str)
                chrom_counts[chrom_str] = 0
                chrom_max_start[chrom_str] = int(chrom_starts.max())
                chrom_min_start[chrom_str] = int(chrom_starts.min())
            chrom_counts[chrom_str] += int(chrom_starts.size)
            chrom_max_start[chrom_str] = max(chrom_max_start[chrom_str], int(chrom_starts.max()))
            chrom_min_start[chrom_str] = min(chrom_min_start[chrom_str], int(chrom_starts.min()))

    if row_count == 0:
        raise ValueError(f"Parquet file contains no rows: {input_file}")
    if dropped_invalid_rows:
        print(f"  Note: dropped {dropped_invalid_rows:,} rows with null chrom/start/end")

    if detected_bin_size is None:
        if delta_counts:
            # Use the most frequent interval width as canonical bin size.
            # Short terminal bins (e.g. at chromosome ends) remain valid.
            detected_bin_size = max(delta_counts.items(), key=lambda kv: (kv[1], kv[0]))[0]
        else:
            detected_bin_size = int(bin_size)

    if detected_bin_size != bin_size:
        print(f"  WARNING: Config bin_size={bin_size} but data has bin_size={detected_bin_size}")
        print(f"  Using detected bin_size={detected_bin_size} for cache build")

    if delta_counts:
        unusual = sorted(d for d in delta_counts.keys() if d != detected_bin_size)
        if unusual:
            preview = unusual[:5]
            print(f"  Note: found non-canonical interval widths {preview}"
                  f" (likely terminal bins at chromosome ends)")

    if detected_bin_size <= 0:
        raise ValueError(f"Invalid bin_size detected: {detected_bin_size}")

    print(f"  Each chromosome should have bins from start=0 to end with stride={detected_bin_size}")

    chrom_num_bins: Dict[str, int] = {}
    for chrom in chrom_order:
        if chrom_min_start[chrom] != 0:
            raise ValueError(
                f"SPARSE MATRIX DETECTED for {chrom}!\n"
                f"  First bin starts at {chrom_min_start[chrom]}, expected 0.\n"
                f"  Data matrix MUST be DENSE (full genome coverage)."
            )
        max_start = chrom_max_start[chrom]
        if max_start % detected_bin_size != 0:
            raise ValueError(
                f"Non-aligned start coordinate in {chrom}: max_start={max_start} "
                f"is not divisible by bin_size={detected_bin_size}"
            )
        chrom_num_bins[chrom] = max_start // detected_bin_size + 1

    mmap_arrays: Dict[str, np.memmap] = {}
    seen_bins: Dict[str, np.ndarray] = {}
    for chrom in chrom_order:
        chrom_file = out_path / f"{chrom}.npy"
        mmap_arrays[chrom] = np.lib.format.open_memmap(
            str(chrom_file),
            mode="w+",
            dtype=np.float32,
            shape=(chrom_num_bins[chrom], len(track_names)),
        )
        seen_bins[chrom] = np.zeros(chrom_num_bins[chrom], dtype=bool)

    second_pass_cols = ["chrom", "start"] + (["end"] if has_end else []) + track_names

    print("  Streaming Parquet write pass...")
    for batch in tqdm(
        pf.iter_batches(columns=second_pass_cols, batch_size=batch_size),
        desc="  [pass2] writing",
        ncols=100,
    ):
        df = batch.to_pandas()
        if df.empty:
            continue

        valid_mask = pd.notna(df["chrom"]) & pd.notna(df["start"])
        if has_end:
            valid_mask &= pd.notna(df["end"])
        df = df.loc[valid_mask]
        if df.empty:
            continue

        chroms = df["chrom"].to_numpy(copy=False)
        starts = df["start"].to_numpy(copy=False)
        if np.any(starts % detected_bin_size != 0):
            raise ValueError(
                f"Found start coordinates not divisible by detected bin_size={detected_bin_size}"
            )

        if has_end:
            ends = df["end"].to_numpy(copy=False)
            deltas = ends - starts
            if np.any(deltas <= 0):
                raise ValueError("Found non-positive interval lengths in data matrix")
            if np.any(deltas > detected_bin_size + 1e-6):
                bad = np.unique(deltas[deltas > detected_bin_size + 1e-6])
                raise ValueError(
                    f"Found interval widths larger than canonical bin_size={detected_bin_size}: {bad[:10]}"
                )

        track_values = df[track_names].to_numpy(dtype=np.float32, copy=False)

        unique_chroms = pd.unique(chroms)
        for chrom in unique_chroms:
            chrom_str = str(chrom)
            if chrom_str not in mmap_arrays:
                raise ValueError(f"Unexpected chromosome in pass2: {chrom_str}")
            mask = chroms == chrom
            chrom_starts = starts[mask]
            chrom_bins = (chrom_starts // detected_bin_size).astype(np.int64)
            if chrom_bins.size == 0:
                continue
            if np.unique(chrom_bins).size != chrom_bins.size:
                raise ValueError(
                    f"Duplicate bins detected within a batch for {chrom_str}"
                )
            if np.any(chrom_bins < 0) or np.any(chrom_bins >= chrom_num_bins[chrom_str]):
                raise ValueError(
                    f"Out-of-range bins detected for {chrom_str} "
                    f"(max valid index {chrom_num_bins[chrom_str] - 1})"
                )
            if np.any(seen_bins[chrom_str][chrom_bins]):
                raise ValueError(f"Duplicate bins detected across batches for {chrom_str}")
            mmap_arrays[chrom_str][chrom_bins, :] = track_values[mask]
            seen_bins[chrom_str][chrom_bins] = True

    for chrom in chrom_order:
        expected = chrom_num_bins[chrom]
        observed = chrom_counts[chrom]
        if observed != expected:
            raise ValueError(
                f"SPARSE MATRIX DETECTED for {chrom}!\n"
                f"  Expected {expected} bins from 0..{chrom_max_start[chrom]} "
                f"(stride={detected_bin_size}), found {observed} rows."
            )
        if not np.all(seen_bins[chrom]):
            missing = int(np.count_nonzero(~seen_bins[chrom]))
            raise ValueError(
                f"SPARSE MATRIX DETECTED for {chrom}!\n"
                f"  Missing {missing} bins during cache build."
            )
        mmap_arrays[chrom].flush()

    meta = {
        "format": "npy-dir-v1",
        "track_names": track_names,
        "bin_size": int(detected_bin_size),
        "chromosomes": chrom_order,
    }
    with (out_path / "metadata.json").open("w") as f:
        json.dump(meta, f)

    print(f"  [OK] Validated: All chromosomes have DENSE coverage (stride={detected_bin_size}bp)")
    print(f"Saved memory-mappable cache to: {output_dir}")
    return track_names, detected_bin_size


def load_from_data_matrix(
    input_file: str,
    bin_size: int = 32,
    track_cols: Optional[List[str]] = None,
    use_mmap: bool = False,
) -> Tuple[Dict[str, np.ndarray], List[str], int]:
    """Load genomic data from a dense matrix (NPZ, Parquet, CSV, or HDF5).

    Supports memory-mapped loading for large datasets to reduce memory footprint.

    Args:
        input_file: Path to matrix file with chrom/start/end and tracks.
        bin_size: Expected bin size; adjusted if the file indicates otherwise.
        track_cols: Optional subset of track column names.
        use_mmap: If True, use memory-mapped files for reduced memory usage.
                 Only applicable for NPZ and HDF5 formats.

    Returns:
        vector_dict: chrom -> array of shape (bins, tracks)
        track_names: Ordered track names used
        bin_size: Detected bin size
    """
    file_path = Path(input_file)
    if file_path.is_dir():
        return load_from_npy_dir(input_file, track_cols=track_cols, use_mmap=use_mmap)
    suffix = file_path.suffix.lower()

    # Handle HDF5 format with optional memory mapping
    if suffix in ['.h5', '.hdf5']:
        print(f"Loading data matrix from HDF5: {input_file}")
        return load_from_hdf5(input_file, track_cols, use_mmap=use_mmap)

    if suffix == '.npz':
        print(f"Loading data matrix from NPZ: {input_file}")
        if use_mmap:
            return load_binned_bigwigs_mmap(input_file)
        return load_binned_bigwigs(input_file)

    elif suffix == '.parquet':
        print(f"Loading data matrix from Parquet: {input_file}")
        try:
            import pyarrow.parquet as pq
        except ImportError:
            raise ImportError("pyarrow required for Parquet. pip install pyarrow")

        table = pq.read_table(input_file)
        df = table.to_pandas()

    elif suffix in ['.csv', '.gz']:
        print(f"Loading data matrix from CSV: {input_file}")
        if input_file.endswith('.gz'):
            import gzip
            df = pd.read_csv(gzip.open(input_file, 'rt'))
        else:
            df = pd.read_csv(input_file)
    else:
        raise ValueError(f"Unsupported format: {suffix}")

    exclude_cols = {'chrom', 'start', 'end', 'region', 'chr', 'chromosome', 'intervalID'}
    available_tracks = [col for col in df.columns if col not in exclude_cols]

    if track_cols is None:
        track_cols = available_tracks
        print(f"  Using all {len(track_cols)} tracks: {', '.join(track_cols[:5])}" +
              (f" ... (+{len(track_cols)-5} more)" if len(track_cols) > 5 else ""))
    else:
        missing_tracks = set(track_cols) - set(available_tracks)
        if missing_tracks:
            raise ValueError(
                f"Selected tracks not found in data: {missing_tracks}\n"
                f"Available tracks: {available_tracks}"
            )
        print(f"  Selected {len(track_cols)} tracks (from {len(available_tracks)} available):")
        print(f"    {', '.join(track_cols)}")

    vector_dict = {}
    chromosomes = df['chrom'].unique()

    print(f"  Converting {len(chromosomes)} chromosomes to arrays...")

    first_chrom_data = df[df['chrom'] == chromosomes[0]].head(10).sort_values('start')
    if 'end' in first_chrom_data.columns:
        detected_bin_size = int(first_chrom_data.iloc[0]['end'] - first_chrom_data.iloc[0]['start'])
        print(f"  Detected bin_size from data: {detected_bin_size}bp")
        if detected_bin_size != bin_size:
            print(f"  WARNING: Config bin_size={bin_size} but data has bin_size={detected_bin_size}")
            print(f"  Using detected bin_size={detected_bin_size} for validation")
            bin_size = detected_bin_size

    print(f"  Each chromosome should have bins from start=0 to end with stride={bin_size}")

    for chrom, chr_data in tqdm(df.groupby('chrom', sort=False, observed=False), desc="  Processing chromosomes", ncols=100):
        chr_data = chr_data.sort_values('start')

        starts = chr_data['start'].values
        if len(starts) > 1:
            if starts[0] != 0:
                raise ValueError(
                    f"SPARSE MATRIX DETECTED for {chrom}!\n"
                    f"  First bin starts at {starts[0]}, expected 0.\n"
                    f"  Data matrix MUST be DENSE (full genome coverage).\n"
                    f"  Use BigWig files (lazy_loading=true) for sparse data instead."
                )

            if len(starts) >= 2:
                actual_stride = int(starts[1] - starts[0])
                if actual_stride != bin_size:
                    print(f"  Adjusting bin_size for {chrom}: {bin_size} -> {actual_stride}")
                    bin_size = actual_stride

            stride_test = starts[1:] - starts[:-1]

            unique_strides = np.unique(stride_test)

            if len(unique_strides) > 1 or (len(unique_strides) == 1 and unique_strides[0] != bin_size):
                # Found non-uniform strides = gaps in coverage
                gaps = np.where(stride_test != bin_size)[0]

                if len(gaps) > 0:
                    first_gap_idx = gaps[0]
                    gap_start = starts[first_gap_idx]
                    gap_size = stride_test[first_gap_idx]
                    raise ValueError(
                        f"SPARSE MATRIX DETECTED for {chrom}!\n"
                        f"  Found gap in coverage at position {gap_start}.\n"
                        f"  Expected stride {bin_size}bp, found {gap_size}bp.\n"
                        f"  Data matrix MUST be DENSE (continuous genome-wide bins).\n"
                        f"  Export full genome with: gtrack.create('track', ..., iterator={bin_size})\n"
                        f"  Or use BigWig files (lazy_loading=true) for sparse data."
                    )

        track_matrix = chr_data[track_cols].values.astype(np.float32)
        vector_dict[chrom] = track_matrix

    print(f"  [OK] Validated: All chromosomes have DENSE coverage (stride={bin_size}bp)")
    print(f"  Loaded {len(vector_dict)} chromosomes, {len(track_cols)} tracks")
    return vector_dict, track_cols, bin_size


def load_from_npy_dir(
    input_dir: str,
    track_cols: Optional[List[str]] = None,
    use_mmap: bool = False,
) -> Tuple[Dict[str, np.ndarray], List[str], int]:
    """Load genomic data from a directory of per-chromosome .npy files.

    Expected layout:
        input_dir/
            metadata.json
            chr1.npy
            chr2.npy
            ...

    Args:
        input_dir: Directory containing per-chromosome .npy files.
        track_cols: Optional subset of track column names.
        use_mmap: If True, memory-map .npy files to reduce RAM usage.

    Returns:
        vector_dict: chrom -> array (bins, tracks)
        track_names: Ordered track names used
        bin_size: Bin size used when exporting
    """
    from pathlib import Path

    input_path = Path(input_dir)
    meta_path = input_path / "metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata.json in {input_dir}")

    with meta_path.open("r") as f:
        meta = json.load(f)

    track_names = meta.get("track_names", [])
    bin_size = int(meta.get("bin_size", 32))
    chromosomes = meta.get("chromosomes")

    if not chromosomes:
        chromosomes = sorted([p.stem for p in input_path.glob("*.npy")])

    if not chromosomes:
        raise ValueError(f"No .npy chromosome files found in {input_dir}")

    if track_cols is None:
        print(f"  Using all {len(track_names)} tracks from cache")
        track_indices = None
    else:
        missing_tracks = set(track_cols) - set(track_names)
        if missing_tracks:
            raise ValueError(
                f"Selected tracks not found in cache: {missing_tracks}\n"
                f"Available tracks: {track_names}"
            )
        track_indices = [track_names.index(t) for t in track_cols]
        track_names = track_cols
        print(f"  Selected {len(track_cols)} tracks from cache")

    vector_dict = {}
    mmap_mode = 'r' if use_mmap else None
    for chrom in chromosomes:
        chrom_path = input_path / f"{chrom}.npy"
        if not chrom_path.exists():
            raise FileNotFoundError(f"Missing chromosome file: {chrom_path}")
        data = np.load(chrom_path, mmap_mode=mmap_mode)
        if track_indices is not None:
            data = data[:, track_indices]
        vector_dict[chrom] = data

    if use_mmap:
        print(f"  Memory-mapped {len(vector_dict)} chromosomes from cache")
    else:
        print(f"  Loaded {len(vector_dict)} chromosomes from cache")
    return vector_dict, track_names, bin_size


def load_binned_bigwigs(input_file: str) -> Tuple[Dict[str, np.ndarray], List[str], int]:
    """Load pre-computed binned BigWig arrays from file.

    Args:
        input_file: Path to NPZ produced by save_binned_bigwigs.

    Returns:
        vector_dict: chrom -> array (bins, tracks)
        track_names: Track names stored in metadata
        bin_size: Bin size used when exporting
    """
    print(f"Loading pre-computed arrays from: {input_file}")
    data = np.load(input_file, allow_pickle=True)

    metadata = data['_metadata'].item()
    track_names = metadata['track_names']
    bin_size = metadata['bin_size']

    vector_dict = {}
    for key in data.keys():
        if key != '_metadata':
            vector_dict[key] = data[key]

    print(f"  Loaded {len(vector_dict)} chromosomes, {len(track_names)} tracks")
    return vector_dict, track_names, bin_size


def load_binned_bigwigs_mmap(input_file: str) -> Tuple[Dict[str, np.ndarray], List[str], int]:
    """Load pre-computed binned BigWig arrays using memory mapping.

    Memory-mapped arrays reduce RAM usage by loading data on-demand from disk.
    This is particularly useful for large genomic datasets that may not fit in RAM.

    Args:
        input_file: Path to NPZ produced by save_binned_bigwigs.

    Returns:
        vector_dict: chrom -> memory-mapped array (bins, tracks)
        track_names: Track names stored in metadata
        bin_size: Bin size used when exporting
    """
    print(f"Loading pre-computed arrays with memory mapping from: {input_file}")

    # Open with mmap_mode='r' for read-only memory mapping
    data = np.load(input_file, allow_pickle=True, mmap_mode='r')

    metadata = data['_metadata'].item()
    track_names = metadata['track_names']
    bin_size = metadata['bin_size']

    vector_dict = {}
    for key in data.keys():
        if key != '_metadata':
            # Arrays are memory-mapped, not loaded into RAM
            vector_dict[key] = data[key]

    print(f"  Memory-mapped {len(vector_dict)} chromosomes, {len(track_names)} tracks")
    return vector_dict, track_names, bin_size


def load_from_hdf5(
    input_file: str,
    track_cols: Optional[List[str]] = None,
    use_mmap: bool = False,
) -> Tuple[Dict[str, np.ndarray], List[str], int]:
    """Load genomic data from HDF5 format with optional memory mapping.

    HDF5 format is ideal for large datasets as it supports:
    - Chunked storage for efficient random access
    - Compression (lzf, gzip)
    - Memory-mapped reading for reduced RAM usage

    Expected HDF5 structure:
        /metadata/track_names  - 1D string dataset
        /metadata/bin_size     - scalar dataset
        /chromosomes/<chrom>   - 2D dataset (bins, tracks)

    Args:
        input_file: Path to HDF5 file.
        track_cols: Optional subset of track column names.
        use_mmap: If True, use memory-mapped reading (data loaded on-demand).

    Returns:
        vector_dict: chrom -> array of shape (bins, tracks)
        track_names: Ordered track names used
        bin_size: Detected bin size
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py required for HDF5 support. Install with: pip install h5py")

    # Open HDF5 file
    # Note: For true memory mapping, we keep the file handle open
    # and return dataset references instead of arrays
    if use_mmap:
        print(f"  Opening HDF5 with memory-mapped access...")
        h5_file = h5py.File(input_file, 'r', driver='sec2')

        # Read metadata
        track_names = list(h5_file['metadata/track_names'][()])
        if isinstance(track_names[0], bytes):
            track_names = [t.decode('utf-8') for t in track_names]
        bin_size = int(h5_file['metadata/bin_size'][()])

        # Handle track selection
        if track_cols is not None:
            track_indices = [track_names.index(t) for t in track_cols if t in track_names]
            track_names = track_cols
        else:
            track_indices = None

        # Return dataset references (lazy loading)
        vector_dict = {}
        for chrom in h5_file['chromosomes'].keys():
            dataset = h5_file['chromosomes'][chrom]
            if track_indices is not None:
                # For selected tracks, we need to load and slice
                # (HDF5 doesn't support arbitrary column indexing efficiently)
                vector_dict[chrom] = dataset[:][:, track_indices]
            else:
                # Return as numpy array (memory mapped internally by h5py)
                vector_dict[chrom] = dataset[:]

        print(f"  Loaded {len(vector_dict)} chromosomes, {len(track_names)} tracks (memory-mapped)")
        return vector_dict, track_names, bin_size

    else:
        # Standard loading - read everything into memory
        with h5py.File(input_file, 'r') as h5_file:
            # Read metadata
            track_names = list(h5_file['metadata/track_names'][()])
            if isinstance(track_names[0], bytes):
                track_names = [t.decode('utf-8') for t in track_names]
            bin_size = int(h5_file['metadata/bin_size'][()])

            # Handle track selection
            if track_cols is not None:
                track_indices = [track_names.index(t) for t in track_cols if t in track_names]
                track_names = track_cols
            else:
                track_indices = None

            # Load chromosome data
            vector_dict = {}
            for chrom in h5_file['chromosomes'].keys():
                data = h5_file['chromosomes'][chrom][:]
                if track_indices is not None:
                    data = data[:, track_indices]
                vector_dict[chrom] = data.astype(np.float32)

        print(f"  Loaded {len(vector_dict)} chromosomes, {len(track_names)} tracks")
        return vector_dict, track_names, bin_size


def save_to_hdf5(
    vector_dict: Dict[str, np.ndarray],
    output_file: str,
    track_names: List[str],
    bin_size: int = 32,
    compression: str = 'lzf',
    chunks: bool = True,
) -> None:
    """Save genomic data to HDF5 format with compression.

    Args:
        vector_dict: chrom -> array (bins, tracks).
        output_file: HDF5 file path.
        track_names: Names to store in metadata.
        bin_size: Bin size used to build the arrays.
        compression: Compression algorithm ('lzf', 'gzip', or None).
        chunks: If True, enable chunked storage for efficient access.
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py required for HDF5 support. Install with: pip install h5py")

    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)

    with h5py.File(output_file, 'w') as h5_file:
        # Create metadata group
        meta_grp = h5_file.create_group('metadata')
        # Store track names as variable-length strings
        dt = h5py.special_dtype(vlen=str)
        meta_grp.create_dataset('track_names', data=track_names, dtype=dt)
        meta_grp.create_dataset('bin_size', data=bin_size)

        # Create chromosomes group
        chrom_grp = h5_file.create_group('chromosomes')

        for chrom, data in tqdm(vector_dict.items(), desc="Saving chromosomes"):
            # Determine optimal chunk shape
            if chunks:
                # Chunk by ~1000 bins and all tracks
                chunk_shape = (min(1000, data.shape[0]), data.shape[1])
            else:
                chunk_shape = None

            chrom_grp.create_dataset(
                chrom,
                data=data.astype(np.float32),
                compression=compression,
                chunks=chunk_shape,
            )

    print(f"Saved HDF5 to: {output_file}")


def save_binned_bigwigs(
    vector_dict: Dict[str, np.ndarray],
    output_file: str,
    track_names: List[str],
    bin_size: int = 32
) -> None:
    """Save pre-computed binned BigWig arrays for fast loading.

    Args:
        vector_dict: chrom -> array (bins, tracks).
        output_file: NPZ file path.
        track_names: Names to store in metadata.
        bin_size: Bin size used to build the arrays.
    """
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)

    arrays_to_save = dict(vector_dict)
    arrays_to_save['_metadata'] = np.array({
        'track_names': track_names,
        'bin_size': bin_size,
        'num_tracks': vector_dict[next(iter(vector_dict))].shape[1]
    }, dtype=object)

    np.savez_compressed(output_file, **arrays_to_save)
    print(f"Saved pre-computed arrays to: {output_file}")


def save_to_npy_dir(
    vector_dict: Dict[str, np.ndarray],
    output_dir: str,
    track_names: List[str],
    bin_size: int = 32
) -> None:
    """Save genomic data to a directory of per-chromosome .npy files.

    Args:
        vector_dict: chrom -> array (bins, tracks).
        output_dir: Output directory.
        track_names: Names to store in metadata.
        bin_size: Bin size used when exporting.
    """
    output_path = Path(output_dir)
    os.makedirs(output_path, exist_ok=True)

    chromosomes = list(vector_dict.keys())
    for chrom in tqdm(chromosomes, desc="Saving chromosomes", ncols=100):
        chrom_path = output_path / f"{chrom}.npy"
        np.save(chrom_path, vector_dict[chrom].astype(np.float32), allow_pickle=False)

    meta = {
        'format': 'npy-dir-v1',
        'track_names': track_names,
        'bin_size': int(bin_size),
        'chromosomes': chromosomes,
    }
    with (output_path / "metadata.json").open("w") as f:
        json.dump(meta, f)

    print(f"Saved memory-mappable cache to: {output_dir}")


def load_multiple_bigwigs(
    bigwig_files: List[str],
    track_names: Optional[List[str]] = None,
    chrom_sizes: Optional[Dict[str, int]] = None,
    bin_size: int = 32,
    lazy: bool = False,
    cache_file: Optional[str] = None
) -> Union[Dict[str, np.ndarray], Dict]:
    """Load multiple BigWig files into stacked arrays or lazy handles.

    Args:
        bigwig_files: Paths to .bw/.bigWig files.
        track_names: Optional names matching the files.
        chrom_sizes: Chromosome lengths; inferred from first BigWig if None.
        bin_size: Bin size for aggregation.
        lazy: If True, return metadata for on-the-fly loading.
        cache_file: Optional NPZ cache path when not lazy.

    Returns:
        Dict of chrom -> array (bins, tracks) or a lazy metadata dict.
    """
    import pyBigWig
    from pathlib import Path

    if not bigwig_files:
        raise ValueError("bigwig_files cannot be empty")

    num_tracks = len(bigwig_files)

    if track_names is None:
        track_names = [os.path.basename(f) for f in bigwig_files]
    elif len(track_names) != num_tracks:
        raise ValueError(f"track_names length ({len(track_names)}) != bigwig_files ({num_tracks})")

    if cache_file and not lazy and Path(cache_file).exists():
        print(f"Loading from cache: {cache_file}")
        vector_dict, cached_track_names, cached_bin_size = load_binned_bigwigs(cache_file)
        if cached_track_names == track_names and cached_bin_size == bin_size:
            print(f"  Cache valid - loaded {len(vector_dict)} chromosomes instantly!")
            return vector_dict
        else:
            print(f"  Cache mismatch, reloading from BigWigs and refreshing cache...")
            try:
                Path(cache_file).unlink()
            except Exception as e:
                print(f"  Warning: could not remove stale cache {cache_file}: {e}")

    if lazy:
        print(f"Prepared {num_tracks} BigWig files for lazy loading")
        return {
            'bigwig_files': bigwig_files,
            'track_names': track_names,
            'num_tracks': num_tracks,
            'lazy': True
        }

    print(f"Loading {num_tracks} BigWig files into memory...")

    if chrom_sizes is None:
        print("  Inferring chromosome sizes from first BigWig...")
        bw = pyBigWig.open(bigwig_files[0])
        chrom_sizes = dict(bw.chroms())
        bw.close()

    vector_dict = {}
    for chrom, size in chrom_sizes.items():
        num_bins = size // bin_size
        vector_dict[chrom] = np.zeros((num_bins, num_tracks), dtype=np.float32)

    for track_idx, (bigwig_file, track_name) in enumerate(zip(bigwig_files, track_names)):
        bw = pyBigWig.open(bigwig_file)
        for chrom in tqdm(chrom_sizes.keys(),
                         desc=f"  [{track_idx+1}/{num_tracks}] {track_name:20s}",
                         leave=False, ncols=100):
            size = chrom_sizes[chrom]
            if chrom not in bw.chroms():
                continue
            num_bins = size // bin_size
            try:
                stats = bw.stats(chrom, 0, num_bins * bin_size, type="mean", nBins=num_bins)
                # Preserve NaNs for masking (don't use nan_to_num)
                vector_dict[chrom][:, track_idx] = np.array(stats, dtype=np.float32)
            except Exception as e:
                tqdm.write(f"  Warning: Could not load {chrom}: {e}")
        bw.close()

    print(f"\n  Successfully loaded {num_tracks} tracks")

    if cache_file:
        save_binned_bigwigs(vector_dict, cache_file, track_names, bin_size)

    return vector_dict


def split_regions_by_chromosome(
    regions_df: pd.DataFrame,
    val_chroms: Optional[List[str]] = None,
    test_chroms: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame]]:
    """Split regions into train/val/test by chromosome.

    Args:
        regions_df: DataFrame with chrom/start/end columns.
        val_chroms: Chromosomes reserved for validation.
        test_chroms: Chromosomes reserved for test (optional).

    Returns:
        train_df, val_df, test_df (test_df may be None).
    """
    if val_chroms is None:
        val_chroms = ['chr8', 'chr10']
    if test_chroms is None:
        test_chroms = ['chr9', 'chr18']

    if isinstance(val_chroms, str):
        val_chroms = [val_chroms]
    if isinstance(test_chroms, str):
        test_chroms = [test_chroms]

    val_mask = regions_df['chrom'].isin(val_chroms)
    test_mask = regions_df['chrom'].isin(test_chroms) if test_chroms else pd.Series([False] * len(regions_df))
    train_mask = ~(val_mask | test_mask)

    train_df = regions_df[train_mask].reset_index(drop=True)
    val_df = regions_df[val_mask].reset_index(drop=True)
    test_df = regions_df[test_mask].reset_index(drop=True) if test_chroms else None

    print(f"\nChromosome-based split:")
    print(f"  Training: {len(train_df):,} regions")
    print(f"  Validation: {len(val_df):,} regions ({', '.join(val_chroms)})")
    if test_df is not None:
        print(f"  Test: {len(test_df):,} regions ({', '.join(test_chroms)})")

    return train_df, val_df, test_df


class BorzoiDataset(Dataset):
    """Dataset producing (sequence, target) pairs for Borzoi models."""

    def __init__(
        self,
        regions_df: pd.DataFrame,
        vector_dict: Union[Dict[str, np.ndarray], Dict],
        genome_fasta: str,
        chrom_sizes: Dict[str, int],
        seq_len: int = 524288,
        pred_len: int = 196608,
        bin_size: int = 32,
        augment_rc: bool = True,
        augment_shift: int = 0,
        squash_targets: bool = True,
        tile_stride: int = None,  # Stride for tiling regions (if None, use pred_len)
        use_random_crop: bool = False,  # NEW: Random cropping like reference implementation
        samples_per_region: int = 1,  # Number of random crops per region per epoch
        sample_regions_fraction: float = 1.0,  # Fraction of regions to randomly sample
        data_bin_size: int = None,  # Resolution of stored data (if coarser than bin_size, targets are interpolated)
        additional_data: list = None,  # Additional data sources: list of (vector_dict, track_names, bin_size)
        track_names: Optional[List[str]] = None,  # Optional full track-name list (main + additional)
        genome_name: Optional[str] = None,  # Optional genome identifier for multi-genome training
    ):
        """Store regions and loading options.

        Args:
            regions_df: BED-like DataFrame with chrom/start/end.
            vector_dict: Pre-loaded arrays or lazy BigWig metadata.
            genome_fasta: Reference FASTA path.
            chrom_sizes: Chromosome sizes.
            seq_len: Input sequence length.
            pred_len: Center prediction window length.
            bin_size: Target bin size in bp.
            augment_rc: Random reverse complement.
            augment_shift: Random center shift in bp.
            squash_targets: Apply squash transform to targets.
            tile_stride: Stride for fixed tiling (defaults to pred_len).
            use_random_crop: Sample random windows within regions instead of fixed tiling.
            samples_per_region: Number of random crops to take from each region per epoch.
                Only applies when use_random_crop=True. Default 1 means each region
                yields one random crop per epoch. Set higher (e.g., 5) to increase
                epoch length by sampling multiple different crops from each region.
            sample_regions_fraction: Fraction of regions to randomly sample (0.0-1.0).
                Default 1.0 uses all regions. Set lower (e.g., 0.1) to use only 10%
                of randomly sampled 1MB regions for faster experimentation.
            data_bin_size: Resolution of stored data in bp. If coarser than bin_size,
                targets are upsampled via nearest-neighbor interpolation.
                If None, assumed equal to bin_size.
            additional_data: List of (vector_dict, track_names, bin_size) tuples for
                extra data at different resolutions. Tracks are appended after main tracks.
            track_names: Optional ordered names for all output tracks.
        """
        # Sample regions if requested
        if sample_regions_fraction < 1.0:
            n_regions = len(regions_df)
            n_sample = max(1, int(n_regions * sample_regions_fraction))
            sampled_regions_df = regions_df.sample(n=n_sample, random_state=42).reset_index(drop=True)
            print(f"  Sampled {n_sample:,}/{n_regions:,} regions ({sample_regions_fraction*100:.1f}%)")
            self.regions_df = sampled_regions_df
        else:
            self.regions_df = regions_df.reset_index(drop=True)
        self.genome_fasta = genome_fasta
        self.genome_name = genome_name
        self.chrom_sizes = chrom_sizes
        self.seq_len = seq_len
        self.pred_len = pred_len
        self.bin_size = bin_size
        self.augment_rc = augment_rc
        self.augment_shift = augment_shift
        self.squash_targets = squash_targets
        self.tile_stride = tile_stride if tile_stride is not None else pred_len
        self.use_random_crop = use_random_crop
        self.samples_per_region = max(1, samples_per_region) if use_random_crop else 1

        # Data resolution: if stored data is coarser than output bin_size, interpolate
        self.data_bin_size = data_bin_size if data_bin_size is not None else bin_size
        self.upsample_targets = (self.bin_size < self.data_bin_size)
        if self.upsample_targets:
            self.upsample_factor = self.data_bin_size // self.bin_size
            print(f"  Target upsampling: {self.data_bin_size}bp data -> {self.bin_size}bp "
                  f"(factor {self.upsample_factor}x nearest-neighbor)")

        # DataLoader uses multiprocessing, not threading
        # Each worker process will lazily initialize its own handles
        self._genome = None
        self._bigwig_handles = None

        self.lazy_loading = isinstance(vector_dict, dict) and vector_dict.get('lazy', False)
        provided_track_names = list(track_names) if track_names is not None else None

        if self.lazy_loading:
            self.bigwig_files = vector_dict['bigwig_files']
            self.track_names = provided_track_names or list(vector_dict['track_names'])
            self.num_tracks = vector_dict['num_tracks']
            self.vector_dict = None
        else:
            self.vector_dict = vector_dict
            self.bigwig_files = None
            first_chrom = next(iter(vector_dict.keys()))
            self.num_tracks = vector_dict[first_chrom].shape[1]
            self.track_names = provided_track_names

        # Additional data sources at different resolutions
        # Store as list of (vector_dict, track_names, bin_size, num_tracks) tuples
        self.additional_data = []
        self.main_num_tracks = self.num_tracks
        for add_vd, add_tn, add_bs in (additional_data or []):
            first_chrom = next(iter(add_vd.keys()))
            n_add = add_vd[first_chrom].shape[1]
            assert n_add == len(add_tn), (
                f"Additional data track count mismatch: {len(add_tn)} names but {n_add} columns"
            )
            self.additional_data.append((add_vd, add_tn, add_bs, n_add))
            self.num_tracks += n_add
            print(f"  Additional tracks: {n_add} at {add_bs}bp ({', '.join(add_tn)})")

        if self.track_names is not None and len(self.track_names) != self.num_tracks:
            raise ValueError(
                f"track_names length ({len(self.track_names)}) must match num_tracks ({self.num_tracks})"
            )

        self.target_bins = calculate_target_bins(self.pred_len, self.bin_size)

        self._generate_tiled_windows()

    def _generate_tiled_windows(self):
        """Generate tiled windows or store regions for random cropping."""
        if self.use_random_crop:
            valid_mask = (self.regions_df['end'] - self.regions_df['start']) >= self.pred_len
            valid_df = self.regions_df[valid_mask]

            self._window_chroms = valid_df['chrom'].values
            self._window_starts = valid_df['start'].values.astype(np.int64)
            self._window_ends = valid_df['end'].values.astype(np.int64)
            self._window_region_idx = valid_df.index.values.astype(np.int64)
            self._num_windows = len(valid_df)
            self.windows = list(range(self._num_windows))
        else:
            all_region_idx = []
            all_chroms = []
            all_pred_starts = []
            all_pred_ends = []

            for region_idx, row in self.regions_df.iterrows():
                chrom = row['chrom']
                region_start = row['start']
                region_end = row['end']
                region_len = region_end - region_start

                if region_len < self.pred_len:
                    continue

                offsets = np.arange(0, region_len - self.pred_len + 1, self.tile_stride, dtype=np.int64)
                pred_starts = region_start + offsets
                pred_ends = pred_starts + self.pred_len

                valid_mask = pred_ends <= region_end
                n_valid = valid_mask.sum()
                if n_valid > 0:
                    all_region_idx.extend([region_idx] * n_valid)
                    all_chroms.extend([chrom] * n_valid)
                    all_pred_starts.extend(pred_starts[valid_mask].tolist())
                    all_pred_ends.extend(pred_ends[valid_mask].tolist())

            self._window_region_idx = np.array(all_region_idx, dtype=np.int64)
            self._window_chroms = np.array(all_chroms)
            self._window_starts = np.array(all_pred_starts, dtype=np.int64)
            self._window_ends = np.array(all_pred_ends, dtype=np.int64)
            self._num_windows = len(all_region_idx)
            self.windows = list(range(self._num_windows))

    @property
    def genome(self):
        """Process-safe lazy loading of genome file handle.

        DataLoader uses multiprocessing, so each worker process will
        have its own instance of this file handle initialized on first access.
        """
        if self._genome is None:
            self._genome = pysam.FastaFile(self.genome_fasta)
        return self._genome

    @property
    def bigwig_handles(self):
        """Process-safe lazy loading of BigWig file handles.

        Each worker process will open its own BigWig handles on first access.
        """
        if self.lazy_loading:
            if self._bigwig_handles is None:
                import pyBigWig
                self._bigwig_handles = [pyBigWig.open(f) for f in self.bigwig_files]
            return self._bigwig_handles
        return None

    def close(self):
        """Explicitly close file handles.

        Call this method when done with the dataset to release resources.
        Prefer using the context manager protocol (with statement) for automatic cleanup.
        """
        if self._genome is not None:
            try:
                self._genome.close()
            except OSError:
                # Expected when file is already closed or handle is invalid
                pass
            except Exception as e:
                warnings.warn(f"Unexpected error closing genome file: {e}", RuntimeWarning)
            self._genome = None

        if self._bigwig_handles is not None:
            for bw in self._bigwig_handles:
                try:
                    bw.close()
                except OSError:
                    # Expected when file is already closed or handle is invalid
                    pass
                except Exception as e:
                    warnings.warn(f"Unexpected error closing BigWig file: {e}", RuntimeWarning)
            self._bigwig_handles = None

    def __enter__(self):
        """Context manager entry - returns self for use in with statements."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - ensures file handles are closed."""
        self.close()
        return False  # Don't suppress exceptions

    def __del__(self):
        """Clean up file handles when the dataset is garbage collected.

        Note: This is a fallback. Prefer using context manager or explicit close().
        """
        self.close()

    def __len__(self) -> int:
        # When using random crops with samples_per_region > 1, multiply the dataset size
        # This allows taking multiple different random crops from each region per epoch
        return len(self.windows) * self.samples_per_region

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, Dict]:
        """Return a single sample.

        Returns:
            seq_tensor: (seq_len, 4) one-hot sequence.
            target_tensor: (target_bins, tracks) coverage.
            metadata: Dict with genomic coordinates and bin size.
        """
        # Map index back to region index when using samples_per_region > 1
        # This allows multiple random crops from the same region per epoch
        region_idx = idx % len(self.windows)

        if self.use_random_crop:
            chrom = self._window_chroms[region_idx]
            region_start = self._window_starts[region_idx]
            region_end = self._window_ends[region_idx]
            region_len = region_end - region_start

            max_offset = region_len - self.pred_len
            if max_offset > 0:
                random_offset = np.random.randint(0, max_offset + 1)
                pred_start = region_start + random_offset
            else:
                pred_start = region_start

            pred_end = pred_start + self.pred_len

            seq, target = self._get_sample(chrom, pred_start, pred_end,
                                          region_start=region_start, region_end=region_end)
        else:
            chrom = self._window_chroms[region_idx]
            pred_start = self._window_starts[region_idx]
            pred_end = self._window_ends[region_idx]
            seq, target = self._get_sample(chrom, pred_start, pred_end,
                                          region_start=pred_start, region_end=pred_end)

        if self.augment_rc and np.random.random() > 0.5:
            seq = reverse_complement_onehot(seq)
            target = target[::-1].copy()

        seq_tensor = torch.from_numpy(seq)
        target_tensor = torch.from_numpy(target)

        # Create metadata with genomic coordinates
        metadata = {
            'chrom': chrom,
            'start': pred_start,
            'end': pred_end,
            'bin_size': self.bin_size,
            'genome_name': self.genome_name or "",
        }

        return seq_tensor, target_tensor, metadata

    def _get_sample(
        self,
        chrom: str,
        pred_start: int,
        pred_end: int,
        region_start: Optional[int] = None,
        region_end: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Fetch sequence and target for a window.

        Returns:
            sequence: (seq_len, 4) one-hot array.
            target: (target_bins, tracks) coverage array.
        """
        center = (pred_start + pred_end) // 2
        half_seq = self.seq_len // 2
        half_pred = self.pred_len // 2

        shift = 0
        if self.augment_shift > 0:
            shift = np.random.randint(-self.augment_shift, self.augment_shift + 1)

            if region_start is not None and region_end is not None:
                # Keep the shifted target fully inside the region bounds
                shifted_pred_start = center + shift - half_pred
                shifted_pred_end = center + shift + half_pred

                if shifted_pred_start < region_start:
                    shift = region_start - (center - half_pred)
                if shifted_pred_end > region_end:
                    shift = region_end - (center + half_pred)

        center_shifted = center + shift
        ctx_start = center_shifted - half_seq
        ctx_end = center_shifted + half_seq

        sequence = self._fetch_sequence(chrom, ctx_start, ctx_end)
        target = self._fetch_target(chrom, center_shifted)

        return sequence, target

    def _fetch_sequence(self, chrom: str, start: int, end: int) -> np.ndarray:
        """Fetch and one-hot encode a sequence."""
        chrom_size = self.chrom_sizes.get(chrom, float('inf'))

        result = np.zeros((self.seq_len, 4), dtype=np.float32)

        left_pad = max(0, -start)
        fetch_start = max(0, start)
        fetch_end = min(chrom_size, end)

        if fetch_end > fetch_start:
            seq_str = self.genome.fetch(chrom, fetch_start, fetch_end).upper()
            seq_onehot = one_hot_encode(seq_str)
            seq_len = len(seq_onehot)
            end_pos = min(left_pad + seq_len, self.seq_len)
            result[left_pad:end_pos] = seq_onehot[:end_pos - left_pad]

        return result

    def _fetch_dense_target(
        self,
        chrom: str,
        pred_start: int,
        pred_end: int,
        vector_dict: Dict[str, np.ndarray],
        source_bin_size: int,
        num_tracks: int,
    ) -> np.ndarray:
        """Fetch target from a single dense data source.

        Crops the relevant region and upsamples to self.bin_size if needed.

        Returns:
            Array of shape (target_bins, num_tracks).
        """
        start_bin = pred_start // source_bin_size
        end_bin = pred_end // source_bin_size
        data_bins = end_bin - start_bin

        if chrom not in vector_dict:
            return np.full((self.target_bins, num_tracks), np.nan, dtype=np.float32)

        chrom_vectors = vector_dict[chrom]
        max_bin = chrom_vectors.shape[0]

        # Initialize with NaN to support masking of missing data
        target = np.full((data_bins, num_tracks), np.nan, dtype=np.float32)

        valid_start_bin = max(0, start_bin)
        valid_end_bin = min(max_bin, end_bin)

        target_offset_start = valid_start_bin - start_bin
        target_offset_end = target_offset_start + (valid_end_bin - valid_start_bin)

        if valid_end_bin > valid_start_bin:
            target[target_offset_start:target_offset_end] = chrom_vectors[valid_start_bin:valid_end_bin]

        # Upsample to target resolution if source is coarser
        if source_bin_size > self.bin_size:
            upsample_factor = source_bin_size // self.bin_size
            target = np.repeat(target, upsample_factor, axis=0)

        return target

    def _fetch_target(self, chrom: str, center: int) -> np.ndarray:
        """Fetch target values for a centered region.

        When data is stored at coarser resolution than bin_size (hires mode),
        data is first loaded at data_bin_size resolution, then upsampled via
        nearest-neighbor interpolation to the target bin_size.

        When additional_data sources exist at different resolutions, each is
        fetched at its native resolution, upsampled to bin_size, and
        concatenated along the track dimension.
        """
        pred_half = self.pred_len // 2
        pred_start = center - pred_half
        pred_end = pred_start + self.pred_len

        if self.lazy_loading:
            # BigWig lazy loading: request bins at output resolution directly
            # Initialize with NaN to support masking
            target = np.full((self.target_bins, self.num_tracks), np.nan, dtype=np.float32)

            chrom_size = self.chrom_sizes.get(chrom)
            bounded_start = max(0, pred_start)
            bounded_end = pred_end if chrom_size is None else min(pred_end, chrom_size)

            if bounded_end <= bounded_start:
                return target

            for track_idx, bw_handle in enumerate(self.bigwig_handles):
                if chrom not in bw_handle.chroms():
                    continue
                try:
                    stats = bw_handle.stats(
                        chrom,
                        bounded_start,
                        bounded_end,
                        type="mean",
                        nBins=self.target_bins,
                    )
                    # Preserve NaNs (don't use nan_to_num)
                    target[:, track_idx] = np.array(stats, dtype=np.float32)
                except RuntimeError:
                    pass
                except Exception as e:
                    warnings.warn(
                        f"Error loading BigWig data for {chrom}:{bounded_start}-{bounded_end}: {e}",
                        RuntimeWarning
                    )
        else:
            # Dense matrix loading: fetch main data
            target = self._fetch_dense_target(
                chrom, pred_start, pred_end,
                self.vector_dict, self.data_bin_size, self.main_num_tracks,
            )

            # Fetch and concatenate additional data sources
            for add_vd, _add_tn, add_bs, add_ntracks in self.additional_data:
                add_target = self._fetch_dense_target(
                    chrom, pred_start, pred_end,
                    add_vd, add_bs, add_ntracks,
                )
                target = np.concatenate([target, add_target], axis=1)

        if self.squash_targets:
            target = squash_coverage(target)

        return target


def _worker_init_fn(worker_id: int) -> None:
    """Initialize per-worker random state for reproducibility.

    Each worker gets a unique seed derived from the base seed and worker_id.
    This ensures reproducible random crops and augmentations across runs.
    """
    worker_info = get_worker_info()
    if worker_info is None:
        return

    # Get base seed from dataset if available, otherwise use worker_id
    dataset = worker_info.dataset
    base_seed = getattr(dataset, 'seed', 0) or 0

    # Create unique seed for this worker
    worker_seed = base_seed + worker_id

    # Seed NumPy's global RNG for this worker
    np.random.seed(worker_seed)

    # Also seed Python's random module if used
    import random
    random.seed(worker_seed)


def create_data_loaders_from_datasets(
    train_dataset: Dataset,
    val_dataset: Dataset,
    batch_size: int,
    num_workers: int,
    accelerator,
    verbose: bool = False,
    seed: Optional[int] = None
) -> tuple:
    """Create DataLoader instances from existing datasets.

    Args:
        train_dataset: Training dataset.
        val_dataset: Validation dataset.
        batch_size: Batch size for both loaders.
        num_workers: Number of worker processes.
        accelerator: Accelerator instance for multi-GPU adjustment.
        verbose: If True, print loader statistics.
        seed: Random seed for worker initialization (for reproducibility).

    Returns:
        Tuple of (train_loader, val_loader).
    """
    # Store seed on datasets for worker_init_fn to access
    if seed is not None:
        train_dataset.seed = seed
        val_dataset.seed = seed + 1000000  # Different seed for validation

    # Adjust num_workers for multi-GPU
    if accelerator.num_processes > 1 and num_workers > 0:
        num_workers = max(1, num_workers // accelerator.num_processes)

    loader_kwargs = {
        'pin_memory': True,
        'num_workers': num_workers,
        'worker_init_fn': _worker_init_fn if num_workers > 0 else None
    }
    if num_workers > 0:
        loader_kwargs['persistent_workers'] = True
        loader_kwargs['prefetch_factor'] = 2

    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True, drop_last=True, **loader_kwargs
    )
    val_loader = DataLoader(
        val_dataset, batch_size=batch_size, shuffle=False, drop_last=False, **loader_kwargs
    )

    if verbose:
        accelerator.print(f"  Batch size: {batch_size}")
        accelerator.print(f"  Training batches: {len(train_loader)}")
        accelerator.print(f"  Validation batches: {len(val_loader)}")

    return train_loader, val_loader


def load_vector_data(config, accelerator, chrom_sizes: Dict[str, int]) -> Tuple[Union[Dict[str, np.ndarray], Dict], List[str], int, list]:
    """Load vector data from matrix or BigWig files.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance for printing.
        chrom_sizes: Chromosome sizes dictionary.

    Returns:
        Tuple of (vector_dict, track_names, data_bin_size, additional_data).
        data_bin_size is the actual resolution of the stored data.
        additional_data is a list of (vector_dict, track_names, bin_size) for extra sources.
    """
    import time

    if hasattr(config.data, 'data_matrix') and config.data.data_matrix:
        accelerator.print(f"\nLoading data from pre-computed matrix: {config.data.data_matrix}")
        start_time = time.time()

        selected_tracks = config.get('data.selected_tracks', None)
        cache_dir = config.get('performance.data_matrix_cache', None)
        use_mmap = bool(config.get('performance.data_matrix_mmap', False))

        if cache_dir:
            cache_path = Path(cache_dir)
            if cache_path.exists():
                accelerator.print(f"  Loading data matrix cache from: {cache_path}")
                vector_dict, track_names, actual_bin_size = load_from_data_matrix(
                    str(cache_path),
                    bin_size=config.model.bin_size,
                    track_cols=selected_tracks,
                    use_mmap=use_mmap
                )
            else:
                if accelerator.is_local_main_process:
                    accelerator.print(f"  Building data matrix cache at: {cache_path}")
                    vector_dict = None
                    temp_dir = cache_path.with_name(f"{cache_path.name}.tmp")
                    if temp_dir.exists():
                        if temp_dir.is_dir():
                            shutil.rmtree(temp_dir)
                        else:
                            temp_dir.unlink()
                    data_matrix_path = Path(config.data.data_matrix)
                    if data_matrix_path.suffix.lower() == ".parquet":
                        batch_size = int(config.get("performance.parquet_streaming_batch_size", 200000))
                        track_names, actual_bin_size = build_parquet_cache_streaming(
                            config.data.data_matrix,
                            str(temp_dir),
                            bin_size=config.model.bin_size,
                            track_cols=selected_tracks,
                            batch_size=batch_size,
                        )
                    else:
                        vector_dict, track_names, actual_bin_size = load_from_data_matrix(
                            config.data.data_matrix,
                            bin_size=config.model.bin_size,
                            track_cols=selected_tracks
                        )
                        save_to_npy_dir(vector_dict, str(temp_dir), track_names, actual_bin_size)
                    if cache_path.exists():
                        if cache_path.is_dir():
                            shutil.rmtree(cache_path)
                        else:
                            cache_path.unlink()
                    os.replace(temp_dir, cache_path)
                    if vector_dict is not None:
                        del vector_dict
                        import gc
                        gc.collect()
                accelerator.wait_for_everyone()
                vector_dict, track_names, actual_bin_size = load_from_data_matrix(
                    str(cache_path),
                    bin_size=config.model.bin_size,
                    track_cols=selected_tracks,
                    use_mmap=use_mmap
                )
        else:
            vector_dict, track_names, actual_bin_size = load_from_data_matrix(
                config.data.data_matrix,
                bin_size=config.model.bin_size,
                track_cols=selected_tracks
            )

        hires_resolution = config.get('model.hires_resolution', None)
        if actual_bin_size != config.model.bin_size and hires_resolution is None:
            # Standard mode: update config bin_size to match data
            if config.model.pred_len % actual_bin_size != 0:
                raise ValueError(
                    f"Data matrix bin_size={actual_bin_size} does not divide pred_len={config.model.pred_len}"
                )
            accelerator.print(f"  Detected bin_size {actual_bin_size} differs from config; updating")
            config.set('model.bin_size', actual_bin_size)
        elif actual_bin_size != config.model.bin_size and hires_resolution is not None:
            # HiRes mode: data is coarser than model output, will interpolate
            accelerator.print(f"  HiRes mode: data at {actual_bin_size}bp, "
                            f"model output at {hires_resolution}bp (will interpolate)")

        elapsed = time.time() - start_time
        accelerator.print(f"  Data matrix loaded in {elapsed:.1f} seconds")

        # Load additional data matrices at different resolutions
        additional_configs = config.get('data.additional_data_matrices', [])
        additional_data = []

        # Remove tracks from main data that are superseded by additional sources
        if additional_configs:
            superseded_tracks = set()
            for add_cfg in additional_configs:
                tc = add_cfg.get('track_col', None) if isinstance(add_cfg, dict) else None
                if tc:
                    superseded_tracks.add(tc)
            overlap = superseded_tracks & set(track_names)
            if overlap:
                accelerator.print(f"  Replacing main tracks with higher-resolution versions: {overlap}")
                keep_indices = [i for i, tn in enumerate(track_names) if tn not in overlap]
                track_names = [track_names[i] for i in keep_indices]
                for chrom in vector_dict:
                    vector_dict[chrom] = vector_dict[chrom][:, keep_indices]

        for add_cfg in additional_configs:
            add_path = add_cfg['path'] if isinstance(add_cfg, dict) else add_cfg
            add_track_col = add_cfg.get('track_col', None) if isinstance(add_cfg, dict) else None
            add_track_cols = [add_track_col] if add_track_col else None
            accelerator.print(f"\nLoading additional data matrix: {add_path}")
            add_vd, add_tn, add_bs = load_from_data_matrix(
                add_path,
                bin_size=config.model.bin_size,
                track_cols=add_track_cols,
                use_mmap=use_mmap,
            )
            accelerator.print(f"  Additional tracks ({add_bs}bp): {', '.join(add_tn)}")
            additional_data.append((add_vd, add_tn, add_bs))

        # Validate total track count (main + additional)
        all_track_names = list(track_names)
        for _, add_tn, _ in additional_data:
            all_track_names.extend(add_tn)
        validate_track_count(config, accelerator, len(all_track_names), all_track_names)

        return vector_dict, track_names, actual_bin_size, additional_data

    elif config.data.bigwig_folder:
        bigwig_folder = Path(config.data.bigwig_folder)
        bigwig_files = sorted(list(bigwig_folder.glob("*.bw")) + list(bigwig_folder.glob("*.bigWig")))

        if not bigwig_files:
            raise ValueError(f"No BigWig files found in {bigwig_folder}")

        bigwig_paths = [str(f) for f in bigwig_files]
        track_names = [f.stem for f in bigwig_files]
        accelerator.print(f"\nAuto-discovered {len(bigwig_paths)} BigWig files")

        validate_track_count(config, accelerator, len(track_names), track_names)

        # For BigWig lazy loading, bin_size is the model output resolution
        bw_bin_size = config.model.bin_size
        vector_dict = load_multiple_bigwigs(
            bigwig_files=bigwig_paths,
            track_names=track_names,
            chrom_sizes=chrom_sizes,
            bin_size=bw_bin_size,
            lazy=config.performance.lazy_loading,
            cache_file=getattr(config.performance, 'cache_file', None)
        )
        return vector_dict, track_names, bw_bin_size, []

    elif config.data.bigwig_files:
        bigwig_paths = [track['path'] for track in config.data.bigwig_files]
        track_names = [track['track'] for track in config.data.bigwig_files]
        accelerator.print(f"\nLoading {len(bigwig_paths)} BigWig files from config")

        validate_track_count(config, accelerator, len(track_names), track_names)

        bw_bin_size = config.model.bin_size
        vector_dict = load_multiple_bigwigs(
            bigwig_files=bigwig_paths,
            track_names=track_names,
            chrom_sizes=chrom_sizes,
            bin_size=bw_bin_size,
            lazy=config.performance.lazy_loading,
            cache_file=getattr(config.performance, 'cache_file', None)
        )
        return vector_dict, track_names, bw_bin_size, []

    else:
        raise ValueError("Must provide one of: data.data_matrix, data.bigwig_folder, or data.bigwig_files")


def validate_track_count(config, accelerator, detected_tracks: int, track_names: List[str]) -> None:
    """Validate or auto-set num_output_tasks based on detected tracks.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance for printing.
        detected_tracks: Number of detected tracks.
        track_names: List of track names.
    """
    if config.model.num_output_tasks is None or config.model.num_output_tasks <= 0:
        config.set('model.num_output_tasks', detected_tracks)
        accelerator.print(f"  Auto-detected num_output_tasks={detected_tracks}")
    elif config.model.num_output_tasks != detected_tracks:
        raise ValueError(
            f"Configuration mismatch: config specifies num_output_tasks={config.model.num_output_tasks}, "
            f"but data has {detected_tracks} tracks."
        )


def load_or_generate_regions(config, accelerator, chrom_sizes: Dict[str, int]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load or generate regions for training and validation.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance for printing.
        chrom_sizes: Chromosome sizes dictionary.

    Returns:
        Tuple of (train_regions, val_regions).
    """
    regions_file = getattr(config.data, 'regions_file', None)

    if regions_file:
        if config.data.split_by_chrom:
            accelerator.print(f"\nLoading sparse regions from: {regions_file}")
            all_regions = load_regions(regions_file)
            train_regions, val_regions, _ = split_regions_by_chromosome(
                all_regions,
                val_chroms=config.data.val_chroms,
                test_chroms=config.data.test_chroms
            )
        else:
            accelerator.print(f"\nLoading pre-split region files")
            train_regions = load_regions(config.data.train_regions)
            val_regions = load_regions(config.data.val_regions)
    else:
        accelerator.print(f"\nGenerating genome-wide regions (no regions_file specified)")
        tile_size = getattr(config.data, 'genome_tile_size', 1000000)
        exclude_chroms = getattr(config.data, 'exclude_chroms', ['chrM'])

        all_regions = generate_genome_wide_regions(
            chrom_sizes=chrom_sizes,
            tile_size=tile_size,
            exclude_chroms=exclude_chroms,
            pred_len=config.model.pred_len
        )
        accelerator.print(f"  Generated {len(all_regions):,} genome-wide tiles (tile_size={tile_size:,}bp)")

        train_regions, val_regions, _ = split_regions_by_chromosome(
            all_regions,
            val_chroms=config.data.val_chroms,
            test_chroms=config.data.test_chroms
        )

    return train_regions, val_regions


def _build_multi_genome_dataset(
    genomes_config: list,
    regions_df: pd.DataFrame,
    vector_dict: Dict,
    chrom_sizes: Dict[str, int],
    seq_len: int,
    pred_len: int,
    bin_size: int,
    augment_rc: bool,
    augment_shift: int,
    squash_targets: bool,
    tile_stride: int,
    use_random_crop: bool,
    default_samples_per_region: int,
    sample_regions_fraction: float,
    data_bin_size: int,
    additional_data: list,
    track_names: List[str],
) -> ConcatDataset:
    """Build a ConcatDataset with one BorzoiDataset per genome.

    All sub-datasets share the same vector_dict, regions_df, and chrom_sizes
    by reference (no memory duplication). Each independently manages its own
    pysam.FastaFile handle.

    Args:
        genomes_config: List of genome dicts with 'fasta', optional 'name' and
            'samples_per_region' keys.
        default_samples_per_region: Fallback samples_per_region when not specified
            per genome.
        Other args: Passed through to BorzoiDataset.

    Returns:
        ConcatDataset wrapping one BorzoiDataset per genome, with a .track_names
        attribute copied from the sub-datasets.
    """
    datasets = []
    for g in genomes_config:
        genome_name = g.get('name', Path(g['fasta']).stem)
        spr = g.get('samples_per_region', default_samples_per_region)

        ds = BorzoiDataset(
            regions_df=regions_df,
            vector_dict=vector_dict,
            genome_fasta=g['fasta'],
            chrom_sizes=chrom_sizes,
            seq_len=seq_len,
            pred_len=pred_len,
            bin_size=bin_size,
            augment_rc=augment_rc,
            augment_shift=augment_shift,
            squash_targets=squash_targets,
            tile_stride=tile_stride,
            use_random_crop=use_random_crop,
            samples_per_region=spr,
            sample_regions_fraction=sample_regions_fraction,
            data_bin_size=data_bin_size,
            additional_data=additional_data,
            track_names=track_names,
            genome_name=genome_name,
        )
        datasets.append(ds)

    concat_ds = ConcatDataset(datasets)
    concat_ds.track_names = track_names
    return concat_ds


def log_dataset_config(config, accelerator, train_regions: pd.DataFrame, train_dataset: Dataset,
                       val_regions: pd.DataFrame, val_dataset: Dataset, samples_per_region: int) -> None:
    """Log dataset configuration details.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance for printing.
        train_regions: Training regions DataFrame.
        train_dataset: Training dataset.
        val_regions: Validation regions DataFrame.
        val_dataset: Validation dataset.
        samples_per_region: Number of samples per region.
    """
    accelerator.print(f"\nDataset configuration:")
    if isinstance(train_dataset, ConcatDataset):
        accelerator.print(f"  Multi-genome training:")
        for i, sub_ds in enumerate(train_dataset.datasets):
            gname = getattr(sub_ds, 'genome_name', f'genome_{i}')
            accelerator.print(f"    [{i}] {gname}: {len(sub_ds):,} samples")
        accelerator.print(f"    Total: {len(train_dataset):,} samples")
    if config.augmentation.use_random_crop:
        accelerator.print(f"  Training mode: Random cropping (like reference implementation)")
        if samples_per_region > 1:
            accelerator.print(f"  Samples per region: {samples_per_region} (×{samples_per_region} epoch length)")
        accelerator.print(f"  Training: {len(train_regions):,} regions → {len(train_dataset):,} samples (randomly cropped each epoch)")
    else:
        accelerator.print(f"  Training mode: Fixed tiling")
        accelerator.print(f"  Tile stride: {config.augmentation.tile_stride or config.model.pred_len} bp")
        accelerator.print(f"  Training: {len(train_regions):,} regions → {len(train_dataset):,} windows")
    val_pred_len = config.get('validation.pred_len', config.model.pred_len)
    if val_pred_len != config.model.pred_len:
        accelerator.print(f"  Validation: {len(val_regions):,} regions → {len(val_dataset):,} windows "
                         f"(pred_len={val_pred_len} bp, center-cropped from {config.model.pred_len} bp)")
    else:
        accelerator.print(f"  Validation: {len(val_regions):,} regions → {len(val_dataset):,} windows (pred_len={val_pred_len} bp)")


def create_data_loaders(config, accelerator, batch_size: Optional[int] = None) -> tuple:
    """Create training and validation data loaders.

    Args:
        config: Configuration object.
        accelerator: Accelerator instance.
        batch_size: Batch size (uses config.training.batch_size if None).

    Returns:
        Tuple of (train_loader, val_loader, train_dataset, val_dataset).
    """
    if batch_size is None:
        batch_size = config.training.batch_size

    accelerator.print("\n" + "=" * 70)
    accelerator.print("Loading Data")
    accelerator.print("=" * 70)

    # Load chromosome sizes
    accelerator.print(f"Loading chromosome sizes from: {config.data.chrom_sizes}")
    chrom_sizes = load_chrom_sizes(config.data.chrom_sizes)
    accelerator.print(f"  Loaded {len(chrom_sizes)} chromosomes")

    # Load data - either from pre-computed matrix or BigWig files
    vector_dict, track_names, data_bin_size, additional_data = load_vector_data(config, accelerator, chrom_sizes)

    # Load or generate regions
    train_regions, val_regions = load_or_generate_regions(config, accelerator, chrom_sizes)

    accelerator.print(f"\nDataset sizes:")
    accelerator.print(f"  Training regions: {len(train_regions):,}")
    accelerator.print(f"  Validation regions: {len(val_regions):,}")

    # Create datasets
    # For hires mode, output bin_size is the target resolution (e.g., 1bp),
    # but data may be stored at coarser resolution (e.g., 32bp).
    hires_resolution = config.get('model.hires_resolution', None)
    output_bin_size = hires_resolution if hires_resolution is not None else config.model.bin_size
    all_track_names = list(track_names)
    for _, add_tn, _ in additional_data:
        all_track_names.extend(add_tn)

    samples_per_region = config.get('augmentation.samples_per_region', 1)
    sample_regions_fraction = config.get('augmentation.sample_regions_fraction', 1.0)

    # Validation pred_len: allows comparing models with different native pred_len
    # using the same effective pred_len for metric computation.
    # - val_pred_len: the bp to evaluate (center-cropped from model output)
    # - stride = val_pred_len for non-overlapping evaluation
    # Default to model.pred_len for backward compatibility
    val_pred_len = config.get('validation.pred_len', config.model.pred_len)

    genomes_config = config.get('data.genomes')

    if genomes_config:
        # Multi-genome: build ConcatDataset for training, primary genome for validation
        accelerator.print(f"\nMulti-genome mode: {len(genomes_config)} genomes")
        for i, g in enumerate(genomes_config):
            spr = g.get('samples_per_region', samples_per_region)
            gname = g.get('name', Path(g['fasta']).stem)
            accelerator.print(f"  [{i}] {gname}: {g['fasta']} (samples_per_region={spr})")

        train_dataset = _build_multi_genome_dataset(
            genomes_config=genomes_config,
            regions_df=train_regions,
            vector_dict=vector_dict,
            chrom_sizes=chrom_sizes,
            seq_len=config.model.seq_len,
            pred_len=config.model.pred_len,
            bin_size=output_bin_size,
            augment_rc=config.augmentation.reverse_complement,
            augment_shift=config.augmentation.shift_max,
            squash_targets=config.augmentation.squash_targets,
            tile_stride=config.augmentation.tile_stride,
            use_random_crop=config.augmentation.use_random_crop,
            default_samples_per_region=samples_per_region,
            sample_regions_fraction=sample_regions_fraction,
            data_bin_size=data_bin_size,
            additional_data=additional_data,
            track_names=all_track_names,
        )

        # Validation uses primary (first) genome only
        primary = genomes_config[0]
        val_dataset = BorzoiDataset(
            regions_df=val_regions,
            vector_dict=vector_dict,
            genome_fasta=primary['fasta'],
            chrom_sizes=chrom_sizes,
            seq_len=config.model.seq_len,
            pred_len=config.model.pred_len,
            bin_size=output_bin_size,
            augment_rc=False,
            augment_shift=0,
            squash_targets=config.augmentation.squash_targets,
            tile_stride=val_pred_len,
            use_random_crop=False,
            sample_regions_fraction=1.0,
            data_bin_size=data_bin_size,
            additional_data=additional_data,
            track_names=all_track_names,
            genome_name=primary.get('name', Path(primary['fasta']).stem),
        )
    else:
        # Single-genome mode (backward compatible)
        train_dataset = BorzoiDataset(
            regions_df=train_regions,
            vector_dict=vector_dict,
            genome_fasta=config.data.genome_fasta,
            chrom_sizes=chrom_sizes,
            seq_len=config.model.seq_len,
            pred_len=config.model.pred_len,
            bin_size=output_bin_size,
            augment_rc=config.augmentation.reverse_complement,
            augment_shift=config.augmentation.shift_max,
            squash_targets=config.augmentation.squash_targets,
            tile_stride=config.augmentation.tile_stride,
            use_random_crop=config.augmentation.use_random_crop,
            samples_per_region=samples_per_region,
            sample_regions_fraction=sample_regions_fraction,
            data_bin_size=data_bin_size,
            additional_data=additional_data,
            track_names=all_track_names,
        )

        val_dataset = BorzoiDataset(
            regions_df=val_regions,
            vector_dict=vector_dict,
            genome_fasta=config.data.genome_fasta,
            chrom_sizes=chrom_sizes,
            seq_len=config.model.seq_len,
            pred_len=config.model.pred_len,
            bin_size=output_bin_size,
            augment_rc=False,
            augment_shift=0,
            squash_targets=config.augmentation.squash_targets,
            tile_stride=val_pred_len,
            use_random_crop=False,
            sample_regions_fraction=1.0,
            data_bin_size=data_bin_size,
            additional_data=additional_data,
            track_names=all_track_names,
        )

    log_dataset_config(config, accelerator, train_regions, train_dataset, val_regions, val_dataset, samples_per_region)

    # Create data loaders
    train_loader, val_loader = create_data_loaders_from_datasets(
        train_dataset, val_dataset, batch_size, config.training.num_workers, accelerator, verbose=True
    )

    accelerator.print("=" * 70)

    return train_loader, val_loader, train_dataset, val_dataset
