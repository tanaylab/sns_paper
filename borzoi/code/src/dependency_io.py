"""
I/O utilities for compute_dependency_matrix.

This module contains file loading and saving functions extracted from
compute_dependency_matrix.py for better code organization.
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml


# =============================================================================
# Configuration Loaders
# =============================================================================

def load_config(config_path: str) -> dict:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_chrom_sizes(path: str) -> dict:
    """Load chromosome sizes from file."""
    chrom_sizes = {}
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def load_track_names(config: dict) -> list:
    """Load track names from config or data."""
    if 'track_names' in config.get('data', {}):
        return config['data']['track_names']
    if 'track_names_file' in config.get('data', {}):
        with open(config['data']['track_names_file'], 'r') as f:
            return [line.strip() for line in f]
    data_matrix = config.get('data', {}).get('data_matrix')
    if data_matrix:
        try:
            suffix = Path(data_matrix).suffix.lower()
            if suffix == '.parquet':
                schema = pq.read_schema(data_matrix)
                exclude = {'chrom', 'chr', 'chromosome', 'start', 'end', 'region'}
                return [name for name in schema.names if name not in exclude]
        except Exception as e:
            print(f"Warning: Could not infer track names from {data_matrix}: {e}")
    return []


def load_operations_tsv(path: str, analysis_window_bp: int, analysis_resolution: int) -> List[Dict]:
    """
    Load batch operations from TSV and validate dimensions.

    All operations must use the same analysis_window_bp and analysis_resolution.

    Args:
        path: Path to TSV file with columns: chrom, center, donor_chrom, donor_start, donor_end
        analysis_window_bp: Analysis window size in bp (must be same for all operations)
        analysis_resolution: Analysis resolution in bp (must be same for all operations)

    Returns:
        List of operation dictionaries with keys: chrom, center, donor_chrom, donor_start, donor_end, operation_id, idx

    Raises:
        ValueError: If TSV is missing required columns
    """
    # Read TSV file
    df = pd.read_csv(path, sep='\t' if path.endswith('.tsv') else ',')

    # Validate required columns
    required = ['chrom', 'center', 'donor_chrom', 'donor_start', 'donor_end']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"TSV missing required columns: {missing}")

    # Convert to list of dictionaries (using itertuples for better performance)
    operations = []
    for idx, row in enumerate(df.itertuples(index=False)):
        op_id = f"{row.chrom}_{row.center}_donor_{row.donor_chrom}_{row.donor_start}-{row.donor_end}"
        operations.append({
            'chrom': str(row.chrom),
            'center': int(row.center),
            'donor_chrom': str(row.donor_chrom),
            'donor_start': int(row.donor_start),
            'donor_end': int(row.donor_end),
            'operation_id': op_id,
            'idx': idx,
        })

    print(f"Loaded {len(operations)} operations from {path}")
    print(f"  All operations will use:")
    print(f"    analysis_window_bp: {analysis_window_bp}")
    print(f"    analysis_resolution: {analysis_resolution}")

    return operations


# =============================================================================
# Observed Data Fetching
# =============================================================================

def fetch_observed_values(
    data_matrix_path: str,
    chrom: str,
    analysis_start: int,
    analysis_end: int,
    track_name: str,
    model_bin_size: int = 32,
    analysis_resolution: int = 192,
) -> np.ndarray:
    """
    Fetch observed values from data_matrix parquet for an analysis window.

    This function fetches observed values at model resolution (e.g., 32bp) and then
    aggregates them using the same mean-pooling logic as predictions, ensuring
    consistent comparison between predicted and observed values.

    Args:
        data_matrix_path: Path to the parquet file containing observed data.
        chrom: Chromosome name (e.g., 'chr5').
        analysis_start: Start genomic position of the analysis window.
        analysis_end: End genomic position of the analysis window.
        track_name: Name of the track column in the parquet file.
        model_bin_size: Bin size in the data_matrix (default: 32bp).
        analysis_resolution: Output resolution (default: 192bp).

    Returns:
        1D numpy array of observed values, aggregated to analysis_resolution
        using the same mean-pooling as predictions.
    """
    if data_matrix_path is None:
        return None

    # Read only the necessary columns for the chromosome
    try:
        # Use predicate pushdown for efficient filtering
        table = pq.read_table(
            data_matrix_path,
            columns=['chrom', 'start', track_name],
            filters=[('chrom', '=', chrom)]
        )
    except Exception as e:
        print(f"  Warning: Failed to load observed values for {chrom}: {e}")
        return None

    if len(table) == 0:
        print(f"  Warning: No data found for {chrom} in {data_matrix_path}")
        return None

    # Convert to pandas for easier indexing
    df = table.to_pandas()

    # Calculate bin indices
    start_bin = analysis_start // model_bin_size
    end_bin = analysis_end // model_bin_size
    num_model_bins = end_bin - start_bin

    # Build observed values array at model resolution
    observed_model_res = np.full(num_model_bins, np.nan, dtype=np.float32)

    # Extract data
    start_positions = df['start'].values
    track_values = df[track_name].values

    # Find indices of bins within the analysis window
    data_bin_indices = (start_positions // model_bin_size).astype(np.int64)
    mask = (data_bin_indices >= start_bin) & (data_bin_indices < end_bin)
    valid_data_bins = data_bin_indices[mask]
    valid_values = track_values[mask]

    # Vectorized placement into output array (much faster than Python loop)
    output_indices = valid_data_bins - start_bin
    valid_mask = (output_indices >= 0) & (output_indices < num_model_bins)
    observed_model_res[output_indices[valid_mask]] = valid_values[valid_mask]

    # Check for missing data
    nan_count = np.sum(np.isnan(observed_model_res))
    if nan_count > 0:
        print(f"  Warning: {nan_count}/{num_model_bins} bins missing observed data")

    # Apply the same aggregation as predictions: mean-pooling from model_bin_size to analysis_resolution
    if analysis_resolution != model_bin_size:
        factor = analysis_resolution // model_bin_size
        n_target_bins = num_model_bins // factor
        trimmed_len = n_target_bins * factor

        # Trim to multiple of factor
        observed_trimmed = observed_model_res[:trimmed_len]

        # Reshape and mean (same as _aggregate_to_resolution)
        observed_reshaped = observed_trimmed.reshape(n_target_bins, factor)
        # Use nanmean to handle any NaN values gracefully
        observed_aggregated = np.nanmean(observed_reshaped, axis=1).astype(np.float32)

        return observed_aggregated

    return observed_model_res


def add_observed_to_result(
    result: Dict,
    data_matrix_path: Optional[str],
    track_name: Optional[str],
    chrom: str,
    model_bin_size: int = 32,
) -> None:
    """
    Add observed values to a result dictionary in-place.

    Args:
        result: Result dictionary with 'readout_positions' and 'metadata' keys
        data_matrix_path: Path to parquet data matrix (can be None)
        track_name: Track name to fetch (can be None)
        chrom: Chromosome name
        model_bin_size: Model bin size in bp
    """
    if not data_matrix_path or not track_name:
        return

    analysis_res = result['metadata']['analysis_resolution']
    readout_positions = result['readout_positions']

    # Compute analysis window from readout positions
    obs_analysis_start = int(readout_positions[0] - analysis_res // 2)
    obs_analysis_end = int(readout_positions[-1] + analysis_res // 2 + analysis_res)

    observed_reference = fetch_observed_values(
        data_matrix_path=data_matrix_path,
        chrom=chrom,
        analysis_start=obs_analysis_start,
        analysis_end=obs_analysis_end,
        track_name=track_name,
        model_bin_size=model_bin_size,
        analysis_resolution=analysis_res,
    )

    if observed_reference is not None:
        result['observed_reference'] = observed_reference
        result['metadata']['observed_data_source'] = data_matrix_path
        result['metadata']['observed_track_name'] = track_name


# =============================================================================
# Save Functions
# =============================================================================

def save_comparison_result(result: Dict, path: str, format: str = "parquet", mode: str = "compare") -> None:
    """
    Save comparison result (1D arrays, not matrix).

    Works for both swap_cg_global and compare modes.

    Args:
        result: Dictionary from compute_global_perturbation() or compare_with_replacement()
        path: Output file path
        format: Output format (parquet, csv, npz)
        mode: Mode name for output file suffix ('global' or 'compare')
    """
    path = Path(path)
    base_path = path.parent / path.stem
    suffix = "_global" if mode == "swap_cg_global" else "_compare"

    # Handle different key names between modes
    mod_key = 'perturbed_prediction' if 'perturbed_prediction' in result else 'modified_prediction'
    has_observed = 'observed_reference' in result

    if format == "npz":
        save_dict = {
            'reference_prediction': result['reference_prediction'],
            'modified_prediction': result[mod_key],
            'difference': result['difference'],
            'readout_positions': result['readout_positions'],
            'metadata': result['metadata'],
        }
        if has_observed:
            save_dict['observed_reference'] = result['observed_reference']
        np.savez_compressed(str(base_path) + f"{suffix}.npz", **save_dict)
    elif format in ("csv", "rds", "parquet"):
        # Create DataFrame with all data
        df_data = {
            'chrom': result['metadata']['chrom'],
            'readout_position': result['readout_positions'],
            'reference_prediction': result['reference_prediction'],
            'modified_prediction': result[mod_key],
            'difference': result['difference'],
        }
        if has_observed:
            df_data['observed_reference'] = result['observed_reference']
        df = pd.DataFrame(df_data)

        if format == "parquet":
            df.to_parquet(str(base_path) + f"{suffix}.parquet", index=False)
        else:
            df.to_csv(str(base_path) + f"{suffix}.csv", index=False)

        # Save metadata
        metadata = result['metadata'].copy()
        for key, value in metadata.items():
            if isinstance(value, np.integer):
                metadata[key] = int(value)
            elif isinstance(value, np.floating):
                metadata[key] = float(value)
        with open(str(base_path) + f"{suffix}_metadata.json", 'w') as f:
            json.dump(metadata, f, indent=2)
    else:
        raise ValueError(f"Unknown format: {format}")

    # Print summary
    ext = format if format == 'parquet' else 'csv'
    print(f"\nComparison result saved:")
    if mode == "swap_cg_global" and 'swap_count' in result['metadata']:
        print(f"  Swap count: {result['metadata']['swap_count']} CG dinucleotides → GC")
    elif mode == "compare":
        print(f"  Replacement: {result['metadata'].get('replacement_start', '?')}-{result['metadata'].get('replacement_end', '?')}")
        print(f"  Replacement length: {result['metadata'].get('replacement_length', '?')}bp")
    print(f"  Output: {base_path}{suffix}.{ext}")


def save_global_perturbation(result: Dict, path: str, format: str = "parquet") -> None:
    """Alias for save_comparison_result with mode='swap_cg_global'."""
    save_comparison_result(result, path, format, mode="swap_cg_global")


def print_output_summary(result: Dict, output_path: str, format: str) -> None:
    """Print detailed summary of saved output files."""
    base_path = Path(output_path).parent / Path(output_path).stem

    print("\n" + "=" * 80)
    print("OUTPUT SUMMARY")
    print("=" * 80)

    # Matrix dimensions
    n_pert, n_readout = result['matrix'].shape
    has_predictions = 'predictions_matrix' in result
    has_observed = 'observed_reference' in result

    print(f"\nRegion: {result['metadata']['chrom']}:{result['metadata']['center']}")
    print(f"Mode: {result['metadata'].get('mode', 'scramble')}")
    if 'donor_info' in result:
        donor = result['donor_info']
        if 'chrom' in donor:
            print(f"Donor: {donor['chrom']}:{donor['start']}-{donor['end']} ({donor['length']} bp)")
        else:
            print(f"Donor: {donor['length']} bp sequence")

    print(f"\nMatrix dimensions:")
    print(f"  Perturbation positions: {n_pert}")
    print(f"  Readout positions: {n_readout}")
    print(f"  Resolution: {result['metadata']['analysis_resolution']} bp")

    # Value ranges
    print(f"\nValue ranges:")
    ref_pred = result['reference_prediction']
    print(f"  Reference prediction: {ref_pred.min():.2f} to {ref_pred.max():.2f} (mean: {ref_pred.mean():.2f})")

    diff_matrix = result['matrix']
    print(f"  Differences: {diff_matrix.min():.2f} to {diff_matrix.max():.2f} (mean: {diff_matrix.mean():.2f})")

    if has_predictions:
        pred_matrix = result['predictions_matrix']
        print(f"  Perturbed predictions: {pred_matrix.min():.2f} to {pred_matrix.max():.2f} (mean: {pred_matrix.mean():.2f})")

    if has_observed:
        obs_ref = result['observed_reference']
        # Handle potential NaN values
        valid_obs = obs_ref[~np.isnan(obs_ref)] if np.any(np.isnan(obs_ref)) else obs_ref
        if len(valid_obs) > 0:
            print(f"  Observed reference: {valid_obs.min():.2f} to {valid_obs.max():.2f} (mean: {valid_obs.mean():.2f})")
            # Compute correlation with predicted reference
            valid_mask = ~np.isnan(obs_ref)
            if np.sum(valid_mask) > 1:
                corr = np.corrcoef(ref_pred[valid_mask], obs_ref[valid_mask])[0, 1]
                print(f"  Pred-Obs correlation: {corr:.4f}")

    # Files created
    print(f"\nFiles created ({format} format):")

    if format == "npz":
        files = [f"{output_path}.npz"]
    elif format in ("csv", "rds"):
        files = [
            f"{base_path}_diff_matrix.csv",
            f"{base_path}_reference.csv",
            f"{base_path}_metadata.json",
            f"{base_path}_load.R"
        ]
        if has_predictions:
            files.insert(1, f"{base_path}_pred_matrix.csv")
    elif format == "parquet":
        files = [
            f"{base_path}.parquet",
            f"{base_path}_diff_wide.parquet"
        ]
        if has_predictions:
            files.insert(1, f"{base_path}_pred_wide.parquet")
    elif format == "hdf5":
        files = [f"{base_path}.h5"]
    else:
        files = []

    for i, filepath in enumerate(files, 1):
        if os.path.exists(filepath):
            size_mb = os.path.getsize(filepath) / (1024 * 1024)
            print(f"  {i}. {os.path.basename(filepath)} ({size_mb:.2f} MB)")

    print("=" * 80)
