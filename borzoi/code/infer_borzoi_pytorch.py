#!/usr/bin/env python3
"""
Genome-wide inference script for fine-tuned Borzoi models.

Supports three input modes:
1. Whole genome: Slide across all chromosomes
2. Selected chromosomes: Slide across specified chromosomes
3. BED regions: Predict on specific regions

Sliding strategy:
- The model takes seq_len bp as input but only predicts the center pred_len bp
- We slide by pred_len to see each bin exactly once (non-overlapping predictions)
- Example for 64k config: seq_len=65536, pred_len=24576, stride=24576

Multi-GPU support via torch.distributed for fast inference.

Output formats:
- BigWig: One file per track (recommended for genome-wide)
- HDF5: Single file with all predictions (recommended for BED regions)
- NPZ: NumPy compressed archive

OOM Safety:
- Automatic batch size reduction on out-of-memory errors
- Configurable minimum batch size and reduction factor
- Automatic retry with smaller batches to prevent crashes
- Each GPU handles OOM independently (no cross-GPU synchronization needed)

Usage:
    # Minimal usage - all parameters from config
    python infer_borzoi_pytorch.py --config configs/pcg_64k.yaml

    # Single GPU - whole genome (with explicit parameters)
    python infer_borzoi_pytorch.py --config configs/pcg_64k.yaml \
        --checkpoint checkpoints/best_model.pt \
        --mode genome \
        --output_dir predictions/

    # Single GPU - selected chromosomes
    python infer_borzoi_pytorch.py --config configs/pcg_64k.yaml \
        --checkpoint checkpoints/best_model.pt \
        --mode chromosomes --chromosomes chr1 chr2 chr3 \
        --output_dir predictions/

    # Single GPU - BED regions
    python infer_borzoi_pytorch.py --config configs/pcg_64k.yaml \
        --checkpoint checkpoints/best_model.pt \
        --mode bed --bed_file regions.bed \
        --output_dir predictions/

    # Multi-GPU (torchrun)
    torchrun --nproc_per_node=4 infer_borzoi_pytorch.py \
        --config configs/pcg_64k.yaml

    # With custom OOM safety settings
    python infer_borzoi_pytorch.py --config configs/pcg_64k.yaml \
        --batch_size 16 \
        --min_batch_size 2 \
        --batch_reduction_factor 0.5

Note: Most parameters can be set in the config file under the 'inference' section.
      Command-line arguments override config values when provided.
"""

import os
import re
import sys
import argparse
import glob
import yaml
import time
import gc
import multiprocessing
from datetime import timedelta
import json
import gzip
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Iterable
from pathlib import Path
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.distributed as dist
from torch.utils.data import Dataset, DataLoader
import pysam

# Match training defaults for quieter TF logs and faster kernels
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
os.environ.setdefault('TF_ENABLE_ONEDNN_OPTS', '0')
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.benchmark = True

# Import from src
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from src.data_loaders import one_hot_encode, load_chrom_sizes
from src.utils import reverse_complement_torch, DotDict, load_yaml_config
from src.perturbation_simple import load_splice_pert, load_seq_pert, SplicePert, SeqPert
from src.training_utils import (
    get_metric_direction,
    is_better,
    get_initial_best_value,
    load_metrics_json,
    is_checkpoint_complete,
)


def _summarize_track_names(track_names: List[str], max_items: int = 12) -> str:
    if len(track_names) <= max_items:
        return ", ".join(track_names)
    head = ", ".join(track_names[:max_items])
    return f"{head}, ... (+{len(track_names) - max_items} more)"


def _load_track_names_from_data_matrix(
    data_matrix_path: str,
    num_tracks: int,
    rank: int = 0,
    cache_dir: Optional[str] = None,
    selected_tracks: Optional[List[str]] = None,
) -> Optional[List[str]]:
    """Load track names from a Parquet data matrix without materializing the full table.

    When ``selected_tracks`` is provided (a subset of columns used during
    training), the function filters the available track names to that subset
    instead of falling back to generic names when the full column count doesn't
    match ``num_tracks``.
    """
    # If selected_tracks is provided and already matches model output count,
    # use it directly — no need to read parquet/cache at all.
    if selected_tracks and len(selected_tracks) == num_tracks:
        if rank == 0:
            print(f"Using selected_tracks from config ({len(selected_tracks)} tracks)")
        return list(selected_tracks)

    if cache_dir:
        meta_path = Path(cache_dir) / "metadata.json"
        if meta_path.exists():
            try:
                with meta_path.open("r", encoding="utf-8") as handle:
                    meta = json.load(handle)
                track_names = meta.get("track_names", [])
                if len(track_names) == num_tracks:
                    if rank == 0:
                        print(f"Using track names from cache metadata: {meta_path}")
                    return track_names
                # Cache has more tracks than model — try filtering by selected_tracks
                if selected_tracks:
                    available = set(track_names)
                    filtered = [t for t in selected_tracks if t in available]
                    if len(filtered) == num_tracks:
                        if rank == 0:
                            print(
                                f"Filtered cache metadata track names by selected_tracks "
                                f"({len(track_names)} -> {len(filtered)})"
                            )
                        return filtered
                if rank == 0:
                    print(
                        f"Warning: Cache metadata has {len(track_names)} track names "
                        f"but model has {num_tracks} output tracks"
                    )
            except Exception as e:
                if rank == 0:
                    print(f"Could not read track names from cache metadata {meta_path}: {e}")

    try:
        import pyarrow.dataset as ds
        import pyarrow.parquet as pq
    except Exception as e:
        if rank == 0:
            print(f"Could not import pyarrow to read track names: {e}")
        return None

    try:
        if os.path.isdir(data_matrix_path):
            dataset = ds.dataset(data_matrix_path, format="parquet")
            columns = dataset.schema.names
        else:
            parquet_file = pq.ParquetFile(data_matrix_path)
            columns = parquet_file.schema.names
    except Exception as e:
        if rank == 0:
            print(f"Could not read Parquet schema for track names: {e}")
        return None

    track_cols = [c for c in columns if c not in ['chrom', 'start', 'end']]
    if len(track_cols) == num_tracks:
        return track_cols

    # Parquet has more tracks than model — try filtering by selected_tracks
    if selected_tracks:
        available = set(track_cols)
        filtered = [t for t in selected_tracks if t in available]
        if len(filtered) == num_tracks:
            if rank == 0:
                print(
                    f"Filtered parquet track columns by selected_tracks "
                    f"({len(track_cols)} -> {len(filtered)})"
                )
            return filtered

    if rank == 0:
        print(f"Warning: Found {len(track_cols)} track columns but model has {num_tracks} output tracks")
    return None


# =============================================================================
# Configuration
# =============================================================================

def load_config(config_path: str) -> DotDict:
    """Load configuration from YAML file."""
    return load_yaml_config(config_path)


def _resolve_inference_configs(
    config: DotDict,
    inference_name: Optional[str],
    run_all: bool,
    rank: int
) -> List[dict]:
    """Resolve the inference config(s) to use, supporting both dict and list formats.

    Args:
        config: Full configuration object.
        inference_name: Name of specific inference config to use (from --inference_name).
        run_all: If True, return all configs with run_after_training=True.
        rank: Process rank for logging.

    Returns:
        List of inference config dicts to run.

    Raises:
        ValueError: If inference_name is specified but not found in config.
    """
    # Get raw inference config (could be dict or list)
    inference_raw = config.get('inference', {})

    # Handle list format (multiple inference configs)
    if isinstance(inference_raw, list):
        if not inference_raw:
            return []

        if inference_name is not None:
            # Find config by name
            for cfg in inference_raw:
                if cfg.get('name') == inference_name:
                    if rank == 0:
                        print(f"Using inference config: {inference_name}")
                    return [cfg]

            # Not found - list available names
            available_names = [cfg.get('name', f'<unnamed at index {i}>') for i, cfg in enumerate(inference_raw)]
            raise ValueError(
                f"Inference config with name '{inference_name}' not found. "
                f"Available configs: {available_names}"
            )

        if run_all:
            # Return all configs with run_after_training: true (or no flag set, for backward compat)
            configs_to_run = []
            for i, cfg in enumerate(inference_raw):
                # Default to True if run_after_training is not specified
                if cfg.get('run_after_training', True):
                    configs_to_run.append(cfg)

            if rank == 0:
                print(f"Running {len(configs_to_run)} inference config(s) with run_after_training=true")
            return configs_to_run
        else:
            # No name specified and not run_all, use first config (backward compat)
            first_cfg = inference_raw[0]
            if rank == 0:
                cfg_name = first_cfg.get('name', 'unnamed')
                print(f"Using first inference config: {cfg_name} (use --run_all to run all configs)")
            return [first_cfg]

    # Handle dict format (single inference config - legacy)
    elif isinstance(inference_raw, dict):
        if inference_name is not None:
            # Check if the dict has a matching name
            if inference_raw.get('name') == inference_name:
                if rank == 0:
                    print(f"Using inference config: {inference_name}")
                return [inference_raw]
            else:
                raise ValueError(
                    f"Inference config with name '{inference_name}' not found. "
                    f"Config has name: {inference_raw.get('name', '<no name>')}"
                )
        if inference_raw:
            return [inference_raw]

    return []


def _resolve_inference_config(config: DotDict, inference_name: Optional[str], rank: int) -> Optional[dict]:
    """Resolve the inference config to use (backward compatibility wrapper).

    Args:
        config: Full configuration object.
        inference_name: Name of specific inference config to use (from --inference_name).
        rank: Process rank for logging.

    Returns:
        The selected inference config dict, or None if not found.

    Raises:
        ValueError: If inference_name is specified but not found in config.
    """
    configs = _resolve_inference_configs(config, inference_name, run_all=False, rank=rank)
    return configs[0] if configs else None


def _resolve_checkpoint_path(
    checkpoint_val: str,
    checkpoint_dir: str,
    rank: int
) -> str:
    """Resolve a checkpoint value to an actual path.

    Args:
        checkpoint_val: Checkpoint value from config (e.g., 'best_model', 'last',
                        'best_model_val_iou', or a direct path).
        checkpoint_dir: Parent checkpoint directory.
        rank: Process rank for logging.

    Returns:
        Resolved checkpoint path.

    Raises:
        ValueError: If checkpoint cannot be found.
    """
    # Check for metric-based checkpoints (e.g., best_model_val_iou)
    if checkpoint_val.startswith('best_model_val_'):
        metric_checkpoint = find_checkpoint_by_metric(checkpoint_dir, checkpoint_val.replace('best_model_', ''), rank)
        if metric_checkpoint is not None:
            return metric_checkpoint
        # Fall through to try as direct path

    if checkpoint_val == 'last':
        # Find the most recently modified checkpoint (step_*, epoch_*, interrupted_epoch_*)
        # This matches the logic in training_utils.find_latest_checkpoint
        checkpoint_dir_path = Path(checkpoint_dir)
        all_checkpoints = []

        for pattern in ["step_*", "epoch_*", "interrupted_epoch_*"]:
            for p in checkpoint_dir_path.glob(pattern):
                if is_checkpoint_complete(p):
                    try:
                        mtime = p.stat().st_mtime
                        all_checkpoints.append((mtime, p))
                    except OSError:
                        pass

        if all_checkpoints:
            all_checkpoints.sort(reverse=True)
            latest_checkpoint = all_checkpoints[0][1]
            if rank == 0:
                print(f"Resolved 'last' checkpoint to: {latest_checkpoint}")
            return str(latest_checkpoint)

        # Fallback: try explicit "last" file/directory names
        for ext in ['.pt', '', '.pth']:
            candidate = os.path.join(checkpoint_dir, f'last{ext}')
            if os.path.exists(candidate):
                return candidate

        raise ValueError(
            f"Checkpoint 'last' not found in {checkpoint_dir}. "
            f"No step_*, epoch_*, or interrupted_epoch_* checkpoints found. "
            f"Please specify --checkpoint or ensure checkpoint exists."
        )

    if checkpoint_val == 'best_model':
        # Resolve relative to checkpoint_dir
        # Try .pt file first
        checkpoint_path = os.path.join(checkpoint_dir, f'{checkpoint_val}.pt')
        if os.path.exists(checkpoint_path):
            return checkpoint_path

        # Try directory (for accelerate checkpoints)
        checkpoint_dir_path = os.path.join(checkpoint_dir, checkpoint_val)
        if os.path.exists(checkpoint_dir_path):
            return checkpoint_dir_path

        # Try .pth extension
        checkpoint_path_pth = os.path.join(checkpoint_dir, f'{checkpoint_val}.pth')
        if os.path.exists(checkpoint_path_pth):
            return checkpoint_path_pth

        raise ValueError(
            f"Checkpoint '{checkpoint_val}' not found in {checkpoint_dir}. "
            f"Please specify --checkpoint or ensure checkpoint exists."
        )

    # Direct path provided
    if os.path.exists(checkpoint_val):
        return checkpoint_val

    # Try relative to checkpoint_dir
    relative_path = os.path.join(checkpoint_dir, checkpoint_val)
    if os.path.exists(relative_path):
        return relative_path

    raise ValueError(
        f"Checkpoint '{checkpoint_val}' not found. "
        f"Please specify --checkpoint or ensure checkpoint exists."
    )


def _read_wandb_run_id(checkpoint_dir: Path) -> Optional[str]:
    run_id_file = checkpoint_dir / "wandb_run_id.txt"
    if not run_id_file.exists():
        return None
    run_id = run_id_file.read_text().strip()
    return run_id or None


def _find_wandb_run_dir(checkpoint_dir: Path) -> Optional[Path]:
    wandb_root = os.environ.get("WANDB_DIR")
    wandb_root = Path(wandb_root) if wandb_root else Path.cwd() / "wandb"
    if not wandb_root.exists():
        return None

    run_id = _read_wandb_run_id(checkpoint_dir)
    candidates: List[Path] = []

    if run_id:
        candidates.extend(sorted(wandb_root.glob(f"run-*-{run_id}")))
        candidates.extend(sorted(wandb_root.glob(f"run-{run_id}")))
        if not candidates:
            candidates.extend(sorted(wandb_root.glob(f"run-*{run_id}*")))

    if candidates:
        return candidates[0]

    latest_run = wandb_root / "latest-run"
    if latest_run.exists():
        return latest_run

    return None


def _read_wandb_project_from_config(checkpoint_dir: Path, run_dir: Optional[Path]) -> Optional[str]:
    config_candidates = [checkpoint_dir / "config.yaml"]
    if run_dir is not None:
        config_candidates.append(run_dir / "files" / "config.yaml")

    for config_path in config_candidates:
        if not config_path.exists():
            continue
        try:
            with open(config_path, "r") as handle:
                config_data = yaml.safe_load(handle)
        except (OSError, yaml.YAMLError):
            continue

        if not isinstance(config_data, dict):
            continue

        logging_cfg = config_data.get("logging")
        if isinstance(logging_cfg, dict):
            if "wandb_project" in logging_cfg:
                return logging_cfg.get("wandb_project")
            nested_value = logging_cfg.get("value")
            if isinstance(nested_value, dict) and "wandb_project" in nested_value:
                return nested_value.get("wandb_project")

        wandb_project = config_data.get("wandb_project")
        if isinstance(wandb_project, dict) and "value" in wandb_project:
            return wandb_project.get("value")
        if isinstance(wandb_project, str):
            return wandb_project

    return None


def _iter_wandb_history(run_dir: Path) -> Optional[Iterable[dict]]:
    history_candidates = [
        run_dir / "files" / "wandb-history.jsonl",
        run_dir / "files" / "wandb-history.jsonl.gz",
        run_dir / "wandb-history.jsonl",
        run_dir / "wandb-history.jsonl.gz",
    ]

    history_path = next((path for path in history_candidates if path.exists()), None)
    if history_path is None:
        return None

    def _reader() -> Iterable[dict]:
        opener = gzip.open if history_path.suffix == ".gz" else open
        with opener(history_path, "rt") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue

    return _reader()


def _parse_wandb_value(value_json: str) -> object:
    try:
        return json.loads(value_json)
    except json.JSONDecodeError:
        try:
            return float(value_json)
        except ValueError:
            return value_json


def _iter_wandb_history_wandb_file(run_dir: Path, rank: int) -> Optional[Iterable[dict]]:
    wandb_files = sorted(run_dir.glob("run-*.wandb"))
    if not wandb_files:
        return None

    try:
        from wandb.sdk.internal.datastore import DataStore
        from wandb.proto import wandb_internal_pb2
    except ImportError:
        if rank == 0:
            print("Warning: wandb internals not available; cannot parse .wandb file.")
        return None

    wandb_file = wandb_files[0]

    def _reader() -> Iterable[dict]:
        store = DataStore()
        try:
            store.open_for_scan(str(wandb_file))
        except Exception as exc:
            if rank == 0:
                print(f"Warning: Failed to open {wandb_file}: {exc}")
            return

        while True:
            try:
                data = store.scan_data()
            except AssertionError:
                break
            if data is None:
                break
            pb = wandb_internal_pb2.Record()
            pb.ParseFromString(data)
            if pb.WhichOneof("record_type") != "history":
                continue
            record: Dict[str, object] = {}
            for item in pb.history.item:
                if not item.nested_key:
                    continue
                if len(item.nested_key) == 1:
                    key = item.nested_key[0]
                else:
                    key = "/".join(item.nested_key)
                record[key] = _parse_wandb_value(item.value_json)
            if pb.history.step and pb.history.step.num:
                record.setdefault("_step", pb.history.step.num)
            if record:
                yield record

    return _reader()


def _iter_wandb_history_api(
    checkpoint_dir: Path,
    run_dir: Optional[Path],
    normalized_metric: str,
    rank: int,
) -> Optional[Iterable[dict]]:
    run_id = _read_wandb_run_id(checkpoint_dir) or os.environ.get("WANDB_RUN_ID")
    if not run_id:
        if rank == 0:
            print("Warning: W&B run id not found; cannot query W&B API.")
        return None

    wandb_project = os.environ.get("WANDB_PROJECT") or _read_wandb_project_from_config(
        checkpoint_dir, run_dir
    )
    if not wandb_project:
        if rank == 0:
            print("Warning: W&B project not found; cannot query W&B API.")
        return None

    try:
        import wandb
    except ImportError:
        if rank == 0:
            print("Warning: wandb not installed; cannot query W&B API.")
        return None

    api = wandb.Api()
    entity = os.environ.get("WANDB_ENTITY") or os.environ.get("WANDB_USERNAME")
    if not entity:
        entity = getattr(api, "default_entity", None)
    if not entity:
        if rank == 0:
            print("Warning: W&B entity not found; set WANDB_ENTITY to query W&B API.")
        return None

    run_path = f"{entity}/{wandb_project}/{run_id}"
    if rank == 0:
        print(f"Querying W&B API: {run_path} for metric '{normalized_metric}'")
    try:
        run = api.run(run_path)
    except Exception as exc:
        if rank == 0:
            print(f"Warning: Failed to load W&B run {run_path}: {exc}")
        return None

    def _reader() -> Iterable[dict]:
        # Include both _step (normal logging) and backfill_step (backfill script)
        for row in run.scan_history(keys=[normalized_metric, "_step", "backfill_step"]):
            # Normalize step: prefer backfill_step if present, otherwise use _step
            step = row.get("backfill_step") or row.get("_step")
            if step is not None:
                row["_step"] = step
            yield row

    return _reader()


def _build_checkpoint_step_map(checkpoint_dir: Path) -> Dict[int, Path]:
    step_map: Dict[int, Path] = {}
    for checkpoint_path in checkpoint_dir.iterdir():
        if not checkpoint_path.is_dir():
            continue
        match = re.match(r"step_(\d+)$", checkpoint_path.name)
        if match:
            step_map[int(match.group(1))] = checkpoint_path
    return step_map


def _find_checkpoint_from_wandb(
    checkpoint_dir: Path,
    normalized_metric: str,
    rank: int,
) -> Optional[str]:
    run_dir = _find_wandb_run_dir(checkpoint_dir)
    if run_dir is None:
        if rank == 0:
            print("Warning: W&B run directory not found; cannot use W&B metrics.")
        return None

    # Try sources in order: jsonl file, binary .wandb file, then API
    # For each source, check if the metric is found; if not, try next source
    # (backfilled metrics only exist in API, not in local files)
    history_sources = [
        ("wandb-history.jsonl", lambda: _iter_wandb_history(run_dir)),
        (".wandb file", lambda: _iter_wandb_history_wandb_file(run_dir, rank)),
        ("W&B API", lambda: _iter_wandb_history_api(checkpoint_dir, run_dir, normalized_metric, rank)),
    ]

    mode = get_metric_direction(normalized_metric)
    best_value = get_initial_best_value(mode)
    best_step = None

    for source_name, get_iter in history_sources:
        history_iter = get_iter()
        if history_iter is None:
            if rank == 0:
                print(f"Warning: {source_name} not available.")
            continue

        # Scan this source for the metric
        for entry in history_iter:
            if normalized_metric not in entry:
                continue
            value = entry.get(normalized_metric)
            if value is None:
                continue
            step = entry.get("_step")
            if step is None:
                continue
            if is_better(value, best_value, mode):
                best_value = value
                best_step = int(step)

        # If we found the metric in this source, stop searching
        if best_step is not None:
            if rank == 0:
                print(f"Found {normalized_metric} in {source_name}.")
            break
        elif rank == 0:
            print(f"Metric {normalized_metric} not found in {source_name}, trying next source...")

    if best_step is None:
        if rank == 0:
            print(f"Warning: Metric {normalized_metric} not found in any W&B source.")
        return None

    step_map = _build_checkpoint_step_map(checkpoint_dir)
    if not step_map:
        if rank == 0:
            print("Warning: No step_* checkpoints found to map W&B metrics.")
        return None

    if best_step in step_map:
        if rank == 0:
            print(
                f"Using checkpoint from W&B {normalized_metric}={best_value:.4f}: {step_map[best_step]}"
            )
        return str(step_map[best_step])

    nearest_step = min(step_map.keys(), key=lambda step_val: abs(step_val - best_step))
    if rank == 0:
        print(
            f"Warning: No checkpoint for W&B step {best_step}; "
            f"using nearest step {nearest_step} for {normalized_metric}={best_value:.4f}."
        )
    return str(step_map[nearest_step])


def find_checkpoint_by_metric(checkpoint_dir: str, metric_name: str, rank: int = 0) -> Optional[str]:
    """Find the best checkpoint for a given metric.

    First tries to find a dedicated best_model_{metric}/ directory.
    If not found, scans all checkpoints' metrics.json files to find the best.

    Args:
        checkpoint_dir: Parent checkpoint directory.
        metric_name: Name of the metric (with or without 'val_' prefix).
        rank: Process rank for logging (only rank 0 prints).

    Returns:
        Path to the best checkpoint for the metric, or None if not found.
    """
    checkpoint_dir = Path(checkpoint_dir)

    # Normalize metric name
    normalized_metric = metric_name if metric_name.startswith('val_') else f'val_{metric_name}'

    # Try dedicated metric checkpoint first
    metric_checkpoint = checkpoint_dir / f'best_model_{normalized_metric}'
    if metric_checkpoint.exists():
        if rank == 0:
            print(f"Using checkpoint for best {normalized_metric}: {metric_checkpoint}")
        return str(metric_checkpoint)

    # Prefer W&B history when available (covers all step_* checkpoints)
    if rank == 0:
        print(f"No best_model_{normalized_metric}/ found, trying W&B history...")
    wandb_checkpoint = _find_checkpoint_from_wandb(checkpoint_dir, normalized_metric, rank)
    if wandb_checkpoint is not None:
        return wandb_checkpoint

    # Fallback: scan all checkpoints for metrics.json
    if rank == 0:
        print(f"W&B history not available for {normalized_metric}, scanning checkpoints for metrics.json...")

    mode = get_metric_direction(normalized_metric)
    best_value = get_initial_best_value(mode)
    best_checkpoint = None

    for checkpoint_path in checkpoint_dir.iterdir():
        if not checkpoint_path.is_dir():
            continue

        metrics_data = load_metrics_json(checkpoint_path)
        if metrics_data is None:
            continue

        metrics = metrics_data.get('metrics', {})
        if normalized_metric in metrics:
            value = metrics[normalized_metric]
            if is_better(value, best_value, mode):
                best_value = value
                best_checkpoint = str(checkpoint_path)

    if best_checkpoint is not None:
        if rank == 0:
            print(f"Found best checkpoint for {normalized_metric}={best_value:.4f}: {best_checkpoint}")
        return best_checkpoint

    # Not found
    if rank == 0:
        print(f"Warning: No checkpoint with {normalized_metric} found in metrics.json or W&B logs")
    return None


# Import unified inference engine (includes AdaptiveBatchManager)
from src.inference_engine import InferenceEngine, AdaptiveBatchManager


# =============================================================================
# Inference Dataset
# =============================================================================

class GenomeInferenceDataset(Dataset):
    """
    Dataset for genome-wide inference.

    Generates windows that tile the genome/regions such that each bin
    is predicted exactly once (using stride = pred_len).
    """

    def __init__(
        self,
        regions: List[Tuple[str, int, int]],  # List of (chrom, start, end) to predict
        genome_fasta: str,
        chrom_sizes: Dict[str, int],
        seq_len: int,
        pred_len: int,
        bin_size: int = 32,
        stride: Optional[int] = None,
        splice_perts: Optional[List[SplicePert]] = None,
        seq_perts: Optional[List[SeqPert]] = None,
    ):
        """
        Initialize the inference dataset.

        Args:
            regions: List of (chrom, start, end) tuples defining prediction regions
            genome_fasta: Path to reference genome FASTA
            chrom_sizes: Dictionary of chromosome sizes
            seq_len: Total input sequence length (e.g., 65536)
            pred_len: Prediction window length (e.g., 24576)
            bin_size: Bin size in bp (default 32)
            stride: Slide stride; defaults to pred_len for non-overlapping windows
            splice_perts: Optional list of splice perturbations to apply
            seq_perts: Optional list of sequence perturbations to apply
        """
        self.genome_fasta = genome_fasta
        self.chrom_sizes = chrom_sizes
        self.seq_len = seq_len
        self.pred_len = pred_len
        self.bin_size = bin_size
        self.stride = stride or pred_len
        self.splice_perts = splice_perts or []
        self.seq_perts = seq_perts or []

        # Calculate context/crop on each side
        self.crop = (seq_len - pred_len) // 2

        # Generate all windows (sliding by pred_len for non-overlapping predictions)
        self.windows = self._generate_windows(regions)

        # Thread-local genome handle
        import threading
        self._local = threading.local()

    @property
    def genome(self):
        """Thread-safe lazy loading of genome file handle."""
        if not hasattr(self._local, 'genome') or self._local.genome is None:
            self._local.genome = pysam.FastaFile(self.genome_fasta)
        return self._local.genome

    def _generate_windows(
        self,
        regions: List[Tuple[str, int, int]]
    ) -> List[Tuple[str, int, int, int, int]]:
        """
        Generate windows that tile regions with stride = pred_len.

        Returns list of (chrom, seq_start, seq_end, pred_start, pred_end)
        """
        windows = []

        for chrom, region_start, region_end in regions:
            # Align to bin boundaries
            region_start = (region_start // self.bin_size) * self.bin_size
            region_end = (region_end // self.bin_size) * self.bin_size

            # Slide with configurable stride (vectorized offsets)
            offsets = np.arange(region_start, region_end, self.stride, dtype=np.int64)
            pred_starts = offsets
            pred_ends = np.minimum(pred_starts + self.pred_len, region_end)

            # Drop windows that would be empty
            valid = pred_starts < pred_ends
            pred_starts = pred_starts[valid]
            pred_ends = pred_ends[valid]

            if pred_starts.size == 0:
                continue

            pred_centers = (pred_starts + pred_ends) // 2
            seq_starts = pred_centers - self.seq_len // 2
            seq_ends = pred_centers + self.seq_len // 2

            windows.extend(
                zip(
                    np.full(pred_starts.size, chrom),
                    seq_starts,
                    seq_ends,
                    pred_starts,
                    pred_ends,
                )
            )

        return windows

    def __len__(self) -> int:
        return len(self.windows)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, Tuple[str, int, int]]:
        """
        Get a single window.

        Returns:
            sequence: One-hot encoded sequence (seq_len, 4)
            metadata: (chrom, pred_start, pred_end)
        """
        chrom, seq_start, seq_end, pred_start, pred_end = self.windows[idx]

        # Fetch and encode sequence
        sequence = self._fetch_sequence(chrom, seq_start, seq_end)
        seq_tensor = torch.from_numpy(sequence)

        return seq_tensor, (chrom, pred_start, pred_end)

    def _fetch_sequence(self, chrom: str, start: int, end: int) -> np.ndarray:
        """Fetch and one-hot encode a sequence with padding for boundaries, applying perturbations."""
        chrom_size = self.chrom_sizes.get(chrom, float('inf'))

        # Handle left boundary
        left_pad_len = max(0, -start)
        fetch_start = max(0, start)
        fetch_end = min(int(chrom_size), end)

        # Fetch base sequence as string
        seq_str = ""
        if fetch_end > fetch_start:
            try:
                seq_str = self.genome.fetch(chrom, fetch_start, fetch_end).upper()
            except Exception as e:
                print(f"Warning: Could not fetch {chrom}:{fetch_start}-{fetch_end}: {e}")

        # Convert to NumPy uint8 array for fast manipulation
        seq_arr = np.frombuffer(seq_str.encode('ascii'), dtype=np.uint8).copy() if seq_str else np.array([], dtype=np.uint8)

        # Pad left if needed
        if left_pad_len > 0:
            left_pad = np.full(left_pad_len, ord('N'), dtype=np.uint8)
            seq_arr = np.concatenate([left_pad, seq_arr])

        # Pad right if needed
        expected_len = self.seq_len
        if len(seq_arr) < expected_len:
            right_pad_len = expected_len - len(seq_arr)
            right_pad = np.full(right_pad_len, ord('N'), dtype=np.uint8)
            seq_arr = np.concatenate([seq_arr, right_pad])
        elif len(seq_arr) > expected_len:
            seq_arr = seq_arr[:expected_len]

        # Apply perturbations if any
        if self.splice_perts or self.seq_perts:
            seq_arr = self._apply_perturbations(seq_arr, chrom, start, end)

        # One-hot encode using the lookup table
        from src.data_loaders import ONE_HOT_LOOKUP
        result = ONE_HOT_LOOKUP[seq_arr]

        return result

    def _apply_perturbations(
        self,
        seq_arr: np.ndarray,
        chrom: str,
        seq_start: int,
        seq_end: int
    ) -> np.ndarray:
        """
        Apply splice and sequence perturbations to a sequence array.

        Args:
            seq_arr: Sequence as uint8 NumPy array (ASCII codes)
            chrom: Chromosome of the window
            seq_start: Genomic start position of the window
            seq_end: Genomic end position of the window

        Returns:
            Modified sequence array
        """
        # Apply splice perturbations first
        for splice in self.splice_perts:
            # Check if perturbation overlaps with this window
            if splice.chrom != chrom:
                continue
            if splice.end <= seq_start or splice.start >= seq_end:
                continue

            # Get donor sequence
            if splice.seq is not None:
                donor_str = splice.seq.upper()
            else:
                try:
                    donor_str = self.genome.fetch(splice.chrom_src, splice.start_src, splice.end_src).upper()
                except Exception as e:
                    print(f"Warning: Could not fetch donor sequence for splice {splice.id}: {e}")
                    continue

            donor_arr = np.frombuffer(donor_str.encode('ascii'), dtype=np.uint8)

            # Convert genomic coords to sequence-relative coords
            rel_start = splice.start - seq_start
            rel_end = splice.end - seq_start

            # Calculate donor offset (when splice starts before window)
            donor_offset = max(0, seq_start - splice.start)

            # Clip to window boundaries
            rel_start = max(0, rel_start)
            rel_end = min(len(seq_arr), rel_end)

            if rel_start >= rel_end:
                continue

            # Apply splice (replace target with donor, accounting for offset)
            target_len = rel_end - rel_start
            seq_arr[rel_start:rel_end] = donor_arr[donor_offset:donor_offset + target_len]

        # Apply sequence perturbations second
        for seq_pert in self.seq_perts:
            # Check if perturbation overlaps with this window
            if seq_pert.chrom != chrom:
                continue
            if seq_pert.end <= seq_start or seq_pert.start >= seq_end:
                continue

            # Get donor sequence
            if seq_pert.seq is not None:
                donor_str = seq_pert.seq.upper()
            else:
                try:
                    donor_str = self.genome.fetch(seq_pert.chrom_src, seq_pert.start_src, seq_pert.end_src).upper()
                except Exception as e:
                    print(f"Warning: Could not fetch donor sequence for seq_pert {seq_pert.id}: {e}")
                    continue

            donor_arr = np.frombuffer(donor_str.encode('ascii'), dtype=np.uint8)

            # Convert genomic coords to sequence-relative coords
            rel_start = seq_pert.start - seq_start
            rel_end = seq_pert.end - seq_start

            # Calculate donor offset (when pert starts before window)
            donor_offset = max(0, seq_start - seq_pert.start)

            # Clip to window boundaries
            rel_start = max(0, rel_start)
            rel_end = min(len(seq_arr), rel_end)

            if rel_start >= rel_end:
                continue

            # Apply perturbation (accounting for offset)
            target_len = rel_end - rel_start
            seq_arr[rel_start:rel_end] = donor_arr[donor_offset:donor_offset + target_len]

        return seq_arr


def collate_with_metadata(batch):
    """Custom collate function that handles metadata."""
    sequences = torch.stack([item[0] for item in batch])
    metadata = [item[1] for item in batch]
    return sequences, metadata


# =============================================================================
# Region Generation
# =============================================================================

def generate_regions_genome(
    chrom_sizes: Dict[str, int],
    exclude_chroms: Optional[List[str]] = None
) -> List[Tuple[str, int, int]]:
    """Generate regions covering the entire genome."""
    if exclude_chroms is None:
        exclude_chroms = []

    # Filter to standard chromosomes and exclude specified ones
    regions = []
    for chrom, size in sorted(chrom_sizes.items()):
        # Skip non-standard chromosomes (scaffolds, etc.)
        if '_' in chrom or 'random' in chrom.lower() or 'un' in chrom.lower():
            continue
        if chrom in exclude_chroms:
            continue
        regions.append((chrom, 0, size))

    return regions


def generate_regions_chromosomes(
    chrom_sizes: Dict[str, int],
    chromosomes: List[str]
) -> List[Tuple[str, int, int]]:
    """Generate regions for specified chromosomes."""
    regions = []
    for chrom in chromosomes:
        if chrom in chrom_sizes:
            regions.append((chrom, 0, chrom_sizes[chrom]))
        else:
            print(f"Warning: Chromosome {chrom} not found in chrom_sizes")
    return regions


def generate_regions_bed(bed_file: str) -> List[Tuple[str, int, int]]:
    """Load regions from a BED file."""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            regions.append((chrom, start, end))
    return regions


# =============================================================================
# Model Loading
# =============================================================================

def load_model(
    checkpoint_path: str,
    config: DotDict,
    device: torch.device
) -> nn.Module:
    """
    Load a fine-tuned Borzoi model from checkpoint.

    Supports multiple checkpoint formats:
    - Accelerate checkpoint directory (with model.safetensors)
    - Standard PyTorch checkpoint (.pt/.pth)
    - State dict only
    """
    from borzoi_pytorch import Borzoi
    from borzoi_pytorch.config_borzoi import BorzoiConfig

    # Mirror training script mapping and head handling
    model_mapping = {
        'borzoi_mouse_rep0': 'johahi/borzoi-replicate-0',
        'borzoi_human_rep0': 'johahi/borzoi-replicate-0',
        'flashzoi': 'johahi/flashzoi-replicate-0',
    }

    model_name = config.model.name.lower()
    hf_model_id = None
    for key, model_id in model_mapping.items():
        if key in model_name:
            hf_model_id = model_id
            break
    if hf_model_id is None:
        hf_model_id = config.model.name

    print(f"Loading model from {hf_model_id}...")

    # Determine species flags with priority: explicit config > path detection
    # Priority 1: Check explicit species field in config.model
    if hasattr(config.model, 'species') and config.model.species:
        species_explicit = str(config.model.species).lower()
        is_mouse_data = species_explicit == 'mouse'
        is_human_data = species_explicit == 'human'
        print(f"Species explicitly set in config.model.species: {config.model.species}")
    # Priority 2: Check explicit species field in config.data
    elif hasattr(config.data, 'species') and config.data.species:
        species_explicit = str(config.data.species).lower()
        is_mouse_data = species_explicit == 'mouse'
        is_human_data = species_explicit == 'human'
        print(f"Species explicitly set in config.data.species: {config.data.species}")
    # Priority 3: Fall back to path detection (legacy behavior)
    else:
        genome_path = str(config.data.genome_fasta).lower() if hasattr(config.data, 'genome_fasta') and config.data.genome_fasta else ""
        is_mouse_data = "mm10" in genome_path or "mouse" in genome_path
        is_human_data = not is_mouse_data
        if any(tag in genome_path for tag in ["hg19", "hg38", "human"]):
            is_human_data = True
        print(f"Species auto-detected from genome path: {'mouse' if is_mouse_data else 'human'}")

    # Calculate target bins
    if config.model.pred_len % config.model.bin_size != 0:
        raise ValueError(f"pred_len ({config.model.pred_len}) must be divisible by bin_size ({config.model.bin_size})")
    target_bins = config.model.pred_len // config.model.bin_size

    # Load configuration
    try:
        borzoi_config = BorzoiConfig.from_pretrained(hf_model_id)
    except Exception:
        borzoi_config = BorzoiConfig()

    # Enable mouse head for mouse data
    if is_mouse_data:
        print("Mouse data detected - enabling mouse_head")
        borzoi_config.enable_mouse_head = True

    borzoi_config.return_center_bins_only = True
    borzoi_config.bins_to_return = target_bins

    # Initialize base model with config
    # Use ignore_mismatched_sizes=True to handle mouse_head being different size
    model = Borzoi.from_pretrained(hf_model_id, config=borzoi_config, ignore_mismatched_sizes=True)

    # Ensure the appropriate head exists (like training)
    if is_human_data and hasattr(model, 'enable_human_head') and not hasattr(model, 'human_head'):
        try:
            model.enable_human_head()
        except Exception as e:
            print(f"Warning: enable_human_head() failed: {e}")
    if is_mouse_data and hasattr(model, 'enable_mouse_head') and not hasattr(model, 'mouse_head'):
        try:
            model.enable_mouse_head()
        except Exception as e:
            print(f"Warning: enable_mouse_head() failed: {e}")

    # First, load checkpoint to determine number of output tasks
    print(f"Loading checkpoint from {checkpoint_path}...")

    checkpoint_path = Path(checkpoint_path)

    # Case 1: Accelerate checkpoint directory (contains model.safetensors)
    if checkpoint_path.is_dir():
        safetensors_path = checkpoint_path / "model.safetensors"
        pytorch_path = checkpoint_path / "pytorch_model.bin"

        if safetensors_path.exists():
            print(f"  Loading from safetensors: {safetensors_path}")
            from safetensors.torch import load_file
            state_dict = load_file(str(safetensors_path))
        elif pytorch_path.exists():
            print(f"  Loading from pytorch_model.bin: {pytorch_path}")
            state_dict = torch.load(pytorch_path, map_location='cpu', weights_only=False)
        else:
            raise FileNotFoundError(
                f"No model.safetensors or pytorch_model.bin found in {checkpoint_path}"
            )
    # Case 2: Direct safetensors file
    elif str(checkpoint_path).endswith('.safetensors'):
        print(f"  Loading from safetensors file")
        from safetensors.torch import load_file
        state_dict = load_file(str(checkpoint_path))
    # Case 3: Standard PyTorch checkpoint
    else:
        checkpoint = torch.load(str(checkpoint_path), map_location='cpu', weights_only=False)

        # Handle different checkpoint formats
        if isinstance(checkpoint, dict):
            if 'model_state_dict' in checkpoint:
                state_dict = checkpoint['model_state_dict']
            elif 'state_dict' in checkpoint:
                state_dict = checkpoint['state_dict']
            else:
                state_dict = checkpoint
        else:
            state_dict = checkpoint

    # Remove 'module.' prefix if present (from DataParallel/DDP)
    if any(k.startswith('module.') for k in state_dict.keys()):
        state_dict = {k.replace('module.', ''): v for k, v in state_dict.items()}

    # Determine which head and number of output tasks from checkpoint
    # Track mouse/human heads independently because some checkpoints include both.
    num_tasks_mouse = None
    num_tasks_human = None
    num_tasks_generic = None
    checkpoint_has_mouse_head = False
    checkpoint_has_human_head = False

    for key in state_dict.keys():
        if 'mouse_head' in key and 'weight' in key and num_tasks_mouse is None:
            num_tasks_mouse = state_dict[key].shape[0]
            checkpoint_has_mouse_head = True
            print(f"  Detected mouse_head with {num_tasks_mouse} output tasks from checkpoint")
        elif 'human_head' in key and 'weight' in key and num_tasks_human is None:
            num_tasks_human = state_dict[key].shape[0]
            checkpoint_has_human_head = True
            print(f"  Detected human_head with {num_tasks_human} output tasks from checkpoint")
        elif 'head' in key and 'weight' in key and num_tasks_generic is None:
            num_tasks_generic = state_dict[key].shape[0]
            print(f"  Detected {num_tasks_generic} output tasks from checkpoint")

    def _replace_head(module: nn.Module, attr: str, num_tasks: int) -> None:
        current_head = getattr(module, attr)
        if hasattr(current_head, 'in_features'):
            in_features = current_head.in_features
        elif hasattr(current_head, 'in_channels'):
            in_features = current_head.in_channels
        elif hasattr(current_head, 'weight'):
            in_features = current_head.weight.shape[1]
        else:
            in_features = 1920
        new_head = nn.Conv1d(in_features, num_tasks, kernel_size=1)
        setattr(module, attr, new_head)
        print(f"  Replaced {attr}: {in_features} -> {num_tasks}")

    # Replace heads with correct size BEFORE loading state dict
    # Choose head based on checkpoint contents, falling back to generic head size.
    if checkpoint_has_mouse_head and hasattr(model, 'mouse_head') and num_tasks_mouse is not None:
        _replace_head(model, 'mouse_head', num_tasks_mouse)
    elif is_mouse_data and hasattr(model, 'mouse_head') and num_tasks_generic is not None:
        _replace_head(model, 'mouse_head', num_tasks_generic)

    if checkpoint_has_human_head and hasattr(model, 'human_head') and num_tasks_human is not None:
        _replace_head(model, 'human_head', num_tasks_human)
    elif (not is_mouse_data) and hasattr(model, 'human_head') and num_tasks_generic is not None:
        _replace_head(model, 'human_head', num_tasks_generic)
    elif hasattr(model, 'head') and num_tasks_generic is not None:
        _replace_head(model, 'head', num_tasks_generic)

    # Load state dict
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    if missing:
        print(f"  Warning: Missing keys: {missing[:5]}{'...' if len(missing) > 5 else ''}")
    if unexpected:
        print(f"  Warning: Unexpected keys: {unexpected[:5]}{'...' if len(unexpected) > 5 else ''}")

    model = model.to(device)
    model.eval()

    # Patch MaxPool1d to avoid INT32 overflow with large batch × channel × length
    from src.model import _patch_maxpool_for_large_tensors

    class _PrintAdapter:
        """Minimal adapter so _patch_maxpool_for_large_tensors can use print()."""
        @staticmethod
        def print(msg: str) -> None:
            print(msg)

    _patch_maxpool_for_large_tensors(model, _PrintAdapter())

    print(f"  Model loaded successfully!")

    # Return both model and is_human flag
    # Prefer mouse head if mouse data and available, otherwise use human head
    is_human = not is_mouse_data
    if is_mouse_data and hasattr(model, 'mouse_head'):
        is_human = False
    elif hasattr(model, 'human_head'):
        is_human = True
    elif hasattr(model, 'head'):
        # Generic head; follow human/default path
        is_human = True

    return model, is_human


# =============================================================================
# Inference Functions
# =============================================================================

def run_inference(
    model: nn.Module,
    dataloader: DataLoader,
    device: torch.device,
    use_rc_average: bool = True,
    mixed_precision: bool = True,
    rank: int = 0,
    world_size: int = 1,
    is_human: bool = True,
    inference_engine: Optional[InferenceEngine] = None,
    use_compile: bool = False,
    compile_mode: str = "reduce-overhead",
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Run inference on the dataloader.

    Args:
        model: Trained model
        dataloader: DataLoader with inference data
        device: Device to run on
        use_rc_average: Whether to average forward and reverse complement predictions
        mixed_precision: Whether to use mixed precision
        rank: Current process rank (for distributed)
        world_size: Total number of processes
        is_human: Whether to use human_head (True) or mouse_head (False)
        inference_engine: Optional pre-configured InferenceEngine. If not provided,
            one will be created internally.
        use_compile: Whether to use torch.compile() for acceleration
        compile_mode: torch.compile mode

    Returns:
        Dictionary mapping chrom -> {positions: array, predictions: array}
    """
    model.eval()

    if inference_engine is None:
        inference_engine = InferenceEngine(
            model=model,
            device=device,
            use_rc_average=use_rc_average,
            mixed_precision=mixed_precision,
            is_human=is_human,
            use_compile=use_compile,
            compile_mode=compile_mode,
        )
        if rank == 0 and inference_engine._compiled:
            print(f"  Using torch.compile() with mode='{compile_mode}'")

    # Collect predictions per chromosome
    results = {}  # chrom -> list of (start, end, predictions)

    # Progress bar only on rank 0
    iterator = tqdm(dataloader, desc=f"Inference (GPU {rank})", disable=(rank != 0))

    for batch_seqs, batch_metadata in iterator:
        try:
            if isinstance(batch_seqs, torch.Tensor):
                batch_seqs_np = batch_seqs.numpy()
            else:
                batch_seqs_np = batch_seqs

            predictions = inference_engine.predict_batch(batch_seqs_np)

            # Store predictions with their genomic coordinates
            for i, (chrom, pred_start, pred_end) in enumerate(batch_metadata):
                if chrom not in results:
                    results[chrom] = []

                # Calculate actual number of bins for this window
                n_bins = (pred_end - pred_start) // dataloader.dataset.bin_size

                results[chrom].append({
                    'start': pred_start,
                    'end': pred_end,
                    'predictions': predictions[i, :n_bins, :].copy()
                })

        except RuntimeError as e:
            if 'out of memory' in str(e).lower():
                # Surface OOM so the outer handler can manage retries or logging
                # In distributed mode, each rank handles OOM independently since
                # they process different amounts of data
                raise RuntimeError("CUDA out of memory")
            else:
                # Propagate other errors to preserve the original stack trace
                raise

    return results


def merge_predictions(
    results: Dict[str, List],
    bin_size: int
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Merge overlapping/adjacent predictions per chromosome.

    Returns:
        Dictionary mapping chrom -> {'starts': array, 'predictions': array}
    """
    merged = {}

    for chrom, windows in results.items():
        # Sort by start position
        windows = sorted(windows, key=lambda x: x['start'])

        # Collect all predictions
        starts = []
        all_preds = []

        for w in windows:
            n_bins = len(w['predictions'])
            bin_starts = np.arange(w['start'], w['end'], bin_size)[:n_bins]
            starts.extend(bin_starts)
            all_preds.append(w['predictions'])

        merged[chrom] = {
            'starts': np.array(starts),
            'predictions': np.vstack(all_preds)
        }

    return merged


# =============================================================================
# Output Writers
# =============================================================================

def _sanitize_filename(name: str) -> str:
    """Sanitize a string to be safe for use as a filename.

    Removes or replaces characters that could cause path traversal or
    other filesystem issues.

    Args:
        name: The input string to sanitize.

    Returns:
        A sanitized string safe for use as a filename.
    """
    # Remove path traversal sequences first
    name = name.replace('..', '_')
    # Replace any character that's not alphanumeric, dash, underscore, or dot
    safe_name = re.sub(r'[^\w\-.]', '_', name)
    # Collapse multiple underscores
    safe_name = re.sub(r'_+', '_', safe_name)
    # Remove leading/trailing underscores and dots (prevent hidden files)
    safe_name = safe_name.strip('_.')
    # Ensure non-empty
    if not safe_name:
        safe_name = 'unnamed'
    return safe_name


def _wait_for_paths(
    paths: List[str],
    poll_interval: float = 5.0,
    label: str = "files",
    min_mtime: Optional[float] = None,
) -> None:
    """Block until all paths exist (optionally requiring a minimum mtime)."""
    def _is_ready(path: str) -> bool:
        if not os.path.exists(path):
            return False
        if min_mtime is None:
            return True
        try:
            return os.path.getmtime(path) >= min_mtime
        except OSError:
            return False

    missing = [p for p in paths if not _is_ready(p)]
    while missing:
        time.sleep(poll_interval)
        missing = [p for p in paths if not _is_ready(p)]


def _get_sync_root(config_path: str, inference_idx: int, sync_dir: Optional[str] = None) -> str:
    """Return a shared sync root path for the current distributed run."""
    master_port = os.environ.get("MASTER_PORT")
    if master_port is not None and not master_port.isdigit():
        master_port = None
    run_id = (
        os.environ.get("TORCHELASTIC_RUN_ID")
        or master_port
        or os.environ.get("SLURM_JOB_ID")
        or os.environ.get("SLURM_STEP_ID")
        or os.environ.get("LSB_JOBID")
        or os.environ.get("PBS_JOBID")
        or os.environ.get("JOB_ID")
        or "default"
    )
    base_dir = (
        sync_dir
        or os.environ.get("BORZOI_SYNC_DIR")
        or os.path.join(os.path.dirname(os.path.abspath(config_path)), ".dist_sync")
    )
    os.makedirs(base_dir, exist_ok=True)
    tag = _sanitize_filename(f"{Path(config_path).stem}_{run_id}_{inference_idx}")
    return os.path.join(base_dir, f"sync_{tag}")


def _write_single_bigwig_worker(args):
    """Worker function for parallel BigWig writing.

    Must be at module level for pickling in multiprocessing.
    """
    import pyBigWig

    track_idx, track_name, output_dir, prefix, chroms_header, predictions, bin_size = args

    safe_name = _sanitize_filename(track_name)
    output_path = os.path.join(output_dir, f"{prefix}{safe_name}.bw")

    bw = pyBigWig.open(output_path, "w")
    bw.addHeader(chroms_header)

    # Iterate in the SAME order as header
    for chrom, chrom_size in chroms_header:
        if chrom not in predictions:
            continue

        data = predictions[chrom]
        starts = data['starts'].astype(np.int64)
        ends = (starts + bin_size).astype(np.int64)
        values = data['predictions'][:, track_idx].astype(np.float64)

        # Filter out bins that exceed chromosome boundaries
        valid_bounds = starts < chrom_size
        starts = starts[valid_bounds]
        ends = ends[valid_bounds]
        values = values[valid_bounds]

        # Clip ends to chromosome boundary
        ends = np.minimum(ends, chrom_size)

        # Sort by position (required for BigWig)
        sort_idx = np.argsort(starts)
        starts = starts[sort_idx]
        ends = ends[sort_idx]
        values = values[sort_idx]

        # Remove duplicates (keep first occurrence)
        _, unique_idx = np.unique(starts, return_index=True)
        unique_idx = np.sort(unique_idx)  # Maintain order
        starts = starts[unique_idx]
        ends = ends[unique_idx]
        values = values[unique_idx]

        # Filter out invalid values
        valid = ~np.isnan(values) & ~np.isinf(values)
        if valid.sum() > 0:
            bw.addEntries(
                [chrom] * valid.sum(),
                starts[valid].tolist(),
                ends=ends[valid].tolist(),
                values=values[valid].tolist()
            )

    bw.close()
    return track_idx


def write_bigwig(
    predictions: Dict[str, Dict],
    output_dir: str,
    chrom_sizes: Dict[str, int],
    track_names: List[str],
    bin_size: int = 32,
    prefix: str = "",
    n_jobs: int = 4
):
    """Write predictions to BigWig files (one per track).

    Args:
        n_jobs: Number of parallel jobs for writing BigWig files (default: 4)
    """
    from multiprocessing import Pool

    os.makedirs(output_dir, exist_ok=True)

    # Sort chromosomes naturally (chr1, chr2, ..., chr10, etc.) for BigWig header
    def chrom_sort_key(chrom):
        # Extract number from chromosome name for natural sorting
        import re
        match = re.match(r'chr(\d+|[XYM])', chrom)
        if match:
            val = match.group(1)
            if val == 'X':
                return (100, chrom)
            elif val == 'Y':
                return (101, chrom)
            elif val == 'M':
                return (102, chrom)
            else:
                return (int(val), chrom)
        return (999, chrom)

    # Get chromosomes that have predictions, sorted properly
    pred_chroms = sorted([c for c in predictions.keys()], key=chrom_sort_key)

    # Build header with only chromosomes that have predictions, in sorted order
    chroms_header = [(c, chrom_sizes[c]) for c in pred_chroms if c in chrom_sizes]

    # Prepare arguments for parallel processing
    worker_args = [
        (i, track_names[i], output_dir, prefix, chroms_header, predictions, bin_size)
        for i in range(len(track_names))
    ]

    # Write BigWig files in parallel
    print(f"Writing {len(track_names)} BigWig files using {n_jobs} parallel jobs...")

    with Pool(n_jobs) as pool:
        list(tqdm(
            pool.imap(_write_single_bigwig_worker, worker_args),
            total=len(track_names),
            desc="Writing BigWig files"
        ))

    print(f"Wrote {len(track_names)} BigWig files to {output_dir}")


def _write_single_misha_track(args):
    """Worker function for parallel Misha track writing via pymisha.

    Accesses shared data (predictions, intervals_df, chrom_order) via
    module-level _misha_shared dict, inherited through fork COW — avoids
    pickling large arrays through the Pool queue.
    """
    (track_idx, full_track_name, misha_binsize, groot, description, overwrite) = args

    import pymisha as pm
    pm.gdb_init(groot)

    shared = _misha_shared
    predictions = shared['predictions']
    intervals_df = shared['intervals_df']
    chrom_order = shared['chrom_order']

    # Extract values per-chrom and concatenate (one track at a time for memory)
    value_chunks = []
    for chrom in chrom_order:
        if chrom not in predictions:
            continue
        data = predictions[chrom]
        starts = data['starts'].astype(np.int64)
        vals = data['predictions'][:, track_idx].astype(np.float64)

        sort_idx = np.argsort(starts)
        starts = starts[sort_idx]
        vals = vals[sort_idx]

        _, unique_idx = np.unique(starts, return_index=True)
        unique_idx = np.sort(unique_idx)
        vals = vals[unique_idx]

        # Replace Inf with NaN (Misha uses NaN as default missing value)
        inf_mask = np.isinf(vals)
        if inf_mask.any():
            vals = vals.copy()
            vals[inf_mask] = np.nan

        value_chunks.append(vals)

    if not value_chunks:
        return track_idx

    values = np.concatenate(value_chunks)

    if pm.gtrack_exists(full_track_name):
        if overwrite:
            pm.gtrack_rm(full_track_name, force=True)
        else:
            return track_idx

    pm.gtrack_create_dense(full_track_name, description, intervals_df, values, misha_binsize)
    return track_idx


# Module-level dict for sharing large data with forked children (COW).
_misha_shared: Dict = {}


def write_misha_tracks(
    predictions: Dict[str, Dict],
    track_names: List[str],
    bin_size: int,
    groot: str,
    track_prefix: str,
    track_suffix: str = "",
    misha_binsize: int = 32,
    description: str = "",
    overwrite: bool = True,
    n_jobs: int = 1,
    direct: bool = True,
) -> None:
    """Write predictions directly to Misha tracks.

    Args:
        predictions: Dict mapping chrom -> {starts, predictions} arrays.
        track_names: List of track names.
        bin_size: Prediction bin size in bp.
        groot: Misha genome database root path.
        track_prefix: Prefix for track names (e.g. "seq.IQ.MMGastru2FM.").
        track_suffix: Suffix for track names.
        misha_binsize: Misha track bin size (default 32).
        description: Track description.
        overwrite: Overwrite existing tracks.
        n_jobs: Number of parallel workers (fork-based).
        direct: Use fast direct binary writer (default True). Set to False
                to use the legacy pymisha gtrack_create_dense path.
    """
    if direct:
        write_misha_tracks_direct(
            predictions=predictions,
            track_names=track_names,
            bin_size=bin_size,
            groot=groot,
            track_prefix=track_prefix,
            track_suffix=track_suffix,
            misha_binsize=misha_binsize,
            description=description,
            overwrite=overwrite,
            n_jobs=n_jobs,
        )
        return

    global _misha_shared

    # Pre-compute shared intervals DataFrame: sort+dedup per-chrom
    interval_rows = []
    chrom_order = []
    for chrom in sorted(predictions.keys()):
        data = predictions[chrom]
        starts = data['starts'].astype(np.int64)
        sort_idx = np.argsort(starts)
        starts = starts[sort_idx]
        _, unique_idx = np.unique(starts, return_index=True)
        unique_idx = np.sort(unique_idx)
        starts = starts[unique_idx]
        ends = starts + bin_size
        n = len(starts)
        interval_rows.append(pd.DataFrame({
            'chrom': np.full(n, chrom, dtype=object),
            'start': starts,
            'end': ends,
        }))
        chrom_order.append(chrom)

    if not interval_rows:
        print("No predictions to write as Misha tracks.")
        return

    intervals_df = pd.concat(interval_rows, ignore_index=True)

    print(f"Writing {len(track_names)} Misha tracks to {groot} using {n_jobs} job(s)...")

    # Store large data in module globals so forked children inherit via COW
    # instead of pickling through the Pool queue.
    _misha_shared = {
        'predictions': predictions,
        'intervals_df': intervals_df,
        'chrom_order': chrom_order,
    }

    # Worker args: only small scalars — safe to pickle
    worker_args = [
        (i, f"{track_prefix}{track_names[i]}{track_suffix}",
         misha_binsize, groot, description, overwrite)
        for i in range(len(track_names))
    ]

    try:
        if n_jobs <= 1:
            import pymisha as pm
            pm.gdb_init(groot)
            for wa in tqdm(worker_args, desc="Writing Misha tracks"):
                _write_single_misha_track(wa)
        else:
            ctx = multiprocessing.get_context('fork')
            with ctx.Pool(n_jobs) as pool:
                list(tqdm(
                    pool.imap(_write_single_misha_track, worker_args),
                    total=len(track_names),
                    desc="Writing Misha tracks",
                ))
    finally:
        _misha_shared = {}

    print(f"Wrote {len(track_names)} Misha tracks.")


# ---------------------------------------------------------------------------
# Direct binary Misha track writer (bypasses slow pymisha gtrack_create_dense)
# ---------------------------------------------------------------------------

def _read_misha_chrom_sizes(groot: str) -> Dict[str, int]:
    """Read chromosome sizes from a Misha genome database root."""
    cs_path = os.path.join(groot, "chrom_sizes.txt")
    chrom_sizes: Dict[str, int] = {}
    with open(cs_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def _track_name_to_dir(groot: str, track_name: str) -> str:
    """Convert a dotted Misha track name to a filesystem directory path.

    Example: "seq.IQ.MMGastru2FM.flashzoi.mm10.metrics.loss.rf1k_atac_MM1"
    becomes: "{groot}/tracks/seq/IQ/MMGastru2FM/flashzoi/mm10/metrics/loss/rf1k_atac_MM1.track"
    """
    parts = track_name.split(".")
    path_parts = parts[:-1] + [parts[-1] + ".track"]
    return os.path.join(groot, "tracks", *path_parts)


def _write_single_track_direct(args):
    """Worker: write one Misha dense track as raw binary files.

    Reads shared data from module-level ``_misha_shared`` (fork COW).
    Does NOT call any pymisha functions — only writes binary files.
    """
    import struct
    import math
    import getpass
    from datetime import datetime

    (track_idx, full_track_name, track_dir, misha_binsize,
     chrom_sizes, description, overwrite) = args

    shared = _misha_shared
    predictions = shared["predictions"]

    if os.path.isdir(track_dir) and not overwrite:
        return track_idx

    os.makedirs(track_dir, exist_ok=True)
    vars_dir = os.path.join(track_dir, "vars")
    os.makedirs(vars_dir, exist_ok=True)

    # Write per-chromosome binary files
    for raw_chrom, chrom_size in chrom_sizes.items():
        chrom_key = f"chr{raw_chrom}" if not raw_chrom.startswith("chr") else raw_chrom
        chrom_file_name = chrom_key  # e.g. "chr1"

        num_bins = math.ceil(chrom_size / misha_binsize)
        # Allocate dense array filled with NaN
        dense = np.full(num_bins, np.nan, dtype=np.float32)

        if chrom_key in predictions:
            data = predictions[chrom_key]
            starts = data["starts"].astype(np.int64)
            vals = data["predictions"][:, track_idx].astype(np.float32)

            # Sort by position
            sort_idx = np.argsort(starts)
            starts = starts[sort_idx]
            vals = vals[sort_idx]

            # Deduplicate
            _, unique_idx = np.unique(starts, return_index=True)
            unique_idx = np.sort(unique_idx)
            starts = starts[unique_idx]
            vals = vals[unique_idx]

            # Replace Inf with NaN
            inf_mask = np.isinf(vals)
            if inf_mask.any():
                vals = vals.copy()
                vals[inf_mask] = np.nan

            # Fill dense array at bin positions
            bin_indices = starts // misha_binsize
            valid = (bin_indices >= 0) & (bin_indices < num_bins)
            dense[bin_indices[valid]] = vals[valid]

        # Write binary: [uint32 binsize][float32 × num_bins]
        out_path = os.path.join(track_dir, chrom_file_name)
        with open(out_path, "wb") as fout:
            fout.write(struct.pack("<I", misha_binsize))
            dense.tofile(fout)

    # Write .attributes
    now_str = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    user = getpass.getuser()
    attrs_path = os.path.join(track_dir, ".attributes")
    with open(attrs_path, "w") as fa:
        fa.write(f"created.by\twrite_misha_tracks_direct(\"{full_track_name}\")\n")
        fa.write(f"created.date\t{now_str}\n")
        fa.write(f"created.user\t{user}\n")
        fa.write(f"description\t{description}\n")

    return track_idx


def write_misha_tracks_direct(
    predictions: Dict[str, Dict],
    track_names: List[str],
    bin_size: int,
    groot: str,
    track_prefix: str,
    track_suffix: str = "",
    misha_binsize: int = 32,
    description: str = "",
    overwrite: bool = True,
    n_jobs: int = 4,
) -> None:
    """Write Misha dense tracks by writing binary files directly.

    ~100x faster than pymisha's ``gtrack_create_dense`` because it
    avoids Python/C bridge overhead and per-track ``gdb_reload()`` calls.
    Calls ``gdb_reload()`` once after all tracks are written.
    """
    global _misha_shared

    chrom_sizes = _read_misha_chrom_sizes(groot)

    print(f"Writing {len(track_names)} Misha tracks (direct) to {groot} "
          f"using {n_jobs} job(s)...")

    # Pre-delete existing tracks if overwrite requested
    if overwrite:
        import shutil
        for tname in track_names:
            full_name = f"{track_prefix}{tname}{track_suffix}"
            track_dir = _track_name_to_dir(groot, full_name)
            if os.path.isdir(track_dir):
                shutil.rmtree(track_dir)

    # Store shared data for fork COW
    _misha_shared = {
        "predictions": predictions,
        "chrom_order": sorted(predictions.keys()),
    }

    worker_args = []
    for i, tname in enumerate(track_names):
        full_name = f"{track_prefix}{tname}{track_suffix}"
        track_dir = _track_name_to_dir(groot, full_name)
        worker_args.append((
            i, full_name, track_dir, misha_binsize,
            chrom_sizes, description, overwrite,
        ))

    try:
        if n_jobs <= 1:
            for wa in tqdm(worker_args, desc="Writing Misha tracks (direct)"):
                _write_single_track_direct(wa)
        else:
            ctx = multiprocessing.get_context("fork")
            with ctx.Pool(n_jobs) as pool:
                list(tqdm(
                    pool.imap(_write_single_track_direct, worker_args),
                    total=len(track_names),
                    desc="Writing Misha tracks (direct)",
                ))
    finally:
        _misha_shared = {}

    # Reload Misha database once so it picks up all new tracks
    try:
        import pymisha as pm
        pm.gdb_init(groot)
        pm.gdb_reload()
        print("Misha database reloaded.")
    except Exception as exc:
        print(f"Warning: could not reload Misha database: {exc}")

    print(f"Wrote {len(track_names)} Misha tracks (direct).")


def write_bigwig_from_rank_files(
    rank_files: List[str],
    output_dir: str,
    chrom_sizes: Dict[str, int],
    track_names: List[str],
    bin_size: int = 32,
    prefix: str = "",
):
    """Write BigWig files by streaming per-chromosome data from rank NPZ files.

    This avoids loading the full genome predictions into memory.
    """
    import pyBigWig

    os.makedirs(output_dir, exist_ok=True)

    def chrom_sort_key(chrom):
        match = re.match(r'chr(\d+|[XYM])', chrom)
        if match:
            val = match.group(1)
            if val == 'X':
                return (100, chrom)
            elif val == 'Y':
                return (101, chrom)
            elif val == 'M':
                return (102, chrom)
            else:
                return (int(val), chrom)
        return (999, chrom)

    chrom_set = set()
    rank_data = [np.load(rank_file, allow_pickle=True) for rank_file in rank_files]
    try:
        for data in rank_data:
            if '_metadata' in data:
                metadata = data['_metadata'].item()
                chrom_set.update(metadata.get('chromosomes', []))
            else:
                for key in data.files:
                    if key.endswith('_starts'):
                        chrom_set.add(key[:-7])

        pred_chroms = sorted([c for c in chrom_set if c in chrom_sizes], key=chrom_sort_key)
        chroms_header = [(c, chrom_sizes[c]) for c in pred_chroms]

        bw_handles = []
        for name in track_names:
            safe_name = _sanitize_filename(name)
            output_path = os.path.join(output_dir, f"{prefix}{safe_name}.bw")
            bw = pyBigWig.open(output_path, "w")
            bw.addHeader(chroms_header)
            bw_handles.append(bw)

        for chrom, chrom_size in tqdm(chroms_header, desc="Writing BigWig files (by chrom)"):
            starts_list = []
            preds_list = []
            s_key = f"{chrom}_starts"
            p_key = f"{chrom}_predictions"

            for data in rank_data:
                if s_key in data and p_key in data:
                    starts_list.append(data[s_key])
                    preds_list.append(data[p_key])

            if not starts_list:
                continue

            starts = np.concatenate(starts_list)
            preds = np.vstack(preds_list)

            valid_bounds = starts < chrom_size
            starts = starts[valid_bounds]
            preds = preds[valid_bounds]

            sort_idx = np.argsort(starts)
            starts = starts[sort_idx]
            preds = preds[sort_idx]

            _, unique_idx = np.unique(starts, return_index=True)
            unique_idx = np.sort(unique_idx)
            starts = starts[unique_idx].astype(np.int64)
            preds = preds[unique_idx]

            # Fast path: dense, fixed-step bins (common for genome-wide inference)
            fixed_step = False
            if starts.size > 1:
                diffs = np.diff(starts)
                fixed_step = np.all(diffs == bin_size)

            all_finite = np.isfinite(preds).all()
            if fixed_step and all_finite:
                start0 = int(starts[0])
                for track_idx, bw in enumerate(bw_handles):
                    values = preds[:, track_idx]
                    bw.addEntries(
                        chrom,
                        start0,
                        values=values,
                        span=bin_size,
                        step=bin_size,
                    )
            else:
                ends = np.minimum(starts + bin_size, chrom_size).astype(np.int64)
                if all_finite:
                    for track_idx, bw in enumerate(bw_handles):
                        values = preds[:, track_idx]
                        bw.addEntries(
                            chrom,
                            starts,
                            ends=ends,
                            values=values,
                        )
                else:
                    for track_idx, bw in enumerate(bw_handles):
                        values = preds[:, track_idx]
                        valid = np.isfinite(values)
                        if valid.any():
                            bw.addEntries(
                                chrom,
                                starts[valid],
                                ends=ends[valid],
                                values=values[valid],
                            )

        for bw in bw_handles:
            bw.close()
    finally:
        for data in rank_data:
            data.close()

    print(f"Wrote {len(track_names)} BigWig files to {output_dir}")


def write_hdf5(
    predictions: Dict[str, Dict],
    output_path: str,
    track_names: List[str],
    bin_size: int = 32
):
    """Write predictions to HDF5 file."""
    import h5py

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    with h5py.File(output_path, 'w') as f:
        # Store metadata
        f.attrs['bin_size'] = bin_size
        f.attrs['num_tracks'] = len(track_names)
        f.create_dataset('track_names', data=[t.encode() for t in track_names])

        # Store per-chromosome data
        for chrom, data in predictions.items():
            grp = f.create_group(chrom)
            grp.create_dataset('starts', data=data['starts'], compression='gzip')
            grp.create_dataset('predictions', data=data['predictions'], compression='gzip')

    print(f"Wrote predictions to {output_path}")


def write_npz(
    predictions: Dict[str, Dict],
    output_path: str,
    track_names: List[str],
    bin_size: int = 32,
    compressed: bool = True,
):
    """Write predictions to NPZ file.

    Args:
        compressed: If True use np.savez_compressed (slower, smaller).
                    If False use np.savez (faster, larger).
    """
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    # Flatten for NPZ format
    arrays_to_save = {
        '_metadata': {
            'bin_size': bin_size,
            'track_names': track_names,
            'chromosomes': list(predictions.keys())
        }
    }

    for chrom, data in predictions.items():
        arrays_to_save[f'{chrom}_starts'] = data['starts']
        arrays_to_save[f'{chrom}_predictions'] = data['predictions']

    save_fn = np.savez_compressed if compressed else np.savez
    save_fn(output_path, **arrays_to_save)
    print(f"Wrote predictions to {output_path}")


def write_parquet(
    predictions: Dict[str, Dict],
    output_dir: str,
    track_names: List[str],
    bin_size: int = 32,
    prefix: str = "",
):
    """
    Write predictions to Parquet shards (one file per chromosome).

    Columns: chrom, start, end, track_0, track_1, ...
    """
    os.makedirs(output_dir, exist_ok=True)

    for chrom, data in predictions.items():
        starts = data['starts'].astype(np.int64)
        ends = (starts + bin_size).astype(np.int64)
        preds = data['predictions']  # shape: (bins, tracks)

        df = pd.DataFrame({
            'chrom': chrom,
            'start': starts,
            'end': ends,
        })

        for i, name in enumerate(track_names):
            df[name] = preds[:, i]

        output_path = os.path.join(output_dir, f"{prefix}{chrom}.parquet")
        df.to_parquet(output_path, index=False)

    print(f"Wrote Parquet shards to {output_dir}")


def _import_bigwigs_to_misha(
    bw_dir: str,
    groot: str,
    track_prefix: str,
    track_suffix: str = "",
    track_name_prefix: str = "",
    track_name_suffix: str = "",
    misha_binsize: int = 32,
    description: str = "",
    overwrite: bool = True,
) -> None:
    """Import BigWig files from a directory into Misha tracks via pymisha.

    Args:
        bw_dir: Directory containing .bw BigWig files.
        groot: Misha genome database root path.
        track_prefix: Prefix for track names (e.g. "seq.IQ.borzoi.epiblast_524k.").
        track_suffix: Suffix appended after the track name.
        track_name_prefix: Prefix for individual track names (before BigWig stem).
        track_name_suffix: Suffix for individual track names (after BigWig stem).
        misha_binsize: Misha track bin size.
        description: Track description string.
        overwrite: Overwrite existing tracks.
    """
    import pymisha as pm

    pm.gdb_init(groot)

    # Create parent directory hierarchy for the track prefix
    pm.gtrack_create_dirs(track_prefix)

    # Find all BigWig files
    bw_files = sorted(Path(bw_dir).glob("*.bw"))
    if not bw_files:
        raise FileNotFoundError(f"No BigWig files found in {bw_dir}")

    print(f"Found {len(bw_files)} BigWig files")

    success_count = 0
    error_count = 0

    for bw_path in bw_files:
        base_name = bw_path.stem

        # Skip files with "intervalID" in the name
        if "intervalid" in base_name.lower():
            print(f"\nSkipping: {bw_path.name} (contains 'intervalID')")
            continue

        # Construct full track name
        track_name = f"{track_name_prefix}{base_name}{track_name_suffix}"
        full_track_name = f"{track_prefix}{track_name}{track_suffix}"

        print(f"\nProcessing: {bw_path.name} -> {full_track_name}")
        pm.gtrack_create_dirs(full_track_name)

        if pm.gtrack_exists(full_track_name):
            if overwrite:
                pm.gtrack_rm(full_track_name, force=True)
            else:
                print(f"  SKIPPING: Track {full_track_name} already exists")
                continue

        try:
            pm.gtrack_import(
                track=full_track_name,
                description=description,
                file=str(bw_path),
                binsize=misha_binsize,
            )
            print(f"  SUCCESS: Created track {full_track_name}")
            success_count += 1
        except Exception as e:
            print(f"  ERROR: Failed to create track: {e}")
            error_count += 1

    # Reload database once after all tracks are written
    try:
        pm.gdb_reload()
        print("\nMisha database reloaded.")
    except Exception as e:
        print(f"\nWarning: Failed to reload Misha database: {e}")

    print(f"\n=== Summary ===")
    print(f"Total files: {len(bw_files)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {error_count}")

    if error_count > 0:
        raise RuntimeError(f"{error_count} track(s) failed to import")


def create_misha_tracks(
    config,
    output_dir: str,
    config_path: str,
    rank: int = 0,
    predictions: Optional[Dict[str, Dict]] = None,
    track_names: Optional[List[str]] = None,
    bin_size: int = 32,
    n_jobs: int = 1,
) -> None:
    """
    Create Misha tracks via pymisha (rank 0 only).

    When predictions are provided, writes tracks directly using pymisha's
    direct binary writer. Falls back to importing BigWig files via pymisha's
    gtrack_import() if predictions are not available.

    Args:
        config: Configuration object with misha_track section.
        output_dir: Directory containing BigWig files (used for BigWig import fallback).
        config_path: Path to YAML config file.
        rank: Process rank (only rank 0 executes).
        predictions: Dict mapping chrom -> {starts, predictions} arrays.
        track_names: List of track names.
        bin_size: Prediction bin size in bp.
        n_jobs: Number of parallel jobs for track writing.
    """
    # Only run on rank 0
    if rank != 0:
        return

    # Check if Misha track creation is enabled
    if not hasattr(config, 'misha_track'):
        return

    misha_cfg = config.misha_track

    if not misha_cfg.get('enabled', False):
        return

    # Validate required fields
    groot = misha_cfg.get('groot')
    track_prefix = misha_cfg.get('track_prefix')

    if not groot or not track_prefix:
        print("\nWarning: Misha track creation enabled but missing required fields (groot, track_prefix). Skipping.")
        return

    print("\n" + "=" * 70)
    print("Creating Misha Tracks")
    print("=" * 70)

    track_suffix = misha_cfg.get('suffix') or misha_cfg.get('track_suffix') or ""
    misha_binsize = int(misha_cfg.get('binsize', 32))
    description = str(misha_cfg.get('description') or '')

    # Try pymisha direct path
    if predictions is not None and track_names is not None:
        try:
            import pymisha  # noqa: F401
        except ImportError:
            print("Warning: pymisha not installed. Cannot create Misha tracks.")
            print("=" * 70)
            return

    if predictions is not None and track_names is not None:
        print(f"Misha genome root: {groot}")
        print(f"Track prefix: {track_prefix}")
        print(f"Track suffix: {track_suffix}")
        print(f"Misha binsize: {misha_binsize}")
        print(f"Parallel jobs: {n_jobs}")
        print()

        try:
            write_misha_tracks(
                predictions=predictions,
                track_names=track_names,
                bin_size=bin_size,
                groot=groot,
                track_prefix=track_prefix,
                track_suffix=track_suffix,
                misha_binsize=misha_binsize,
                description=description,
                overwrite=True,
                n_jobs=n_jobs,
            )
            print("\nMisha tracks created successfully!")
        except Exception as e:
            print(f"\nError: Misha track creation failed: {e}")
            import traceback
            traceback.print_exc()

        print("=" * 70)
        return

    # Fallback: import BigWig files via pymisha (no in-memory predictions)
    print("No in-memory predictions available. Importing BigWig files via pymisha.")

    bw_dir = misha_cfg.get('bw_dir') or output_dir
    track_name_prefix = str(misha_cfg.get('track_name_prefix') or '')
    track_name_suffix = str(misha_cfg.get('track_name_suffix') or '')

    try:
        import pymisha as pm
    except ImportError:
        print("Error: pymisha not installed. Cannot create Misha tracks.")
        print("=" * 70)
        return

    print(f"Misha genome root: {groot}")
    print(f"Track prefix: {track_prefix}")
    print(f"BigWig directory: {bw_dir}")
    print(f"Misha binsize: {misha_binsize}")
    print()

    try:
        _import_bigwigs_to_misha(
            bw_dir=bw_dir,
            groot=groot,
            track_prefix=track_prefix,
            track_suffix=track_suffix,
            track_name_prefix=track_name_prefix,
            track_name_suffix=track_name_suffix,
            misha_binsize=misha_binsize,
            description=description,
            overwrite=True,
        )
        print("\nMisha tracks created successfully!")
    except Exception as e:
        print(f"\nError: Misha track creation failed: {e}")
        import traceback
        traceback.print_exc()

    print("=" * 70)


# =============================================================================
# Deferred Postprocessing
# =============================================================================

def generate_finalize_script(
    output_dir: str,
    config_path: str,
    working_config,
    chrom_sizes_path: str,
    track_names: List[str],
    bin_size: int,
    bigwig_jobs: int,
    prefix: str,
    npz_paths: "str | List[str]",
    output_format: Optional[List[str]] = None,
    misha_cfg: Optional[dict] = None,
) -> None:
    """Generate finalize_inference.sh and finalize_inference.py for deferred CPU postprocessing.

    Saves manifest files (.track_names.json, .effective_config.yaml) and creates
    executable scripts that convert NPZ predictions to BigWig and/or Misha tracks.

    When misha is enabled and bigwig is the only output format, BigWig writing
    is skipped (it was only an intermediary for the R-based Misha import).

    Args:
        output_dir: Directory containing NPZ predictions and where scripts are written.
        config_path: Path to the original YAML config file.
        working_config: Fully-resolved config object (DotDict or dict).
        chrom_sizes_path: Path to chromosome sizes file.
        track_names: List of track names.
        bin_size: Bin size in bp.
        bigwig_jobs: Number of parallel jobs for BigWig/Misha writing.
        prefix: Prefix for output files.
        npz_paths: Path(s) to NPZ prediction file(s). A single string for
                   single-GPU runs or a list of rank file paths for multi-GPU.
        output_format: List of output format strings (e.g. ["bigwig"]).
        misha_cfg: Resolved misha_track config dict (or None if disabled).
    """
    output_dir = os.path.abspath(output_dir)
    if isinstance(npz_paths, str):
        npz_paths = [npz_paths]
    npz_paths = [os.path.abspath(p) for p in npz_paths]
    chrom_sizes_path = os.path.abspath(chrom_sizes_path)
    infer_script_dir = os.path.abspath(os.path.dirname(__file__))

    if output_format is None:
        output_format = ['bigwig']

    # Determine what to write
    misha_enabled = misha_cfg is not None
    # Skip BigWig when misha is the sole goal (bigwig was only an intermediary)
    non_bigwig_formats = [f for f in output_format if f != 'bigwig' and f != 'npz']
    write_bigwig_flag = 'bigwig' in output_format and (not misha_enabled or len(non_bigwig_formats) > 0)

    # Save .track_names.json (prefix makes names unique when multiple runs share a directory)
    track_names_path = os.path.join(output_dir, f'.{prefix}track_names.json')
    with open(track_names_path, 'w', encoding='utf-8') as f:
        json.dump(track_names, f)
    print(f"Saved track names to {track_names_path}")

    # Save .effective_config.yaml
    effective_config_path = os.path.join(output_dir, f'.{prefix}effective_config.yaml')

    def _to_plain_dict(value):
        if isinstance(value, dict):
            return {k: _to_plain_dict(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [_to_plain_dict(v) for v in value]
        return value

    with open(effective_config_path, 'w', encoding='utf-8') as f:
        yaml.safe_dump(_to_plain_dict(working_config), f, sort_keys=False)
    print(f"Saved effective config to {effective_config_path}")

    # Resolve misha variables for the generated script
    misha_groot = ''
    misha_track_prefix = ''
    misha_track_suffix = ''
    misha_binsize = 32
    misha_description = ''
    if misha_enabled:
        misha_groot = misha_cfg.get('groot', '')
        misha_track_prefix = misha_cfg.get('track_prefix', '')
        misha_track_suffix = misha_cfg.get('suffix') or misha_cfg.get('track_suffix') or ''
        misha_binsize = int(misha_cfg.get('binsize', 32))
        misha_description = str(misha_cfg.get('description') or '')

    # --- Generate finalize_inference.py ---
    py_script_path = os.path.join(output_dir, f'{prefix}finalize_inference.py')
    py_content = f'''\
#!/usr/bin/env python3
"""Auto-generated script for deferred postprocessing: NPZ -> BigWig/Misha."""

import os
import sys
import json
import numpy as np

# --- Editable variables ---
INFER_SCRIPT_DIR = {infer_script_dir!r}
CHROM_SIZES = {chrom_sizes_path!r}
OUTPUT_DIR = {output_dir!r}
NPZ_PATHS = {npz_paths!r}
TRACK_NAMES_PATH = {track_names_path!r}
BIN_SIZE = {bin_size!r}
BIGWIG_JOBS = {bigwig_jobs!r}
PREFIX = {prefix!r}
WRITE_BIGWIG = {write_bigwig_flag!r}
MISHA_ENABLED = {misha_enabled!r}
MISHA_GROOT = {misha_groot!r}
MISHA_TRACK_PREFIX = {misha_track_prefix!r}
MISHA_TRACK_SUFFIX = {misha_track_suffix!r}
MISHA_BINSIZE = {misha_binsize!r}
MISHA_DESCRIPTION = {misha_description!r}
N_JOBS = {bigwig_jobs!r}
# --- End editable variables ---

# Allow importing from the inference script directory
sys.path.insert(0, INFER_SCRIPT_DIR)
from src.data_loaders import load_chrom_sizes


def load_and_merge_npz(paths):
    """Load one or more NPZ files and merge into a single predictions dict."""
    all_predictions = {{}}

    for path in paths:
        print(f"Loading {{path}}...")
        data = np.load(path, allow_pickle=True)
        metadata = data["_metadata"].item()

        for chrom in metadata["chromosomes"]:
            if chrom not in all_predictions:
                all_predictions[chrom] = {{"starts": [], "predictions": []}}
            all_predictions[chrom]["starts"].append(data[f"{{chrom}}_starts"])
            all_predictions[chrom]["predictions"].append(data[f"{{chrom}}_predictions"])

    # Concatenate and sort per chromosome
    predictions = {{}}
    for chrom, arrays in all_predictions.items():
        starts = np.concatenate(arrays["starts"])
        preds = np.vstack(arrays["predictions"])
        sort_idx = np.argsort(starts)
        predictions[chrom] = {{"starts": starts[sort_idx], "predictions": preds[sort_idx]}}

    return predictions


def main():
    print(f"Loading chrom sizes from {{CHROM_SIZES}}...")
    chrom_sizes = load_chrom_sizes(CHROM_SIZES)

    print(f"Loading track names from {{TRACK_NAMES_PATH}}...")
    with open(TRACK_NAMES_PATH, "r", encoding="utf-8") as f:
        track_names = json.load(f)

    predictions = load_and_merge_npz(NPZ_PATHS)

    if WRITE_BIGWIG:
        from infer_borzoi_pytorch import write_bigwig
        print(f"Writing BigWig files to {{OUTPUT_DIR}}...")
        write_bigwig(
            predictions=predictions,
            output_dir=OUTPUT_DIR,
            chrom_sizes=chrom_sizes,
            track_names=track_names,
            bin_size=BIN_SIZE,
            prefix=PREFIX,
            n_jobs=BIGWIG_JOBS,
        )
        print("BigWig conversion complete.")

    if MISHA_ENABLED:
        from infer_borzoi_pytorch import write_misha_tracks
        print(f"Writing Misha tracks to {{MISHA_GROOT}}...")
        write_misha_tracks(
            predictions=predictions,
            track_names=track_names,
            bin_size=BIN_SIZE,
            groot=MISHA_GROOT,
            track_prefix=MISHA_TRACK_PREFIX,
            track_suffix=MISHA_TRACK_SUFFIX,
            misha_binsize=MISHA_BINSIZE,
            description=MISHA_DESCRIPTION,
            overwrite=True,
            n_jobs=N_JOBS,
        )
        print("Misha track creation complete.")

    print("Finalize complete.")

if __name__ == "__main__":
    main()
'''
    with open(py_script_path, 'w', encoding='utf-8') as f:
        f.write(py_content)
    os.chmod(py_script_path, 0o755)
    print(f"Generated {py_script_path}")

    # --- Generate finalize_inference.sh ---
    sh_script_path = os.path.join(output_dir, f'{prefix}finalize_inference.sh')

    sh_content = f'''\
#!/usr/bin/env bash
# Auto-generated script for deferred postprocessing.
# Run this on a CPU node after GPU inference has completed.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"

echo "Running finalize_inference.py..."
python "${{SCRIPT_DIR}}/{prefix}finalize_inference.py"
echo "Finalize complete."
'''
    with open(sh_script_path, 'w', encoding='utf-8') as f:
        f.write(sh_content)
    os.chmod(sh_script_path, 0o755)
    print(f"Generated {sh_script_path}")
    if misha_enabled and not write_bigwig_flag:
        print(f"  (BigWig skipped — Misha tracks will be written directly from NPZ)")
    print(f"\nTo complete postprocessing, run:\n  bash {sh_script_path}")


# =============================================================================
# Single Inference Config Runner
# =============================================================================

def run_single_inference(
    inf_cfg: dict,
    config: DotDict,
    args,
    device: torch.device,
    rank: int,
    world_size: int,
    config_path: str,
    inference_idx: int,
    total_inferences: int,
) -> None:
    """Run inference for a single inference config.

    Args:
        inf_cfg: The inference configuration dict.
        config: Full configuration object.
        args: Parsed command-line arguments.
        device: PyTorch device.
        rank: Process rank.
        world_size: Total number of processes.
        config_path: Path to the YAML config file.
        inference_idx: Index of this inference config (for logging).
        total_inferences: Total number of inference configs (for logging).
    """
    checkpoint_dir = getattr(config.logging, 'checkpoint_dir', './checkpoints')

    # Build config-specific header
    inf_name = inf_cfg.get('name', f'config_{inference_idx}')
    checkpoint_val = inf_cfg.get('checkpoint', 'best_model')
    metric_val = inf_cfg.get('metric')

    if rank == 0:
        print("\n" + "=" * 70)
        print(f"INFERENCE {inference_idx + 1}/{total_inferences}: {inf_name}")
        print(f"Checkpoint: {checkpoint_val}")
        print("=" * 70)

    # Resolve checkpoint path - command line --checkpoint overrides config
    if args.checkpoint is not None:
        checkpoint_path = args.checkpoint
    elif args.metric is not None:
        # User requested checkpoint by metric via CLI
        metric_checkpoint = find_checkpoint_by_metric(checkpoint_dir, args.metric, rank)
        if metric_checkpoint is not None:
            checkpoint_path = metric_checkpoint
        else:
            raise ValueError(
                f"Metric '{args.metric}' not found in W&B history or checkpoint metrics.json files. "
                f"Cannot proceed with metric-based checkpoint selection. "
                f"Either run backfill_wandb_metrics.py with --metrics {args.metric}, "
                f"or specify an explicit --checkpoint path."
            )
    elif metric_val is not None:
        # Per-inference metric-based selection (from config)
        metric_checkpoint = find_checkpoint_by_metric(checkpoint_dir, metric_val, rank)
        if metric_checkpoint is not None:
            checkpoint_path = metric_checkpoint
        else:
            raise ValueError(
                f"Metric '{metric_val}' not found in W&B history or checkpoint metrics.json files. "
                f"Cannot proceed with metric-based checkpoint selection. "
                f"Either run backfill_wandb_metrics.py with --metrics {metric_val}, "
                f"or specify an explicit checkpoint in config (remove 'metric' field and add 'checkpoint' field)."
            )
    else:
        checkpoint_path = _resolve_checkpoint_path(checkpoint_val, checkpoint_dir, rank)

    if rank == 0:
        print(f"Using checkpoint: {checkpoint_path}")

    # Resolve output directory - command line --output_dir overrides config
    if args.output_dir is not None:
        output_dir = args.output_dir
    else:
        output_dir = inf_cfg.get('output_dir')
        if output_dir is None:
            output_dir = os.path.join(checkpoint_dir, 'predictions')

    if rank == 0:
        print(f"Output directory: {output_dir}")

    sync_dir = inf_cfg.get('sync_dir')
    sync_root = _get_sync_root(config_path, inference_idx, sync_dir=sync_dir) if world_size > 1 else None

    # Use a shared sync_start across all ranks to avoid race conditions.
    # Without this, ranks that finish the previous inference earlier enter
    # the merge_done poll loop sooner and wake up later (due to the 5-second
    # poll interval), giving them a sync_start AFTER rank 0 has already
    # written files for the current inference — causing those files to be
    # rejected as stale and the rank to wait forever.
    if world_size > 1 and sync_root is not None:
        sync_start_path = f"{sync_root}.sync_start"
        if rank == 0:
            # Clean stale sync files from previous runs before writing new ones
            for stale in glob.glob(f"{sync_root}.*"):
                try:
                    os.remove(stale)
                except OSError:
                    pass
            sync_start = time.time()
            with open(sync_start_path, "w") as f:
                f.write(str(sync_start))
        else:
            _wait_for_paths([sync_start_path], label="sync_start")
            with open(sync_start_path, "r") as f:
                sync_start = float(f.read().strip())
    else:
        sync_start = time.time()

    # Resolve mode
    mode = args.mode if args.mode is not None else inf_cfg.get('mode', 'genome')

    # Resolve chromosomes (for mode=chromosomes)
    chromosomes = args.chromosomes
    if chromosomes is None and mode == 'chromosomes':
        chromosomes = inf_cfg.get('chromosomes')

    # Resolve BED file (for mode=bed)
    bed_file = args.bed_file
    if bed_file is None and mode == 'bed':
        bed_file = inf_cfg.get('bed_file')

    # Resolve exclude_chroms (for mode=genome)
    exclude_chroms = args.exclude_chroms
    if exclude_chroms is None and mode == 'genome':
        exclude_chroms = inf_cfg.get('exclude_chroms')
        if exclude_chroms is None:
            exclude_chroms = getattr(config.data, 'exclude_chroms', None)

    # Resolve batch size
    batch_size = args.batch_size
    if batch_size is None:
        batch_size = inf_cfg.get('batch_size')
        if batch_size is None:
            batch_size = getattr(config.training, 'batch_size', 8)

    # Resolve num_workers
    num_workers = args.num_workers
    if num_workers is None:
        num_workers = inf_cfg.get('num_workers')
        if num_workers is None:
            num_workers = getattr(config.training, 'num_workers', 4)

    # Resolve stride
    stride = args.stride
    if stride is None:
        stride = inf_cfg.get('stride')

    # Resolve output format
    output_format = args.output_format
    if output_format is None:
        output_format = inf_cfg.get('output_format', ['bigwig'])

    # Resolve defer_postprocess
    defer_postprocess = getattr(args, 'defer_postprocess', False) or inf_cfg.get('defer_postprocess', False)

    if defer_postprocess:
        if 'npz' not in output_format:
            output_format = list(output_format)
            output_format.append('npz')
        if rank == 0:
            print("[defer_postprocess] Forcing NPZ output; BigWig writing and Misha track creation will be deferred.")

    # Resolve track_names
    track_names_arg = args.track_names
    if track_names_arg is None:
        track_names_arg = inf_cfg.get('track_names')

    # Resolve bigwig_jobs
    bigwig_jobs = args.bigwig_jobs
    if bigwig_jobs is None:
        bigwig_jobs = inf_cfg.get('bigwig_jobs', 4)

    # Resolve RC average and mixed precision
    if args.no_rc_average:
        use_rc_average = False
    else:
        use_rc_average = inf_cfg.get('use_rc_average', True)

    if args.no_mixed_precision:
        use_mixed_precision = False
    else:
        use_mixed_precision = inf_cfg.get('mixed_precision', True)

    # torch.compile()
    use_compile = args.use_compile or inf_cfg.get('use_compile', False)

    # Resolve genome_fasta - per-inference config can override
    genome_fasta = args.genome_fasta
    if genome_fasta is None:
        genome_fasta = inf_cfg.get('genome_fasta')
        if genome_fasta is None:
            genome_fasta = config.data.genome_fasta

    # Create a working config copy with potentially overridden genome_fasta
    working_config = DotDict(dict(config))
    working_config['data'] = DotDict(dict(config.data))
    working_config['data']['genome_fasta'] = genome_fasta

    # Override species if provided via command line
    if args.species is not None:
        if 'model' not in working_config:
            working_config['model'] = {}
        working_config['model']['species'] = args.species

    # Merge per-inference misha_track with top-level misha_track
    if 'misha_track' in inf_cfg:
        per_inference_misha = inf_cfg.get('misha_track', {})
        if per_inference_misha:
            top_level_misha = config.get('misha_track', {})
            if isinstance(top_level_misha, dict):
                merged_misha = dict(top_level_misha)
            else:
                merged_misha = {}
            if isinstance(per_inference_misha, dict):
                merged_misha.update(per_inference_misha)
            working_config['misha_track'] = merged_misha
            if rank == 0:
                print(f"Using per-inference misha_track settings (merged with top-level)")

    # Override Misha track parameters if provided via command line
    if any([args.misha_groot, args.misha_track_prefix, args.misha_suffix, args.misha_binsize]):
        if 'misha_track' not in working_config:
            working_config['misha_track'] = {}

        if args.misha_groot is not None:
            working_config['misha_track']['groot'] = args.misha_groot
        if args.misha_track_prefix is not None:
            working_config['misha_track']['track_prefix'] = args.misha_track_prefix
        if args.misha_suffix is not None:
            working_config['misha_track']['suffix'] = args.misha_suffix
        if args.misha_binsize is not None:
            working_config['misha_track']['binsize'] = args.misha_binsize

    # Load chromosome sizes
    chrom_sizes = load_chrom_sizes(config.data.chrom_sizes)

    # Generate regions based on mode
    if mode == 'genome':
        if rank == 0:
            print(f"Mode: Whole genome inference")
        regions = generate_regions_genome(chrom_sizes, exclude_chroms)
    elif mode == 'chromosomes':
        if not chromosomes:
            raise ValueError("--chromosomes required for mode=chromosomes")
        if rank == 0:
            print(f"Mode: Selected chromosomes: {chromosomes}")
        regions = generate_regions_chromosomes(chrom_sizes, chromosomes)
    elif mode == 'bed':
        if not bed_file:
            raise ValueError("--bed_file required for mode=bed")
        if rank == 0:
            print(f"Mode: BED regions from {bed_file}")
        regions = generate_regions_bed(bed_file)

    if rank == 0:
        total_bp = sum(r[2] - r[1] for r in regions)
        print(f"Total regions: {len(regions)}, Total bp: {total_bp:,}")

    # Split regions for distributed
    if world_size > 1:
        regions = split_regions_for_rank(regions, rank, world_size)
        if rank == 0:
            print(f"Each GPU handles ~{len(regions)} regions")

    # Load perturbations if provided
    splice_perts = None
    seq_perts = None

    if args.splice_pert:
        if rank == 0:
            print(f"\nLoading splice perturbations from {args.splice_pert}")
        splice_perts = load_splice_pert(args.splice_pert)
        if rank == 0:
            print(f"  Loaded {len(splice_perts)} splice operations")

    if args.seq_pert:
        if rank == 0:
            print(f"Loading sequence perturbations from {args.seq_pert}")
        seq_perts = load_seq_pert(args.seq_pert)
        if rank == 0:
            print(f"  Loaded {len(seq_perts)} sequence perturbations")

    # Create dataset
    genome_fasta_path = working_config['data']['genome_fasta']
    if rank == 0:
        print(f"\nUsing genome FASTA: {genome_fasta_path}")

    # Ensure FASTA index exists before all ranks try to access it
    fai_path = genome_fasta_path + '.fai'
    if rank == 0:
        if not os.path.exists(fai_path):
            print(f"Creating FASTA index: {fai_path}")
            _fasta = pysam.FastaFile(genome_fasta_path)
            _fasta.close()
            print(f"FASTA index created successfully")

    # Synchronize all ranks before proceeding (file-based to avoid NCCL timeouts)
    if world_size > 1:
        _wait_for_paths([fai_path], label="FASTA index")

    dataset = GenomeInferenceDataset(
        regions=regions,
        genome_fasta=genome_fasta_path,
        chrom_sizes=chrom_sizes,
        seq_len=config.model.seq_len,
        pred_len=config.model.pred_len,
        bin_size=config.model.bin_size,
        stride=stride,
        splice_perts=splice_perts,
        seq_perts=seq_perts,
    )

    if rank == 0:
        print(f"Total windows: {len(dataset)}")
        print(f"Window size: {config.model.seq_len} bp")
        print(f"Prediction size: {config.model.pred_len} bp ({config.model.pred_len // config.model.bin_size} bins)")

    # Load model
    model, is_human = load_model(checkpoint_path, working_config, device)

    if rank == 0:
        print(f"Using {'human' if is_human else 'mouse'} head for inference")

    # Get number of output tracks from the correct head
    num_tracks = None

    if is_human:
        if hasattr(model, 'human_head') and isinstance(model.human_head, nn.Conv1d):
            num_tracks = model.human_head.out_channels
        elif hasattr(model, 'head') and isinstance(model.head, nn.Conv1d):
            num_tracks = model.head.out_channels
    else:
        if hasattr(model, 'mouse_head') and isinstance(model.mouse_head, nn.Conv1d):
            num_tracks = model.mouse_head.out_channels

    if num_tracks is None:
        with torch.no_grad():
            test_input = torch.zeros(1, 4, config.model.seq_len, device=device)
            test_output = model(test_input, is_human=is_human)
            if isinstance(test_output, tuple):
                test_output = test_output[0]
            num_tracks = test_output.shape[1]

    if rank == 0:
        print(f"Number of output tracks: {num_tracks}")

    # Set track names
    track_names = None
    if track_names_arg:
        track_names = track_names_arg
    elif rank == 0:
        if hasattr(config.data, 'data_matrix') and config.data.data_matrix:
            if rank == 0:
                print(f"Loading track names from {config.data.data_matrix}...")
            selected_tracks = config.get('data.selected_tracks', None)
            track_names = _load_track_names_from_data_matrix(
                config.data.data_matrix,
                num_tracks,
                rank=rank,
                cache_dir=config.get('performance.data_matrix_cache', None),
                selected_tracks=selected_tracks,
            )
            if track_names and rank == 0:
                print(f"Loaded track names: {_summarize_track_names(track_names)}")
        if track_names is None:
            track_names = [f'track_{i}' for i in range(num_tracks)]
            if rank == 0:
                print(f"Using generic track names: {_summarize_track_names(track_names)}")
        if sync_root is not None:
            track_names_path = f"{sync_root}.track_names.json"
            with open(track_names_path, "w", encoding="utf-8") as handle:
                json.dump(track_names, handle)
    elif sync_root is not None:
        track_names_path = f"{sync_root}.track_names.json"
        _wait_for_paths([track_names_path], label="track names", min_mtime=sync_start)
        with open(track_names_path, "r", encoding="utf-8") as handle:
            track_names = json.load(handle)

    if track_names is None:
        track_names = [f'track_{i}' for i in range(num_tracks)]

    # Initialize batch manager for OOM handling
    batch_manager = AdaptiveBatchManager(
        initial_batch_size=batch_size,
        min_batch_size=args.min_batch_size,
        reduction_factor=args.batch_reduction_factor,
        recovery_delay=args.oom_recovery_delay,
    )

    if rank == 0:
        print(f"\nOOM protection enabled:")
        print(f"  Initial batch size: {batch_manager.current_batch_size}")
        print(f"  Minimum batch size: {batch_manager.min_batch_size}")

    # Run inference with OOM recovery
    results = None
    max_retries = 10

    for attempt in range(max_retries):
        try:
            if rank == 0:
                if attempt == 0:
                    print(f"\nStarting inference with batch_size={batch_manager.current_batch_size}...")
                else:
                    print(f"\nRetrying inference (attempt {attempt + 1}/{max_retries}) with batch_size={batch_manager.current_batch_size}...")

            loader_kwargs = {
                "dataset": dataset,
                "batch_size": batch_manager.current_batch_size,
                "shuffle": False,
                "num_workers": num_workers,
                "pin_memory": True,
                "collate_fn": collate_with_metadata,
            }
            if num_workers > 0:
                loader_kwargs["prefetch_factor"] = 2
                loader_kwargs["persistent_workers"] = True

            dataloader = DataLoader(**loader_kwargs)

            results = run_inference(
                model=model,
                dataloader=dataloader,
                device=device,
                use_rc_average=use_rc_average,
                mixed_precision=use_mixed_precision,
                rank=rank,
                world_size=world_size,
                is_human=is_human,
                use_compile=use_compile,
                compile_mode=args.compile_mode,
            )

            if batch_manager.oom_count > 0:
                print(f"\n[Rank {rank}] Inference completed successfully after {batch_manager.oom_count} OOM recovery attempt(s)")
            break

        except RuntimeError as e:
            if 'out of memory' in str(e).lower():
                print(f"\n[Rank {rank}] OOM Error detected during inference!")

                del dataloader
                torch.cuda.empty_cache()
                time.sleep(batch_manager.recovery_delay)

                if not batch_manager.reduce_batch_size():
                    error_msg = f"Cannot reduce batch size below {batch_manager.min_batch_size}. OOM cannot be recovered."
                    print(f"\n[Rank {rank}] ERROR: {error_msg}")
                    raise RuntimeError(error_msg) from e

                print(f"[Rank {rank}] Reducing batch size to {batch_manager.current_batch_size}")
            else:
                raise

    if results is None:
        raise RuntimeError(f"Inference failed after {max_retries} attempts")

    # Release GPU memory
    if rank == 0:
        print(f"\nReleasing GPU memory...")

    del model
    torch.cuda.empty_cache()

    if rank == 0:
        print(f"GPU memory released. Continuing with CPU-only operations...")

    # Merge predictions
    predictions = merge_predictions(results, config.model.bin_size)

    # Gather predictions from all ranks (if distributed)
    bigwig_written = False
    merge_done_path = None
    rank_npz_paths = None
    if world_size > 1:
        if sync_root is None:
            sync_root = _get_sync_root(config_path, inference_idx, sync_dir=sync_dir)
        if rank == 0:
            output_dir_path = f"{sync_root}.output_dir"
            with open(output_dir_path, "w", encoding="utf-8") as handle:
                handle.write(output_dir)
        else:
            output_dir_path = f"{sync_root}.output_dir"
            _wait_for_paths([output_dir_path], label="output_dir", min_mtime=sync_start)
            with open(output_dir_path, "r", encoding="utf-8") as handle:
                output_dir = handle.read().strip()
        rank_output = os.path.join(output_dir, f'.rank_{rank}_predictions.npz')
        os.makedirs(output_dir, exist_ok=True)
        write_npz(predictions, rank_output, track_names, config.model.bin_size,
                  compressed=not defer_postprocess)
        rank_done = f"{sync_root}.rank_{rank}.done"
        Path(rank_done).touch()
        del predictions
        gc.collect()

        if rank == 0:
            all_rank_done = [f"{sync_root}.rank_{r}.done" for r in range(world_size)]
            _wait_for_paths(all_rank_done, label="rank prediction files", min_mtime=sync_start)
            # Stream directly to BigWig only when misha is disabled (otherwise
            # we need the merged predictions dict for misha track creation).
            _misha_active = (
                hasattr(working_config, 'misha_track')
                and isinstance(working_config.get('misha_track', None), dict)
                and working_config.misha_track.get('enabled', False)
            )
            stream_bigwig_only = set(output_format) == {'bigwig'} and not _misha_active
            rank_files = [
                os.path.join(output_dir, f'.rank_{r}_predictions.npz')
                for r in range(world_size)
            ]

            if stream_bigwig_only and not defer_postprocess:
                print("Writing BigWig files from rank predictions (streaming by chromosome)...")
                write_bigwig_from_rank_files(
                    rank_files=rank_files,
                    output_dir=output_dir,
                    chrom_sizes=chrom_sizes,
                    track_names=track_names,
                    bin_size=config.model.bin_size,
                    prefix=args.prefix,
                )

                for rank_file in rank_files:
                    if os.path.exists(rank_file):
                        os.remove(rank_file)

                bigwig_written = True
            elif defer_postprocess:
                # Keep rank files for the finalize script to merge
                print("Keeping rank prediction files for deferred postprocessing...")
                rank_npz_paths = []
                for rank_file in rank_files:
                    # Rename from hidden (.rank_N_...) to visible (rank_N_...)
                    final_name = os.path.join(
                        os.path.dirname(rank_file),
                        os.path.basename(rank_file).lstrip('.')
                    )
                    os.rename(rank_file, final_name)
                    rank_npz_paths.append(final_name)
                predictions = None  # No in-memory merge needed
            else:
                print("Merging predictions from all GPUs...")
                all_predictions = {}

                for rank_file in rank_files:
                    data = np.load(rank_file, allow_pickle=True)
                    metadata = data['_metadata'].item()

                    for chrom in metadata['chromosomes']:
                        if chrom not in all_predictions:
                            all_predictions[chrom] = {
                                'starts': [],
                                'predictions': []
                            }
                        all_predictions[chrom]['starts'].append(data[f'{chrom}_starts'])
                        all_predictions[chrom]['predictions'].append(data[f'{chrom}_predictions'])

                    os.remove(rank_file)

                for chrom in all_predictions:
                    starts = np.concatenate(all_predictions[chrom]['starts'])
                    preds = np.vstack(all_predictions[chrom]['predictions'])
                    sort_idx = np.argsort(starts)
                    all_predictions[chrom] = {
                        'starts': starts[sort_idx],
                        'predictions': preds[sort_idx]
                    }

                predictions = all_predictions
            merge_done_path = f"{sync_root}.merge_done"
        else:
            merge_done_path = f"{sync_root}.merge_done"

    # Write output (only rank 0)
    npz_path = None
    if rank == 0:
        os.makedirs(output_dir, exist_ok=True)

        for fmt in output_format:
            if fmt == 'bigwig':
                if defer_postprocess:
                    continue
                if bigwig_written:
                    continue
                write_bigwig(
                    predictions=predictions,
                    output_dir=output_dir,
                    chrom_sizes=chrom_sizes,
                    track_names=track_names,
                    bin_size=config.model.bin_size,
                    prefix=args.prefix,
                    n_jobs=bigwig_jobs,
                )
            elif fmt == 'hdf5':
                if defer_postprocess and predictions is None:
                    continue
                write_hdf5(
                    predictions=predictions,
                    output_path=os.path.join(output_dir, f'{args.prefix}predictions.h5'),
                    track_names=track_names,
                    bin_size=config.model.bin_size,
                )
            elif fmt == 'npz':
                if defer_postprocess and predictions is None:
                    continue
                npz_path = os.path.join(output_dir, f'{args.prefix}predictions.npz')
                write_npz(
                    predictions=predictions,
                    output_path=npz_path,
                    track_names=track_names,
                    bin_size=config.model.bin_size,
                )
            elif fmt == 'parquet':
                if defer_postprocess and predictions is None:
                    continue
                write_parquet(
                    predictions=predictions,
                    output_dir=output_dir,
                    track_names=track_names,
                    bin_size=config.model.bin_size,
                    prefix=args.prefix,
                )

        # Misha tracks — independent of output_format
        if not defer_postprocess and predictions is not None:
            create_misha_tracks(
                config=working_config,
                output_dir=output_dir,
                config_path=config_path,
                rank=rank,
                predictions=predictions,
                track_names=track_names,
                bin_size=config.model.bin_size,
                n_jobs=bigwig_jobs,
            )

        if defer_postprocess:
            # Resolve misha config for finalize script
            resolved_misha_cfg = None
            if hasattr(working_config, 'misha_track'):
                mc = working_config.get('misha_track', {}) if hasattr(working_config, 'get') else working_config.get('misha_track', {})
                if isinstance(mc, dict) and mc.get('enabled', False):
                    resolved_misha_cfg = mc

            # Use rank file paths if available (multi-GPU), otherwise single NPZ
            if rank_npz_paths:
                finalize_npz_paths = rank_npz_paths
            elif npz_path is not None:
                finalize_npz_paths = npz_path
            else:
                finalize_npz_paths = os.path.join(output_dir, f'{args.prefix}predictions.npz')
            generate_finalize_script(
                output_dir=output_dir,
                config_path=config_path,
                working_config=working_config,
                chrom_sizes_path=config.data.chrom_sizes,
                track_names=track_names,
                bin_size=config.model.bin_size,
                bigwig_jobs=bigwig_jobs,
                prefix=args.prefix,
                npz_paths=finalize_npz_paths,
                output_format=output_format,
                misha_cfg=resolved_misha_cfg,
            )

    if world_size > 1 and merge_done_path is not None:
        if rank == 0:
            Path(merge_done_path).touch()
        else:
            _wait_for_paths([merge_done_path], label="merge completion", min_mtime=sync_start)

    print(f"Inference {inference_idx + 1}/{total_inferences} complete!")


# =============================================================================
# Distributed Utils
# =============================================================================

def setup_distributed():
    """Initialize distributed training if available."""
    if 'RANK' in os.environ and 'WORLD_SIZE' in os.environ:
        rank = int(os.environ['RANK'])
        world_size = int(os.environ['WORLD_SIZE'])
        local_rank = int(os.environ.get('LOCAL_RANK', 0))

        # Ensure the correct device is set before NCCL init to avoid rank/GPU mismatches.
        if torch.cuda.is_available():
            torch.cuda.set_device(local_rank)
        dist.init_process_group(backend='nccl', timeout=timedelta(hours=1))

        return rank, world_size, local_rank
    else:
        return 0, 1, 0


def cleanup_distributed():
    """Clean up distributed resources."""
    if dist.is_initialized():
        dist.destroy_process_group()


def split_regions_for_rank(
    regions: List[Tuple[str, int, int]],
    rank: int,
    world_size: int
) -> List[Tuple[str, int, int]]:
    """Split regions among distributed ranks."""
    # Sort regions by size for better load balancing
    regions_with_size = [(r, r[2] - r[1]) for r in regions]
    regions_with_size.sort(key=lambda x: -x[1])  # Sort by size descending

    # Distribute to ranks in round-robin fashion
    rank_regions = [[] for _ in range(world_size)]
    for i, (region, _) in enumerate(regions_with_size):
        rank_regions[i % world_size].append(region)

    return rank_regions[rank]


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Genome-wide inference for fine-tuned Borzoi models',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument('--config', type=str, required=True,
                        help='Path to configuration YAML file')
    
    # Optional arguments (will be filled from config if not provided)
    parser.add_argument('--checkpoint', type=str, default=None,
                        help='Path to model checkpoint (or use inference.checkpoint from config)')
    parser.add_argument('--metric', type=str, default=None,
                        help='Select checkpoint by best metric (e.g., val_pearson, val_genome_wide_pearson). '
                             'Looks for best_model_{metric}/ directory first, then scans metrics.json files.')
    parser.add_argument('--output_dir', type=str, default=None,
                        help='Output directory for predictions (or use inference.output_dir from config)')
    parser.add_argument('--genome_fasta', type=str, default=None,
                        help='Path to genome FASTA file (overrides data.genome_fasta from config)')

    # Species selection (for model head)
    parser.add_argument('--species', type=str, choices=['mouse', 'human'], default=None,
                        help='Species for model head selection (mouse or human). Overrides auto-detection from genome path.')

    # Input mode
    parser.add_argument('--mode', type=str, choices=['genome', 'chromosomes', 'bed'], default=None,
                        help='Inference mode: genome (all chroms), chromosomes (selected), bed (regions)')
    parser.add_argument('--chromosomes', type=str, nargs='+', default=None,
                        help='Chromosomes to predict (for mode=chromosomes)')
    parser.add_argument('--bed_file', type=str, default=None,
                        help='BED file with regions (for mode=bed)')
    parser.add_argument('--exclude_chroms', type=str, nargs='+', default=None,
                        help='Chromosomes to exclude (for mode=genome)')

    # Inference options
    parser.add_argument('--batch_size', type=int, default=None,
                        help='Batch size for inference (or use inference.batch_size from config)')
    parser.add_argument('--num_workers', type=int, default=None,
                        help='Number of data loading workers (or use inference.num_workers from config)')
    parser.add_argument('--no_rc_average', action='store_true',
                        help='Disable reverse complement averaging (overrides inference.use_rc_average)')
    parser.add_argument('--no_mixed_precision', action='store_true',
                        help='Disable mixed precision (overrides inference.mixed_precision)')
    parser.add_argument('--use_compile', action='store_true',
                        help='Use torch.compile() for 10-30%% faster inference (requires PyTorch 2.0+)')
    parser.add_argument('--compile_mode', type=str, default='reduce-overhead',
                        choices=['default', 'reduce-overhead', 'max-autotune'],
                        help='torch.compile mode (default: reduce-overhead)')
    parser.add_argument('--stride', type=int, default=None,
                        help='Stride for sliding windows (defaults to pred_len, or use inference.stride from config)')

    # OOM safety options
    parser.add_argument('--min_batch_size', type=int, default=1,
                        help='Minimum batch size for OOM recovery (default: 1)')
    parser.add_argument('--batch_reduction_factor', type=float, default=0.5,
                        help='Factor to reduce batch size on OOM (default: 0.5)')
    parser.add_argument('--oom_recovery_delay', type=float, default=2.0,
                        help='Delay in seconds after OOM before retry (default: 2.0)')

    # Perturbation options
    parser.add_argument('--splice_pert', type=str, default=None,
                        help='Path to splice perturbation table (TSV) - applies region transplantations during inference')
    parser.add_argument('--seq_pert', type=str, default=None,
                        help='Path to sequence perturbation table (TSV) - applies sequence edits during inference')

    # Output options
    parser.add_argument('--output_format', type=str, nargs='+', default=None,
                        choices=['bigwig', 'hdf5', 'npz', 'parquet'],
                        help='Output format(s) (or use inference.output_format from config)')
    parser.add_argument('--track_names', type=str, nargs='+', default=None,
                        help='Track names (or use inference.track_names from config)')
    parser.add_argument('--prefix', type=str, default='',
                        help='Prefix for output files')
    parser.add_argument('--bigwig_jobs', type=int, default=None,
                        help='Number of parallel jobs for writing BigWig files (or use inference.bigwig_jobs from config)')

    # Misha track options
    parser.add_argument('--misha_groot', type=str, default=None,
                        help='Misha genome root directory (overrides misha_track.groot from config)')
    parser.add_argument('--misha_track_prefix', type=str, default=None,
                        help='Misha track name prefix (overrides misha_track.track_prefix from config)')
    parser.add_argument('--misha_suffix', type=str, default=None,
                        help='Misha track name suffix (overrides misha_track.suffix from config)')
    parser.add_argument('--misha_binsize', type=int, default=None,
                        help='Misha track binsize (overrides misha_track.binsize from config)')

    # Deferred postprocessing
    parser.add_argument('--defer_postprocess', action='store_true',
                        help='Save NPZ predictions only; skip BigWig writing and Misha track creation. '
                             'Generates a finalize_inference.sh script for later CPU postprocessing.')

    # Multiple inference config support
    parser.add_argument('--inference_name', type=str, default=None,
                        help='Name of specific inference config to run (when config has multiple inference configs as a list)')
    parser.add_argument('--run_all', action='store_true', default=True,
                        help='Run all inference configs with run_after_training=true (default: True)')
    parser.add_argument('--no_run_all', action='store_true',
                        help='Only run the first inference config (disables --run_all)')

    args = parser.parse_args()

    # Handle --no_run_all flag
    run_all = args.run_all and not args.no_run_all

    # Setup distributed
    rank, world_size, local_rank = setup_distributed()
    device = torch.device(f'cuda:{local_rank}' if torch.cuda.is_available() else 'cpu')

    if rank == 0:
        print(f"Running inference with {world_size} GPU(s)")
        print(f"Device: {device}")

    try:
        # Load config first
        config = load_config(args.config)

        # Resolve inference configs to run
        inference_configs = _resolve_inference_configs(
            config=config,
            inference_name=args.inference_name,
            run_all=run_all,
            rank=rank
        )

        if not inference_configs:
            if rank == 0:
                print("No inference configs found to run.")
                print("Please add an 'inference' section to your config file or specify --checkpoint and --output_dir.")
            return

        if rank == 0:
            print(f"\n{'=' * 70}")
            print(f"TOTAL INFERENCE CONFIGS TO RUN: {len(inference_configs)}")
            print(f"{'=' * 70}")

        # Run each inference config
        for idx, inf_cfg in enumerate(inference_configs):
            run_single_inference(
                inf_cfg=inf_cfg,
                config=config,
                args=args,
                device=device,
                rank=rank,
                world_size=world_size,
                config_path=args.config,
                inference_idx=idx,
                total_inferences=len(inference_configs),
            )

        if rank == 0:
            print("\n" + "=" * 70)
            print(f"ALL {len(inference_configs)} INFERENCE(S) COMPLETE!")
            print("=" * 70)

    finally:
        cleanup_distributed()


if __name__ == '__main__':
    main()
