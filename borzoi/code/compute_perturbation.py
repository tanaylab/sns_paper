#!/usr/bin/env python3
"""
Simplified perturbation analysis for Borzoi models.

This script uses a clean table-based API with two-stage perturbations:
1. Splice perturbations (transplantations)
2. Sequence perturbations (fine-grained edits)

Example usage (single GPU):
    python compute_perturbation.py \
        --config configs/config.yaml \
        --checkpoint checkpoints/best_model.pt \
        --regions regions.bed \
        --splice_pert splice.pert \
        --seq_pert seq.pert \
        --track_name EB4_cnt \
        --output results/perturbations.parquet

Example usage (multi-GPU via torchrun):
    CUDA_VISIBLE_DEVICES=0,1,2,3 torchrun --nproc_per_node=4 compute_perturbation.py \
        --config configs/config.yaml \
        --checkpoint checkpoints/best_model.pt \
        --regions regions.bed \
        --splice_pert splice.pert \
        --seq_pert seq.pert \
        --track_name all \
        --output results/perturbations.parquet \
        --distributed

Example usage (baseline + seq.pert only, no splice):
    python compute_perturbation.py \
        --config configs/config.yaml \
        --checkpoint checkpoints/best_model.pt \
        --regions regions.bed \
        --seq_pert seq.pert \
        --track_name EB4_cnt,EB5_cnt \
        --output results/perturbations.parquet
"""

import argparse
import os
import sys
import time
from pathlib import Path
from typing import List, Dict, Tuple

import pandas as pd
import torch
import yaml

# Quieter TF logs
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
os.environ.setdefault('TF_ENABLE_ONEDNN_OPTS', '0')

# Import local modules
from src.dependency_io import load_config, load_chrom_sizes, load_track_names
from src.dependency_distributed import setup_distributed, cleanup_distributed, print_rank0
from src.perturbation_simple import (
    load_regions,
    load_splice_pert,
    load_seq_pert,
    GenomeCache,
    SimplePerturbationAnalyzer,
)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Simplified perturbation analysis for Borzoi models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Required arguments
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config YAML file (contains genome_fasta, chrom_sizes, etc.)",
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        default=None,
        help="Path to model checkpoint. If not provided, uses inference.checkpoint from config or best_model from checkpoint_dir.",
    )
    parser.add_argument(
        "--metric",
        type=str,
        default=None,
        help="Select checkpoint by best metric (e.g., val_pearson). Looks for best_model_{metric}/ directory.",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Path to region table (BED or CSV with chrom, start, end, region_id)",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output parquet file path",
    )

    # Perturbation tables (at least one required)
    parser.add_argument(
        "--splice_pert",
        type=str,
        help="Path to splice perturbation table (TSV with splice operations)",
    )
    parser.add_argument(
        "--seq_pert",
        type=str,
        help="Path to sequence perturbation table (TSV with sequence edits)",
    )

    # Track selection (optional - defaults to "all" if not provided)
    track_group = parser.add_mutually_exclusive_group()
    track_group.add_argument(
        "--track_name",
        type=str,
        default=None,
        help='Track name(s) to analyze. Can be: single name, comma-separated list, or "all". Defaults to "all" if neither track_name nor track_idx provided.',
    )
    track_group.add_argument(
        "--track_idx",
        type=str,
        default=None,
        help="Track index/indices to analyze (comma-separated). Alternative to --track_name",
    )

    # Optional arguments
    parser.add_argument(
        "--normalize_regions",
        type=int,
        help="Normalize all regions to this size (e.g., 524288). Overridden by per-region settings in BED file.",
    )
    parser.add_argument(
        "--data_matrix",
        type=str,
        default=None,
        help="Path to observed data parquet file (for adding observed columns). If not provided, uses data.data_matrix from config.",
    )
    parser.add_argument(
        "--inference_batch_size",
        type=int,
        default=64,
        help="Batch size for inference (default: 64)",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=16,
        help="Number of workers for parallel sequence loading (default: 16)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=60427,
        help="Random seed (default: 60427)",
    )

    # Inference options
    parser.add_argument(
        "--no_rc_average",
        action="store_true",
        help="Disable reverse complement averaging",
    )
    parser.add_argument(
        "--no_mixed_precision",
        action="store_true",
        help="Disable mixed precision inference",
    )
    parser.add_argument(
        "--use_compile",
        action="store_true",
        help="Enable torch.compile() for faster inference",
    )
    parser.add_argument(
        "--compile_mode",
        type=str,
        default="default",
        choices=["default", "reduce-overhead", "max-autotune"],
        help="Compile mode if --use_compile is enabled (default: default)",
    )

    # Distributed training
    parser.add_argument(
        "--distributed",
        action="store_true",
        help="Enable distributed training (use with torchrun)",
    )

    # Verbose output
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress bars and verbose output",
    )
    
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Measure and print timing statistics",
    )

    args = parser.parse_args()

    # Validation
    if not args.splice_pert and not args.seq_pert:
        parser.error("At least one of --splice_pert or --seq_pert is required")

    # Default track_name to "all" if neither provided
    if not args.track_name and not args.track_idx:
        args.track_name = "all"

    return args


def resolve_tracks(
    args,
    config: dict,
) -> Tuple[List[int], List[str]]:
    """
    Resolve track indices and names from arguments and config.

    Returns:
        (track_indices, track_names) tuple
    """
    # Load all available track names from config
    all_track_names = load_track_names(config)

    if args.track_name:
        # Parse track_name argument
        if args.track_name.lower() == "all":
            # Use all tracks
            track_names = all_track_names
            track_indices = list(range(len(track_names)))
        else:
            # Parse comma-separated list
            track_names = [name.strip() for name in args.track_name.split(',')]
            track_indices = []
            for name in track_names:
                try:
                    idx = all_track_names.index(name)
                    track_indices.append(idx)
                except ValueError:
                    raise ValueError(
                        f"Track '{name}' not found in config. "
                        f"Available tracks: {all_track_names[:10]}..."
                    )
    else:
        # Parse track_idx argument
        idx_strs = args.track_idx.split(',')
        track_indices = [int(s.strip()) for s in idx_strs]

        # Validate indices
        max_idx = len(all_track_names) - 1
        for idx in track_indices:
            if idx < 0 or idx > max_idx:
                raise ValueError(
                    f"Track index {idx} out of range [0, {max_idx}]. "
                    f"Total tracks: {len(all_track_names)}"
                )

        # Get names
        track_names = [all_track_names[idx] for idx in track_indices]

    return track_indices, track_names


def load_model_for_inference(
    checkpoint_path: str,
    config: dict,
    device: torch.device,
) -> Tuple:
    """
    Load Borzoi model for inference.

    Returns:
        (model, is_human) tuple
    """
    from infer_borzoi_pytorch import load_model
    from src.utils import DotDict

    # Resolve checkpoint path (handle Accelerate directories)
    ckpt_path = Path(checkpoint_path)
    if not ckpt_path.exists():
        if ckpt_path.suffix:
            candidate = ckpt_path.with_suffix('')
            if candidate.exists():
                ckpt_path = candidate
        if not ckpt_path.exists():
            best_dir = ckpt_path.parent / "best_model"
            if best_dir.exists():
                ckpt_path = best_dir

    if not ckpt_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

    model, is_human = load_model(
        str(ckpt_path),
        DotDict(config),
        device,
    )

    return model, is_human


def main():
    """Main entry point."""
    t_start = time.time()
    args = parse_args()

    # Setup distributed (if enabled)
    if args.distributed:
        rank, world_size, device = setup_distributed()
    else:
        rank = 0
        world_size = 1
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # Load config
    print_rank0(f"Loading config from {args.config}", rank)
    config = load_config(args.config)

    genome_fasta = config['data']['genome_fasta']
    chrom_sizes_file = config['data']['chrom_sizes']

    print_rank0(f"  Genome: {genome_fasta}", rank)
    print_rank0(f"  Chrom sizes: {chrom_sizes_file}", rank)

    chrom_sizes = load_chrom_sizes(chrom_sizes_file)

    # Resolve checkpoint path from config if not provided
    checkpoint_path = args.checkpoint
    if checkpoint_path is None:
        # Try to resolve from config
        checkpoint_dir = config.get('logging', {}).get('checkpoint_dir', './checkpoints')

        if args.metric is not None:
            # User requested checkpoint by metric
            from infer_borzoi_pytorch import find_checkpoint_by_metric
            metric_checkpoint = find_checkpoint_by_metric(checkpoint_dir, args.metric, rank)
            if metric_checkpoint is not None:
                checkpoint_path = metric_checkpoint
            else:
                # Fallback to best_model
                print_rank0(f"Metric {args.metric} not found, falling back to best_model", rank)
                best_model_path = Path(checkpoint_dir) / 'best_model'
                if best_model_path.exists():
                    checkpoint_path = str(best_model_path)

        if checkpoint_path is None:
            # Try inference.checkpoint from config
            inf_cfg = config.get('inference', {})
            checkpoint_val = inf_cfg.get('checkpoint', 'best_model')

            if checkpoint_val in ['best_model', 'last']:
                # Resolve relative to checkpoint_dir
                for ext in ['', '.pt', '.pth']:
                    candidate = Path(checkpoint_dir) / f'{checkpoint_val}{ext}'
                    if candidate.exists():
                        checkpoint_path = str(candidate)
                        break
            else:
                # Direct path in config
                checkpoint_path = checkpoint_val

        if checkpoint_path is None:
            raise ValueError(
                "No checkpoint found. Provide --checkpoint, --metric, or set inference.checkpoint in config."
            )

        print_rank0(f"  Resolved checkpoint from config: {checkpoint_path}", rank)

    # Resolve tracks
    track_indices, track_names = resolve_tracks(args, config)
    print_rank0(f"  Analyzing {len(track_names)} track(s): {track_names[:5]}{'...' if len(track_names) > 5 else ''}", rank)

    # Load model
    print_rank0(f"\nLoading model from {checkpoint_path}", rank)
    model, is_human = load_model_for_inference(checkpoint_path, config, device)
    model.eval()

    # Get model params from config
    seq_len = config.get('model', {}).get('seq_len', 524288)
    pred_len = config.get('model', {}).get('pred_len', 196608)
    bin_size = config.get('model', {}).get('bin_size', 32)

    print_rank0(f"  Model params: seq_len={seq_len}, pred_len={pred_len}, bin_size={bin_size}", rank)

    # Load tables
    print_rank0(f"\nLoading input tables:", rank)
    print_rank0(f"  Regions: {args.regions}", rank)
    regions = load_regions(args.regions, normalize_to=args.normalize_regions)
    print_rank0(f"    Loaded {len(regions)} regions", rank)

    splice_perts = []
    if args.splice_pert:
        print_rank0(f"  Splice perturbations: {args.splice_pert}", rank)
        splice_perts = load_splice_pert(args.splice_pert)
        print_rank0(f"    Loaded {len(splice_perts)} splice operations", rank)

    seq_perts = []
    if args.seq_pert:
        print_rank0(f"  Sequence perturbations: {args.seq_pert}", rank)
        seq_perts = load_seq_pert(args.seq_pert)
        print_rank0(f"    Loaded {len(seq_perts)} sequence perturbations", rank)

    # Create genome cache
    print_rank0(f"\nInitializing genome cache with {args.num_workers * 2} handles", rank)
    genome_cache = GenomeCache(
        genome_fasta=genome_fasta,
        num_handles=args.num_workers * 2,
    )

    # Create analyzer
    print_rank0(f"\nCreating analyzer", rank)
    analyzer = SimplePerturbationAnalyzer(
        model=model,
        genome_cache=genome_cache,
        chrom_sizes=chrom_sizes,
        device=device,
        seq_len=seq_len,
        pred_len=pred_len,
        model_bin_size=bin_size,
        is_human=is_human,
        use_rc_average=not args.no_rc_average,
        mixed_precision=not args.no_mixed_precision,
        use_compile=args.use_compile,
        compile_mode=args.compile_mode,
        seed=args.seed,
    )
    
    t_init = time.time()

    # Resolve data_matrix from config if not provided
    data_matrix_path = args.data_matrix
    if data_matrix_path is None:
        data_matrix_path = config.get('data', {}).get('data_matrix')
        if data_matrix_path:
            print_rank0(f"  Using data_matrix from config: {data_matrix_path}", rank)

    # Compute
    print_rank0(f"\nComputing perturbation effects", rank)
    df = analyzer.compute(
        regions=regions,
        splice_perts=splice_perts,
        seq_perts=seq_perts,
        track_indices=track_indices,
        track_names=track_names,
        data_matrix_path=data_matrix_path,
        inference_batch_size=args.inference_batch_size,
        num_workers=args.num_workers,
        show_progress=not args.quiet,
        profile=args.profile,
    )
    
    t_compute = time.time()

    # Save results (only rank 0 writes to avoid corruption)
    if rank == 0:
        print_rank0(f"\nSaving results to {args.output}", rank)
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(args.output, index=False)
    
    # Synchronize before timing/printing (so rank 0 waits if it finished early)
    if args.distributed:
        torch.distributed.barrier()

    t_write = time.time()

    if rank == 0:
        print(f"  Shape: {df.shape}")
        print(f"  Columns ({len(df.columns)}): {list(df.columns)[:10]}{'...' if len(df.columns) > 10 else ''}")

    if args.profile and rank == 0:
        print(f"\n--- Top-level Timing ---")
        print(f"Initialization: {t_init - t_start:.2f}s")
        print(f"Output writing: {t_write - t_compute:.2f}s")
        print(f"Total Wall Time: {t_write - t_start:.2f}s")
        print("------------------------")
        
    print_rank0(f"\nDone!", rank)

    # Cleanup
    genome_cache.close()

    if args.distributed:
        cleanup_distributed()


if __name__ == "__main__":
    main()
