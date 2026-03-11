"""Configuration utilities for Borzoi fine-tuning."""

import yaml
import argparse
from typing import Any, Dict, List, Optional
from pathlib import Path
import copy


class DotDict(dict):
    """Dictionary with dot notation for nested keys."""

    def __init__(self, d: Dict[str, Any]):
        super().__init__()
        for key, value in d.items():
            if isinstance(value, dict):
                value = DotDict(value)
            self[key] = value

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"Config has no attribute '{key}'")

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"Config has no attribute '{key}'")


class Config:
    """Configuration management with YAML loading and CLI overrides.

    This class provides a hierarchical configuration system for Borzoi fine-tuning.
    Configurations can be loaded from YAML files and overridden via command-line arguments.

    Example YAML configuration:
        data:
          genome_fasta: /path/to/genome.fa
          chrom_sizes: /path/to/chrom.sizes
          data_matrix: /path/to/coverage.parquet
        model:
          name: borzoi_mouse_rep0
          freeze_trunk: false
        training:
          batch_size: 2
          epochs: 50
    """

    # Default configuration with comprehensive documentation
    DEFAULT_CONFIG = {
        # =========================================================================
        # DATA CONFIGURATION
        # =========================================================================
        'data': {
            'genome_fasta': None,  # REQUIRED: Path to reference genome FASTA file
            'genomes': None,  # Multi-genome list: [{'fasta': path, 'name': str, 'samples_per_region': int}]
            'chrom_sizes': None,   # REQUIRED: Path to chromosome sizes file (tab-separated: chrom\tsize)
            # Target data source (exactly one must be specified):
            'data_matrix': None,   # Path to pre-computed matrix (NPZ, Parquet, or CSV)
            'bigwig_files': [],    # List of BigWig files with track info: [{'path': ..., 'track': ...}]
            'bigwig_folder': None, # Directory containing *.bw or *.bigWig files
            'additional_data_matrices': [],  # Additional data at different resolutions: [{'path': ..., 'track_col': ...}]
            # Region specification:
            'regions_file': None,  # BED file with training regions (optional, defaults to genome-wide)
            'train_regions': None, # Pre-split training regions BED (when split_by_chrom=False)
            'val_regions': None,   # Pre-split validation regions BED (when split_by_chrom=False)
            'test_regions': None,  # Pre-split test regions BED (optional)
            'val_chroms': ['chr8', 'chr10'],   # Validation chromosomes (when split_by_chrom=True)
            'test_chroms': ['chr9', 'chr18'],  # Test chromosomes (when split_by_chrom=True)
            'split_by_chrom': True,            # If True, split data by chromosome; if False, use pre-split files
        },
        # =========================================================================
        # MODEL CONFIGURATION
        # =========================================================================
        'model': {
            'name': 'Borzoi_mouse_rep0',  # Model identifier (HuggingFace model ID or local path)
            'species': None,     # Explicit species: 'mouse' or 'human'. If None, auto-detected from genome path.
                                 # Set explicitly to avoid relying on path-based detection which can be fragile.
            'seq_len': 524288,   # Input sequence length in bp (must be power of 2 for efficiency)
            'pred_len': None,    # Prediction window in bp (None = auto: seq_len * 3/8)
            'bin_size': 32,      # Bin size in bp for model output
            'freeze_trunk': False,      # If True, only train the head (linear probe)
            'train_from_scratch': False, # If True, initialize model with random weights instead of loading pretrained
            'num_output_tasks': None,   # Number of output tracks (auto-detected from data if None)
            # High-resolution decoder (U-Net extension beyond 32bp):
            'hires_resolution': None,          # Target resolution in bp: 1, 2, 4, 8, or 16. None = standard 32bp Borzoi.
            'hires_decoder_channels': None,    # Custom decoder channel widths (list). None = auto [512, 256, 128, 64, 32].
        },
        # =========================================================================
        # TRAINING CONFIGURATION
        # =========================================================================
        'training': {
            'batch_size': 2,      # Batch size per GPU
            'epochs': 50,         # Total training epochs
            'learning_rate': 6e-5,      # Base learning rate
            'optimizer': 'adamw',       # Optimizer: 'adamw' (recommended for transformers)
            # Weight decay settings (different for trunk vs transformer for stability):
            'weight_decay': 1e-8,              # Weight decay for convolutional trunk layers
            'weight_decay_transformer': 2e-8,  # Higher decay for transformer layers (more prone to overfitting)
            # Loss configuration:
            'loss': 'poisson_multinomial',    # Loss function: 'poisson', 'poisson_multinomial', 'mse', 'mse_log', 'cosine_mse', 'crested'
            'loss_center_crop_bp': None,      # Center crop length in bp for loss computation (reduces edge effects)
            'loss_center_crop_bins': None,    # Alternative: specify crop in bins directly
            # Mixed precision and memory:
            'mixed_precision': True,          # Use FP16 mixed precision (recommended for memory efficiency)
            'gradient_checkpointing': False,  # Trade compute for memory (slower but fits larger batches)
            'compile_model': False,           # Use torch.compile (experimental, may speed up training)
            'clip_grad_norm': 0.2,            # Gradient clipping threshold (prevents exploding gradients)
            'warmup_steps': 1000,             # Learning rate warmup steps
            'num_workers': 4,                 # DataLoader worker processes
            'seed': 42,                       # Random seed for reproducibility
            # Two-stage training (linear probe + fine-tune):
            'head_learning_rate': 0.0001,         # Learning rate for head during linear probe stage
            'finetune_learning_rate': 0.00005,    # Learning rate for full model during fine-tune stage
            'linear_probe_steps': None,           # Steps for linear probe (if None, uses linear_probe_epochs)
            'linear_probe_epochs': 1,             # Epochs for linear probe stage before unfreezing trunk
            'linear_probe_early_stopping_patience': None,  # Early stop linear probe to jump to fine-tuning
            'force_finetune_on_resume': False,    # If True, override checkpoint stage and resume in finetune
            # Early stopping:
            'early_stopping_patience': None,  # Evaluations without improvement before stopping (None=disabled)
            'early_stopping_metric': 'val_loss',  # Metric key to monitor (e.g., 'val_genome_wide_pearson')
            'early_stopping_mode': 'min',     # 'min' for loss, 'max' for correlation metrics
            'min_improvement': 1e-4,          # Minimum metric improvement to reset patience counter
            # Best model selection:
            'best_model_metric': 'val_loss',  # Metric for primary best_model checkpoint
            'best_model_mode': None,          # 'min' or 'max', auto-detected from metric if None
            # Learning rate scheduling:
            'reduce_lr_on_plateau': False,    # Reduce LR when validation plateaus
            'lr_patience': 5,                 # Evaluations to wait before reducing LR
            'lr_factor': 0.5,                 # Factor to multiply LR by on reduction
            'lr_min': 1e-7,                   # Minimum learning rate floor
            # OOM (Out-Of-Memory) handling:
            'auto_batch_size': True,              # Enable automatic batch size reduction on OOM
            'min_batch_size': 1,                  # Minimum batch size before giving up
            'batch_size_reduction_factor': 0.5,   # Factor to reduce batch size by on OOM
            'use_gradient_accumulation': True,    # Maintain effective batch size via accumulation
            'max_gradient_accumulation_steps': 16, # Maximum accumulation steps
            'oom_recovery_delay': 2.0,            # Seconds to wait after OOM before retrying
            'max_batch_size': None,               # Upper bound for auto-scaling (defaults to initial batch_size)
            'enable_auto_batch_growth': True,     # Try to grow batch size back up after stable steps
            'batch_size_growth_factor': 1.25,     # Growth factor when scaling up
            'batch_growth_interval': 400,         # Steps without OOM before attempting to grow
        },
        # =========================================================================
        # DATA AUGMENTATION
        # =========================================================================
        'augmentation': {
            'reverse_complement': True,   # Random reverse complement augmentation
            'shift_max': 3,               # Maximum random shift in bins (for jittering)
            'squash_targets': True,       # Apply Borzoi's squash transformation to targets
            'tile_stride': None,          # Stride for tiling regions (None = pred_len, non-overlapping)
            'use_random_crop': False,     # If True, randomly crop within regions instead of fixed tiling
            # 'samples_per_region': 1,    # Number of random crops per region per epoch (only with use_random_crop)
        },
        # =========================================================================
        # PERFORMANCE SETTINGS
        # =========================================================================
        'performance': {
            'lazy_loading': True,      # Load BigWig data on-the-fly (saves memory, slightly slower)
            'cache_file': None,        # Path to cache pre-computed BigWig arrays (speeds up repeated runs)
            'data_matrix_cache': None, # Directory for memory-mapped cache of dense matrices
            'data_matrix_mmap': False, # If True, memory-map cached .npy arrays
            'parquet_streaming_batch_size': 1000000, # Rows per Parquet batch when building streaming cache
        },
        # =========================================================================
        # DISTRIBUTED TRAINING
        # =========================================================================
        'distributed': {
            'enabled': True,   # Enable distributed training via Accelerate
            'num_gpus': 1,     # Number of GPUs to use (auto-detected by Accelerate)
        },
        # =========================================================================
        # LOGGING AND CHECKPOINTING
        # =========================================================================
        'logging': {
            'use_wandb': False,               # Enable Weights & Biases logging
            'wandb_project': 'borzoi-finetune', # W&B project name
            'wandb_run_name': None,           # W&B run name (auto-generated if None)
            'checkpoint_dir': './checkpoints', # Directory for saving checkpoints
            'save_best_only': True,           # Only save best model (by validation loss)
            'eval_every_n': 256,              # Evaluate every N steps
            'save_every_n': 5000,             # Save checkpoint every N steps
            'log_every_n': 10,                # Log metrics every N steps
            # Multi-metric best model tracking:
            'track_best_metrics': ['val_loss'],  # Metrics to track best checkpoints for (creates best_model_{metric}/)
            'save_metrics_json': True,           # Save metrics.json with each best checkpoint
        },
        # =========================================================================
        # VALIDATION SETTINGS
        # =========================================================================
        'validation': {
            'region_metrics': {
                'enabled': False,        # Enable region-based R² metrics
                'bed_file': None,        # BED file with regions for aggregation (e.g., gene bodies, enhancers)
                'aggregation_mode': 'both',  # 'per_track' (R² per track), 'aggregate' (all tracks), or 'both'
                'overlap_mode': 'strict',    # 'strict' (bin must be fully contained in region) or 'partial' (weighted by overlap)
            }
        },
        # =========================================================================
        # POST-TRAINING INFERENCE
        # =========================================================================
        'inference': {
            'run_after_training': False,    # Run genome-wide inference after training completes
            'mode': 'genome',               # 'genome' (all chromosomes), 'chromosomes' (specific list), or 'bed' (regions file)
            'chromosomes': None,            # List of chromosomes for mode='chromosomes'
            'bed_file': None,               # BED file path for mode='bed'
            'exclude_chroms': ['chrM'],     # Chromosomes to exclude for mode='genome'
            'batch_size': None,             # Inference batch size (defaults to training.batch_size)
            'output_dir': None,             # Output directory (defaults to checkpoint_dir/predictions)
            'output_format': ['bigwig'],    # Output formats: 'bigwig', 'hdf5', 'npz', 'parquet'
            'bigwig_jobs': 4,               # Parallel jobs for BigWig writing
            'track_names': None,            # Track names for output files (auto-detected if None)
            'num_workers': None,            # DataLoader workers (defaults to training.num_workers)
            'stride': None,                 # Inference stride (defaults to pred_len for non-overlapping)
            'use_rc_average': True,         # Average forward and reverse complement predictions
            'mixed_precision': True,        # Use AMP for faster inference
            'checkpoint': 'best_model',     # Which checkpoint to use: 'best_model', 'last', or explicit path
            'defer_postprocess': False,     # If True, save NPZ only and generate finalize_inference.sh for later CPU postprocessing
        },
        # =========================================================================
        # MISHA TRACK CREATION (Weizmann-specific genomic database)
        # =========================================================================
        'misha_track': {
            'enabled': False,    # Enable Misha track creation from BigWig predictions
            'groot': None,       # REQUIRED: Misha genome root directory (e.g., '/data/genomes/mm10')
                                 # This directory should contain the Misha genome database
            'track_prefix': None, # REQUIRED: Full track path using dots as directory separators
                                 # Example: 'seq.predictions.borzoi.pcg_524k' creates:
                                 # {groot}/tracks/seq/predictions/borzoi/pcg_524k/
            'track_suffix': '',  # Optional suffix appended to track_prefix
            'track_name_prefix': '',  # Prefix for individual track names (e.g., 'pred_')
            'track_name_suffix': '',  # Suffix for individual track names (e.g., '_524k')
            'description': None,      # Track description (auto-generated from config if None)
            'bw_dir': None,           # BigWig source directory (defaults to inference.output_dir)
            'track_names': None,      # Track names (auto-detected from .bw files if None)
            'binsize': 20,            # Misha track bin size for intervals
        },
    }

    def __init__(self, config_dict: Optional[Dict[str, Any]] = None, source_path: Optional[Path] = None):
        self._config = copy.deepcopy(self.DEFAULT_CONFIG)
        if config_dict is not None:
            self._deep_merge(self._config, config_dict)
        self._dotdict = DotDict(self._config)
        self._source_path = source_path  # Original config file path, if loaded from yaml

    @property
    def source_path(self) -> Optional[Path]:
        """Return the original config file path, if loaded from yaml."""
        return self._source_path

    @property
    def source_name(self) -> Optional[str]:
        """Return just the filename of the original config file."""
        return self._source_path.name if self._source_path else None

    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'Config':
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Config file not found: {yaml_path}")

        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)

        return cls(config_dict, source_path=yaml_path)

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'Config':
        if hasattr(args, 'config') and args.config is not None:
            config = cls.from_yaml(args.config)
        else:
            config = cls()

        overrides = {
            'batch_size': 'training.batch_size',
            'learning_rate': 'training.learning_rate',
            'epochs': 'training.epochs',
            'seq_len': 'model.seq_len',
            'pred_len': 'model.pred_len',
            'freeze_trunk': 'model.freeze_trunk',
            'train_from_scratch': 'model.train_from_scratch',
            'use_wandb': 'logging.use_wandb',
            'wandb_project': 'logging.wandb_project',
            'checkpoint_dir': 'logging.checkpoint_dir',
            'num_workers': 'training.num_workers',
            'force_finetune_on_resume': 'training.force_finetune_on_resume',
        }

        for arg_name, config_path in overrides.items():
            if hasattr(args, arg_name) and getattr(args, arg_name) is not None:
                config.set(config_path, getattr(args, arg_name))

        return config

    def get(self, key_path: str, default: Any = None) -> Any:
        keys = key_path.split('.')
        value = self._config

        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default

        return value

    def set(self, key_path: str, value: Any) -> None:
        keys = key_path.split('.')
        target = self._config

        for key in keys[:-1]:
            if key not in target:
                target[key] = {}
            target = target[key]

        target[keys[-1]] = value
        self._dotdict = DotDict(self._config)

    def __getattr__(self, name: str) -> Any:
        if name.startswith('_'):
            return super().__getattribute__(name)
        return getattr(self._dotdict, name)

    def to_dict(self) -> Dict[str, Any]:
        return copy.deepcopy(self._config)

    def get_inference_configs(self) -> List[Dict[str, Any]]:
        """Return inference configs as a list (handles both dict and list formats).

        Supports backward compatibility:
        - If inference is a dict, wraps it in a list
        - If inference is a list, returns as-is (filtering out None entries)
        - If inference is empty/missing, returns empty list

        Returns:
            List of inference configuration dicts.
        """
        inference_cfg = self.get('inference', {})
        if isinstance(inference_cfg, list):
            # Filter out None entries that could come from YAML nulls
            return [cfg for cfg in inference_cfg if cfg is not None and isinstance(cfg, dict)]
        elif isinstance(inference_cfg, dict):
            return [inference_cfg] if inference_cfg else []
        return []

    def save(self, yaml_path: str) -> None:
        yaml_path = Path(yaml_path)
        yaml_path.parent.mkdir(parents=True, exist_ok=True)

        with open(yaml_path, 'w') as f:
            yaml.dump(self._config, f, default_flow_style=False, sort_keys=False)

    def validate(self, check_file_exists: bool = True) -> None:
        """Validate required fields, source settings, and optionally file existence.

        Args:
            check_file_exists: If True, verify that required files exist on disk.
                              Set to False for unit tests or when paths are placeholders.
        """
        genomes = self.get('data.genomes')
        if genomes is not None:
            if not isinstance(genomes, list) or len(genomes) == 0:
                raise ValueError("data.genomes must be a non-empty list")
            for i, g in enumerate(genomes):
                if not isinstance(g, dict) or 'fasta' not in g:
                    raise ValueError(f"data.genomes[{i}] must be a dict with at least a 'fasta' key")
            # Auto-populate genome_fasta from first genome for downstream compatibility
            if self.get('data.genome_fasta') is None:
                self.set('data.genome_fasta', genomes[0]['fasta'])
            required = ['data.chrom_sizes']
        else:
            required = ['data.genome_fasta', 'data.chrom_sizes']

        for field in required:
            if self.get(field) is None:
                raise ValueError(f"Required field '{field}' is not set")

        # Validate required files exist if requested
        if check_file_exists:
            self._validate_file_paths()

        if self.data.split_by_chrom:
            if not self.data.val_chroms:
                raise ValueError("When split_by_chrom=True, must provide data.val_chroms")
        else:
            if not self.data.train_regions:
                raise ValueError("When split_by_chrom=False, must provide data.train_regions")
            if not self.data.val_regions:
                raise ValueError("When split_by_chrom=False, must provide data.val_regions")

        has_data_matrix = bool(hasattr(self.data, 'data_matrix') and self.data.data_matrix)
        has_bigwig_files = bool(self.data.bigwig_files and len(self.data.bigwig_files) > 0)
        has_bigwig_folder = bool(self.data.bigwig_folder)

        num_sources = sum([has_data_matrix, has_bigwig_files, has_bigwig_folder])

        if num_sources == 0:
            raise ValueError("Must provide one of: data.data_matrix, data.bigwig_files, or data.bigwig_folder")

        if num_sources > 1:
            raise ValueError("Cannot specify multiple data sources")

        if self.model.num_output_tasks is None:
            if has_data_matrix:
                print("  Note: num_output_tasks will be auto-detected from data_matrix")
                self.set('model.num_output_tasks', -1)
            elif has_bigwig_files:
                self.set('model.num_output_tasks', len(self.data.bigwig_files))
            elif has_bigwig_folder:
                from pathlib import Path
                folder = Path(self.data.bigwig_folder)
                if not folder.exists():
                    raise ValueError(f"BigWig folder not found: {folder}")
                bw_files = list(folder.glob("*.bw")) + list(folder.glob("*.bigWig"))
                if not bw_files:
                    raise ValueError(f"No BigWig files found in {folder}")
                self.set('model.num_output_tasks', len(bw_files))

        # Validate seq_len is compatible with 128x encoder downsampling
        if self.model.seq_len % 128 != 0:
            raise ValueError(
                f"seq_len ({self.model.seq_len}) must be divisible by 128 "
                f"(required by the encoder's 128x downsampling)"
            )

        # Auto-scale pred_len if not set
        if self.model.pred_len is None:
            auto_pred_len = self.model.seq_len * 3 // 8
            self.set('model.pred_len', auto_pred_len)

        # Validate sequence lengths
        if self.model.pred_len >= self.model.seq_len:
            raise ValueError(f"pred_len ({self.model.pred_len}) must be < seq_len ({self.model.seq_len})")

        if (self.model.seq_len - self.model.pred_len) % 2 != 0:
            raise ValueError("(seq_len - pred_len) must be even")

        if self.model.pred_len % self.model.bin_size != 0:
            raise ValueError(f"pred_len must be divisible by bin_size")

        # Validate optional loss cropping
        crop_bp = self.get('training.loss_center_crop_bp', None)
        crop_bins = self.get('training.loss_center_crop_bins', None)

        if crop_bp is not None and crop_bp <= 0:
            raise ValueError("training.loss_center_crop_bp must be positive when set")

        if crop_bins is not None and crop_bins <= 0:
            raise ValueError("training.loss_center_crop_bins must be positive when set")

        if crop_bp is not None and crop_bp % self.model.bin_size != 0:
            raise ValueError(
                f"loss_center_crop_bp ({crop_bp}) must be divisible by bin_size ({self.model.bin_size})"
            )

        resolved_crop_bp = None
        if crop_bp is not None:
            resolved_crop_bp = crop_bp
        if crop_bins is not None:
            resolved_from_bins = crop_bins * self.model.bin_size
            if resolved_crop_bp is not None and resolved_crop_bp != resolved_from_bins:
                raise ValueError(
                    "loss_center_crop_bp and loss_center_crop_bins are both set but do not match "
                    f"({resolved_crop_bp} bp vs {resolved_from_bins} bp)"
                )
            resolved_crop_bp = resolved_from_bins

        if resolved_crop_bp is not None:
            if resolved_crop_bp > self.model.pred_len:
                raise ValueError(
                    f"loss_center_crop ({resolved_crop_bp} bp) cannot exceed pred_len ({self.model.pred_len} bp)"
                )

    def _validate_file_paths(self) -> None:
        """Validate that required file paths exist on disk."""
        # Required files that must exist
        required_files = [
            ('data.genome_fasta', self.get('data.genome_fasta')),
            ('data.chrom_sizes', self.get('data.chrom_sizes')),
        ]

        # Multi-genome FASTA files
        genomes = self.get('data.genomes')
        if genomes:
            for i, g in enumerate(genomes):
                required_files.append((f'data.genomes[{i}].fasta', g['fasta']))

        # Optional files that should exist if specified
        optional_files = [
            ('data.data_matrix', self.get('data.data_matrix')),
            ('data.regions_file', self.get('data.regions_file')),
            ('data.train_regions', self.get('data.train_regions')),
            ('data.val_regions', self.get('data.val_regions')),
            ('data.test_regions', self.get('data.test_regions')),
        ]

        missing_files = []

        for field_name, file_path in required_files:
            if file_path is not None:
                path = Path(file_path)
                if not path.exists():
                    missing_files.append((field_name, file_path))

        for field_name, file_path in optional_files:
            if file_path is not None:
                path = Path(file_path)
                if not path.exists():
                    missing_files.append((field_name, file_path))

        # Check bigwig folder if specified
        bw_folder = self.get('data.bigwig_folder')
        if bw_folder is not None:
            folder = Path(bw_folder)
            if not folder.exists():
                missing_files.append(('data.bigwig_folder', bw_folder))
            elif not folder.is_dir():
                raise ValueError(f"data.bigwig_folder is not a directory: {bw_folder}")

        # Check individual bigwig files if specified
        bw_files = self.get('data.bigwig_files', [])
        if bw_files:
            for i, bw_entry in enumerate(bw_files):
                if isinstance(bw_entry, dict) and 'path' in bw_entry:
                    bw_path = bw_entry['path']
                    if not Path(bw_path).exists():
                        missing_files.append((f'data.bigwig_files[{i}].path', bw_path))

        if missing_files:
            error_lines = ["The following required files do not exist:"]
            for field_name, file_path in missing_files:
                error_lines.append(f"  - {field_name}: {file_path}")
            error_lines.append("\nPlease verify file paths in your configuration.")
            raise FileNotFoundError('\n'.join(error_lines))

    @staticmethod
    def _deep_merge(base: Dict, override: Dict) -> None:
        for key, value in override.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                Config._deep_merge(base[key], value)
            else:
                base[key] = value

    def __repr__(self) -> str:
        return f"Config({yaml.dump(self._config, default_flow_style=False)})"


def get_default_config() -> Config:
    return Config()


def create_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Fine-tune Borzoi model on custom genomic data (PyTorch + Accelerate)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--config', type=str, default=None,
                        help='Path to YAML configuration file')

    parser.add_argument('--genome_fasta', type=str, help='Path to genome FASTA file')
    parser.add_argument('--chrom_sizes', type=str, help='Path to chromosome sizes file')
    parser.add_argument('--train_regions', type=str, help='Path to training regions BED')
    parser.add_argument('--val_regions', type=str, help='Path to validation regions BED')

    parser.add_argument('--seq_len', type=int, help='Input sequence length')
    parser.add_argument('--pred_len', type=int, help='Prediction window length')
    parser.add_argument('--freeze_trunk', action='store_true', default=None,
                        help='Freeze trunk, train only head')
    parser.add_argument('--train_from_scratch', action='store_true', default=None,
                        help='Initialize model with random weights instead of loading pretrained')

    parser.add_argument('--batch_size', type=int, help='Batch size per GPU')
    parser.add_argument('--epochs', type=int, help='Number of training epochs')
    parser.add_argument('--learning_rate', type=float, help='Learning rate')
    parser.add_argument('--num_workers', type=int, help='DataLoader workers')

    parser.add_argument('--use_wandb', action='store_true', default=None,
                        help='Enable Weights & Biases logging')
    parser.add_argument('--wandb_project', type=str, help='W&B project name')
    parser.add_argument('--checkpoint_dir', type=str, help='Directory for checkpoints')

    parser.add_argument('--resume', action='store_true', help='Resume from latest checkpoint')
    parser.add_argument('--checkpoint_path', type=str, help='Explicit checkpoint path to resume from')
    parser.add_argument('--force_finetune_on_resume', action='store_true', default=None,
                        help='Force full fine-tuning when resuming (override checkpoint stage)')

    return parser
