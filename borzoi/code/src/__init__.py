"""
Borzoi Fine-Tuning Library (PyTorch + Accelerate version)

Modules:
- config: Configuration management
- data_loaders: Data loading and preprocessing
- losses: Custom loss functions
- metrics: Custom metrics
"""

__version__ = "0.1.0"

from .config import Config, get_default_config, create_arg_parser
from .data_loaders import (
    load_chrom_sizes,
    load_regions,
    load_multiple_bigwigs,
    save_binned_bigwigs,
    load_binned_bigwigs,
    load_from_data_matrix,
    split_regions_by_chromosome,
    BorzoiDataset,
)
from .losses import (
    get_loss,
    poisson_loss,
    poisson_multinomial_torch,
    mse_loss,
    mse_log_loss,
    cosine_mse_log_loss,
    cosine_mse_raw_loss,
    pearson_loss,
    PoissonMultinomialLoss,
    PoissonLoss,
    MSELogLoss,
    CosineMSELogLoss,
    set_nan_checking_enabled,
    is_nan_checking_enabled,
)
from .metrics import (
    get_metrics,
    PearsonCorrelation,
    SpearmanCorrelation,
    R2Score,
)
from .utils import (
    one_hot_encode_torch,
    reverse_complement_torch,
    squash_coverage_torch,
    unsquash_coverage_torch,
    count_parameters,
    freeze_layers,
    set_seed,
    predict_batch,
    predict_with_rc_average,
)

__all__ = [
    # Config
    "Config",
    "get_default_config",
    "create_arg_parser",
    # Data loaders
    "load_chrom_sizes",
    "load_regions",
    "load_multiple_bigwigs",
    "split_regions_by_chromosome",
    "BorzoiDataset",
    # Losses
    "get_loss",
    "poisson_loss",
    "poisson_multinomial_torch",
    "mse_loss",
    "mse_log_loss",
    "cosine_mse_log_loss",
    "PoissonMultinomialLoss",
    "PoissonLoss",
    "MSELogLoss",
    "CosineMSELogLoss",
    "set_nan_checking_enabled",
    "is_nan_checking_enabled",
    # Metrics
    "get_metrics",
    "PearsonCorrelation",
    "SpearmanCorrelation",
    "R2Score",
    # Utils
    "one_hot_encode_torch",
    "reverse_complement_torch",
    "squash_coverage_torch",
    "unsquash_coverage_torch",
    "count_parameters",
    "freeze_layers",
    "set_seed",
    "predict_batch",
    "predict_with_rc_average",
]
