"""
Utility functions for Borzoi fine-tuning with PyTorch.
Includes functions adapted from the official training-borzoi repository.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import yaml
from typing import Tuple, Optional, Dict, List, Any


# =============================================================================
# Configuration Utilities
# =============================================================================

class DotDict(dict):
    """
    Dictionary with dot notation access.

    Enables accessing nested dictionary values using dot notation:
        config = DotDict({'model': {'name': 'borzoi'}})
        print(config.model.name)  # 'borzoi'

    Also supports nested get() with dot-separated keys:
        config.get('model.name', 'default')
    """

    def __getattr__(self, key: str) -> Any:
        try:
            value = self[key]
            if isinstance(value, dict) and not isinstance(value, DotDict):
                return DotDict(value)
            return value
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{key}'")

    def __setattr__(self, key: str, value: Any) -> None:
        self[key] = value

    def __delattr__(self, key: str) -> None:
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{key}'")

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a value with support for dot-separated keys.

        Args:
            key: Key to look up. Can be dot-separated for nested access (e.g., 'model.name').
            default: Default value if key is not found.

        Returns:
            The value if found, otherwise default.
        """
        if '.' not in key:
            return super().get(key, default)

        # Handle dot-separated keys
        keys = key.split('.')
        current = self
        for k in keys:
            if isinstance(current, dict):
                current = current.get(k)
                if current is None:
                    return default
            else:
                return default
        return current


def load_yaml_config(config_path: str) -> DotDict:
    """
    Load configuration from YAML file.

    Args:
        config_path: Path to YAML configuration file.

    Returns:
        Configuration as DotDict for dot notation access.
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return DotDict(config if config else {})

# Numerical stability constant for power operations on near-zero values
SQUASH_EPSILON = 1e-6

# Threshold for squash transformation where sqrt scaling kicks in.
# Borzoi's squash applies y^0.75, then sqrt(excess) for values above threshold.
# 384.0 corresponds to raw coverage of ~8000, above which additional compression helps.
SQUASH_THRESHOLD = 384.0

# Precompute ASCII -> base index mapping once, cache per device for reuse
_BYTE_TO_INDEX_CPU = torch.full((256,), -1, dtype=torch.long)
for base, idx in (('A', 0), ('C', 1), ('G', 2), ('T', 3)):
    _BYTE_TO_INDEX_CPU[ord(base)] = idx
    _BYTE_TO_INDEX_CPU[ord(base.lower())] = idx

_BYTE_TO_INDEX_CACHE: Dict[Optional[torch.device], torch.Tensor] = {None: _BYTE_TO_INDEX_CPU}


def _get_byte_mapping(device: Optional[torch.device]) -> torch.Tensor:
    if device in _BYTE_TO_INDEX_CACHE:
        return _BYTE_TO_INDEX_CACHE[device]

    mapping = _BYTE_TO_INDEX_CPU.to(device)
    _BYTE_TO_INDEX_CACHE[device] = mapping
    return mapping


# =============================================================================
# Sequence Utilities
# =============================================================================

def one_hot_encode_torch(sequence: str, device: torch.device = None) -> torch.Tensor:
    """
    Convert a DNA sequence to one-hot encoding as a PyTorch tensor.
    
    Args:
        sequence: DNA string (A, C, G, T, N)
        device: Target device for the tensor
        
    Returns:
        Tensor of shape (len(sequence), 4)
    """
    mapping = _get_byte_mapping(device)

    seq_bytes = torch.tensor(np.frombuffer(sequence.encode('ascii'), dtype=np.uint8), dtype=torch.long, device=device)
    idxs = mapping[seq_bytes]

    valid_mask = idxs >= 0
    encoded = torch.zeros(seq_bytes.numel(), 4, dtype=torch.float32, device=device)
    if valid_mask.any():
        encoded[valid_mask] = F.one_hot(idxs[valid_mask], num_classes=4).float()

    return encoded


def reverse_complement_torch(seq_tensor: torch.Tensor) -> torch.Tensor:
    """
    Reverse complement a one-hot encoded sequence tensor.
    
    Args:
        seq_tensor: Tensor of shape (..., seq_len, 4)
        
    Returns:
        Reverse complemented tensor of same shape
    """
    # Reverse along sequence dimension and swap A<->T, C<->G
    return seq_tensor.flip(-2)[..., [3, 2, 1, 0]]


# =============================================================================
# Coverage Transformations
# =============================================================================

def squash_coverage_torch(y: torch.Tensor, threshold: float = SQUASH_THRESHOLD) -> torch.Tensor:
    """
    Apply Borzoi's squash transformation to coverage values.

    y_squashed = y^(3/4) if y^(3/4) <= threshold
                 else threshold + sqrt(y^(3/4) - threshold)

    Args:
        y: Raw coverage values
        threshold: Threshold for additional sqrt transform (see SQUASH_THRESHOLD)

    Returns:
        Transformed coverage values
    """
    y_exp = torch.pow(y + SQUASH_EPSILON, 0.75)
    
    mask = y_exp > threshold
    result = y_exp.clone()
    result[mask] = threshold + torch.sqrt(y_exp[mask] - threshold)
    
    return result


def unsquash_coverage_torch(y_squashed: torch.Tensor, threshold: float = SQUASH_THRESHOLD) -> torch.Tensor:
    """
    Inverse of squash_coverage - convert squashed values back to raw coverage.
    """
    mask = y_squashed > threshold
    y_exp = y_squashed.clone()
    
    residual = y_squashed[mask] - threshold
    y_exp[mask] = threshold + residual ** 2
    
    y_raw = torch.pow(y_exp, 4.0 / 3.0)
    
    return y_raw


# =============================================================================
# Model Utilities (from training-borzoi)
# =============================================================================

def count_parameters(model: nn.Module, trainable_only: bool = True) -> int:
    """Count model parameters."""
    if trainable_only:
        return sum(p.numel() for p in model.parameters() if p.requires_grad)
    return sum(p.numel() for p in model.parameters())


def print_model_summary(model: nn.Module, accelerator=None):
    """Print a summary of model parameters."""
    total_params = count_parameters(model, trainable_only=False)
    trainable_params = count_parameters(model, trainable_only=True)
    
    msg = f"Model Summary:\n"
    msg += f"  Total parameters: {total_params:,}\n"
    msg += f"  Trainable parameters: {trainable_params:,}\n"
    msg += f"  Non-trainable parameters: {total_params - trainable_params:,}\n"
    
    if accelerator is not None:
        accelerator.print(msg)
    else:
        print(msg)


def freeze_layers(model: nn.Module, layer_names: List[str], freeze: bool = True):
    """
    Freeze or unfreeze specific layers by name pattern.
    
    Args:
        model: PyTorch model
        layer_names: List of name patterns to match
        freeze: If True, freeze; if False, unfreeze
    """
    for name, param in model.named_parameters():
        if any(pattern in name for pattern in layer_names):
            param.requires_grad = not freeze


def get_layer_groups(model: nn.Module) -> Dict[str, List[str]]:
    """
    Group model parameters by layer type for differential learning rates.
    
    Returns dict with keys: 'trunk', 'transformer', 'head', 'other'
    """
    groups = {
        'trunk': [],
        'transformer': [],
        'head': [],
        'other': []
    }
    
    for name, param in model.named_parameters():
        if not param.requires_grad:
            continue
            
        name_lower = name.lower()
        if 'head' in name_lower:
            groups['head'].append(name)
        elif 'transformer' in name_lower or 'attention' in name_lower:
            groups['transformer'].append(name)
        elif 'conv' in name_lower or 'stem' in name_lower:
            groups['trunk'].append(name)
        else:
            groups['other'].append(name)
    
    return groups


# =============================================================================
# Genomic Interval Utilities
# =============================================================================

def expand_region(
    chrom: str,
    start: int,
    end: int,
    seq_len: int,
    chrom_sizes: Dict[str, int]
) -> Tuple[int, int, int, int]:
    """
    Expand a genomic region to the required sequence length, centered on the region.
    
    Args:
        chrom: Chromosome name
        start: Region start position
        end: Region end position
        seq_len: Required total sequence length
        chrom_sizes: Dictionary of chromosome sizes
        
    Returns:
        Tuple of (expanded_start, expanded_end, left_pad, right_pad)
    """
    region_len = end - start
    expand_each_side = (seq_len - region_len) // 2
    
    expanded_start = start - expand_each_side
    expanded_end = end + expand_each_side
    
    if (expanded_end - expanded_start) != seq_len:
        expanded_end = expanded_start + seq_len
    
    left_pad = 0
    right_pad = 0
    
    if expanded_start < 0:
        left_pad = -expanded_start
        expanded_start = 0
        
    chrom_size = chrom_sizes.get(chrom, float('inf'))
    if expanded_end > chrom_size:
        right_pad = expanded_end - chrom_size
        expanded_end = chrom_size
        
    return expanded_start, expanded_end, left_pad, right_pad


def tile_region(
    region_start: int,
    region_end: int,
    pred_len: int,
    padding: int
) -> List[Tuple[int, int, int, int]]:
    """
    Tile a region into overlapping windows for inference.
    
    Args:
        region_start: Start of region to tile
        region_end: End of region to tile
        pred_len: Length of prediction window (stride)
        padding: Context padding on each side
        
    Returns:
        List of (window_start, window_end, pred_start, pred_end) tuples
    """
    windows = []
    current_pred_start = region_start
    
    while current_pred_start < region_end:
        current_pred_end = min(current_pred_start + pred_len, region_end)
        
        window_start = current_pred_start - padding
        window_end = current_pred_end + padding
        
        windows.append((window_start, window_end, current_pred_start, current_pred_end))
        current_pred_start = current_pred_end
        
    return windows


# =============================================================================
# Training Utilities
# =============================================================================

def set_seed(seed: int):
    """Set random seeds for reproducibility."""
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# =============================================================================
# Inference Utilities
# =============================================================================

@torch.no_grad()
def predict_batch(
    model: nn.Module,
    sequences: torch.Tensor,
    device: torch.device = None,
) -> torch.Tensor:
    """
    Run inference on a batch of sequences.
    
    Args:
        model: Trained model
        sequences: One-hot encoded sequences of shape (batch, seq_len, 4)
        device: Device to run on
        
    Returns:
        Predictions of shape (batch, bins, tracks)
    """
    model.eval()
    
    if device is not None:
        sequences = sequences.to(device)
    
    outputs = model(sequences)
    
    if isinstance(outputs, tuple):
        outputs = outputs[0]
    
    return outputs


@torch.no_grad()
def predict_with_rc_average(
    model: nn.Module,
    sequences: torch.Tensor,
    device: torch.device = None,
) -> torch.Tensor:
    """
    Predict with reverse complement averaging for more robust predictions.
    
    Args:
        model: Trained model
        sequences: One-hot encoded sequences of shape (batch, seq_len, 4)
        device: Device to run on
        
    Returns:
        Averaged predictions of shape (batch, bins, tracks)
    """
    # Forward prediction
    pred_fwd = predict_batch(model, sequences, device)
    
    # Reverse complement prediction
    seq_rc = reverse_complement_torch(sequences)
    pred_rc = predict_batch(model, seq_rc, device)
    
    # Reverse the RC predictions to align with forward
    pred_rc = pred_rc.flip(1)
    
    # Average
    pred_avg = (pred_fwd + pred_rc) / 2.0
    
    return pred_avg
