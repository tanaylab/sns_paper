"""Borzoi model loading utilities."""

import math
from typing import Tuple

import torch
import torch.nn as nn
from accelerate import Accelerator

from .config import Config

# Default number of input channels for Borzoi head layers.
# This is the output dimension of the Borzoi trunk's final convolutional layer
# before the species-specific heads. This value is architecture-dependent and
# matches the official borzoi-pytorch implementation.
BORZOI_HEAD_IN_FEATURES = 1920


def _lecun_normal_init(module: nn.Module) -> None:
    """Apply LeCun normal initialization to a Conv1d module.

    This replicates the Keras/TensorFlow LeCun normal initialization used in the
    original Borzoi training. The initialization uses a truncated normal distribution
    with stddev = sqrt(1/fan_in) / 0.87962566103423978.

    The magic constant 0.87962566103423978 is the stddev of a standard normal
    distribution truncated to [-2, 2], used to scale the truncated normal to
    have unit variance.

    Reference: original Borzoi training code (borzoi_model.py)

    Args:
        module: A Conv1d module to initialize.
    """
    with torch.no_grad():
        fan_in = module.weight.size(1) * module.weight.size(2)
        scale = 1.0 / fan_in
        stddev = math.sqrt(scale) / 0.87962566103423978
        bound = 2.0 * stddev
        nn.init.trunc_normal_(module.weight, std=stddev, a=-bound, b=bound)


def _init_weights_keras(model: nn.Module, accelerator: Accelerator) -> None:
    """Apply Keras-compatible weight initialization for training from scratch.

    This replicates the initialization scheme used in the original Borzoi training
    from the original Borzoi training code:
    - Conv1d layers: LeCun normal initialization
    - Attention Q/K/V projections: Kaiming normal (nonlinearity='relu')
    - Attention output projection: Zeros
    - FFN linear layers: Kaiming normal (nonlinearity='relu')
    - Relative position biases: Kaiming normal
    - All biases: Zeros
    - BatchNorm: weight=1.0, bias=0.0
    - LayerNorm: weight=1.0, bias=0.0

    Args:
        model: The Borzoi model to initialize.
        accelerator: Accelerator instance for logging.
    """
    initialized_counts = {
        'conv1d_lecun': 0,
        'linear_kaiming': 0,
        'linear_zeros': 0,
        'norm_layers': 0,
        'attention_biases': 0,
    }

    for name, module in model.named_modules():
        # Conv1d layers: LeCun normal (Keras default for conv layers)
        if isinstance(module, nn.Conv1d):
            _lecun_normal_init(module)
            if module.bias is not None:
                nn.init.zeros_(module.bias)
            initialized_counts['conv1d_lecun'] += 1

        # Linear layers: context-dependent initialization
        elif isinstance(module, nn.Linear):
            # Attention output projections (to_out): initialize to zeros for residual stability
            if 'to_out' in name or 'out_proj' in name:
                nn.init.zeros_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
                initialized_counts['linear_zeros'] += 1
            else:
                # All other linear layers (Q/K/V projections, FFN): Kaiming normal
                nn.init.kaiming_normal_(module.weight, nonlinearity='relu')
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
                initialized_counts['linear_kaiming'] += 1

        # Normalization layers
        elif isinstance(module, (nn.LayerNorm, nn.BatchNorm1d)):
            if module.weight is not None:
                module.weight.data.fill_(1.0)
            if module.bias is not None:
                module.bias.data.zero_()
            initialized_counts['norm_layers'] += 1

    # Initialize relative position biases in attention layers (if they exist as Parameters)
    for name, param in model.named_parameters():
        if 'rel_content_bias' in name or 'rel_pos_bias' in name:
            nn.init.kaiming_normal_(param.data, nonlinearity='relu')
            initialized_counts['attention_biases'] += 1

    accelerator.print(f"Keras-compatible initialization applied:")
    accelerator.print(f"  Conv1d (LeCun normal): {initialized_counts['conv1d_lecun']} layers")
    accelerator.print(f"  Linear (Kaiming normal): {initialized_counts['linear_kaiming']} layers")
    accelerator.print(f"  Linear (zeros, attn out): {initialized_counts['linear_zeros']} layers")
    accelerator.print(f"  Norm layers: {initialized_counts['norm_layers']} layers")
    accelerator.print(f"  Attention biases (Kaiming): {initialized_counts['attention_biases']} params")


def _patch_maxpool_for_large_tensors(model: nn.Module, accelerator: Accelerator) -> None:
    """Replace MaxPool1d(kernel_size=2) with a reshape-based implementation.

    PyTorch's CUDA max_pool1d kernel uses 32-bit integer indexing internally.
    With large batch sizes and long sequences (e.g. batch=24, channels=512,
    length=524288 → 6.4B elements), the linear index overflows INT32_MAX
    (2,147,483,647), causing ``RuntimeError: integer out of range``.

    We use ``torch.maximum(x[:, :, 0::2], x[:, :, 1::2])`` which operates on
    strided views (no copies) and avoids allocating an int64 indices tensor
    that ``Tensor.max(dim=...)`` would create.

    Args:
        model: Borzoi or BorzoiHiRes model.
        accelerator: Accelerator instance for logging.
    """
    import types

    def _safe_forward(self, x: torch.Tensor) -> torch.Tensor:
        return torch.maximum(x[:, :, 0::2], x[:, :, 1::2])

    borzoi = model.encoder if hasattr(model, 'encoder') else model

    patched_count = 0
    for module in borzoi.modules():
        if isinstance(module, nn.MaxPool1d) and module.kernel_size == 2 and module.stride == 2:
            module.forward = types.MethodType(_safe_forward, module)
            patched_count += 1

    if patched_count:
        accelerator.print(
            f"Patched {patched_count} MaxPool1d layers with reshape-based pooling "
            f"(avoids INT32 overflow in CUDA kernel for large tensors)"
        )


def _patch_position_embeddings(model: nn.Module, seq_len: int, accelerator: Accelerator) -> None:
    """Monkey-patch Attention position buffers for sequence lengths beyond 524kb.

    The borzoi-pytorch Attention class hardcodes a position embedding buffer sized
    for 4096 transformer positions (seq_len=524288 / 128). For longer sequences,
    we replace the buffer with one sized for the actual transformer sequence length.

    Args:
        model: Borzoi or BorzoiHiRes model.
        seq_len: Input sequence length in bp.
        accelerator: Accelerator instance for logging.
    """
    from borzoi_pytorch.pytorch_borzoi_transformer import Attention, get_positional_embed

    transformer_seq_len = seq_len // 128
    if transformer_seq_len <= 4096:
        return  # Standard Borzoi, no patching needed

    borzoi = model.encoder if hasattr(model, 'encoder') else model

    patched_count = 0
    for module in borzoi.modules():
        if isinstance(module, Attention):
            new_positions = get_positional_embed(
                transformer_seq_len,
                module.num_rel_pos_features,
                module.positions.device,
            )
            module.register_buffer("positions", new_positions, persistent=False)
            patched_count += 1

    accelerator.print(
        f"Patched {patched_count} attention layers: position buffer "
        f"resized from 4096 to {transformer_seq_len} tokens"
    )


def load_borzoi_model(config: Config, accelerator: Accelerator) -> nn.Module:
    """Load Borzoi model using borzoi-pytorch.

    Supports loading from HuggingFace Hub or local weights.

    Args:
        config: Configuration object with model settings.
        accelerator: Accelerator instance for distributed training.

    Returns:
        Loaded Borzoi model.
    """
    from borzoi_pytorch import Borzoi
    from borzoi_pytorch.config_borzoi import BorzoiConfig

    accelerator.print("\n" + "=" * 70)
    accelerator.print("Loading Borzoi Model")
    accelerator.print("=" * 70)

    model_name = config.model.name.lower()

    # Map common names to HuggingFace model IDs
    model_mapping = {
        'borzoi_mouse_rep0': 'johahi/borzoi-replicate-0',
        'borzoi_human_rep0': 'johahi/borzoi-replicate-0',
        'flashzoi': 'johahi/flashzoi-replicate-0',
    }

    # Try to find model ID
    hf_model_id = None
    for key, model_id in model_mapping.items():
        if key in model_name:
            hf_model_id = model_id
            break

    if hf_model_id is None:
        hf_model_id = config.model.name

    accelerator.print(f"Loading model from: {hf_model_id}")

    # Determine if we need mouse head (config.model.species takes priority, then path detection)
    is_mouse_data, is_human_data = detect_species(config)
    explicit_species = config.get('model.species', None)
    if explicit_species:
        accelerator.print(f"Using explicit species from config: {explicit_species}")
    else:
        accelerator.print(f"Species auto-detected from genome path")

    # Load configuration
    try:
        borzoi_config = BorzoiConfig.from_pretrained(hf_model_id)
    except Exception as e:
        accelerator.print(f"Could not load config from {hf_model_id}, using defaults")
        borzoi_config = BorzoiConfig()

    # Enable mouse head for mouse data
    if is_mouse_data:
        accelerator.print("Mouse data detected - enabling mouse_head")
        borzoi_config.enable_mouse_head = True

    # Align crop/output length with requested prediction window
    if config.model.pred_len % config.model.bin_size != 0:
        raise ValueError(
            f"pred_len ({config.model.pred_len}) must be divisible by bin_size ({config.model.bin_size})"
        )

    target_bins = config.model.pred_len // config.model.bin_size
    max_bins = config.model.seq_len // config.model.bin_size
    if target_bins > max_bins:
        accelerator.print(
            f"Warning: pred_len {config.model.pred_len} exceeds available bins for seq_len {config.model.seq_len}. "
            f"Clipping to {max_bins * config.model.bin_size} bp ({max_bins} bins)."
        )
        target_bins = max_bins

    borzoi_config.return_center_bins_only = True
    borzoi_config.bins_to_return = target_bins

    # Load model (with optional random initialization)
    train_from_scratch = config.get('model.train_from_scratch', False)

    if train_from_scratch:
        accelerator.print("=" * 70)
        accelerator.print("TRAINING FROM SCRATCH - Initializing with random weights")
        accelerator.print("=" * 70)
        model = Borzoi(borzoi_config)

        # Apply Keras-compatible initialization (LeCun normal for conv, Kaiming for linear)
        # This replicates the original Borzoi training initialization from:
        # original Borzoi training code (borzoi_model.py)
        accelerator.print("Applying Keras-compatible initialization...")
        _init_weights_keras(model, accelerator)
        accelerator.print("Random initialization complete")
    else:
        # Load pretrained weights (existing behavior)
        try:
            model = Borzoi.from_pretrained(hf_model_id, config=borzoi_config, ignore_mismatched_sizes=True)
            accelerator.print(f"Loaded pretrained weights from {hf_model_id}")
            if is_mouse_data:
                accelerator.print("  Note: Using flashzoi config with borzoi weights - mouse_head initialized randomly")
        except Exception as e:
            accelerator.print(f"Warning: Could not load pretrained weights: {e}")
            accelerator.print("Falling back to random initialization")
            model = Borzoi(borzoi_config)
            model.apply(model._init_weights)

    # Patch MaxPool1d to avoid INT32 overflow with large batch × channel × length
    _patch_maxpool_for_large_tensors(model, accelerator)

    # Patch position embeddings for extended sequence lengths (>524kb)
    _patch_position_embeddings(model, config.model.seq_len, accelerator)

    # Determine if we need to modify the head
    num_output_tasks = config.model.num_output_tasks

    if num_output_tasks is not None and num_output_tasks > 0:
        _setup_model_head(
            model,
            accelerator,
            num_output_tasks,
            is_mouse_data,
            is_human_data,
            train_from_scratch,
        )

    # Enable gradient checkpointing if requested
    if config.training.gradient_checkpointing:
        _enable_gradient_checkpointing(model, accelerator)

    # Freeze trunk if requested
    if config.model.freeze_trunk:
        _freeze_trunk(model, accelerator)
    else:
        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        accelerator.print(f"\nTotal parameters: {total_params:,}")
        accelerator.print(f"Trainable parameters: {trainable_params:,}")
        init_status = "Random (Keras-compatible: LeCun/Kaiming)" if train_from_scratch else "Pretrained (HuggingFace)"
        accelerator.print(f"Initialization: {init_status}")

    accelerator.print("=" * 70)

    return model


def load_borzoi_hires_model(config: Config, accelerator: Accelerator) -> nn.Module:
    """Load BorzoiHiRes model: pretrained Borzoi encoder + U-Net decoder.

    Wraps a standard Borzoi model with additional decoder stages to achieve
    sub-32bp resolution output. The pretrained encoder is frozen by default.

    Args:
        config: Configuration with model.hires_resolution set.
        accelerator: Accelerator instance for distributed training.

    Returns:
        BorzoiHiRes model.
    """
    from .model_hires import BorzoiHiRes

    accelerator.print("\n" + "=" * 70)
    accelerator.print("Loading High-Resolution Borzoi Model")
    accelerator.print("=" * 70)

    target_res = config.model.hires_resolution
    accelerator.print(f"Target resolution: {target_res}bp (native Borzoi: 32bp)")

    # Load the base Borzoi encoder
    base_model = load_borzoi_model(config, accelerator)

    # Get decoder channel config
    decoder_channels = config.get('model.hires_decoder_channels', None)

    # Build HiRes wrapper
    from borzoi_pytorch.config_borzoi import BorzoiConfig
    borzoi_config = BorzoiConfig(dim=1536)  # for channel dimension reference

    model = BorzoiHiRes(
        borzoi_config=borzoi_config,
        pretrained_model=base_model,
        num_output_tracks=config.model.num_output_tasks,
        target_resolution=target_res,
        decoder_channels=decoder_channels,
        freeze_encoder=config.model.freeze_trunk,
    )

    if config.training.gradient_checkpointing:
        model.enable_gradient_checkpointing()
        accelerator.print("Decoder gradient checkpointing enabled")

    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    encoder_params = sum(p.numel() for p in model.encoder.parameters())
    decoder_params = total_params - encoder_params

    accelerator.print(f"\nBorzoiHiRes ({target_res}bp):")
    accelerator.print(f"  Encoder parameters: {encoder_params:,}")
    accelerator.print(f"  Decoder parameters: {decoder_params:,}")
    accelerator.print(f"  Total parameters: {total_params:,}")
    accelerator.print(f"  Trainable parameters: {trainable_params:,}")
    accelerator.print("=" * 70)

    return model


def _replace_head(
    module: nn.Module,
    attr: str,
    num_tasks: int,
    accelerator: Accelerator,
    train_from_scratch: bool = False,
) -> None:
    """Replace a Conv1d head with a new one for the specified number of tasks.

    The new head is automatically moved to the same device and dtype as the
    original head to ensure consistency.

    Args:
        module: The model containing the head to replace.
        attr: Name of the head attribute (e.g., 'mouse_head', 'human_head').
        num_tasks: Number of output tasks/tracks for the new head.
        accelerator: Accelerator instance for logging.
    """
    import torch

    current_head = getattr(module, attr)

    # Determine device and dtype from existing head
    # This ensures the new head is on the same device as the rest of the model
    device = None
    dtype = None
    try:
        param = next(current_head.parameters())
        device = param.device
        dtype = param.dtype
    except StopIteration:
        # No parameters in current head; try to get device from module
        try:
            param = next(module.parameters())
            device = param.device
            dtype = param.dtype
        except StopIteration:
            pass

    # Try to infer input features from existing head
    if hasattr(current_head, 'in_features'):
        in_features = current_head.in_features
    elif hasattr(current_head, 'in_channels'):
        in_features = current_head.in_channels
    elif hasattr(current_head, 'weight'):
        in_features = current_head.weight.shape[1]
    else:
        # Fall back to Borzoi's default trunk output dimension
        in_features = BORZOI_HEAD_IN_FEATURES

    new_head = nn.Conv1d(in_features, num_tasks, kernel_size=1)
    if train_from_scratch:
        _lecun_normal_init(new_head)
        if new_head.bias is not None:
            nn.init.zeros_(new_head.bias)

    # Move new head to the same device/dtype as the original
    if device is not None:
        new_head = new_head.to(device=device, dtype=dtype)

    setattr(module, attr, new_head)

    device_str = f" on {device}" if device is not None else ""
    accelerator.print(f"Replaced {attr} with {num_tasks} output tasks (in_features={in_features}){device_str}")


def _setup_model_head(
    model: nn.Module,
    accelerator: Accelerator,
    num_output_tasks: int,
    is_mouse_data: bool,
    is_human_data: bool,
    train_from_scratch: bool,
) -> None:
    """Setup the appropriate head for the model based on data type."""
    if is_mouse_data and not is_human_data:
        if not hasattr(model, 'mouse_head') and hasattr(model, 'enable_mouse_head'):
            accelerator.print("Mouse data detected; enabling mouse_head")
            try:
                model.enable_mouse_head()
            except Exception as e:
                accelerator.print(f"Warning: enable_mouse_head() failed: {e}")
        if hasattr(model, 'mouse_head'):
            accelerator.print("Replacing mouse_head for mouse genome")
            _replace_head(model, 'mouse_head', num_output_tasks, accelerator, train_from_scratch=train_from_scratch)
        elif hasattr(model, 'human_head'):
            accelerator.print("Falling back to human_head replacement (mouse_head unavailable)")
            _replace_head(model, 'human_head', num_output_tasks, accelerator, train_from_scratch=train_from_scratch)
        # Freeze the unused human head if present
        if hasattr(model, 'human_head'):
            for p in model.human_head.parameters():
                p.requires_grad = False
    else:
        if not hasattr(model, 'human_head') and hasattr(model, 'enable_human_head'):
            accelerator.print("Human data detected; enabling human_head")
            try:
                model.enable_human_head()
            except Exception as e:
                accelerator.print(f"Warning: enable_human_head() failed: {e}")
        if hasattr(model, 'human_head'):
            accelerator.print("Replacing human_head")
            _replace_head(model, 'human_head', num_output_tasks, accelerator, train_from_scratch=train_from_scratch)
        elif hasattr(model, 'head'):
            _replace_head(model, 'head', num_output_tasks, accelerator, train_from_scratch=train_from_scratch)
        if hasattr(model, 'mouse_head'):
            for p in model.mouse_head.parameters():
                p.requires_grad = False


def _enable_gradient_checkpointing(model: nn.Module, accelerator: Accelerator) -> None:
    """Enable gradient checkpointing for memory efficiency."""
    if hasattr(model, 'gradient_checkpointing_enable'):
        # Force support flag for models that have the method but don't set the flag (e.g. Borzoi)
        try:
            model.supports_gradient_checkpointing = True
            model.gradient_checkpointing = True
            if hasattr(type(model), 'supports_gradient_checkpointing'):
                type(model).supports_gradient_checkpointing = True
            
            # Add a manual _set_gradient_checkpointing if it's missing or failing
            if not hasattr(model, '_set_gradient_checkpointing'):
                def _set_gradient_checkpointing(module, value=False):
                    if hasattr(module, "gradient_checkpointing"):
                        module.gradient_checkpointing = value
                model._set_gradient_checkpointing = _set_gradient_checkpointing
        except Exception:
            pass
        
        try:
            model.gradient_checkpointing_enable()
            accelerator.print("\nGradient checkpointing enabled (memory-efficient mode)")
        except Exception as e:
            # If standard enable fails, try setting the attribute directly on all submodules
            accelerator.print(f"\nWarning: Standard gradient checkpointing enable failed: {e}")
            accelerator.print("Attempting manual submodule checkpointing activation...")
            count = 0
            for m in model.modules():
                if hasattr(m, "gradient_checkpointing"):
                    m.gradient_checkpointing = True
                    count += 1
            if count > 0:
                accelerator.print(f"  Manual activation applied to {count} submodules")
            else:
                accelerator.print("  No submodules with 'gradient_checkpointing' attribute found.")
    elif hasattr(model, 'set_gradient_checkpointing'):
        model.set_gradient_checkpointing(True)
        accelerator.print("\nGradient checkpointing enabled (memory-efficient mode)")
    else:
        accelerator.print("\nWarning: Model does not support gradient checkpointing")


def _freeze_trunk(model: nn.Module, accelerator: Accelerator) -> None:
    """Freeze all layers except head layers."""
    accelerator.print("\nFreezing trunk layers...")
    trainable_count = 0
    total_count = 0

    for name, param in model.named_parameters():
        total_count += 1
        if 'head' in name.lower():
            param.requires_grad = True
            trainable_count += 1
        else:
            param.requires_grad = False

    accelerator.print(f"Trainable parameters: {trainable_count}/{total_count}")


def _detect_species_from_path(genome_path: str, warn_on_ambiguity: bool = False) -> Tuple[bool, bool]:
    """Detect species from a genome path string.

    Species detection logic:
    - If path contains "mm10" or "mouse" -> mouse data
    - If path contains "hg19", "hg38", or "human" -> human data
    - If both detected -> mouse takes precedence (with optional warning)
    - If neither detected -> defaults to human (with optional warning)

    Note: This is a fallback method. Prefer setting model.species explicitly in config.

    Args:
        genome_path: Path to genome FASTA file.
        warn_on_ambiguity: If True, emit warnings for ambiguous or undetected species.

    Returns:
        Tuple of (is_mouse, is_human).
    """
    import warnings

    path_lower = genome_path.lower() if genome_path else ""
    is_mouse = "mm10" in path_lower or "mouse" in path_lower
    is_human_explicit = any(tag in path_lower for tag in ["hg19", "hg38", "human"])

    if is_mouse and is_human_explicit:
        # Ambiguous - path contains both mouse and human indicators
        if warn_on_ambiguity:
            warnings.warn(
                f"Genome path contains both mouse and human indicators: {genome_path}. "
                "Defaulting to mouse based on 'mm10'/'mouse' taking precedence. "
                "Consider setting model.species explicitly in config.",
                UserWarning
            )
        is_human = False
    elif is_mouse:
        is_human = False
    elif is_human_explicit:
        is_human = True
    else:
        # Neither detected - default to human
        if warn_on_ambiguity:
            warnings.warn(
                f"Could not detect species from genome path: {genome_path}. "
                "Defaulting to human. Set model.species='mouse' or 'human' in config, "
                "or add 'mm10', 'mouse', 'hg19', 'hg38', or 'human' to the genome path.",
                UserWarning
            )
        is_human = True

    return is_mouse, is_human


def detect_species(config: Config) -> Tuple[bool, bool]:
    """Detect whether data is mouse or human based on config.

    Detection priority:
    1. If model.species is set explicitly ('mouse' or 'human'), use that
    2. Otherwise, fall back to path-based detection from data.genome_fasta

    Args:
        config: Configuration object with model.species and data.genome_fasta.

    Returns:
        Tuple of (is_mouse, is_human).
    """
    # First, check for explicit species setting
    explicit_species = config.get('model.species', None)

    if explicit_species is not None:
        species_lower = str(explicit_species).lower().strip()
        if species_lower in ('mouse', 'mm10', 'mm39'):
            return True, False  # is_mouse=True, is_human=False
        elif species_lower in ('human', 'hg19', 'hg38'):
            return False, True  # is_mouse=False, is_human=True
        else:
            import warnings
            warnings.warn(
                f"Unknown model.species value: '{explicit_species}'. "
                "Expected 'mouse' or 'human'. Falling back to path-based detection.",
                UserWarning
            )

    # Fall back to path-based detection
    genome_path = str(config.data.genome_fasta) if config.data.genome_fasta else ""
    return _detect_species_from_path(genome_path, warn_on_ambiguity=True)
