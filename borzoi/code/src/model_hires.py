"""High-resolution (sub-32bp) Borzoi model with U-Net decoder.

This module provides a wrapper around pretrained Borzoi that adds additional
U-Net decoder stages to achieve resolutions finer than the native 32bp.
Skip connections from the encoder are captured via forward hooks, avoiding
any modification to the pretrained borzoi-pytorch source.
"""

import math
from typing import Dict, List, Optional

import torch
import torch.nn as nn


class SkipConnectionCapture:
    """Captures intermediate encoder activations via forward hooks.

    Registers hooks on specific ConvBlock modules inside the Borzoi encoder
    to capture skip connections at each resolution level for the decoder.

    Captured skips (for input length L):
        conv_dna:     (B, 512,  L/2)   - 2x stride
        res_tower_0:  (B, 608,  L/2)   - 2x stride (before pool)
        res_tower_2:  (B, 736,  L/4)   - 4x stride (before pool)
        res_tower_4:  (B, 896,  L/8)   - 8x stride (before pool)
        res_tower_6:  (B, 1056, L/16)  - 16x stride (before pool)
    """

    def __init__(self, borzoi_model: nn.Module):
        self.model = borzoi_model
        self.skips: Dict[str, torch.Tensor] = {}
        self._hooks: List[torch.utils.hooks.RemovableHook] = []

    def _make_hook(self, name: str):
        def hook_fn(module, input, output):
            self.skips[name] = output
        return hook_fn

    def register_hooks(self):
        """Register forward hooks on encoder submodules."""
        self.remove_hooks()
        self.skips.clear()

        # conv_dna output: (B, 512, L/2) after its internal Conv1d + MaxPool(2)
        h = self.model.conv_dna.register_forward_hook(self._make_hook('conv_dna'))
        self._hooks.append(h)

        # res_tower ConvBlocks at indices 0, 2, 4, 6 (before their subsequent MaxPool).
        # res_tower is Sequential: [ConvBlock, MaxPool, ConvBlock, MaxPool, ..., ConvBlock]
        # The MaxPool instances are all the same object, so we hook the unique ConvBlocks.
        for idx in [0, 2, 4, 6]:
            h = self.model.res_tower[idx].register_forward_hook(
                self._make_hook(f'res_tower_{idx}')
            )
            self._hooks.append(h)

    def remove_hooks(self):
        """Remove all registered hooks."""
        for h in self._hooks:
            h.remove()
        self._hooks.clear()

    def clear(self):
        """Clear captured skip tensors (call between forward passes)."""
        self.skips.clear()


class HiResDecoderStage(nn.Module):
    """Single decoder upsampling stage: upsample 2x -> concat skip -> project -> conv.

    Args:
        in_channels: Channels from previous decoder stage.
        skip_channels: Channels from encoder skip connection.
        out_channels: Output channels for this stage.
        kernel_size: Kernel size for the convolution block.
    """

    def __init__(
        self,
        in_channels: int,
        skip_channels: int,
        out_channels: int,
        kernel_size: int = 5,
    ):
        super().__init__()
        self.upsample = nn.Upsample(scale_factor=2, mode='linear', align_corners=False)
        self.project = nn.Conv1d(in_channels + skip_channels, out_channels, kernel_size=1)
        self.project_no_skip = nn.Conv1d(in_channels, out_channels, kernel_size=1)
        self.norm = nn.BatchNorm1d(out_channels, eps=0.001)
        self.act = nn.GELU(approximate='tanh')
        self.conv = nn.Conv1d(out_channels, out_channels, kernel_size=kernel_size, padding='same')

    def forward(self, x: torch.Tensor, skip: Optional[torch.Tensor] = None) -> torch.Tensor:
        x = self.upsample(x)

        if skip is not None:
            # Handle potential size mismatch from rounding
            if x.shape[-1] != skip.shape[-1]:
                x = nn.functional.interpolate(x, size=skip.shape[-1], mode='linear', align_corners=False)
            x = torch.cat([x, skip], dim=1)
            x = self.project(x)
        else:
            x = self.project_no_skip(x)

        x = self.norm(x)
        x = self.act(x)
        x = self.conv(x)
        return x


class BorzoiHiRes(nn.Module):
    """Borzoi with high-resolution U-Net decoder.

    Wraps a pretrained Borzoi encoder and adds multi-stage U-Net decoder
    stages to upsample from 32bp to the target resolution (1, 2, 4, 8, or 16bp).
    Skip connections are captured from the encoder via forward hooks.

    The pretrained encoder is frozen by default; only the new decoder and head
    are trained. The standard 32bp Borzoi path is completely unchanged.

    Architecture (example for 1bp, 5 new stages):
        Encoder (pretrained, frozen):
            Input 1bp -> conv_dna(2x) -> res_tower(16x) -> unet1(2x) -> transformer(2x)
            -> existing unet decoder (2x, 2x) -> crop -> embeddings at 32bp

        New decoder:
            32bp --[stage0: skip=res_tower_6 @16bp]--> 16bp
            16bp --[stage1: skip=res_tower_4 @8bp]-->  8bp
             8bp --[stage2: skip=res_tower_2 @4bp]-->  4bp
             4bp --[stage3: skip=res_tower_0 @2bp]-->  2bp
             2bp --[stage4: skip=conv_dna    @2bp]-->  1bp (no skip, just upsample)
        -> head -> softplus -> output at 1bp
    """
    supports_gradient_checkpointing = True

    # Skip connection specs: (hook_name, channels, encoder_stride_bp)
    # Ordered from coarsest to finest, matching decoder stage order.
    ENCODER_SKIPS = [
        ('res_tower_6', 1056, 16),  # 16bp
        ('res_tower_4', 896,  8),   # 8bp
        ('res_tower_2', 736,  4),   # 4bp
        ('res_tower_0', 608,  2),   # 2bp
        # conv_dna is at 2bp, so the last stage (2->1) has no matching skip.
        # We use None for the final stage.
    ]

    def __init__(
        self,
        borzoi_config=None,
        num_output_tracks: int = 1,
        target_resolution: int = 1,
        decoder_channels: Optional[List[int]] = None,
        freeze_encoder: bool = True,
        pretrained_model: nn.Module = None,
    ):
        super().__init__()

        assert target_resolution in (1, 2, 4, 8, 16, 32), \
            f"target_resolution must be power of 2 in [1..32], got {target_resolution}"

        self.target_resolution = target_resolution
        self.num_new_stages = int(math.log2(32 // target_resolution))

        if decoder_channels is None:
            all_channels = [512, 256, 128, 64, 32]
            decoder_channels = all_channels[:self.num_new_stages]

        assert len(decoder_channels) == self.num_new_stages, \
            f"Need {self.num_new_stages} decoder channel specs, got {len(decoder_channels)}"

        # Encoder
        if pretrained_model is not None:
            self.encoder = pretrained_model
        else:
            from borzoi_pytorch import Borzoi
            self.encoder = Borzoi(borzoi_config)

        # Skip capture
        self.skip_capture = SkipConnectionCapture(self.encoder)
        self.skip_capture.register_hooks()

        # Build decoder stages
        encoder_out_channels = borzoi_config.dim if borzoi_config else 1536

        self.decoder_stages = nn.ModuleList()
        # Use available skips for the first N-1 stages, None for the last if going to 1bp
        skip_specs = self.ENCODER_SKIPS[:self.num_new_stages]

        self._skip_names = []
        self._skip_strides = []

        in_ch = encoder_out_channels
        for i in range(self.num_new_stages):
            if i < len(skip_specs):
                skip_name, skip_ch, skip_stride = skip_specs[i]
                # The last stage (2bp->1bp) has no 1bp skip; conv_dna is at 2bp
                # so we only use it if we're going from 4->2, not 2->1
                self._skip_names.append(skip_name)
                self._skip_strides.append(skip_stride)
            else:
                # No skip available for this stage
                skip_ch = encoder_out_channels  # placeholder, won't be used
                self._skip_names.append(None)
                self._skip_strides.append(None)

            out_ch = decoder_channels[i]
            self.decoder_stages.append(
                HiResDecoderStage(
                    in_channels=in_ch,
                    skip_channels=skip_ch,
                    out_channels=out_ch,
                    kernel_size=5,
                )
            )
            in_ch = out_ch

        # Output head
        self.head_norm = nn.BatchNorm1d(decoder_channels[-1], eps=0.001)
        self.head_act = nn.GELU(approximate='tanh')
        self.head = nn.Conv1d(decoder_channels[-1], num_output_tracks, kernel_size=1)
        self.final_softplus = nn.Softplus()

        self._gradient_checkpointing = False

        if freeze_encoder:
            for param in self.encoder.parameters():
                param.requires_grad = False
        else:
            # Even if trunk is not frozen, we must freeze the internal Borzoi heads
            # and final join layers because they are not used in the HiRes forward pass
            # (which stops at get_embs_after_crop) and would cause DDP errors.
            unused_patterns = ['head', 'final_joined_convs', 'final_softplus']
            frozen_count = 0
            for name, param in self.encoder.named_parameters():
                if any(pattern in name for pattern in unused_patterns):
                    param.requires_grad = False
                    frozen_count += 1

    def enable_gradient_checkpointing(self):
        """Enable gradient checkpointing for both encoder and decoder stages."""
        self._gradient_checkpointing = True
        if hasattr(self.encoder, 'gradient_checkpointing_enable'):
            # Force support flag if needed (similar to load_borzoi_model)
            try:
                self.encoder.supports_gradient_checkpointing = True
                self.encoder.gradient_checkpointing = True
                if hasattr(type(self.encoder), 'supports_gradient_checkpointing'):
                    type(self.encoder).supports_gradient_checkpointing = True
            except Exception:
                pass
            try:
                self.encoder.gradient_checkpointing_enable()
            except Exception:
                # Manual fallback if encoder's method fails
                for m in self.encoder.modules():
                    if hasattr(m, "gradient_checkpointing"):
                        m.gradient_checkpointing = True

    # Alias for compatibility with _enable_gradient_checkpointing in model.py
    gradient_checkpointing_enable = enable_gradient_checkpointing

    def disable_gradient_checkpointing(self):
        """Disable gradient checkpointing for decoder stages."""
        self._gradient_checkpointing = False

    def forward(self, x: torch.Tensor, is_human: bool = True, **kwargs) -> torch.Tensor:
        """Forward pass.

        Args:
            x: Input one-hot DNA, shape (B, 4, L).
            is_human: Ignored (kept for API compatibility with Borzoi forward).
            **kwargs: Ignored (compatibility with Borzoi's data_parallel_training etc).

        Returns:
            Predictions at target resolution, shape (B, num_tracks, crop_bins * 32/target_res).
        """
        self.skip_capture.clear()

        # Run encoder: returns (B, dim, crop_bins) at 32bp resolution
        with torch.set_grad_enabled(
            any(p.requires_grad for p in self.encoder.parameters())
        ):
            embs_32bp = self.encoder.get_embs_after_crop(x)

        crop_bins = embs_32bp.shape[-1]

        # Run decoder stages
        h = embs_32bp
        for i, stage in enumerate(self.decoder_stages):
            skip_name = self._skip_names[i]
            skip_stride = self._skip_strides[i]

            if skip_name is not None:
                raw_skip = self.skip_capture.skips.get(skip_name)
                if raw_skip is not None:
                    # Compute target length for the skip at this decoder level.
                    # After upsampling h by 2x, we get h_len * 2.
                    # The skip at stride S covers L/S positions of the full input.
                    # But we're working with cropped output, so center-crop skip
                    # to match the upsampled decoder feature map.
                    target_len = h.shape[-1] * 2
                    skip = self._center_crop_skip(raw_skip, target_len)
                else:
                    skip = None
            else:
                skip = None

            if self._gradient_checkpointing and self.training:
                h = torch.utils.checkpoint.checkpoint(stage, h, skip, use_reentrant=False)
            else:
                h = stage(h, skip)

        # Head
        h = self.head_norm(h)
        h = self.head_act(h)
        h = self.final_softplus(self.head(h.float()))

        return h

    @staticmethod
    def _center_crop_skip(skip: torch.Tensor, target_length: int) -> torch.Tensor:
        """Center-crop a skip connection tensor to match decoder spatial size."""
        current_length = skip.shape[-1]
        if current_length == target_length:
            return skip
        if current_length < target_length:
            # Skip is smaller; pad symmetrically with zeros
            pad = target_length - current_length
            return nn.functional.pad(skip, (pad // 2, pad - pad // 2))
        # Center crop
        trim = (current_length - target_length) // 2
        return skip[:, :, trim:trim + target_length]
