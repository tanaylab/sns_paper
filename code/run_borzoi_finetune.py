#!/usr/bin/env python
import os
# IMPORTANT: set backend before importing keras
os.environ["KERAS_BACKEND"] = "torch"

import argparse
from pathlib import Path
import zipfile
import tempfile

import numpy as np
import anndata
import keras

import crested
from crested.utils import read_bigwig_region, one_hot_encode_sequence
from crested.tl.zoo.utils._attention import MultiheadAttention  # force registration


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Finetune Borzoi to predict per-bin profiles on genomic windows (PyTorch backend)"
    )

    # Required
    parser.add_argument("--bigwigs-folder", required=True,
                        help="Path to folder containing bigwig files")
    parser.add_argument("--regions-file", required=True,
                        help="Consensus regions BED file (CGDDs / peaks / etc.)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for results")
    parser.add_argument("--project-name", required=True,
                        help="Project name (used in logging / checkpoints)")

    # Window / prediction config
    parser.add_argument("--seq-len", type=int, default=64000,
                        help="Input sequence length (window size), e.g. 64000 bp")
    parser.add_argument("--prediction-fraction", type=float, default=0.5,
                        help="Fraction of the window to predict in the center (e.g. 0.5 => center 32kb of 64kb)")
    parser.add_argument("--bin-size-bp", type=int, default=32,
                        help="Bin size in bp for bigWig binning (Borzoi uses 32 bp)")

    # Training
    parser.add_argument("--batch-size", type=int, default=4,
                        help="Batch size (each sample is a full window)")
    parser.add_argument("--epochs", type=int, default=10,
                        help="Number of training epochs")
    parser.add_argument("--initial-lr", type=float, default=1e-5,
                        help="Learning rate for Adam")
    parser.add_argument("--logger", choices=["wandb", "tensorboard", "none"],
                        default="none", help="Logging backend (only tensorboard is wired here)")

    # Resume
    parser.add_argument("--resume", action="store_true",
                        help="Resume from latest checkpoint in output_dir/checkpoints")
    parser.add_argument("--checkpoint-path", default=None,
                        help="Explicit .keras checkpoint to resume from (overrides auto-search)")

    # Split strategy
    parser.add_argument("--split-strategy", choices=["chr", "region"],
                        default="chr", help="Train/val/test split strategy")
    parser.add_argument("--val-size", type=float, default=0.1,
                        help="Validation size if split_strategy='region'")
    parser.add_argument("--test-size", type=float, default=0.1,
                        help="Test size if split_strategy='region'")
    parser.add_argument("--val-chroms", nargs="+", default=["chr8", "chr10"],
                        help="Validation chromosomes if split_strategy='chr'")
    parser.add_argument("--test-chroms", nargs="+", default=["chr9", "chr18"],
                        help="Test chromosomes if split_strategy='chr'")

    # Borzoi base model
    parser.add_argument("--borzoi-model", default="Borzoi_mouse_rep0",
                        help="Name for crested.get_model, e.g. Borzoi_mouse_rep0 or Borzoi_human_rep0")

    # Paths
    parser.add_argument("--chromsizes-file",
                        default="/home/aviezerl/proj/motif_reg/comparison/traj-for-crested/mm10.chrom.sizes",
                        help="Chromosome sizes file")
    parser.add_argument("--genome-file",
                        default="/home/michalel/tanaylab/data/mm10/mm10.fa",
                        help="Genome FASTA file")

    return parser.parse_args()


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def write_run_script(args, output_dir: Path):
    """Write a shell script that can reproduce the run with the same parameters."""
    script_path = output_dir / "reproduce_run.sh"
    python_script = os.path.abspath(__file__)

    cmd = [
        "#!/bin/bash",
        "",
        "# This script was automatically generated to reproduce the Borzoi profile run",
        "",
        f"python {python_script} \\",
    ]

    def add_flag(name, value):
        cmd.append(f"    --{name.replace('_','-')} {value} \\")

    # Required
    add_flag("bigwigs-folder", args.bigwigs_folder)
    add_flag("regions-file", args.regions_file)
    add_flag("output-dir", args.output_dir)
    add_flag("project-name", args.project_name)

    # Window / pred
    add_flag("seq-len", args.seq_len)
    add_flag("prediction-fraction", args.prediction_fraction)
    add_flag("bin-size-bp", args.bin_size_bp)

    # Training
    add_flag("batch-size", args.batch_size)
    add_flag("epochs", args.epochs)
    add_flag("initial-lr", args.initial_lr)
    add_flag("logger", args.logger)

    # Split
    add_flag("split-strategy", args.split_strategy)
    if args.split_strategy == "chr":
        cmd.append(f"    --val-chroms {' '.join(args.val_chroms)} \\")
        cmd.append(f"    --test-chroms {' '.join(args.test_chroms)} \\")
    else:
        add_flag("val-size", args.val_size)
        add_flag("test-size", args.test_size)

    # Borzoi / paths
    add_flag("borzoi-model", args.borzoi_model)
    add_flag("chromsizes-file", args.chromsizes_file)
    add_flag("genome-file", args.genome_file)

    # Resume flags only if set
    if args.resume:
        cmd.append("    --resume \\")
    if args.checkpoint_path:
        add_flag("checkpoint-path", args.checkpoint_path)

    # Remove trailing backslash
    cmd[-1] = cmd[-1].rstrip(" \\")
    with open(script_path, "w") as f:
        f.write("\n".join(cmd))
    os.chmod(script_path, 0o755)


def find_latest_checkpoint(output_dir: Path):
    """Find latest .keras checkpoint in output_dir/checkpoints."""
    ckpt_dir = output_dir / "checkpoints"
    if not ckpt_dir.exists():
        return None

    ckpts = list(ckpt_dir.glob("*.keras"))
    if not ckpts:
        return None

    # prefer numbered checkpoints, fall back to others
    numbered = [cp for cp in ckpts if cp.stem.split(".")[0].isdigit()]
    if numbered:
        numbered_sorted = sorted(numbered, key=lambda p: int(p.stem.split(".")[0]))
        return numbered_sorted[-1]

    # else just pick first (e.g. best_model.keras)
    return ckpts[0]


def create_borzoi_profile_model(
    seq_len: int,
    num_tracks: int,
    prediction_fraction: float,
    bin_size_bp: int,
    pretrained_model_path: str = None,
    model_name: str = None,
):
    """
    Build a Borzoi-based profile model that predicts profiles for the CENTER of the window.

    - Base Borzoi: outputs [target_length, pretrained_num_classes]
    - Head:
      - Dense along channels to num_tracks (your cell types / marks)
      - Cropping1D to keep only central fraction of bins
    """
    # decide pretrained_num_classes from model_name (human vs mouse)
    if model_name and "human" in model_name.lower():
        pretrained_num_classes = 7611
    else:
        pretrained_num_classes = 2608

    # Borzoi "native" binning: target_length = seq_len / 32
    target_length = seq_len // bin_size_bp

    base_model = crested.tl.zoo.borzoi(
        seq_len=seq_len,
        target_length=target_length,
        num_classes=pretrained_num_classes,
    )

    if pretrained_model_path:
        with zipfile.ZipFile(pretrained_model_path) as zf, tempfile.TemporaryDirectory() as tmpdir:
            weights_path = zf.extract("model.weights.h5", tmpdir)
            base_model.load_weights(weights_path)

    # Borzoi final conv layer in CREsted is "final_conv_activation"
    x = base_model.get_layer("final_conv_activation").output  # (B, target_length, pretrained_num_classes)

    # Map channel dimension to num_tracks with a position-wise Dense
    x = keras.layers.Dense(num_tracks, activation="softplus", name="dense_profile")(x)

    # Crop to center fraction in the bin dimension
    center_bins = int(round(target_length * prediction_fraction))
    if center_bins <= 0:
        raise ValueError("prediction_fraction too small; no bins left to predict.")
    if center_bins > target_length:
        raise ValueError("prediction_fraction > 1 is invalid.")

    left_crop = (target_length - center_bins) // 2
    right_crop = target_length - center_bins - left_crop

    if left_crop > 0 or right_crop > 0:
        x = keras.layers.Cropping1D(cropping=(left_crop, right_crop),
                                    name="center_crop")(x)
    # shape now: (B, center_bins, num_tracks)

    model = keras.Model(inputs=base_model.inputs,
                        outputs=x,
                        name="Borzoi_profile")

    # Store a few attributes for convenience
    model.seq_len = seq_len
    model.target_length = target_length
    model.center_bins = center_bins
    model.bin_size_bp = bin_size_bp
    return model


# ---------------------------------------------------------------------
# Data: import, split, windowing, generator / Sequence
# ---------------------------------------------------------------------
def load_or_preprocess_anndata(args, output_dir: Path) -> anndata.AnnData:
    """
    Use CREsted to:
      - import bigWigs + regions into AnnData (X = scalar summaries)
      - add train/val/test split
      - change the region width to seq_len (so var[start,end] are full windows)
    We ignore X for labels (we'll read bigWigs per bin ourselves), but
    we reuse var for coordinates + splitting.
    """
    adata_path = output_dir / "preprocessed_data.h5ad"
    if adata_path.exists():
        print(f"Loading preprocessed AnnData from {adata_path}")
        return anndata.read_h5ad(adata_path)

    print("Importing bigWigs via crested.import_bigwigs()...")
    adata = crested.import_bigwigs(
        bigwigs_folder=args.bigwigs_folder,
        regions_file=args.regions_file,
        chromsizes_file=args.chromsizes_file,
        target="mean",
        target_region_width=None,
    )

    # Split
    if args.split_strategy == "chr":
        crested.pp.train_val_test_split(
            adata,
            strategy="chr",
            val_chroms=args.val_chroms,
            test_chroms=args.test_chroms,
        )
    else:
        crested.pp.train_val_test_split(
            adata,
            strategy="region",
            val_size=args.val_size,
            test_size=args.test_size,
            random_state=42,
        )

    print("Split counts:")
    print(adata.var["split"].value_counts())

    # Resize regions to full window size (seq_len)
    print(f"Resizing regions to {args.seq_len} bp windows...")
    crested.pp.change_regions_width(
        adata,
        target_width=args.seq_len,
        chromsizes_file=args.chromsizes_file,
    )

    adata.write_h5ad(adata_path)
    return adata


def make_index_splits(adata: anndata.AnnData):
    """Return numpy index arrays for train/val/test based on var['split']."""
    split = adata.var["split"].values
    idx_train = np.where(split == "train")[0]
    idx_val = np.where(split == "val")[0]
    idx_test = np.where(split == "test")[0]
    return idx_train, idx_val, idx_test


def get_bigwig_paths(adata: anndata.AnnData, bigwigs_folder: str):
    """
    crested.import_bigwigs usually stores a 'file_path' column in obs with full path to bigWigs.
    We simply use those; otherwise fall back to {obs_name}.bw in the given folder.
    """
    if "file_path" in adata.obs:
        paths = list(adata.obs["file_path"].values)
    else:
        folder = Path(bigwigs_folder)
        paths = []
        for name in adata.obs_names:
            for ext in (".bw", ".bigWig", ".bigwig"):
                candidate = folder / f"{name}{ext}"
                if candidate.exists():
                    paths.append(str(candidate))
                    break
            else:
                raise FileNotFoundError(f"Could not find bigwig for {name} in {bigwigs_folder}")
    return paths


def window_generator(
    adata: anndata.AnnData,
    genome: crested.Genome,
    bigwig_paths,
    indices,
    seq_len: int,
    bin_size_bp: int,
    prediction_fraction: float,
):
    """
    Python generator that yields (X, Y) for each region index in `indices`:

    - X: one-hot encoded DNA sequence, shape (seq_len, 4)
    - Y: binned bigWig signal, center bins only, shape (center_bins, n_tracks)
    """
    chr_arr = adata.var["chr"].values
    start_arr = adata.var["start"].values
    end_arr = adata.var["end"].values

    n_tracks = len(bigwig_paths)
    target_length = seq_len // bin_size_bp
    center_bins = int(round(target_length * prediction_fraction))
    left_crop = (target_length - center_bins) // 2

    for idx in indices:
        chrom = chr_arr[idx]
        start = int(start_arr[idx])
        end = int(end_arr[idx])

        # Fetch sequence from genome, handle edge effects
        seq = genome.fetch(chrom, start, end).upper()
        if len(seq) != seq_len:
            if len(seq) < seq_len:
                pad = seq_len - len(seq)
                seq = seq + "N" * pad
            else:
                seq = seq[:seq_len]

        x = one_hot_encode_sequence(seq, expand_dim=False).astype(np.float32)  # (seq_len, 4)

        # Build profiles
        profiles = []
        for bw_path in bigwig_paths:
            values, _positions = read_bigwig_region(
                bigwig_file=bw_path,
                coordinates=(chrom, start, end),
                bin_size=bin_size_bp,
                target="mean",
            )  # values shape: [target_length] ideally

            if len(values) != target_length:
                if len(values) > target_length:
                    values = values[:target_length]
                else:
                    pad = target_length - len(values)
                    values = np.pad(values, (0, pad), constant_values=0.0)

            center_vals = values[left_crop:left_crop + center_bins]
            profiles.append(center_vals)

        y = np.stack(profiles, axis=-1).astype(np.float32)  # (center_bins, n_tracks)
        yield x, y


class WindowSequence(keras.utils.Sequence):
    """
    Keras Sequence wrapper around window_generator for use with any backend (including torch).
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        genome: crested.Genome,
        bigwig_paths,
        indices,
        seq_len: int,
        bin_size_bp: int,
        prediction_fraction: float,
        batch_size: int,
        shuffle: bool = True,
    ):
        self.adata = adata
        self.genome = genome
        self.bigwig_paths = bigwig_paths
        self.indices = np.array(indices, dtype=int)
        self.seq_len = seq_len
        self.bin_size_bp = bin_size_bp
        self.prediction_fraction = prediction_fraction
        self.batch_size = batch_size
        self.shuffle = shuffle

        self.target_length = seq_len // bin_size_bp
        self.center_bins = int(round(self.target_length * prediction_fraction))
        self.n_tracks = len(bigwig_paths)

        self.on_epoch_end()

    def __len__(self):
        # number of batches per epoch
        return int(np.ceil(len(self.indices) / self.batch_size))

    def on_epoch_end(self):
        if self.shuffle:
            np.random.shuffle(self.indices)

    def __getitem__(self, batch_index):
        batch_inds = self.indices[
            batch_index * self.batch_size:(batch_index + 1) * self.batch_size
        ]
        xs = []
        ys = []
        for x, y in window_generator(
            self.adata,
            self.genome,
            self.bigwig_paths,
            batch_inds,
            self.seq_len,
            self.bin_size_bp,
            self.prediction_fraction,
        ):
            xs.append(x)
            ys.append(y)
        xs = np.stack(xs, axis=0)  # (B, seq_len, 4)
        ys = np.stack(ys, axis=0)  # (B, center_bins, n_tracks)
        return xs, ys


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_run_script(args, output_dir)

    # AnnData with regions + splits + resized windows
    adata = load_or_preprocess_anndata(args, output_dir)

    # Determine number of tracks
    n_tracks = adata.n_obs  # obs = bigWig tracks

    # Create Genome object and register (optional)
    genome = crested.Genome(args.genome_file, args.chromsizes_file)
    crested.register_genome(genome)

    # Determine Borzoi target length / center bins consistently
    target_length = args.seq_len // args.bin_size_bp
    center_bins = int(round(target_length * args.prediction_fraction))
    print(f"Target_length (bins): {target_length}, center_bins: {center_bins}")

    # BigWig paths
    bigwig_paths = get_bigwig_paths(adata, args.bigwigs_folder)
    assert len(bigwig_paths) == n_tracks, "Number of bigwigs != n_obs"

    # Indices for splits
    idx_train, idx_val, idx_test = make_index_splits(adata)
    print(f"Train: {len(idx_train)}, Val: {len(idx_val)}, Test: {len(idx_test)}")

    # Model: create or load
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)

    if args.resume and (args.checkpoint_path or find_latest_checkpoint(output_dir) is not None):
        if args.checkpoint_path:
            ckpt = Path(args.checkpoint_path)
        else:
            ckpt = find_latest_checkpoint(output_dir)
        print(f"Resuming from checkpoint {ckpt}")
        custom_objects = {"MultiheadAttention": MultiheadAttention}
        model = keras.models.load_model(ckpt, custom_objects=custom_objects)
        # infer meta (fallback if not stored as attributes)
        if not hasattr(model, "center_bins"):
            model.center_bins = center_bins
        if not hasattr(model, "seq_len"):
            model.seq_len = args.seq_len
    else:
        print("Creating Borzoi profile model from pretrained weights...")
        model_path, _ = crested.get_model(args.borzoi_model)  # zip + output names
        model = create_borzoi_profile_model(
            seq_len=args.seq_len,
            num_tracks=n_tracks,
            prediction_fraction=args.prediction_fraction,
            bin_size_bp=args.bin_size_bp,
            pretrained_model_path=model_path,
            model_name=args.borzoi_model,
        )

    # Loss & metrics (same CosineMSELogLoss as peak regression, but now per-bin)
    optimizer = keras.optimizers.Adam(learning_rate=args.initial_lr)
    loss = crested.tl.losses.CosineMSELogLoss(max_weight=100)
    metrics = [
        keras.metrics.MeanAbsoluteError(name="mae"),
        keras.metrics.MeanSquaredError(name="mse"),
        keras.metrics.CosineSimilarity(axis=-1, name="cosine"),
    ]

    model.compile(optimizer=optimizer, loss=loss, metrics=metrics)
    model.summary()

    # Keras Sequences (backend-agnostic, works with torch backend)
    train_seq = WindowSequence(
        adata, genome, bigwig_paths,
        idx_train, args.seq_len, args.bin_size_bp,
        args.prediction_fraction, args.batch_size,
        shuffle=True,
    )
    val_seq = WindowSequence(
        adata, genome, bigwig_paths,
        idx_val, args.seq_len, args.bin_size_bp,
        args.prediction-fraction if hasattr(args, "prediction-fraction") else args.prediction_fraction,
        args.batch_size,
        shuffle=False,
    )
    test_seq = WindowSequence(
        adata, genome, bigwig_paths,
        idx_test, args.seq_len, args.bin_size_bp,
        args.prediction_fraction, args.batch_size,
        shuffle=False,
    )

    # Callbacks
    ckpt_cb = keras.callbacks.ModelCheckpoint(
        filepath=str(checkpoint_dir / "{epoch:03d}.keras"),
        monitor="val_loss",
        save_best_only=False,
        save_weights_only=False,
        verbose=1,
    )
    es_cb = keras.callbacks.EarlyStopping(
        monitor="val_loss",
        patience=5,
        restore_best_weights=True,
    )
    cbs = [ckpt_cb, es_cb]

    if args.logger == "tensorboard":
        tb_logdir = output_dir / "tensorboard"
        tb_cb = keras.callbacks.TensorBoard(log_dir=str(tb_logdir))
        cbs.append(tb_cb)

    # Train
    print("Starting training (PyTorch backend)...")
    model.fit(
        train_seq,
        validation_data=val_seq,
        epochs=args.epochs,
        callbacks=cbs,
    )

    # Evaluate
    print("Evaluating on test set...")
    test_results = model.evaluate(test_seq, return_dict=True)
    print("Test metrics:", test_results)

    print("Done.")


if __name__ == "__main__":
    main()
