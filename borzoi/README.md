# Borzoi/Flashzoi Training, Inference, and In-Silico Genome Analysis

All deep learning code, training configurations, and inference scripts for Fig 5 and Figs S7-S8 of the paper. Models predict H3K27me3 and H3K4me3 CUT&Tag tracks from DNA sequence at 32bp resolution.

## Overview

```
borzoi/
  code/              Training and inference Python library
    src/             Core modules (data loading, model, losses, metrics, etc.)
    scripts/         Genome perturbation utilities
  configs/           57 YAML training configurations (5 experiment categories)
  inference/         Shell scripts for in-silico genome experiments
```

## Model types

The paper uses two model architectures:

- **Flashzoi**: The Borzoi architecture trained from scratch on our epigenomic data (`train_from_scratch: true`, model name `flashzoi`)
- **Foundation Model (FM)**: The pre-trained Borzoi model (`johahi/borzoi-replicate-0`) fine-tuned on our data via a two-stage process (linear probe then full fine-tuning)

## Hardware requirements

All models were trained on 8 NVIDIA L40S GPUs (48GB VRAM each) using mixed-precision (fp16) training with PyTorch Distributed Data Parallel via `torchrun`.

- The 524kb receptive field models (the largest) use batch size 4 per GPU and require ~40GB VRAM per GPU.
- Smaller receptive fields (1kb–64kb) use proportionally larger batch sizes and train faster.
- Inference requires at least one GPU. Multi-GPU inference is supported and recommended for genome-wide predictions.

## Installation

```bash
cd borzoi/code
pip install -r requirements.txt
```

Key dependencies: PyTorch >= 2.0, Accelerate >= 0.20, Transformers >= 4.30, [borzoi-pytorch](https://pypi.org/project/borzoi-pytorch/) >= 0.1, pysam, pyBigWig.

## Training data

All models are trained on the same target data: normalized CUT&Tag H3K27me3 and H3K4me3 tracks from mouse ESC-derived embryoid bodies, stored as `pcg_norm_1000bp.parquet`. This file contains 1kb-binned signal values genome-wide.

- **Train chromosomes**: all autosomes except chr4, chr8, chr9, chr10, chr14, chr15
- **Validation chromosomes**: chr8, chr9
- **Test chromosomes**: chr4, chr10, chr14, chr15

## Training and inference code

The `code/` directory contains all Python code needed for training, inference, and in-silico perturbation analysis.

### Entry points

| Script | Description |
|--------|-------------|
| `code/train_borzoi_pytorch.py` | Training (PyTorch + Accelerate, multi-GPU distributed) |
| `code/infer_borzoi_pytorch.py` | Genome-wide inference (BigWig/HDF5/Parquet output) |
| `code/compute_perturbation.py` | In-silico perturbation analysis (transfection experiments) |
| `code/scripts/create_perturbed_genome.py` | Create perturbed genome FASTAs (silicus+/mm10- variants) |

### Core modules (`code/src/`)

| Module | Description |
|--------|-------------|
| `config.py` | YAML configuration parsing with defaults |
| `data_loaders.py` | Parquet/BigWig data loading, one-hot encoding, Borzoi squash transform |
| `model.py` | Borzoi model loading (HuggingFace or from-scratch initialization) |
| `model_hires.py` | U-Net decoder for high-resolution predictions |
| `losses.py` | Loss functions (Poisson, MSE, Cosine+MSE, Pearson) |
| `metrics.py` | Evaluation metrics (Pearson, Spearman, R-squared, genome-wide Pearson) |
| `training_utils.py` | Optimizer/scheduler setup, checkpointing, W&B integration |
| `evaluation.py` | Validation loop |
| `inference_engine.py` | Inference API with reverse-complement averaging |
| `inference.py` | Post-training inference orchestration |
| `perturbation_simple.py` | Splice and sequence perturbation operations |
| `batch_manager.py` | Adaptive batch sizing for OOM recovery |
| `utils.py` | Reverse complement, seed setting, general utilities |

## Config directory structure

### `configs/flashzoi_rf/` — Flashzoi from-scratch, receptive field comparison (Fig 5A-B)

10 configs training the Flashzoi architecture from random initialization on mm10, with receptive fields from 1kb to 524kb:

| Config | Receptive field | seq_len |
|--------|----------------|---------|
| `flashzoi_rf1k.yaml` | 1 kb | 1,024 |
| `flashzoi_rf2k.yaml` | 2 kb | 2,048 |
| `flashzoi_rf4k.yaml` | 4 kb | 4,096 |
| `flashzoi_rf8k.yaml` | 8 kb | 8,192 |
| `flashzoi_rf16k.yaml` | 16 kb | 16,384 |
| `flashzoi_rf32k.yaml` | 32 kb | 32,768 |
| `flashzoi_rf64k.yaml` | 64 kb | 65,536 |
| `flashzoi_rf128k.yaml` | 128 kb | 131,072 |
| `flashzoi_rf256k.yaml` | 256 kb | 262,144 |
| `flashzoi_rf524k.yaml` | 524 kb | 524,288 |

The `flashzoi_rf524k` model is also used as the base model for in-silico genome inference experiments (see `inference/`).

### `configs/borzoi_finetuned_rf/` — Foundation model RF comparison (Fig 5A-B)

10 configs fine-tuning the pre-trained Borzoi model on our data, with matching receptive fields:

| Config | Receptive field | Base model |
|--------|----------------|------------|
| `borzoi_finetuned_rf1k.yaml` | 1 kb | `johahi/borzoi-replicate-0` |
| ... | ... | ... |
| `borzoi_finetuned_rf524k.yaml` | 524 kb | `johahi/borzoi-replicate-0` |

Key difference from flashzoi: uses two-stage training (linear probe for 40 epochs, then full fine-tuning) and lower learning rate (5e-5 vs 6e-5).

### `configs/silicus_from_scratch/` — Flashzoi trained on in-silico genomes (Fig 5C)

12 configs training Flashzoi from scratch on synthetic genomes. The silicus genome is a synthetic genome sampled to match the CG and GC content distribution of mm10 but without any biological sequence features. Variants splice specific mm10 features back in:

**Single-feature addition (silicus+X):**
- `silicusPlusCGD.yaml` — silicus + CpG-dense domains
- `silicusPlusCRE.yaml` — silicus + cis-regulatory elements (ATAC peaks)
- `silicusPlusCTCF.yaml` — silicus + CTCF motif sites
- `silicusPlusExon.yaml` — silicus + exonic sequences
- `silicusPlusTE.yaml` — silicus + transposable elements (LINE/LTR)
- `silicusPlusRandom.yaml` — silicus + random 10% of 1kb tiles (control)

**Progressive telescope series (silicus+CGD+CRE+...):**
- `silicusPlusCGDCre.yaml` → `silicusPlusCGDCreCtcf.yaml` → `silicusPlusCGDCreCtcfExon.yaml` → `silicusPlusCGDCreCtcfExonTE.yaml`

**Base silicus variants:**
- `silicus55.yaml` — silicus with merged high-GC and high-CG bins
- `mus_silicus_cg_gc_lower.yaml` — base silicus genome

### `configs/silicus_finetuned/` — Flashzoi fine-tuned from mm10 to synthetic genomes (Fig 5, Fig S7-S8)

15 configs that take the mm10-trained Flashzoi checkpoint (`flashzoi_rf524k` best model) and fine-tune it on synthetic genomes. This tests whether the mm10-learned representations transfer to synthetic genomes.

Includes silicus variants (same as above) plus markovius variants:
- `markovius.yaml` — Markov-chain sampled genome (no CG/GC structure)
- `markoviusPlusCGDCre.yaml` — markovius + CGD + CRE

Each config runs inference on both the training genome and mm10 to compare predictions.

### `configs/fm_silicus/` — Foundation model fine-tuned on synthetic genomes (Fig S7-S8)

10 configs fine-tuning the pre-trained Borzoi on synthetic genomes. The FM telescope series uses a different feature ordering than the flashzoi series (adds GC first, then CGD+CTCF together):

silicus → +GC → +CGD+CTCF → +CGD+CTCF+CRE → +CGD+CTCF+CRE+Exon → +CGD+CTCF+CRE+Exon+TE → +CGD+CTCF+CRE+Exon+TE+UTR3

Also includes `fm_random.yaml`, `fm_markov.yaml`, `fm_markov_no_repeats.yaml`.

## `inference/` — Inference-only experiments (Fig 5C)

Shell scripts that run the mm10-trained model on modified genomes without retraining:

| Script | Model used | Genomes |
|--------|-----------|---------|
| `run_silicus_plus_inference.sh` | Flashzoi rf524k (mm10) | silicus, silicus+CGD, +CRE, +CTCF, +Exon, +TE, +Random |
| `run_silicus_telescopic_inference.sh` | Flashzoi rf524k (mm10) | silicus+CGD+CRE, +CTCF, +Exon, +TE (progressive) |
| `run_mm10_minus_inference.sh` | Flashzoi rf524k (mm10) | mm10 minus CGD, CRE, CTCF, Exon, TE, Random, etc. |
| `run_mm10_minus_inference_fm.sh` | Borzoi FM rf524k (mm10) | Same mm10-minus genomes |

## Running training

```bash
cd borzoi/code

# Train Flashzoi from scratch (8 GPUs)
torchrun --nproc_per_node=8 train_borzoi_pytorch.py \
    --config ../configs/flashzoi_rf/flashzoi_rf524k.yaml

# Fine-tune Borzoi foundation model (8 GPUs)
torchrun --nproc_per_node=8 train_borzoi_pytorch.py \
    --config ../configs/borzoi_finetuned_rf/borzoi_finetuned_rf524k.yaml
```

## Running inference

```bash
cd borzoi/code

# Genome-wide inference on a modified genome
torchrun --nproc_per_node=8 infer_borzoi_pytorch.py \
    --config ../configs/flashzoi_rf/flashzoi_rf524k.yaml \
    --genome_fasta genomes/silicus_plus/silicusPlusCGD/silicusPlusCGD.fa \
    --output_dir predictions/silicusPlusCGD

# Or use the inference shell scripts (which set track naming conventions):
cd ../inference
./run_silicus_plus_inference.sh
./run_mm10_minus_inference.sh
```

## In-silico transfection experiments (Fig 5D-E)

The transfection experiments (Jermann et al. / Schubeler, Fig 5D; Mendenhall et al. / Bernstein, Fig 5E) use the `compute_perturbation.py` script to run in-silico sequence perturbations with the fine-tuned Borzoi FM model.

### How the transfection results were generated

The transfection analysis splices donor CGI or BAC sequences into safe-harbor loci in the mm10 genome and predicts the resulting H3K27me3/H3K4me3 signals:

```bash
cd borzoi/code

# Schubeler transfections (Fig 5D): splice CGIs into safe-harbor loci
torchrun --nproc_per_node=8 compute_perturbation.py \
    --config ../configs/borzoi_finetuned_rf/borzoi_finetuned_rf524k.yaml \
    --checkpoint <path_to_best_model.safetensors> \
    --regions data/files/schubler/schubler_intervals_to_transfect.csv \
    --splice_pert <splice_perturbation_table> \
    --track_name EB4_cnt,EB4_cnt_k4 \
    --output data/files/schubler/results.csv \
    --distributed

# Bernstein transfections (Fig 5E): splice BAC constructs into safe-harbor loci
torchrun --nproc_per_node=8 compute_perturbation.py \
    --config ../configs/borzoi_finetuned_rf/borzoi_finetuned_rf524k.yaml \
    --checkpoint <path_to_best_model.safetensors> \
    --regions <bernstein_regions> \
    --splice_pert <bernstein_splice_table> \
    --track_name EB4_cnt,EB4_cnt_k4 \
    --output data/files/schubler/bernstein/results.csv \
    --distributed
```

### Pre-computed transfection results

The pre-computed results are provided in `data/files/schubler/`:

| File | Description |
|------|-------------|
| `results.csv` | Schubeler CGI transfections — baseline + per-CGI predictions at readout positions |
| `results1.csv` | Extended Schubeler results (additional CGI variants) |
| `bernstein/results.csv` | Bernstein BAC transfections |
| `bernstein_shifted/results.csv` | Bernstein BAC transfections with shifted positions |
| `cg_titration/results.csv` | CpG density titration experiment |
| `cg_addition/results.csv` | CpG addition experiment |
| `schubler_intervals_to_transfect.csv` | CGI definitions (coordinates, CG/GC content) |
| `seq_defs.csv` | CGI and synthetic sequences used as donors |

## Pre-trained model weights

Trained model checkpoints for all configurations are available on [HuggingFace](https://huggingface.co/tanaylab):

### mm10 receptive field series (Fig 5A-B)

| Model series | Receptive fields | HuggingFace |
|---|---|---|
| **Flashzoi** (from scratch) | 1k – 524k | [`sns-paper-flashzoi-rf{1k..524k}`](https://huggingface.co/collections/tanaylab/sns-paper-flashzoi-models-69b14882da177c6ec6623b6c) |
| **Borzoi FM** (fine-tuned) | 1k – 524k | [`sns-paper-borzoi-finetuned-rf{1k..524k}`](https://huggingface.co/collections/tanaylab/sns-paper-borzoi-fine-tuned-models-69b148d031a38501cc77ed44) |

### In-silico genome models (Fig 5C, Fig S7-S8)

| Model series | Description | HuggingFace |
|---|---|---|
| **Flashzoi silicus** (from scratch) | Trained from scratch on synthetic genomes (2 models) | [`sns-paper-flashzoi-*-from-scratch`](https://huggingface.co/collections/tanaylab/sns-paper-flashzoi-silicus-from-scratch-69b14fe4c25e7d1ad53458e6) |
| **Flashzoi silicus** (fine-tuned) | mm10-trained Flashzoi fine-tuned on synthetic genomes (10 models: telescope series + markovius) | [`sns-paper-flashzoi-*-finetuned`](https://huggingface.co/collections/tanaylab/sns-paper-flashzoi-silicus-fine-tuned-69b14fe5e3b1e2000df408e5) |
| **Borzoi FM silicus** (fine-tuned) | Pre-trained Borzoi FM fine-tuned on synthetic genomes (10 models: FM telescope + markov + random) | [`sns-paper-borzoi-fm-*`](https://huggingface.co/collections/tanaylab/sns-paper-borzoi-fm-silicus-fine-tuned-69b14fe746c6e7a9648abac8) |

Each repo contains `model.safetensors` (weights) and `config.yaml` (training configuration). All models are under the [`tanaylab`](https://huggingface.co/tanaylab) organization.

The mm10 RF models can be used directly with `infer_borzoi_pytorch.py` for inference, including on synthetic genomes (see `inference/` scripts). The silicus fine-tuned models were created by taking the mm10-trained checkpoint and fine-tuning on synthetic genomes.

## Data paths

All config paths are relative. Before running, set up the following:
- `mm10/mm10.fa` — mm10 reference genome FASTA (+ .fai index)
- `mm10/mm10.chrom.sizes` — chromosome sizes file
- `data/pcg_norm_1000bp.parquet` — training data (available from the paper's data repository)
- `genomes/` — synthetic genome FASTAs (created by `analysis/create_genomes.ipynb`)

Pre-computed model predictions for all experiments are provided in the misha tracks (included in the `mm10.tar.gz` data archive), so re-running training is not required to reproduce the figures.
