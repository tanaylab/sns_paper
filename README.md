# Seed and Spread (SnS) model of the pluripotent epigenome paper companion repository

Welcome to the companion repository for the 
**Reverse engineering the genomic encodings of the pluripotent epigenome**

## Overview
During early embryogenesis the epigenome is globally erased and rebuilt in pluripotent cells of the **epiblast**.

Because this state is established **de novo**, it must largely be **encoded in DNA sequence**.

This project asks:  
**Can the pluripotent epigenome be predicted directly from sequence?**  
**Can genomic transformers “explain”, not just predict, the epigenome?**  




This repository contains all analysis code for the seed-and-spread (SNS) model of the pluripotent epigenome. The model predicts genome-wide distributions of Polycomb (H3K27me3), Trithorax (H3K4me3) and DNA methylation.
The model follows a **seed-and-spread framework**:

1. **CpG-dense domains (CGDDs)** act as epigenomic seeds  
2. **Motif logic** determines PcG vs TxG identity at seeds  
3. **Polycomb spreading** generates chromosomal PcG domains  
DNA methylation acts as an **antagonist of PcG spreading**.

## Repository structure

```
sns_paper/
  analysis/          Jupyter notebooks reproducing all main and supplementary figures
  borzoi/            Borzoi/Flashzoi training, inference, and configs (Fig 5)
    code/            Training and inference Python code (borzoi-finetune library)
    configs/         YAML training configurations organized by experiment type
    inference/       Shell scripts for in-silico genome inference experiments
  code/              R and Python source code (utility functions, models, data preparation)
  data/              Genomic data, misha track databases, and precomputed intermediates
    mm10/            Mouse genome (misha database with tracks)
    hg19/            Human genome (misha database with tracks)
  tables/            Supplementary tables and large precomputed R objects
  figs/              Figure outputs (PPTX and PNG)
  output/            Analysis output files
```

## Figure-to-notebook mapping

### Main figures

| Figure | Notebook | Description |
|--------|----------|-------------|
| Fig 1 | `analysis/sns_Fig1.ipynb` | CGDD definition, motif enrichment, PcG/TrxG seed classification (mouse + human) |
| Fig 2 | `analysis/sns_Fig2_Fig3EDF.ipynb` | Genome-wide seed-and-spread model, spreading from seeds |
| Fig 3 | `analysis/sns_Fig3.ipynb` | Model errors, chromosomal spreading, imprinting loci |
| Fig 4 | `analysis/sns_Fig4.ipynb` | DNA methylation antagonizes PcG (WT vs DKO analysis) |
| Fig 5 | `analysis/sns_Fig5_borzoi_figure.ipynb` | Borzoi pre-training, receptive fields, in-silico perturbations |
| Fig 5 | `analysis/sns_Fig5_borzoi_AUPRC.ipynb` | Borzoi AUPRC performance metrics |
| Fig 5 | `analysis/sns_Fig5_borzoi_performance.ipynb` | Borzoi R-squared performance metrics |
| Fig 5D | `analysis/sns_Fig5_schubeler_transfections.ipynb` | In-silico reproduction of Jermann et al. integration experiments |
| Fig 5E | `analysis/sns_Fig5_bernstein_transfections.ipynb` | In-silico reproduction of Mendenhall et al. integration experiments |
| Fig 6 | `analysis/sns_Fig6.ipynb` | PcG/TrxG differentiation in mesoderm vs ectoderm |

### Supplementary figures

| Figure | Notebook |
|--------|----------|
| Fig S1 | `analysis/sns_Fig1.ipynb` (RNA, ATAC, CUT&Tag vs ChIP comparisons) |
| Fig S2 | `analysis/sns_Fig1.ipynb` (motif clustering, model performance) |
| Fig S3 | `analysis/sns_Fig1.ipynb` (human CGDD analysis) |
| Fig S4 | `analysis/sns_Fig2_Fig3EDF.ipynb` (human genome-wide predictions) |
| Fig S5 | `analysis/sns_FigS5_and_S6dkoRNA.ipynb` |
| Fig S6 | `analysis/sns_FigS5_and_S6dkoRNA.ipynb` + `analysis/sns_Fig3.ipynb` |
| Fig S7-S8 | `analysis/sns_Fig5_borzoi_AUPRC.ipynb` + `analysis/sns_Fig5_borzoi_figure.ipynb` |
| Fig S9 | `analysis/sns_FigS9_UMAP.ipynb` (multiome UMAP analysis) |

### Key utility scripts

| File | Description |
|------|-------------|
| `code/seq2epi_utils.r` | Core functions for CGDD detection, motif analysis, and the SNS model |
| `code/fig_fun.r` | Figure generation and plotting functions |
| `code/hgpcg_domains_hg19.r` | Human genome CGDD analysis functions |
| `code/borzoi_utils.R` | Borzoi model evaluation and track utilities |
| `code/borzoi-plot.R` | Borzoi figure plotting functions |
| `code/plot_gw_rsqr.R` | Genome-wide R-squared plotting |
| `code/CRE_dev.R` | CRE accessibility analysis (data preparation) |
### Data preparation scripts

These scripts document how intermediate data was generated. They reference local paths and are not meant to be re-run directly:

| File | Description |
|------|-------------|
| `analysis/create_genomes.ipynb` | Create synthetic genomes (silicus, markovius, mm10-minus variants) for Fig 5 |
| `code/prepare_data_borzoi.R` / `.ipynb` | Import Borzoi predictions into misha tracks |
| `code/norm_run_sns_cnt_atac_norm.ipynb` | CUT&Tag normalization using ATAC tracks |
| `code/find_atac_peaks.ipynb` | ATAC peak calling |
| `code/IQ-*.ipynb` | IceQream model training for human CGDD classification |
| `code/build_IQ_responses_hg19.ipynb` | Build IQ response features for human genome |
| `code/Prepare-data-for-IQ-hg.ipynb` | Prepare data for IQ human analysis |

## Data availability

### Raw sequencing data

All newly generated CUT&Tag, ATAC-seq, multiome, and RNA-seq data are available from GEO (accession pending).

### Data downloads

All supplementary data files are hosted on Amazon S3 (`sns-paper` bucket, us-east-1).

| Archive | Source Directory | Compressed Size | URL |
|---------|-----------------|----------------|-----|
| `files.tar.gz` | `data/files/` | 8.4 GiB | https://sns-paper.s3.amazonaws.com/files.tar.gz |
| `mm10.tar.gz` | `data/mm10/` | 7.7 GiB | https://sns-paper.s3.amazonaws.com/mm10.tar.gz |
| `hg19.tar.gz` | `data/hg19/` | 10.3 GiB | https://sns-paper.s3.amazonaws.com/hg19.tar.gz |
| `tables.tar.gz` | `tables/` | 583 MiB | https://sns-paper.s3.amazonaws.com/tables.tar.gz |
| `output.tar.gz` | `output/` | 4.7 MiB | https://sns-paper.s3.amazonaws.com/output.tar.gz |

Download and extract:

```bash
wget https://sns-paper.s3.amazonaws.com/files.tar.gz
wget https://sns-paper.s3.amazonaws.com/mm10.tar.gz
wget https://sns-paper.s3.amazonaws.com/hg19.tar.gz
wget https://sns-paper.s3.amazonaws.com/tables.tar.gz
wget https://sns-paper.s3.amazonaws.com/output.tar.gz

#tar xzf files.tar.gz -C data/
tar -xzf files.tar.gz --strip-components=1 -C /data
tar xzf mm10.tar.gz -C data/
tar xzf hg19.tar.gz -C data/
tar xzf tables.tar.gz
tar xzf output.tar.gz
```

The `data/` directory includes:

- `data/mm10/` — Mouse misha genome database with CUT&Tag, ATAC, and model prediction tracks
- `data/hg19/` — Human misha genome database with CUT&Tag tracks
- `data/` — Intermediate analysis files (CGDD annotations, motif features, model predictions, in-silico transfection results)

### Supplementary tables

| File | Description |
|------|-------------|
| `tables/TableS1.csv` | missclassified CGDDs |
| `tables/tableS2.csv` | Model underestimations of the spreading |
| `tables/tableS3.xlsx` | New mesoderm domains |
| `tables/tableS4.xlsx` | ATAC-seq tracks used for CnT normalization |

## Dependencies

### R packages

The analysis was developed using R 4.4.1. The following packages are required:

**From CRAN:**
`tidyverse`, `data.table`, `glmnet`, `zoo`, `pheatmap`, `ggrepel`, `umap`, `doMC`, `doParallel`, `patchwork`, `officer`, `rvg`, `PRROC`

**From CRAN (Tanay lab):**
- [`misha`](https://github.com/tanaylab/misha) — Genomic database and interval analysis
- [`tgstat`](https://github.com/tanaylab/tgstat) — Statistical utilities
- [`tglkmeans`](https://github.com/tanaylab/tglkmeans) — K-means clustering

**From Bioconductor:**
`Biostrings`

**From GitHub (Tanay lab):**
- [`misha.ext`](https://github.com/tanaylab/misha.ext) — Extended misha utilities
- [`tgutil`](https://github.com/tanaylab/tgutil) — General utilities
- [`prego`](https://github.com/tanaylab/prego) — PWM motif scanning
- [`iceqream`](https://github.com/tanaylab/iceqream) — Quantitative CRE accessibility models
- [`shaman`](https://github.com/tanaylab/shaman) — Hi-C analysis (used for Fig 3D only)

Install GitHub packages with:
```r
# install.packages("remotes")
remotes::install_github("tanaylab/misha.ext")
remotes::install_github("tanaylab/tgutil")
remotes::install_github("tanaylab/prego")
remotes::install_github("tanaylab/iceqream")
remotes::install_github("tanaylab/shaman")
```

## Running the analysis

1. Clone this repository and download the data from the associated data repository.
2. Install the R and Python dependencies listed above.
3. Open the Jupyter notebooks in the `analysis/` directory. Each notebook is self-contained and loads data using `here()` relative paths from the repository root.
4. The notebooks are meant to be run with an R kernel (IRkernel) in Jupyter.

### Borzoi/Flashzoi training and in-silico genome analysis

The `borzoi/` directory contains all training code, configurations, and inference scripts for the deep learning models in Fig 5:

- `borzoi/code/` — Training, inference, and perturbation analysis Python code
- `borzoi/configs/` — 57 YAML training configurations across 5 experiment categories
- `borzoi/inference/` — Shell scripts for in-silico genome perturbation experiments

Pre-computed model predictions are provided in the misha tracks, so re-training is **not required** to reproduce the figures. Trained model checkpoints are available on HuggingFace (see `borzoi/README.md`).

## License

This code is provided for academic research purposes accompanying the manuscript.
