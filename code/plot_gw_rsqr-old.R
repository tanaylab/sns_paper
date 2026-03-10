#!/usr/bin/env Rscript
# Script to plot whole genome R^2 comparison across model configurations

library(here)
# library(misha)
devtools::load_all("/home/aviezerl/src/misha")
library(misha.ext)
library(data.table)
library(tidyverse)
library(ggrepel)
library(tgutil)

options(gmax.data.size = 1e10)
options(gmultitasking = FALSE)

# Setup paths
link_dir <- paste0(here(), "/")
setwd(here())
# gsetroot(paste0(link_dir, "data/mm10/"))
gsetroot("/home/aviezerl/mm10")
gdb.reload()
gdataset.load(paste0(link_dir, "data/mm10/"), force = T)

# Source utilities
source(paste0(link_dir, "code/seq2epi_utils.r"))
source(paste0(link_dir, "code/fig_fun.r"))

# Initialize pipeline and load data
mod <- init_pipe()
mod$gw$cg_trace = readRDS(here('data/cg_trace_mm10.rds'))
mod$gw$feats = readRDS(here('data/feats_mm10.rds'))
mod$gw$feats35 = readRDS(here('data/feats35_mm10.rds'))
mod$gw$feats_iqdn = readRDS(here('data/feats_iqdn_mm10.rds'))
cg_trace <- mod$gw$cg_trace

# Filter data
f_norp <- cg_trace$d_ltr != 0 &
          cg_trace$d_line != 0 &
          rowSums(is.na(cg_trace)) == 0 &
          cg_trace$start > 3e6
cg_trace_f <- cg_trace[f_norp, ]

# Setup train/test split based on Borzoi folds
borz_folds <- gintervals.load("borzoi.folds")
btrain <- borz_folds[borz_folds$type == "train", ]
btrainc <- gintervals.canonic(btrain)
gvtrack.create("borz_tr_d", btrainc, "distance")

test_chroms <- c("chr4", "chr10", "chr14", "chr15")

# Track configuration table with explicit genome information
# train_genome: genome used for model training
# infer_genome: genome used for inference/prediction
# is_k27: TRUE for H3K27me3 tracks (used in R^2 calculation), FALSE for K4 tracks
#
# Genome rules:
# - Receptive field series (brz*) + SNS: train=mm10, infer=mm10
# - Silicus models WITHOUT "finetune" in track name: train=mm10, infer=silicus
# - Silicus models WITH "finetune" in track name:
#   - If "mm10" in col_name: train=silicus, infer=mm10
#   - Otherwise: train=silicus, infer=silicus

track_tbl <- tribble(
    ~track_name, ~col_name, ~model_type, ~features, ~train_genome, ~infer_genome, ~is_k27,

    # ==========================================================================
    # Receptive Field Series: train=mm10, infer=mm10
    # ==========================================================================
    "seq.IQ.pcg.borzoi.finetune_norm_1k_EB4_cnt_cropped",
        "brz1k", "borzoi", "1k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_1k_EB4_cnt_k4_cropped",
        "brz1k_k4", "borzoi", "1k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_2k_EB4_cnt_cropped",
        "brz2k", "borzoi", "2k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_2k_EB4_cnt_k4_cropped",
        "brz2k_k4", "borzoi", "2k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_4k_EB4_cnt_cropped",
        "brz4k", "borzoi", "4k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_4k_EB4_cnt_k4_cropped",
        "brz4k_k4", "borzoi", "4k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_8k_EB4_cnt_cropped",
        "brz8k", "borzoi", "8k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_8k_EB4_cnt_k4_cropped",
        "brz8k_k4", "borzoi", "8k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_16k_EB4_cnt_cropped",
        "brz16k", "borzoi", "16k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_16k_EB4_cnt_k4_cropped",
        "brz16k_k4", "borzoi", "16k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_32k_EB4_cnt_cropped",
        "brz32k", "borzoi", "32k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_32k_EB4_cnt_k4_cropped",
        "brz32k_k4", "borzoi", "32k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_64k_EB4_cnt_cropped",
        "brz64k", "borzoi", "64k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_64k_EB4_cnt_k4_cropped",
        "brz64k_k4", "borzoi", "64k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_128k_EB4_cnt_cropped",
        "brz128k", "borzoi", "128k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_128k_EB4_cnt_k4_cropped",
        "brz128k_k4", "borzoi", "128k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_256k_EB4_cnt_cropped",
        "brz256k", "borzoi", "256k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_256k_EB4_cnt_k4_cropped",
        "brz256k_k4", "borzoi", "256k", "mm10", "mm10", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_EB4_cnt_cropped",
        "brz524k", "borzoi", "524k", "mm10", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_EB4_cnt_k4_cropped",
        "brz524k_k4", "borzoi", "524k", "mm10", "mm10", FALSE,

    # ==========================================================================
    # SNS Models: train=mm10, infer=mm10
    # ==========================================================================
    "jk.epipcg.pred.eb4_xgb_lm_seeds",
        "sns_lm", "sns", "linear", "mm10", "mm10", TRUE,
    "jk.epipcg.pred.eb4_xgb_brz2k_seeds",
        "sns_brz2k", "sns", "borzoi_2k", "mm10", "mm10", TRUE,

    # ==========================================================================
    # Silicus Models (no finetune): train=mm10, infer=silicus
    # ==========================================================================
    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_ctcf.finetune_norm_524k_EB4_cnt_cropped",
        "sil_ctcf", "silicus", "+CTCF", "mm10", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_ctcf.finetune_norm_524k_EB4_cnt_k4_cropped",
        "sil_ctcf_k4", "silicus", "+CTCF", "mm10", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250.finetune_norm_524k_EB4_cnt_cropped",
        "sil_cgd", "silicus", "+CGD", "mm10", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250.finetune_norm_524k_EB4_cnt_k4_cropped",
        "sil_cgd_k4", "silicus", "+CGD", "mm10", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250_with_ctcf.finetune_norm_524k_EB4_cnt_cropped",
        "sil_cgd_ctcf", "silicus", "+CGD+CTCF", "mm10", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250_with_ctcf.finetune_norm_524k_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_k4", "silicus", "+CGD+CTCF", "mm10", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_epi.finetune_norm_524k_EB4_cnt_cropped",
        "sil_cgd_cre_epi", "silicus", "+CGD+CRE(epi)", "mm10", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_epi.finetune_norm_524k_EB4_cnt_k4_cropped",
        "sil_cgd_cre_epi_k4", "silicus", "+CGD+CRE(epi)", "mm10", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_marginal.finetune_norm_524k_EB4_cnt_cropped",
        "sil_cgd_cre_mrg", "silicus", "+CGD+CRE(mrg)", "mm10", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_marginal.finetune_norm_524k_EB4_cnt_k4_cropped",
        "sil_cgd_cre_mrg_k4", "silicus", "+CGD+CRE(mrg)", "mm10", "silicus", FALSE,

    # ==========================================================================
    # Silicus Finetuned Models: train=silicus, infer varies
    # ==========================================================================
    # +CGD+CTCF variants
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_EB4_cnt_cropped",
        "sil_cgd_ctcf_finetune", "silicus_ft", "+CGD+CTCF", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_finetune_k4", "silicus_ft", "+CGD+CTCF", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_mm10_EB4_cnt_cropped",
        "sil_cgd_ctcf_finetune_mm10", "silicus_ft", "+CGD+CTCF", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_mm10_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_finetune_mm10_k4", "silicus_ft", "+CGD+CTCF", "silicus", "mm10", FALSE,

    # +CGD+CTCF+CRE variants
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_mouse_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_epi", "silicus_ft", "+CGD+CTCF+CRE", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_mouse_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_epi_k4", "silicus_ft", "+CGD+CTCF+CRE", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_epi_human", "silicus_ft", "+CGD+CTCF+CRE (human head)", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_epi_human_k4", "silicus_ft", "+CGD+CTCF+CRE (human head)", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_mouse_mm10_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_epi_mm10", "silicus_ft", "+CGD+CTCF+CRE", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_mouse_mm10_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_epi_mm10_k4", "silicus_ft", "+CGD+CTCF+CRE", "silicus", "mm10", FALSE,

    # +CGD+CTCF+CRE+Exons variants
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_mouse_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_exons_mouse", "silicus_ft", "+CGD+CTCF+CRE+Exons", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_mouse_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_exons_mouse_k4", "silicus_ft", "+CGD+CTCF+CRE+Exons", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_mm10_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_exons_mm10", "silicus_ft", "+CGD+CTCF+CRE+Exons", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_mm10_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_exons_mm10_k4", "silicus_ft", "+CGD+CTCF+CRE+Exons", "silicus", "mm10", FALSE,

    # +CGD+CTCF+CRE+Exons+Line+LTR variants
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_mouse_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_mouse", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_mouse_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_mouse_k4", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_mm10_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_mm10", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_mm10_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_mm10_k4", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR", "silicus", "mm10", FALSE,

    # +CGD+CTCF+CRE+Exons+Line+LTR+UTR3 variants
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mouse", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR+UTR3", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_k4", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR+UTR3", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_mm10_EB4_cnt_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mm10", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR+UTR3", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_mm10_EB4_cnt_k4_cropped",
        "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mm10_k4", "silicus_ft", "+CGD+CTCF+CRE+Exons+Line+LTR+UTR3", "silicus", "mm10", FALSE,

    # ==========================================================================
    # Silicus Base Models (no features): train=silicus, infer varies
    # ==========================================================================
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_mouse_EB4_cnt_cropped",
        "sil_base", "silicus_ft", "base", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_mouse_EB4_cnt_k4_cropped",
        "sil_base_k4", "silicus_ft", "base", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_mm10_EB4_cnt_cropped",
        "sil_base_mm10", "silicus_ft", "base", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_mm10_EB4_cnt_k4_cropped",
        "sil_base_mm10_k4", "silicus_ft", "base", "silicus", "mm10", FALSE,

    # ==========================================================================
    # Silicus GC Models: train=silicus, infer varies
    # ==========================================================================
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_gc_mouse_EB4_cnt_cropped",
        "sil_gc", "silicus_ft", "GC", "silicus", "silicus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_gc_mouse_EB4_cnt_k4_cropped",
        "sil_gc_k4", "silicus_ft", "GC", "silicus", "silicus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_gc_mouse_mm10_EB4_cnt_cropped",
        "sil_gc_mm10", "silicus_ft", "GC", "silicus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_silicus_gc_mouse_mm10_EB4_cnt_k4_cropped",
        "sil_gc_mm10_k4", "silicus_ft", "GC", "silicus", "mm10", FALSE,

    # ==========================================================================
    # Markov5 Models: train=markov5, infer varies
    # ==========================================================================
    "seq.IQ.pcg.borzoi.finetune_norm_524k_markov_no_repeats_mouse_EB4_cnt_cropped",
        "markov5", "markov", "Markov5", "markov5", "markov5", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_markov_no_repeats_mouse_EB4_cnt_k4_cropped",
        "markov5_k4", "markov", "Markov5", "markov5", "markov5", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_markov_no_repeats_mm10_EB4_cnt_cropped",
        "markov5_mm10", "markov", "Markov5", "markov5", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_markov_no_repeats_mm10_EB4_cnt_k4_cropped",
        "markov5_mm10_k4", "markov", "Markov5", "markov5", "mm10", FALSE,

    # ==========================================================================
    # Random Models: train=random, infer varies
    # ==========================================================================
    "seq.IQ.pcg.borzoi.finetune_norm_524k_random_mouse_EB4_cnt_cropped",
        "random", "random", "Random", "random", "random", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_random_mouse_EB4_cnt_k4_cropped",
        "random_k4", "random", "Random", "random", "random", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_random_mm10_EB4_cnt_cropped",
        "random_mm10", "random", "Random", "random", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_random_mm10_EB4_cnt_k4_cropped",
        "random_mm10_k4", "random", "Random", "random", "mm10", FALSE,

    # ==========================================================================
    # mm10 Minus Models (mm10 with silicus substitution): train=mm10Minus, infer varies
    # ==========================================================================
    # 10% silicus substitution
    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre10_EB4_cnt_cropped",
        "mm10minus_sil10", "mm10minus", "10% silicus", "mm10Minus", "mm10Minus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre10_EB4_cnt_k4_cropped",
        "mm10minus_sil10_k4", "mm10minus", "10% silicus", "mm10Minus", "mm10Minus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre10_mm10_EB4_cnt_cropped",
        "mm10minus_sil10_mm10", "mm10minus", "10% silicus", "mm10Minus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre10_mm10_EB4_cnt_k4_cropped",
        "mm10minus_sil10_mm10_k4", "mm10minus", "10% silicus", "mm10Minus", "mm10", FALSE,

    # 20% silicus substitution
    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre20_EB4_cnt_cropped",
        "mm10minus_sil20", "mm10minus", "20% silicus", "mm10Minus", "mm10Minus", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre20_EB4_cnt_k4_cropped",
        "mm10minus_sil20_k4", "mm10minus", "20% silicus", "mm10Minus", "mm10Minus", FALSE,

    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre20_mm10_EB4_cnt_cropped",
        "mm10minus_sil20_mm10", "mm10minus", "20% silicus", "mm10Minus", "mm10", TRUE,
    "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre20_mm10_EB4_cnt_k4_cropped",
        "mm10minus_sil20_mm10_k4", "mm10minus", "20% silicus", "mm10Minus", "mm10", FALSE

    # # 30% silicus substitution
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre30_EB4_cnt_cropped",
    #     "mm10minus_sil30", "mm10minus", "30% silicus", "mm10Minus", "mm10Minus", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre30_EB4_cnt_k4_cropped",
    #     "mm10minus_sil30_k4", "mm10minus", "30% silicus", "mm10Minus", "mm10Minus", FALSE,
    #
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre30_mm10_EB4_cnt_cropped",
    #     "mm10minus_sil30_mm10", "mm10minus", "30% silicus", "mm10Minus", "mm10", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre30_mm10_EB4_cnt_k4_cropped",
    #     "mm10minus_sil30_mm10_k4", "mm10minus", "30% silicus", "mm10Minus", "mm10", FALSE,

    # # 40% silicus substitution
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre40_EB4_cnt_cropped",
    #     "mm10minus_sil40", "mm10minus", "40% silicus", "mm10Minus", "mm10Minus", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre40_EB4_cnt_k4_cropped",
    #     "mm10minus_sil40_k4", "mm10minus", "40% silicus", "mm10Minus", "mm10Minus", FALSE,
    #
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre40_mm10_EB4_cnt_cropped",
    #     "mm10minus_sil40_mm10", "mm10minus", "40% silicus", "mm10Minus", "mm10", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre40_mm10_EB4_cnt_k4_cropped",
    #     "mm10minus_sil40_mm10_k4", "mm10minus", "40% silicus", "mm10Minus", "mm10", FALSE,

    # # 50% silicus substitution
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre50_EB4_cnt_cropped",
    #     "mm10minus_sil50", "mm10minus", "50% silicus", "mm10Minus", "mm10Minus", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre50_EB4_cnt_k4_cropped",
    #     "mm10minus_sil50_k4", "mm10minus", "50% silicus", "mm10Minus", "mm10Minus", FALSE,
    #
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre50_mm10_EB4_cnt_cropped",
    #     "mm10minus_sil50_mm10", "mm10minus", "50% silicus", "mm10Minus", "mm10", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre50_mm10_EB4_cnt_k4_cropped",
    #     "mm10minus_sil50_mm10_k4", "mm10minus", "50% silicus", "mm10Minus", "mm10", FALSE,

    # # 60% silicus substitution
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre60_EB4_cnt_cropped",
    #     "mm10minus_sil60", "mm10minus", "60% silicus", "mm10Minus", "mm10Minus", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre60_EB4_cnt_k4_cropped",
    #     "mm10minus_sil60_k4", "mm10minus", "60% silicus", "mm10Minus", "mm10Minus", FALSE,
    #
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre60_mm10_EB4_cnt_cropped",
    #     "mm10minus_sil60_mm10", "mm10minus", "60% silicus", "mm10Minus", "mm10", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre60_mm10_EB4_cnt_k4_cropped",
    #     "mm10minus_sil60_mm10_k4", "mm10minus", "60% silicus", "mm10Minus", "mm10", FALSE,

    # # 70% silicus substitution
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre70_EB4_cnt_cropped",
    #     "mm10minus_sil70", "mm10minus", "70% silicus", "mm10Minus", "mm10Minus", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre70_EB4_cnt_k4_cropped",
    #     "mm10minus_sil70_k4", "mm10minus", "70% silicus", "mm10Minus", "mm10Minus", FALSE,
    #
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre70_mm10_EB4_cnt_cropped",
    #     "mm10minus_sil70_mm10", "mm10minus", "70% silicus", "mm10Minus", "mm10", TRUE,
    # "seq.IQ.pcg.borzoi.finetune_norm_524k_mm10_minus_silicus_cg_ctcf_cre70_mm10_EB4_cnt_k4_cropped",
    #     "mm10minus_sil70_mm10_k4", "mm10minus", "70% silicus", "mm10Minus", "mm10", FALSE
)

# Generate display labels from track metadata
track_tbl <- track_tbl %>%
    mutate(
        genome_label = case_when(
            train_genome == "mm10" & infer_genome == "mm10" ~ "",
            train_genome == "mm10" & infer_genome == "silicus" ~ "[mm10->sil]",
            train_genome == "silicus" & infer_genome == "silicus" ~ "[sil->sil]",
            train_genome == "silicus" & infer_genome == "mm10" ~ "[sil->mm10]",
            train_genome == "markov5" & infer_genome == "markov5" ~ "[mrk->mrk]",
            train_genome == "markov5" & infer_genome == "mm10" ~ "[mrk->mm10]",
            train_genome == "random" & infer_genome == "random" ~ "[rnd->rnd]",
            train_genome == "random" & infer_genome == "mm10" ~ "[rnd->mm10]",
            train_genome == "mm10Minus" & infer_genome == "mm10Minus" ~ "[mm10-->mm10-]",
            train_genome == "mm10Minus" & infer_genome == "mm10" ~ "[mm10-->mm10]",
            TRUE ~ paste0("[", train_genome, "->", infer_genome, "]")
        ),
        label = case_when(
            model_type == "borzoi" ~ features,
            model_type == "sns" ~ paste0("SNS\n(", features, ")"),
            model_type %in% c("silicus", "silicus_ft", "markov", "random", "mm10minus") ~ paste0(features, "\n", genome_label)
        )
    )

# Track names for extraction
track_names <- track_tbl$track_name
col_names <- track_tbl$col_name

stopifnot(all(gtrack.exists(track_names)))
stopifnot(length(track_names) == length(unique(col_names)))

# Extract prediction tracks
message("Extracting prediction tracks...")
gext2 <- gextract(track_names,
    iterator = cg_trace_f, intervals = cg_trace_f,
    colnames = col_names
) %>% arrange(intervalID) %cache_rds% here("data/gext2.rds")
gext2[is.na(gext2)] <- 0

# Extract distance to training regions
cmp_gw <- gextract("borz_tr_d",
    intervals = cg_trace_f, iterator = cg_trace_f,
    colnames = "borz_tr_d"
)
cg_trace_f$borz_tr_d <- cmp_gw$borz_tr_d

# Create train/test indicator
# 0 = train chromosomes, within 524k of training
# 1 = train chromosomes, far from training
# 2 = test chromosomes, within 524k of training
# 3 = test chromosomes, far from training (strictest test)
train_mod_gw <- ifelse(cmp_gw$chrom %in% test_chroms, 2, 0) +
    ifelse(cmp_gw$borz_tr_d > 524000, 1, 0)

message("Train/test split:")
print(table(train_mod_gw))

# ==============================================================================
# Create filters for region-specific analyses
# ==============================================================================

# Non-repeat filter (regions NOT on repeats)
message("Creating non-repeat filter...")
gvtrack.create("rmsk_d", "intervs.global.rmsk", "distance")
rmsk_dist <- gextract("rmsk_d",
    intervals = cg_trace_f, iterator = cg_trace_f,
    colnames = "rmsk_d"
)
f_no_rmsk <- rmsk_dist$rmsk_d != 0
f_rmsk <- rmsk_dist$rmsk_d == 0  # Regions ON repeats
message("Regions not on repeats: ", sum(f_no_rmsk), " / ", length(f_no_rmsk))
message("Regions on repeats: ", sum(f_rmsk), " / ", length(f_rmsk))

# Exon filter
message("Creating exon filter...")
gvtrack.create("exon_d", "intervs.global.exon", "distance")
exon_dist <- gextract("exon_d",
    intervals = cg_trace_f, iterator = cg_trace_f,
    colnames = "exon_d"
)
f_exon <- exon_dist$exon_d == 0
message("Regions on exons: ", sum(f_exon), " / ", length(f_exon))

# CRE filter
message("Creating CRE filter...")
cre_peaks <- readRDS("/home/aviezerl/src/borzoi_finetune/data/peaks_multieb_marginal.rds")$peaks
cre_neighbors <- gintervals.neighbors1(cg_trace_f, cre_peaks)
f_cre <- cre_neighbors$dist == 0
message("Regions on CREs: ", sum(f_cre), " / ", length(f_cre))

# Promoter filter
message("Creating promoter filter...")
promoters <- misha.ext::get_promoters()
prom_neighbors <- gintervals.neighbors1(cg_trace_f, promoters)
f_promoter <- prom_neighbors$dist == 0
message("Regions on promoters: ", sum(f_promoter), " / ", length(f_promoter))

# K27 and K4 columns
k27_cols <- track_tbl %>% filter(is_k27) %>% pull(col_name)
k4_cols <- track_tbl %>% filter(!is_k27) %>% pull(col_name)

# gext2 <- gext2 %>% gintervals.neighbors("borzoi.test_clean") %>% 
    # mutate(type = ifelse(dist == 0 & chrom %in% test_chroms, "test", "train"))

# Calculate R^2 for a given set of columns and response
# Optional extra_filter allows subsetting (e.g., non-repeat regions)
calc_rsqr <- function(pred_cols, resp_vec, gext_data, train_mod, extra_filter = NULL) {
    cor_te <- c()
    cor_tr <- c()
    f_te <- train_mod == 3
    f_tr <- train_mod == 0

    # Apply extra filter if provided
    if (!is.null(extra_filter)) {
        f_te <- f_te & extra_filter
        f_tr <- f_tr & extra_filter
    }

    for (col in pred_cols) {
        pred <- as.numeric(gext_data[, col])
        cor_te <- c(cor_te, cor(pred[f_te], resp_vec[f_te], method = "pearson")^2)
        cor_tr <- c(cor_tr, cor(pred[f_tr], resp_vec[f_tr], method = "pearson")^2)
    }
    names(cor_te) <- pred_cols
    names(cor_tr) <- pred_cols

    list(test = cor_te, train = cor_tr)
}

stopifnot(all(gext2$chrom == cg_trace_f$chrom))
stopifnot(all(gext2$start == cg_trace_f$start))
stopifnot(all(gext2$end == cg_trace_f$end))


# Calculate R^2 for K27 predictions
message("Calculating K27 R^2 values...")
k27_rsqr <- calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw)
cor_k27_gw_te <- k27_rsqr$test
cor_k27_gw_tr <- k27_rsqr$train

# Extract K4 response from track (similar to K27)
message("Extracting K4 response...")
gvtrack.create("obs_k4", "jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3", "sum")
k4_resp <- gextract("obs_k4",
    intervals = cg_trace_f, iterator = cg_trace_f,
    colnames = "k4_1k"
)
k4_resp$lk4_1k <- log2(1 + k4_resp$k4_1k)

# Calculate R^2 for K4 predictions
message("Calculating K4 R^2 values...")
k4_rsqr <- calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw)
cor_k4_gw_te <- k4_rsqr$test
cor_k4_gw_tr <- k4_rsqr$train

# ==============================================================================
# Calculate R^2 for non-repeat regions (variant plot)
# ==============================================================================
message("Calculating K27 R^2 values (non-repeat regions)...")
k27_rsqr_norp <- calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw, f_no_rmsk)
cor_k27_norp_te <- k27_rsqr_norp$test
cor_k27_norp_tr <- k27_rsqr_norp$train

message("Calculating K4 R^2 values (non-repeat regions)...")
k4_rsqr_norp <- calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw, f_no_rmsk)
cor_k4_norp_te <- k4_rsqr_norp$test
cor_k4_norp_tr <- k4_rsqr_norp$train

# ==============================================================================
# Calculate R^2 for repeat regions
# ==============================================================================
message("Calculating K27 R^2 values (repeat regions)...")
k27_rsqr_rmsk <- calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw, f_rmsk)
cor_k27_rmsk_te <- k27_rsqr_rmsk$test
cor_k27_rmsk_tr <- k27_rsqr_rmsk$train

message("Calculating K4 R^2 values (repeat regions)...")
k4_rsqr_rmsk <- calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw, f_rmsk)
cor_k4_rmsk_te <- k4_rsqr_rmsk$test
cor_k4_rmsk_tr <- k4_rsqr_rmsk$train

# ==============================================================================
# Calculate R^2 for exon regions
# ==============================================================================
message("Calculating K27 R^2 values (exon regions)...")
k27_rsqr_exon <- calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw, f_exon)
cor_k27_exon_te <- k27_rsqr_exon$test
cor_k27_exon_tr <- k27_rsqr_exon$train

message("Calculating K4 R^2 values (exon regions)...")
k4_rsqr_exon <- calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw, f_exon)
cor_k4_exon_te <- k4_rsqr_exon$test
cor_k4_exon_tr <- k4_rsqr_exon$train

# ==============================================================================
# Calculate R^2 for CRE regions
# ==============================================================================
message("Calculating K27 R^2 values (CRE regions)...")
k27_rsqr_cre <- calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw, f_cre)
cor_k27_cre_te <- k27_rsqr_cre$test
cor_k27_cre_tr <- k27_rsqr_cre$train

message("Calculating K4 R^2 values (CRE regions)...")
k4_rsqr_cre <- calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw, f_cre)
cor_k4_cre_te <- k4_rsqr_cre$test
cor_k4_cre_tr <- k4_rsqr_cre$train

# ==============================================================================
# Calculate R^2 for promoter regions
# ==============================================================================
message("Calculating K27 R^2 values (promoter regions)...")
k27_rsqr_prom <- calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw, f_promoter)
cor_k27_prom_te <- k27_rsqr_prom$test
cor_k27_prom_tr <- k27_rsqr_prom$train

message("Calculating K4 R^2 values (promoter regions)...")
k4_rsqr_prom <- calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw, f_promoter)
cor_k4_prom_te <- k4_rsqr_prom$test
cor_k4_prom_tr <- k4_rsqr_prom$train

# ==============================================================================
# Plotting function for R^2 comparison
# ==============================================================================
create_rsqr_plot <- function(cor_te, cor_tr, track_tbl, is_k27_filter,
                              title, y_max = 0.88) {

    # Get relevant columns based on filter
    model_cols <- names(cor_te)

    # Prepare data for plotting
    cor_df <- data.frame(
        name = model_cols,
        Test = cor_te,
        Train = cor_tr
    )

    # Add metadata from track_tbl
    cor_df <- cor_df %>%
        left_join(
            track_tbl %>%
                filter(is_k27 == is_k27_filter) %>%
                select(col_name, model_type, features, train_genome, infer_genome, label),
            by = c("name" = "col_name")
        )

    cor_df_gg <- cor_df %>%
        pivot_longer(cols = c(Test, Train), names_to = "Set", values_to = "R2")

    # Define factor order for x-axis (grouped by model type and genome configuration)
    # For K4, we need to use the _k4 suffix versions
    if (is_k27_filter) {
        name_order <- c(
            # mm10 -> mm10 models (Borzoi RF series)
            "brz1k", "brz2k", "brz4k", "brz8k", "brz16k", "brz32k",
            "brz64k", "brz128k", "brz256k", "brz524k",
            # mm10 -> mm10 models (SNS)
            "sns_lm", "sns_brz2k",
            # mm10 -> silicus models
            "sil_ctcf", "sil_cgd", "sil_cgd_ctcf",
            "sil_cgd_cre_epi", "sil_cgd_cre_mrg",
            # silicus/markov/random -> silicus/markov/random models (synthetic -> synthetic)
            "sil_base", "sil_gc", "markov5", "random",
            "sil_cgd_ctcf_finetune", "sil_cgd_ctcf_cre_epi", "sil_cgd_ctcf_cre_epi_human",
            "sil_cgd_ctcf_cre_exons_mouse", "sil_cgd_ctcf_cre_exons_line_ltr_mouse",
            "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mouse",
            # silicus/markov/random -> mm10 models (synthetic -> real)
            "sil_base_mm10", "sil_gc_mm10", "markov5_mm10", "random_mm10",
            "sil_cgd_ctcf_finetune_mm10", "sil_cgd_ctcf_cre_epi_mm10",
            "sil_cgd_ctcf_cre_exons_mm10", "sil_cgd_ctcf_cre_exons_line_ltr_mm10",
            "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mm10",
            # mm10Minus -> mm10Minus models (minus -> minus)
            "mm10minus_sil10", "mm10minus_sil20",
            # mm10Minus -> mm10 models (minus -> real)
            "mm10minus_sil10_mm10", "mm10minus_sil20_mm10"
        )
        # Line series: RF series (solid) and mm10Minus series (dashed)
        line_names_rf <- c(
            "brz1k", "brz2k", "brz4k", "brz8k", "brz16k", "brz32k",
            "brz64k", "brz128k", "brz256k", "brz524k"
        )
        line_names_minus <- c(
            "mm10minus_sil10", "mm10minus_sil20"
        )
        line_names_minus_mm10 <- c(
            "mm10minus_sil10_mm10", "mm10minus_sil20_mm10"
        )
    } else {
        name_order <- c(
            # mm10 -> mm10 models
            "brz1k_k4", "brz2k_k4", "brz4k_k4", "brz8k_k4", "brz16k_k4", "brz32k_k4",
            "brz64k_k4", "brz128k_k4", "brz256k_k4", "brz524k_k4",
            # mm10 -> silicus models
            "sil_ctcf_k4", "sil_cgd_k4", "sil_cgd_ctcf_k4",
            "sil_cgd_cre_epi_k4", "sil_cgd_cre_mrg_k4",
            # silicus/markov/random -> silicus/markov/random models (synthetic -> synthetic)
            "sil_base_k4", "sil_gc_k4", "markov5_k4", "random_k4",
            "sil_cgd_ctcf_finetune_k4", "sil_cgd_ctcf_cre_epi_k4", "sil_cgd_ctcf_cre_epi_human_k4",
            "sil_cgd_ctcf_cre_exons_mouse_k4", "sil_cgd_ctcf_cre_exons_line_ltr_mouse_k4",
            "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_k4",
            # silicus/markov/random -> mm10 models (synthetic -> real)
            "sil_base_mm10_k4", "sil_gc_mm10_k4", "markov5_mm10_k4", "random_mm10_k4",
            "sil_cgd_ctcf_finetune_mm10_k4", "sil_cgd_ctcf_cre_epi_mm10_k4",
            "sil_cgd_ctcf_cre_exons_mm10_k4", "sil_cgd_ctcf_cre_exons_line_ltr_mm10_k4",
            "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mm10_k4",
            # mm10Minus -> mm10Minus models (minus -> minus)
            "mm10minus_sil10_k4", "mm10minus_sil20_k4",
            # mm10Minus -> mm10 models (minus -> real)
            "mm10minus_sil10_mm10_k4", "mm10minus_sil20_mm10_k4"
        )
        # Line series: RF series (solid) and mm10Minus series (dashed)
        line_names_rf <- c(
            "brz1k_k4", "brz2k_k4", "brz4k_k4", "brz8k_k4", "brz16k_k4", "brz32k_k4",
            "brz64k_k4", "brz128k_k4", "brz256k_k4", "brz524k_k4"
        )
        line_names_minus <- c(
            "mm10minus_sil10_k4", "mm10minus_sil20_k4"
        )
        line_names_minus_mm10 <- c(
            "mm10minus_sil10_mm10_k4", "mm10minus_sil20_mm10_k4"
        )
    }

    # Filter to only existing models
    name_order <- name_order[name_order %in% model_cols]

    # Create label map from track_tbl
    label_map <- track_tbl %>%
        filter(is_k27 == is_k27_filter, col_name %in% name_order) %>%
        distinct(col_name, label) %>%
        deframe()

    # Get genome info for each model in the actual display order
    track_genome_info <- track_tbl %>%
        filter(is_k27 == is_k27_filter, col_name %in% name_order) %>%
        select(col_name, train_genome, infer_genome)

    # Create ordered genome info matching the FILTERED name_order
    ordered_genome_info <- data.frame(col_name = name_order) %>%
        left_join(track_genome_info, by = "col_name") %>%
        filter(!is.na(train_genome))  # Remove any models not in track_tbl

    # Assign region category to each model based on genome configuration
    ordered_genome_info <- ordered_genome_info %>%
        mutate(
            pos = row_number(),
            region = case_when(
                train_genome == "mm10" & infer_genome == "mm10" ~ "mm10_mm10",
                train_genome == "mm10" & infer_genome == "silicus" ~ "mm10_sil",
                train_genome %in% c("silicus", "markov5", "random") &
                    infer_genome %in% c("silicus", "markov5", "random") ~ "silicus_silicus",
                train_genome %in% c("silicus", "markov5", "random") &
                    infer_genome == "mm10" ~ "silicus_mm10",
                train_genome == "mm10Minus" & infer_genome == "mm10Minus" ~ "minus_minus",
                train_genome == "mm10Minus" & infer_genome == "mm10" ~ "minus_mm10",
                TRUE ~ "other"
            )
        )

    # Calculate x-axis positions for region boundaries based on WHERE each region ends
    # in the actual ordered sequence (not just counting)
    find_last_pos <- function(region_name) {
        positions <- ordered_genome_info$pos[ordered_genome_info$region == region_name]
        if (length(positions) > 0) max(positions) else NA
    }

    # Get the last position of each region
    last_mm10_mm10 <- find_last_pos("mm10_mm10")
    last_mm10_sil <- find_last_pos("mm10_sil")
    last_silicus_silicus <- find_last_pos("silicus_silicus")
    last_silicus_mm10 <- find_last_pos("silicus_mm10")
    last_minus_minus <- find_last_pos("minus_minus")
    last_minus_mm10 <- find_last_pos("minus_mm10")

    # Calculate boundaries - each boundary is after the last model in that region
    # If a region doesn't exist, use the previous boundary
    x1_end <- if (!is.na(last_mm10_mm10)) last_mm10_mm10 + 0.5 else 0.5
    x2_end <- if (!is.na(last_mm10_sil)) last_mm10_sil + 0.5 else x1_end
    x3_end <- if (!is.na(last_silicus_silicus)) last_silicus_silicus + 0.5 else x2_end
    x4_end <- if (!is.na(last_silicus_mm10)) last_silicus_mm10 + 0.5 else x3_end
    x5_end <- if (!is.na(last_minus_minus)) last_minus_minus + 0.5 else x4_end
    x6_end <- if (!is.na(last_minus_mm10)) last_minus_mm10 + 0.5 else x5_end

    # Create the plot
    p <- cor_df_gg %>%
        mutate(
            name = factor(name, levels = name_order),
            r2_label = round(R2, 3)
        ) %>%
        filter(!is.na(name)) %>%
        ggplot(aes(x = name, y = R2, color = Set, group = Set)) +
        # Background rectangles for genome regions
        annotate("rect", xmin = 0.5, xmax = x1_end, ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = "#4DAF4A") +
        annotate("rect", xmin = x1_end, xmax = x2_end, ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = "#984EA3") +
        annotate("rect", xmin = x2_end, xmax = x3_end, ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = "#FF7F00") +
        annotate("rect", xmin = x3_end, xmax = x4_end, ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = "#E41A1C") +
        annotate("rect", xmin = x4_end, xmax = x5_end, ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = "#377EB8") +
        annotate("rect", xmin = x5_end, xmax = x6_end + 0.5, ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = "#A65628") +
        # Region labels at top
        annotate("text", x = (0.5 + x1_end) / 2, y = y_max * 0.94,
            label = "mm10 -> mm10", size = 3.2, color = "#2E7D32", fontface = "bold") +
        annotate("text", x = (x1_end + x2_end) / 2, y = y_max * 0.94,
            label = "mm10 -> silicus", size = 3.2, color = "#6A1B9A", fontface = "bold") +
        annotate("text", x = (x2_end + x3_end) / 2, y = y_max * 0.94,
            label = "silicus -> silicus", size = 3.2, color = "#E65100", fontface = "bold") +
        annotate("text", x = (x3_end + x4_end) / 2, y = y_max * 0.94,
            label = "silicus -> mm10", size = 3.2, color = "#B71C1C", fontface = "bold") +
        annotate("text", x = (x4_end + x5_end) / 2, y = y_max * 0.94,
            label = "minus -> minus", size = 3.2, color = "#1565C0", fontface = "bold") +
        annotate("text", x = (x5_end + x6_end + 0.5) / 2, y = y_max * 0.94,
            label = "minus -> mm10", size = 3.2, color = "#5D4037", fontface = "bold") +
        # Secondary labels
        annotate("text", x = (0.5 + 10.5) / 2, y = y_max * 0.90,
            label = "(Borzoi RF series)", size = 2.8, color = "gray50", fontface = "italic") +
        # Vertical dashed lines for each x position (from label to points)
        geom_vline(xintercept = 1:length(name_order), linetype = "dashed", color = "gray80", alpha = 0.3) +
        # RF series line (solid)
        geom_line(
            data = function(d) filter(d, name %in% line_names_rf),
            linewidth = 1,
            alpha = 0.7
        ) +
        # mm10Minus -> mm10Minus series line (dashed)
        geom_line(
            data = function(d) filter(d, name %in% line_names_minus),
            linewidth = 1,
            alpha = 0.7,
            linetype = "dashed"
        ) +
        # mm10Minus -> mm10 series line (dashed)
        geom_line(
            data = function(d) filter(d, name %in% line_names_minus_mm10),
            linewidth = 1,
            alpha = 0.7,
            linetype = "dashed"
        ) +
        geom_point(size = 3.5) +
        geom_text_repel(
            aes(label = r2_label),
            size = 2.8,
            max.overlaps = Inf,
            min.segment.length = 0,
            box.padding = 0.25,
            point.padding = 0.2,
            show.legend = FALSE
        ) +
        scale_color_manual(
            values = c("Train" = "#2166AC", "Test" = "#D6604D"),
            labels = c("Train" = "Train", "Test" = "Test")
        ) +
        scale_x_discrete(labels = label_map) +
        scale_y_continuous(
            limits = c(0, y_max),
            breaks = seq(0, floor(y_max * 10) / 10, 0.2),
            expand = c(0.02, 0)
        ) +
        labs(
            x = "Model Configuration",
            y = expression(R^2),
            color = "Dataset",
            title = title,
            subtitle = "Train/Infer genome shown above: mm10 = real mouse, silicus = synthetic genome"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5),
            axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
            axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y = element_text(size = 12),
            legend.position = "top",
            legend.title = element_text(face = "bold"),
            legend.text = element_text(size = 11),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
            plot.margin = margin(20, 20, 20, 20)
        )

    return(p)
}

# ==============================================================================
# Create K27 plot
# ==============================================================================
message("Creating K27 plot...")
p_k27 <- create_rsqr_plot(
    cor_k27_gw_te, cor_k27_gw_tr, track_tbl,
    is_k27_filter = TRUE,
    title = "Whole Genome H3K27me3 Prediction Performance",
    y_max = 0.88
)

# ==============================================================================
# Create K4 plot
# ==============================================================================
message("Creating K4 plot...")
p_k4 <- create_rsqr_plot(
    cor_k4_gw_te, cor_k4_gw_tr, track_tbl,
    is_k27_filter = FALSE,
    title = "Whole Genome H3K4me3 Prediction Performance",
    y_max = 0.88
) 

# ==============================================================================
# Save plots
# ==============================================================================

# Save K27 plot
output_file_k27 <- paste0(link_dir, "figures/gw_rsqr_k27_plot.pdf")
message("Saving K27 plot to: ", output_file_k27)
ggsave(output_file_k27, p_k27, width = 14, height = 10)

output_png_k27 <- paste0(link_dir, "figures/gw_rsqr_k27_plot.png")
ggsave(output_png_k27, p_k27, width = 14, height = 10, dpi = 300)

# Save K4 plot
output_file_k4 <- paste0(link_dir, "figures/gw_rsqr_k4_plot.pdf")
message("Saving K4 plot to: ", output_file_k4)
ggsave(output_file_k4, p_k4, width = 14, height = 10)

output_png_k4 <- paste0(link_dir, "figures/gw_rsqr_k4_plot.png")
ggsave(output_png_k4, p_k4, width = 14, height = 10, dpi = 300)

# Also save legacy filename for backwards compatibility
ggsave(paste0(link_dir, "figures/gw_rsqr_plot.pdf"), p_k27, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_plot.png"), p_k27, width = 14, height = 10, dpi = 300)

# ==============================================================================
# Create non-repeat variant plots
# ==============================================================================
message("Creating K27 non-repeat plot...")
p_k27_norp <- create_rsqr_plot(
    cor_k27_norp_te, cor_k27_norp_tr, track_tbl,
    is_k27_filter = TRUE,
    title = "Whole Genome H3K27me3 Prediction Performance (Non-Repeat Regions)",
    y_max = 0.88
)

message("Creating K4 non-repeat plot...")
p_k4_norp <- create_rsqr_plot(
    cor_k4_norp_te, cor_k4_norp_tr, track_tbl,
    is_k27_filter = FALSE,
    title = "Whole Genome H3K4me3 Prediction Performance (Non-Repeat Regions)",
    y_max = 0.88
)

# Save non-repeat K27 plot
output_file_k27_norp <- paste0(link_dir, "figures/gw_rsqr_k27_norp_plot.pdf")
message("Saving K27 non-repeat plot to: ", output_file_k27_norp)
ggsave(output_file_k27_norp, p_k27_norp, width = 14, height = 10)

output_png_k27_norp <- paste0(link_dir, "figures/gw_rsqr_k27_norp_plot.png")
ggsave(output_png_k27_norp, p_k27_norp, width = 14, height = 10, dpi = 300)

# Save non-repeat K4 plot
output_file_k4_norp <- paste0(link_dir, "figures/gw_rsqr_k4_norp_plot.pdf")
message("Saving K4 non-repeat plot to: ", output_file_k4_norp)
ggsave(output_file_k4_norp, p_k4_norp, width = 14, height = 10)

output_png_k4_norp <- paste0(link_dir, "figures/gw_rsqr_k4_norp_plot.png")
ggsave(output_png_k4_norp, p_k4_norp, width = 14, height = 10, dpi = 300)

# ==============================================================================
# Create repeat region plots
# ==============================================================================
message("Creating K27 repeat plot...")
p_k27_rmsk <- create_rsqr_plot(
    cor_k27_rmsk_te, cor_k27_rmsk_tr, track_tbl,
    is_k27_filter = TRUE,
    title = "H3K27me3 Prediction Performance (Repeat Regions)",
    y_max = 0.88
)

message("Creating K4 repeat plot...")
p_k4_rmsk <- create_rsqr_plot(
    cor_k4_rmsk_te, cor_k4_rmsk_tr, track_tbl,
    is_k27_filter = FALSE,
    title = "H3K4me3 Prediction Performance (Repeat Regions)",
    y_max = 0.88
)

ggsave(paste0(link_dir, "figures/gw_rsqr_k27_rmsk_plot.pdf"), p_k27_rmsk, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k27_rmsk_plot.png"), p_k27_rmsk, width = 14, height = 10, dpi = 300)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_rmsk_plot.pdf"), p_k4_rmsk, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_rmsk_plot.png"), p_k4_rmsk, width = 14, height = 10, dpi = 300)

# ==============================================================================
# Create exon region plots
# ==============================================================================
message("Creating K27 exon plot...")
p_k27_exon <- create_rsqr_plot(
    cor_k27_exon_te, cor_k27_exon_tr, track_tbl,
    is_k27_filter = TRUE,
    title = "H3K27me3 Prediction Performance (Exon Regions)",
    y_max = 0.88
)

message("Creating K4 exon plot...")
p_k4_exon <- create_rsqr_plot(
    cor_k4_exon_te, cor_k4_exon_tr, track_tbl,
    is_k27_filter = FALSE,
    title = "H3K4me3 Prediction Performance (Exon Regions)",
    y_max = 0.88
)

ggsave(paste0(link_dir, "figures/gw_rsqr_k27_exon_plot.pdf"), p_k27_exon, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k27_exon_plot.png"), p_k27_exon, width = 14, height = 10, dpi = 300)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_exon_plot.pdf"), p_k4_exon, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_exon_plot.png"), p_k4_exon, width = 14, height = 10, dpi = 300)

# ==============================================================================
# Create CRE region plots
# ==============================================================================
message("Creating K27 CRE plot...")
p_k27_cre <- create_rsqr_plot(
    cor_k27_cre_te, cor_k27_cre_tr, track_tbl,
    is_k27_filter = TRUE,
    title = "H3K27me3 Prediction Performance (CRE Regions)",
    y_max = 0.88
)

message("Creating K4 CRE plot...")
p_k4_cre <- create_rsqr_plot(
    cor_k4_cre_te, cor_k4_cre_tr, track_tbl,
    is_k27_filter = FALSE,
    title = "H3K4me3 Prediction Performance (CRE Regions)",
    y_max = 0.88
)

ggsave(paste0(link_dir, "figures/gw_rsqr_k27_cre_plot.pdf"), p_k27_cre, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k27_cre_plot.png"), p_k27_cre, width = 14, height = 10, dpi = 300)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_cre_plot.pdf"), p_k4_cre, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_cre_plot.png"), p_k4_cre, width = 14, height = 10, dpi = 300)

# ==============================================================================
# Create promoter region plots
# ==============================================================================
message("Creating K27 promoter plot...")
p_k27_prom <- create_rsqr_plot(
    cor_k27_prom_te, cor_k27_prom_tr, track_tbl,
    is_k27_filter = TRUE,
    title = "H3K27me3 Prediction Performance (Promoter Regions)",
    y_max = 0.88
)

message("Creating K4 promoter plot...")
p_k4_prom <- create_rsqr_plot(
    cor_k4_prom_te, cor_k4_prom_tr, track_tbl,
    is_k27_filter = FALSE,
    title = "H3K4me3 Prediction Performance (Promoter Regions)",
    y_max = 0.88
)

ggsave(paste0(link_dir, "figures/gw_rsqr_k27_prom_plot.pdf"), p_k27_prom, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k27_prom_plot.png"), p_k27_prom, width = 14, height = 10, dpi = 300)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_prom_plot.pdf"), p_k4_prom, width = 14, height = 10)
ggsave(paste0(link_dir, "figures/gw_rsqr_k4_prom_plot.png"), p_k4_prom, width = 14, height = 10, dpi = 300)

message("Done!")

# ==============================================================================
# Print summary statistics
# ==============================================================================
message("\n=== H3K27me3 R^2 Summary ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k27_gw_te, 3))
message("\nTrain set:")
print(round(cor_k27_gw_tr, 3))

message("\n=== H3K4me3 R^2 Summary ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k4_gw_te, 3))
message("\nTrain set:")
print(round(cor_k4_gw_tr, 3))

message("\n=== H3K27me3 R^2 Summary (Non-Repeat Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k27_norp_te, 3))
message("\nTrain set:")
print(round(cor_k27_norp_tr, 3))

message("\n=== H3K4me3 R^2 Summary (Non-Repeat Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k4_norp_te, 3))
message("\nTrain set:")
print(round(cor_k4_norp_tr, 3))

message("\n=== H3K27me3 R^2 Summary (Repeat Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k27_rmsk_te, 3))
message("\nTrain set:")
print(round(cor_k27_rmsk_tr, 3))

message("\n=== H3K4me3 R^2 Summary (Repeat Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k4_rmsk_te, 3))
message("\nTrain set:")
print(round(cor_k4_rmsk_tr, 3))

message("\n=== H3K27me3 R^2 Summary (Exon Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k27_exon_te, 3))
message("\nTrain set:")
print(round(cor_k27_exon_tr, 3))

message("\n=== H3K4me3 R^2 Summary (Exon Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k4_exon_te, 3))
message("\nTrain set:")
print(round(cor_k4_exon_tr, 3))

message("\n=== H3K27me3 R^2 Summary (CRE Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k27_cre_te, 3))
message("\nTrain set:")
print(round(cor_k27_cre_tr, 3))

message("\n=== H3K4me3 R^2 Summary (CRE Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k4_cre_te, 3))
message("\nTrain set:")
print(round(cor_k4_cre_tr, 3))

message("\n=== H3K27me3 R^2 Summary (Promoter Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k27_prom_te, 3))
message("\nTrain set:")
print(round(cor_k27_prom_tr, 3))

message("\n=== H3K4me3 R^2 Summary (Promoter Regions) ===")
message("Test set (held out chromosomes and regions):")
print(round(cor_k4_prom_te, 3))
message("\nTrain set:")
print(round(cor_k4_prom_tr, 3))
