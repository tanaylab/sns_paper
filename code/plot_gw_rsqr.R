#!/usr/bin/env Rscript

# ==============================================================================
# 1. Libraries & Setup
# ==============================================================================

library(here)
library(misha)
library(misha.ext)
library(tidyverse)
library(ggrepel)
library(tgutil)
library(doMC)

options(gmax.data.size = 1e10)
options(gmultitasking = FALSE)

# Source utility functions
source(here("code", "borzoi_utils.R"))

# ==============================================================================
# 2. Track Configuration
# ==============================================================================
# Note: Utility functions (load_genome_variant_index, slugify_label,
# add_nomenclature_columns) are now in borzoi_utils.R

#' Build track table from organized group definitions
#'
#' To add a new track:
#' - For Borzoi RF series: add size to borzoi_rf_sizes
#' - For SNS models: add row to sns_tracks tibble
#' - For other models: add row to the appropriate section
build_track_table <- function() {

    # --- Borzoi RF Series (mm10 -> mm10) ---
    # Just specify the sizes, track names are auto-generated
    borzoi_rf_sizes <- c("1k", "2k", "4k", "8k", "16k", "32k", "64k", "128k", "256k", "524k")

    borzoi_rf <- tibble(
        size = borzoi_rf_sizes,
        config_order = seq_along(borzoi_rf_sizes)
    ) %>%
        mutate(
            # K27 tracks
            track_name_k27 = paste0("seq.IQ.pcg.borzoi.finetune_norm_", size, "_EB4_cnt_cropped"),
            col_name_k27 = paste0("brz", size),
            # K4 tracks
            track_name_k4 = paste0("seq.IQ.pcg.borzoi.finetune_norm_", size, "_EB4_cnt_k4_cropped"),
            col_name_k4 = paste0("brz", size, "_k4"),
            # Metadata
            model_type = "borzoi",
            features = size,
            train_genome = "mm10",
            infer_genome = "mm10"
        )

    # Pivot to long format (one row per track)
    borzoi_rf_long <- bind_rows(
        borzoi_rf %>%
            select(track_name = track_name_k27, col_name = col_name_k27,
                   model_type, features, train_genome, infer_genome, config_order) %>%
            mutate(order_family = "borzoi_rf",
                   train_variant_key = "mm10",
                   infer_variant_key = "mm10") %>%
            mutate(is_k27 = TRUE),
        borzoi_rf %>%
            select(track_name = track_name_k4, col_name = col_name_k4,
                   model_type, features, train_genome, infer_genome, config_order) %>%
            mutate(order_family = "borzoi_rf",
                   train_variant_key = "mm10",
                   infer_variant_key = "mm10") %>%
            mutate(is_k27 = FALSE)
    )

    # --- From Scratch 524k Models (mm10 -> mm10) ---
    from_scratch_524k <- bind_rows(
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_EB4_cnt_cropped",
            col_name = "from_scratch_524k",
            model_type = "from_scratch",
            features = "524k (lr=6e-5)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 1,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = TRUE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_EB4_cnt_k4_cropped",
            col_name = "from_scratch_524k_k4",
            model_type = "from_scratch",
            features = "524k (lr=6e-5)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 1,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = FALSE
        ),
        # Added 4k tracks
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_4k_EB4_cnt_cropped",
            col_name = "from_scratch_4k",
            model_type = "from_scratch",
            features = "4k (lr=6e-5)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 0,
            order_family = "from_scratch_4k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = TRUE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_4k_EB4_cnt_k4_cropped",
            col_name = "from_scratch_4k_k4",
            model_type = "from_scratch",
            features = "4k (lr=6e-5)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 0,
            order_family = "from_scratch_4k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = FALSE
        ),
        # Added 524k lr=1e-4
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus4_EB4_cnt_cropped",
            col_name = "from_scratch_524k_lr1eminus4",
            model_type = "from_scratch",
            features = "524k (lr=1e-4)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 1.5,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = TRUE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus4_EB4_cnt_k4_cropped",
            col_name = "from_scratch_524k_lr1eminus4_k4",
            model_type = "from_scratch",
            features = "524k (lr=1e-4)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 1.5,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = FALSE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus5_EB4_cnt_cropped",
            col_name = "from_scratch_524k_lr1eminus5",
            model_type = "from_scratch",
            features = "524k (lr=1e-5)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 2,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = TRUE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus5_EB4_cnt_k4_cropped",
            col_name = "from_scratch_524k_lr1eminus5_k4",
            model_type = "from_scratch",
            features = "524k (lr=1e-5)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 2,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = FALSE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus5_EB4_cnt_val_genome_wide_pearson",
            col_name = "from_scratch_524k_lr1eminus5_val_genome_wide_pearson",
            model_type = "from_scratch",
            features = "524k (lr=1e-5, gw-pearson)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 3,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = TRUE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus5_EB4_cnt_k4_val_genome_wide_pearson",
            col_name = "from_scratch_524k_lr1eminus5_k4_val_genome_wide_pearson",
            model_type = "from_scratch",
            features = "524k (lr=1e-5, gw-pearson)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 3,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = FALSE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus5_EB4_cnt_val_iou",
            col_name = "from_scratch_524k_lr1eminus5_val_iou",
            model_type = "from_scratch",
            features = "524k (lr=1e-5, val-iou)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 4,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = TRUE
        ),
        tibble(
            track_name = "seq.IQ.pcg.flashzoi.from_scratch_norm_524k_lr_1eminus5_EB4_cnt_k4_val_iou",
            col_name = "from_scratch_524k_lr1eminus5_k4_val_iou",
            model_type = "from_scratch",
            features = "524k (lr=1e-5, val-iou)",
            train_genome = "mm10",
            infer_genome = "mm10",
            config_order = 4,
            order_family = "from_scratch_524k",
            train_variant_key = "mm10",
            infer_variant_key = "mm10",
            is_k27 = FALSE
        )
    )

    # --- SNS Models (mm10 -> mm10, K27 only) ---
    sns_tracks <- tibble(
        track_name = c("jk.epipcg.pred.eb4_xgb_lm_seeds", "jk.epipcg.pred.eb4_xgb_brz2k_seeds"),
        col_name = c("sns_lm", "sns_brz2k"),
        features = c("linear", "borzoi_2k"),
        model_type = "sns",
        train_genome = "mm10",
        infer_genome = "mm10",
        is_k27 = TRUE,
        config_order = seq_len(2),
        order_family = "sns",
        train_variant_key = "mm10",
        infer_variant_key = "mm10"
    )

    # --- Silicus Models (mm10 -> silicus) ---
    # Models trained on mm10, inference on silicus genome
    silicus_mm10_to_sil <- tribble(
        ~col_prefix, ~features, ~track_suffix, ~variant_key,
        "sil_ctcf", "+CTCF", "mus_silicus_cg_gc_lower_with_ctcf", "mus_silicus_cg_gc_lower_with_ctcf",
        "sil_cgd", "+CGD", "mus_silicus_cg_gc_lower_with_CGD_exp250", "mus_silicus_cg_gc_lower_with_CGD_exp250",
        "sil_cgd_ctcf", "+CGD+CTCF", "mus_silicus_cg_gc_lower_with_CGD_exp250_with_ctcf", "mus_silicus_cg_gc_lower_with_CGD_exp250_with_ctcf",
        "sil_cgd_cre_epi", "+CGD+CRE(epi)", "mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_epi", "mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_epi",
        "sil_cgd_cre_mrg", "+CGD+CRE(mrg)", "mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_marginal", "mus_silicus_cg_gc_lower_with_CGD_exp250_with_cre_marginal"
    ) %>%
        mutate(
            config_order = row_number(),
            model_type = "silicus",
            train_genome = "mm10",
            infer_genome = "silicus"
        )

    silicus_mm10_to_sil_long <- bind_rows(
        silicus_mm10_to_sil %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.", track_suffix, ".finetune_norm_524k_EB4_cnt_cropped"),
                col_name = col_prefix,
                is_k27 = TRUE,
                order_family = "mm10_silicus",
                train_variant_key = "mm10",
                infer_variant_key = variant_key
            ) %>% select(-col_prefix, -track_suffix),
        silicus_mm10_to_sil %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.", track_suffix, ".finetune_norm_524k_EB4_cnt_k4_cropped"),
                col_name = paste0(col_prefix, "_k4"),
                is_k27 = FALSE,
                order_family = "mm10_silicus",
                train_variant_key = "mm10",
                infer_variant_key = variant_key
            ) %>% select(-col_prefix, -track_suffix)
    )

    # --- Silicus Finetuned Models ---
    # Define finetuned silicus models with their configurations
    silicus_ft_configs <- tribble(
        ~col_prefix, ~features, ~track_infix, ~infer_genome, ~variant_key,
        # +CGD+CTCF variants
        "sil_cgd_ctcf_finetune", "+CGD+CTCF", "silicus_cgd_ctcf", "silicus", "silicus_cgd_ctcf",
        "sil_cgd_ctcf_finetune_mm10", "+CGD+CTCF", "silicus_cgd_ctcf_mm10", "mm10", "silicus_cgd_ctcf",
        # +CGD+CTCF+CRE variants
        "sil_cgd_ctcf_cre_epi", "+CGD+CTCF+CRE", "silicus_cgd_ctcf_cre_mouse", "silicus", "silicus_cgd_ctcf_cre_mouse",
        "sil_cgd_ctcf_cre_epi_human", "+CGD+CTCF+CRE (human head)", "silicus_cgd_ctcf_cre", "silicus", "silicus_cgd_ctcf_cre",
        "sil_cgd_ctcf_cre_epi_mm10", "+CGD+CTCF+CRE", "silicus_cgd_ctcf_cre_mouse_mm10", "mm10", "silicus_cgd_ctcf_cre_mouse",
        # +CGD+CTCF+CRE+Exons variants
        "sil_cgd_ctcf_cre_exons_mouse", "+CGD+CTCF+CRE+Exons", "silicus_cgd_ctcf_cre_exons_mouse", "silicus", "silicus_cgd_ctcf_cre_exons_mouse",
        "sil_cgd_ctcf_cre_exons_mm10", "+CGD+CTCF+CRE+Exons", "silicus_cgd_ctcf_cre_exons_mm10", "mm10", "silicus_cgd_ctcf_cre_exons_mouse",
        # +CGD+CTCF+CRE+Exons+Line+LTR variants
        "sil_cgd_ctcf_cre_exons_line_ltr_mouse", "+CGD+CTCF+CRE+Exons+Line+LTR", "silicus_cgd_ctcf_cre_exons_line_ltr_mouse", "silicus", "silicus_cgd_ctcf_cre_exons_line_ltr_mouse",
        "sil_cgd_ctcf_cre_exons_line_ltr_mm10", "+CGD+CTCF+CRE+Exons+Line+LTR", "silicus_cgd_ctcf_cre_exons_line_ltr_mm10", "mm10", "silicus_cgd_ctcf_cre_exons_line_ltr_mouse",
        # +CGD+CTCF+CRE+Exons+Line+LTR+UTR3 variants
        "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mouse", "+CGD+CTCF+CRE+Exons+Line+LTR+UTR3", "silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse", "silicus", "silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse",
        "sil_cgd_ctcf_cre_exons_line_ltr_utr3_mm10", "+CGD+CTCF+CRE+Exons+Line+LTR+UTR3", "silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse_mm10", "mm10", "silicus_cgd_ctcf_cre_exons_line_ltr_utr3_mouse"
    ) %>%
        mutate(
            config_order = row_number(),
            model_type = "silicus_ft",
            train_genome = "silicus"
        )

    silicus_ft_long <- bind_rows(
        silicus_ft_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_cropped"),
                col_name = col_prefix,
                is_k27 = TRUE,
                order_family = "silicus_ft",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix),
        silicus_ft_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_k4_cropped"),
                col_name = paste0(col_prefix, "_k4"),
                is_k27 = FALSE,
                order_family = "silicus_ft",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix)
    )

    # --- Silicus Base and GC Models ---
    silicus_base_gc <- tribble(
        ~col_prefix, ~features, ~track_infix, ~infer_genome, ~variant_key,
        "sil_base", "base", "silicus_mouse", "silicus", "silicus_mouse",
        "sil_base_mm10", "base", "silicus_mm10", "mm10", "silicus_mouse",
        "sil_gc", "GC", "silicus_gc_mouse", "silicus", "silicus_gc_mouse",
        "sil_gc_mm10", "GC", "silicus_gc_mouse_mm10", "mm10", "silicus_gc_mouse"
    ) %>%
        mutate(
            config_order = row_number(),
            core_order = if_else(features == "base", 1, 2),
            model_type = "silicus_ft",
            train_genome = "silicus"
        )

    silicus_base_gc_long <- bind_rows(
        silicus_base_gc %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_cropped"),
                col_name = col_prefix,
                is_k27 = TRUE,
                order_family = "silicus_core",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix),
        silicus_base_gc %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_k4_cropped"),
                col_name = paste0(col_prefix, "_k4"),
                is_k27 = FALSE,
                order_family = "silicus_core",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix)
    )

    # --- Markov5 Models ---
    markov_configs <- tribble(
        ~col_prefix, ~track_infix, ~infer_genome, ~variant_key,
        "markov5", "markov_no_repeats_mouse", "markov5", "markov5",
        "markov5_mm10", "markov_no_repeats_mm10", "mm10", "markov5"
    ) %>%
        mutate(
            config_order = row_number(),
            core_order = 3,
            model_type = "markov",
            features = "Markov5",
            train_genome = "markov5"
        )

    markov_long <- bind_rows(
        markov_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_cropped"),
                col_name = col_prefix,
                is_k27 = TRUE,
                order_family = "markov",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix),
        markov_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_k4_cropped"),
                col_name = paste0(col_prefix, "_k4"),
                is_k27 = FALSE,
                order_family = "markov",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix)
    )

    # --- Random Models ---
    random_configs <- tribble(
        ~col_prefix, ~track_infix, ~infer_genome, ~variant_key,
        "random", "random_mouse", "random", "random",
        "random_mm10", "random_mm10", "mm10", "random"
    ) %>%
        mutate(
            config_order = row_number(),
            core_order = 4,
            model_type = "random",
            features = "Random",
            train_genome = "random"
        )

    random_long <- bind_rows(
        random_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_cropped"),
                col_name = col_prefix,
                is_k27 = TRUE,
                order_family = "random",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix),
        random_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_k4_cropped"),
                col_name = paste0(col_prefix, "_k4"),
                is_k27 = FALSE,
                order_family = "random",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix)
    )

    # --- mm10 Minus Models (mm10 with silicus substitution) ---
    mm10minus_configs <- tribble(
        ~col_prefix, ~features, ~track_infix, ~infer_genome, ~variant_key,
        "mm10minus_sil10", "10% silicus", "mm10_minus_silicus_cg_ctcf_cre10", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre10",
        "mm10minus_sil10_mm10", "10% silicus", "mm10_minus_silicus_cg_ctcf_cre10_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre10",
        "mm10minus_sil20", "20% silicus", "mm10_minus_silicus_cg_ctcf_cre20", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre20",
        "mm10minus_sil20_mm10", "20% silicus", "mm10_minus_silicus_cg_ctcf_cre20_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre20",
        "mm10minus_sil30", "30% silicus", "mm10_minus_silicus_cg_ctcf_cre30", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre30",
        "mm10minus_sil30_mm10", "30% silicus", "mm10_minus_silicus_cg_ctcf_cre30_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre30",
        "mm10minus_sil50", "50% silicus", "mm10_minus_silicus_cg_ctcf_cre50", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre50",
        "mm10minus_sil50_mm10", "50% silicus", "mm10_minus_silicus_cg_ctcf_cre50_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre50",
        "mm10minus_sil60", "60% silicus", "mm10_minus_silicus_cg_ctcf_cre60", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre60",
        "mm10minus_sil60_mm10", "60% silicus", "mm10_minus_silicus_cg_ctcf_cre60_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre60",
        "mm10minus_sil90", "90% silicus", "mm10_minus_silicus_cg_ctcf_cre90", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre90",
        "mm10minus_sil90_mm10", "90% silicus", "mm10_minus_silicus_cg_ctcf_cre90_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre90",
        "mm10minus_sil100", "100% silicus", "mm10_minus_silicus_cg_ctcf_cre100", "mm10Minus", "mm10_minus_silicus_cg_ctcf_cre100",
        "mm10minus_sil100_mm10", "100% silicus", "mm10_minus_silicus_cg_ctcf_cre100_mm10", "mm10", "mm10_minus_silicus_cg_ctcf_cre100",
    ) %>%
        mutate(
            config_order = row_number(),
            model_type = "mm10minus",
            train_genome = "mm10Minus"
        )

    mm10minus_long <- bind_rows(
        mm10minus_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_cropped"),
                col_name = col_prefix,
                is_k27 = TRUE,
                order_family = "mm10minus",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix),
        mm10minus_configs %>%
            mutate(
                track_name = paste0("seq.IQ.pcg.borzoi.finetune_norm_524k_", track_infix, "_EB4_cnt_k4_cropped"),
                col_name = paste0(col_prefix, "_k4"),
                is_k27 = FALSE,
                order_family = "mm10minus",
                train_variant_key = variant_key,
                infer_variant_key = if_else(infer_genome == "mm10", "mm10", variant_key)
            ) %>% select(-col_prefix, -track_infix)
    )

    # --- Combine all tracks ---
    track_tbl <- bind_rows(
        borzoi_rf_long,
        from_scratch_524k,
        sns_tracks,
        silicus_mm10_to_sil_long,
        silicus_base_gc_long,
        markov_long,
        random_long,
        silicus_ft_long,
        mm10minus_long
    )

    track_tbl <- track_tbl %>%
        mutate(
            display_region = case_when(
                model_type == "mm10minus" & infer_genome == "mm10Minus" ~ "minus_minus",
                model_type == "mm10minus" & infer_genome == "mm10" ~ "minus_mm10",
                train_genome == "mm10" & infer_genome == "mm10" ~ "mm10_mm10",
                train_genome == "mm10" & infer_genome == "silicus" ~ "mm10_silicus",
                train_genome %in% c("silicus", "markov5", "random") &
                    infer_genome %in% c("silicus", "markov5", "random") ~ "synth_synth",
                train_genome %in% c("silicus", "markov5", "random") &
                    infer_genome == "mm10" ~ "synth_mm10",
                TRUE ~ "other"
            ),
            region_order = case_when(
                display_region == "mm10_mm10" ~ 1,
                display_region == "mm10_silicus" ~ 2,
                display_region == "synth_synth" ~ 3,
                display_region == "synth_mm10" ~ 4,
                display_region == "minus_minus" ~ 5,
                display_region == "minus_mm10" ~ 6,
                TRUE ~ 7
            ),
            order_group = case_when(
                display_region == "mm10_mm10" & order_family == "borzoi_rf" ~ 1,
                display_region == "mm10_mm10" & order_family == "from_scratch_524k" ~ 1.5,
                display_region == "mm10_mm10" & order_family == "from_scratch_4k" ~ 1.5,
                display_region == "mm10_mm10" & order_family == "sns" ~ 2,
                display_region %in% c("synth_synth", "synth_mm10") &
                    order_family %in% c("silicus_core", "markov", "random") ~ 1,
                display_region %in% c("synth_synth", "synth_mm10") &
                    order_family == "silicus_ft" ~ 2,
                TRUE ~ 1
            ),
            order_within_group = coalesce(core_order, config_order)
        )

    # --- Generate display labels ---
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
                train_genome == "mm10Minus" & infer_genome == "mm10Minus" ~ "[mm10Minus->mm10-]",
                train_genome == "mm10Minus" & infer_genome == "mm10" ~ "[mm10Minus->mm10]",
                TRUE ~ paste0("[", train_genome, "->", infer_genome, "]")
            ),
            label = case_when(
                model_type == "borzoi" ~ features,
                model_type == "from_scratch" ~ paste0(features, "(from scratch)"),
                model_type == "sns" ~ paste0("SNS\n(", features, ")"),
                model_type %in% c("silicus", "silicus_ft", "markov", "random", "mm10minus") ~
                    paste0(features, "\n", genome_label)
            )
        )

    variant_index <- load_genome_variant_index()
    track_tbl <- add_nomenclature_columns(track_tbl, variant_index)

    return(track_tbl)
}

# ==============================================================================
# 3. Filter Definitions
# ==============================================================================

#' Create all genomic region filters
#'
#' To add a new filter: add an entry to the returned list with:
#'   - filter: logical vector (TRUE = include region)
#'   - label: display label for plot title
#'   - suffix: file name suffix
create_filters <- function(cg_trace_f) {
    message("Creating genomic region filters...")

    # Create virtual tracks for distance calculations
    gvtrack.create("rmsk_d", "intervs.global.rmsk", "distance")
    gvtrack.create("exon_d", "intervs.global.exon", "distance")

    # Extract distances in batch
    dists <- gextract(c("rmsk_d", "exon_d"),
        intervals = cg_trace_f, iterator = cg_trace_f,
        colnames = c("rmsk_d", "exon_d")
    )

    # Load CRE and promoter data
    cre_peaks <- readRDS(here("data", "files", "peaks_multieb_marginal.rds"))$peaks
    promoters <- misha.ext::get_promoters(upstream = 2000, downstream = 250)

    cre_neighbors <- gintervals.neighbors1(cg_trace_f, cre_peaks)
    prom_neighbors <- gintervals.neighbors1(cg_trace_f, promoters)

    cgd <- tgutil::fread(here("data", "files", "mm10_cgdom.csv")) |>
        mutate(start = start - 250, end = end + 250)
    cgd_neighbors <- gintervals.neighbors1(cg_trace_f, cgd)    

    filters <- list(
        all = list(
            filter = rep(TRUE, nrow(cg_trace_f)),
            label = "",
            suffix = ""
        ),
        norp = list(
            filter = dists$rmsk_d != 0,
            label = "(Non-Repeat Regions)",
            suffix = "_norp"
        ),
        rmsk = list(
            filter = dists$rmsk_d == 0,
            label = "(Repeat Regions)",
            suffix = "_rmsk"
        ),
        exon = list(
            filter = dists$exon_d == 0,
            label = "(Exon Regions)",
            suffix = "_exon"
        ),
        cre = list(
            filter = cre_neighbors$dist == 0,
            label = "(CRE Regions)",
            suffix = "_cre"
        ),
        prom = list(
            filter = prom_neighbors$dist == 0,
            label = "(Promoter Regions)",
            suffix = "_prom"
        ),
        cgd = list(
            filter = cgd_neighbors$dist == 0,
            label = "(CGD Regions)",
            suffix = "_cgd"
        )
    )

    # Log filter statistics
    for (fname in names(filters)) {
        n_pass <- sum(filters[[fname]]$filter)
        message(sprintf("  %s: %d / %d regions", fname, n_pass, length(filters[[fname]]$filter)))
    }

    return(filters)
}

# ==============================================================================
# 4. Data Loading
# ==============================================================================

setup_misha <- function() {
    link_dir <- paste0(here(), "/")
    setwd(here())
    gsetroot(here("data", "mm10"))
    gdb.reload()
    gdataset.load(paste0(link_dir, "data/mm10/"), force = TRUE)

    return(link_dir)
}

# Legacy function for backwards compatibility
# This now only loads what we actually need
load_data <- function(link_dir) {
    cg_trace_f <- load_cg_trace_filtered()
    return(list(mod = NULL, cg_trace_f = cg_trace_f))
}

setup_train_test_split <- function(cg_trace_f) {
    # Setup train/test split based on Borzoi folds
    borz_folds <- gintervals.load("borzoi.folds")
    btrain <- borz_folds[borz_folds$type == "train", ]
    btrainc <- gintervals.canonic(btrain)
    gvtrack.create("borz_tr_d", btrainc, "distance")

    test_chroms <- c("chr4", "chr10", "chr14", "chr15")

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

    return(list(cg_trace_f = cg_trace_f, train_mod_gw = train_mod_gw))
}

# ==============================================================================
# 5. R^2 Calculation
# ==============================================================================

#' Calculate R^2 for a given set of columns and response
calc_rsqr <- function(pred_cols, resp_vec, gext_data, train_mod, extra_filter = NULL) {
    cor_te <- c()
    cor_tr <- c()
    f_te <- train_mod == 3
    f_tr <- train_mod == 0

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

#' Calculate R^2 for all filter x mark combinations
calc_all_rsqr <- function(gext2, cg_trace_f, k4_resp, track_tbl, train_mod_gw, filters) {
    k27_cols <- track_tbl %>% filter(is_k27) %>% pull(col_name)
    k4_cols <- track_tbl %>% filter(!is_k27) %>% pull(col_name)

    results <- list()
    for (fname in names(filters)) {
        message(sprintf("Calculating R^2 for filter: %s", fname))
        f <- filters[[fname]]$filter
        results[[fname]] <- list(
            k27 = calc_rsqr(k27_cols, cg_trace_f$lk27_1k, gext2, train_mod_gw, f),
            k4 = calc_rsqr(k4_cols, k4_resp$lk4_1k, gext2, train_mod_gw, f)
        )
    }

    return(results)
}

# ==============================================================================
# 5b. Peak Overlap Calculation (using pre-extracted data)
# ==============================================================================

#' Calculate peak overlap for all tracks 
#'
#' @param gext_data Extracted data frame with all prediction columns
#' @param obs_k27_vec Observed K27 values vector
#' @param obs_k4_vec Observed K4 values vector
#' @param track_tbl Track table
#' @param filters Filter definitions (same as for R^2)
#' @param T_q Quantile threshold (default 0.98)
calc_all_peak_overlap <- function(gext_data, obs_k27_vec, obs_k4_vec, track_tbl, filters,
                                   T_q = 0.98) {
    require(matrixStats)

    k27_cols <- track_tbl %>% filter(is_k27) %>% pull(col_name)
    k4_cols <- track_tbl %>% filter(!is_k27) %>% pull(col_name)
    chrom_vec <- gext_data$chrom
    test_chroms <- c("chr4", "chr10", "chr14", "chr15")

    # Pre-extract full matrices once (avoid repeated subsetting)
    k27_mat_full <- as.matrix(gext_data[, k27_cols, drop = FALSE])
    k4_mat_full <- as.matrix(gext_data[, k4_cols, drop = FALSE])

    doMC::registerDoMC(cores = 50)
    results <- plyr::llply(names(filters), function(fname) {
        f <- filters[[fname]]$filter
        if (is.null(f)) f <- rep(TRUE, nrow(gext_data))

        message(sprintf("Calculating peak overlap %s...", filters[[fname]]$label))

        # Subset indices
        idx <- which(f)
        chrom_f <- chrom_vec[idx]
        is_test <- chrom_f %in% test_chroms
        is_train <- !is_test

        # --- K27 tracks ---
        obs_f <- obs_k27_vec[idx]
        q_obs <- quantile(obs_f, T_q, na.rm = TRUE)
        top_obs <- obs_f > q_obs

        n_obs_test <- sum(top_obs & is_test)
        n_obs_total <- sum(top_obs)
        on_test1_k27 <- if (n_obs_total > 0) n_obs_test / n_obs_total else NA

        pred_mat <- k27_mat_full[idx, , drop = FALSE]
        q_preds <- colQuantiles(pred_mat, probs = T_q, na.rm = TRUE)
        top_pred_mat <- t(t(pred_mat) > q_preds)  # Fast threshold comparison

        top_both_mat <- top_pred_mat & top_obs

        n_pred_test <- colSums2(top_pred_mat & is_test)
        n_pred_train <- colSums2(top_pred_mat & is_train)
        n_pred_total <- colSums2(top_pred_mat)
        n_both_test <- colSums2(top_both_mat & is_test)
        n_both_train <- colSums2(top_both_mat & is_train)

        k27_df <- data.frame(
            col_name = k27_cols,
            train_p = ifelse(n_pred_train > 0, n_both_train / n_pred_train, NA),
            test_p = ifelse(n_pred_test > 0, n_both_test / n_pred_test, NA),
            on_test1 = on_test1_k27,
            on_test2 = ifelse(n_pred_total > 0, n_pred_test / n_pred_total, NA)
        )

        # --- K4 tracks ---
        obs_f <- obs_k4_vec[idx]
        q_obs <- quantile(obs_f, T_q, na.rm = TRUE)
        top_obs <- obs_f > q_obs

        n_obs_test <- sum(top_obs & is_test)
        n_obs_total <- sum(top_obs)
        on_test1_k4 <- if (n_obs_total > 0) n_obs_test / n_obs_total else NA

        pred_mat <- k4_mat_full[idx, , drop = FALSE]
        q_preds <- colQuantiles(pred_mat, probs = T_q, na.rm = TRUE)
        top_pred_mat <- t(t(pred_mat) > q_preds)

        top_both_mat <- top_pred_mat & top_obs

        n_pred_test <- colSums2(top_pred_mat & is_test)
        n_pred_train <- colSums2(top_pred_mat & is_train)
        n_pred_total <- colSums2(top_pred_mat)
        n_both_test <- colSums2(top_both_mat & is_test)
        n_both_train <- colSums2(top_both_mat & is_train)

        k4_df <- data.frame(
            col_name = k4_cols,
            train_p = ifelse(n_pred_train > 0, n_both_train / n_pred_train, NA),
            test_p = ifelse(n_pred_test > 0, n_both_test / n_pred_test, NA),
            on_test1 = on_test1_k4,
            on_test2 = ifelse(n_pred_total > 0, n_pred_test / n_pred_total, NA)
        )

        list(k27_df = k27_df, k4_df = k4_df)
    }, .parallel = TRUE)
    names(results) <- names(filters)

    return(results)
}

# ==============================================================================
# 6. Plotting
# ==============================================================================

#' Define which model types to exclude from plots
#' To exclude a series: add its model_type to this vector
#' Available model_types: "borzoi", "from_scratch", "sns", "silicus", "silicus_ft", "markov", "random", "mm10minus"
DEFAULT_EXCLUDE_MODELS <- c()  # Empty = include all

#' Get the ordered list of model names for plotting
#'
#' Order: RF series -> SNS (all in mm10->mm10), then mm10->silicus,
#'        synth->synth, synth->mm10, minus->minus, minus->mm10
#'
#' @param track_tbl Track table
#' @param is_k27_filter TRUE for K27 tracks, FALSE for K4
#' @param exclude_models Character vector of model_types to exclude
get_name_order <- function(track_tbl, is_k27_filter, exclude_models = DEFAULT_EXCLUDE_MODELS) {
    # Filter out excluded models
    tbl <- track_tbl %>%
        filter(is_k27 == is_k27_filter, !model_type %in% exclude_models)

    # Debug output
    # message("Display region assignments:")
    # print(tbl %>% distinct(model_type, train_genome, infer_genome, display_region, region_order) %>% arrange(region_order))

    # Sort by numeric values for reliability
    tbl <- tbl %>%
        arrange(region_order, order_group, order_within_group, col_name) %>%
        pull(col_name)

    return(tbl)
}

#' Compute region boundaries for background shading
#' Groups silicus/markov/random into "synth" regions
compute_region_boundaries <- function(track_tbl, name_order, is_k27_filter, exclude_models = DEFAULT_EXCLUDE_MODELS) {
    # Region configuration: colors and labels for DISPLAY regions (not genome pairs)
    region_config <- tribble(
        ~display_region, ~color, ~label, ~text_color,
        "mm10_mm10", "#4DAF4A", "mm10 -> mm10", "#2E7D32",
        "mm10_silicus", "#984EA3", "mm10 -> silicus", "#6A1B9A",
        "synth_synth", "#FF7F00", "synth -> synth", "#E65100",
        "synth_mm10", "#E41A1C", "synth -> mm10", "#B71C1C",
        "minus_minus", "#377EB8", "minus -> minus", "#1565C0",
        "minus_mm10", "#A65628", "minus -> mm10", "#5D4037"
    )

    # Get position and genome info for each model
    ordered_info <- tibble(col_name = name_order, pos = seq_along(name_order)) %>%
        left_join(
            track_tbl %>%
                filter(is_k27 == is_k27_filter, !model_type %in% exclude_models) %>%
                select(col_name, model_type, display_region),
            by = "col_name"
        )

    # Calculate boundaries for each display region
    boundaries <- ordered_info %>%
        filter(!is.na(display_region)) %>%
        group_by(display_region) %>%
        summarise(
            start_pos = min(pos),
            end_pos = max(pos),
            .groups = "drop"
        ) %>%
        left_join(region_config, by = "display_region") %>%
        filter(!is.na(color)) %>%
        # Convert to boundary coordinates (between positions)
        mutate(
            xmin = start_pos - 0.5,
            xmax = end_pos + 0.5
        ) %>%
        arrange(start_pos)

    return(boundaries)
}

#' Get the RF series column names for line connection
get_rf_series_names <- function(track_tbl, is_k27_filter, exclude_models = DEFAULT_EXCLUDE_MODELS) {
    if ("borzoi" %in% exclude_models) return(character(0))
    track_tbl %>%
        filter(is_k27 == is_k27_filter, model_type == "borzoi") %>%
        pull(col_name)
}

#' Get mm10Minus series column names for line connection
get_minus_series_names <- function(track_tbl, is_k27_filter, target_infer_genome, exclude_models = DEFAULT_EXCLUDE_MODELS) {
    if ("mm10minus" %in% exclude_models) return(character(0))
    track_tbl %>%
        filter(is_k27 == is_k27_filter, model_type == "mm10minus",
               infer_genome == target_infer_genome) %>%
        pull(col_name)
}

#' Create R^2 comparison plot
#'
#' @param cor_te Named vector of test R^2 values
#' @param cor_tr Named vector of train R^2 values
#' @param track_tbl Track table
#' @param is_k27_filter TRUE for K27, FALSE for K4
#' @param title Plot title
#' @param y_max Maximum y-axis value
#' @param exclude_models Character vector of model_types to exclude from the plot
create_rsqr_plot <- function(cor_te, cor_tr, track_tbl, is_k27_filter, title, y_max = 0.88,
                              exclude_models = DEFAULT_EXCLUDE_MODELS) {
    model_cols <- names(cor_te)

    # Filter out excluded model types from results
    if (length(exclude_models) > 0) {
        excluded_cols <- track_tbl %>%
            filter(model_type %in% exclude_models) %>%
            pull(col_name)
        model_cols <- setdiff(model_cols, excluded_cols)
        cor_te <- cor_te[model_cols]
        cor_tr <- cor_tr[model_cols]
    }

    # Prepare data
    cor_df <- data.frame(
        name = model_cols,
        Test = cor_te,
        Train = cor_tr
    ) %>%
        left_join(
            track_tbl %>%
                filter(is_k27 == is_k27_filter, !model_type %in% exclude_models) %>%
                select(col_name, model_type, features, train_genome, infer_genome, label),
            by = c("name" = "col_name")
        ) %>%
        pivot_longer(cols = c(Test, Train), names_to = "Set", values_to = "R2")

    # Get ordering and boundaries
    name_order <- get_name_order(track_tbl, is_k27_filter, exclude_models)
    name_order <- name_order[name_order %in% model_cols]
    boundaries <- compute_region_boundaries(track_tbl, name_order, is_k27_filter, exclude_models)

    # Get series names for line connections
    rf_series <- get_rf_series_names(track_tbl, is_k27_filter, exclude_models)
    minus_series <- get_minus_series_names(track_tbl, is_k27_filter, "mm10Minus", exclude_models)
    minus_mm10_series <- get_minus_series_names(track_tbl, is_k27_filter, "mm10", exclude_models)

    # Create label map
    label_map <- track_tbl %>%
        filter(is_k27 == is_k27_filter, col_name %in% name_order) %>%
        distinct(col_name, label) %>%
        deframe()

    # Build plot
    p <- cor_df %>%
        mutate(
            name = factor(name, levels = name_order),
            r2_label = round(R2, 3)
        ) %>%
        filter(!is.na(name)) %>%
        ggplot(aes(x = name, y = R2, color = Set, group = Set))

    # Add background rectangles
    for (i in seq_len(nrow(boundaries))) {
        b <- boundaries[i, ]
        p <- p + annotate("rect",
            xmin = b$xmin, xmax = b$xmax,
            ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = b$color
        )
        # Add region label at top
        p <- p + annotate("text",
            x = (b$xmin + b$xmax) / 2, y = y_max * 0.94,
            label = b$label, size = 3.2,
            color = b$text_color, fontface = "bold"
        )
    }

    # Add Borzoi RF series label
    rf_series_in_plot <- rf_series[rf_series %in% name_order]
    if (length(rf_series_in_plot) > 0) {
        rf_positions <- which(name_order %in% rf_series_in_plot)
        if (length(rf_positions) > 0) {
            p <- p + annotate("text",
                x = mean(rf_positions), y = y_max * 0.90,
                label = "(Borzoi RF series)", size = 2.8,
                color = "gray50", fontface = "italic"
            )
        }
    }

    # Add vertical guide lines
    p <- p + geom_vline(
        xintercept = seq_along(name_order),
        linetype = "dashed", color = "gray80", alpha = 0.3
    )

    # Add line connections for series
    rf_series_in_plot <- rf_series[rf_series %in% name_order]
    if (length(rf_series_in_plot) > 0) {
        p <- p + geom_line(
            data = function(d) filter(d, name %in% rf_series_in_plot),
            linewidth = 1, alpha = 0.7
        )
    }
    minus_series_in_plot <- minus_series[minus_series %in% name_order]
    if (length(minus_series_in_plot) > 0) {
        p <- p + geom_line(
            data = function(d) filter(d, name %in% minus_series_in_plot),
            linewidth = 1, alpha = 0.7
        )
    }
    minus_mm10_in_plot <- minus_mm10_series[minus_mm10_series %in% name_order]
    if (length(minus_mm10_in_plot) > 0) {
        p <- p + geom_line(
            data = function(d) filter(d, name %in% minus_mm10_in_plot),
            linewidth = 1, alpha = 0.7
        )
    }

    # Add points and labels
    p <- p +
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
        scale_x_discrete(limits = name_order, labels = label_map) +
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
            subtitle = "Train/Infer genome shown above: mm10 = real mouse, synth = synthetic genome"
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

#' Create peak overlap comparison plot
#'
#' @param overlap_df Data frame with peak overlap results (from calc_all_peak_overlap)
#' @param track_tbl Track table
#' @param is_k27_filter TRUE for K27, FALSE for K4
#' @param title Plot title
#' @param y_max Maximum y-axis value
#' @param exclude_models Character vector of model_types to exclude from the plot
create_peak_overlap_plot <- function(overlap_df, track_tbl, is_k27_filter, title, y_max = 1.0,
                                      exclude_models = DEFAULT_EXCLUDE_MODELS) {
    # Filter out excluded model types
    if (length(exclude_models) > 0) {
        excluded_cols <- track_tbl %>%
            filter(model_type %in% exclude_models) %>%
            pull(col_name)
        overlap_df <- overlap_df %>% filter(!col_name %in% excluded_cols)
    }

    # Join with track metadata
    plot_df <- overlap_df %>%
        left_join(
            track_tbl %>%
                filter(is_k27 == is_k27_filter, !model_type %in% exclude_models) %>%
                select(col_name, model_type, features, train_genome, infer_genome, label),
            by = "col_name"
        ) %>%
        filter(!is.na(model_type)) %>%
        pivot_longer(
            cols = c(train_p, test_p),
            names_to = "Set",
            values_to = "Overlap"
        ) %>%
        mutate(Set = case_when(
            Set == "train_p" ~ "Train",
            Set == "test_p" ~ "Test"
        ))

    # Get ordering and boundaries
    name_order <- get_name_order(track_tbl, is_k27_filter, exclude_models)
    name_order <- name_order[name_order %in% overlap_df$col_name]
    boundaries <- compute_region_boundaries(track_tbl, name_order, is_k27_filter, exclude_models)

    # Get series names for line connections
    rf_series <- get_rf_series_names(track_tbl, is_k27_filter, exclude_models)
    minus_series <- get_minus_series_names(track_tbl, is_k27_filter, "mm10Minus", exclude_models)
    minus_mm10_series <- get_minus_series_names(track_tbl, is_k27_filter, "mm10", exclude_models)

    # Create label map
    label_map <- track_tbl %>%
        filter(is_k27 == is_k27_filter, col_name %in% name_order) %>%
        distinct(col_name, label) %>%
        deframe()

    # Build plot
    p <- plot_df %>%
        mutate(
            col_name = factor(col_name, levels = name_order),
            overlap_label = round(Overlap, 3)
        ) %>%
        filter(!is.na(col_name), !is.na(Overlap)) %>%
        ggplot(aes(x = col_name, y = Overlap, color = Set, group = Set))

    # Add background rectangles
    for (i in seq_len(nrow(boundaries))) {
        b <- boundaries[i, ]
        p <- p + annotate("rect",
            xmin = b$xmin, xmax = b$xmax,
            ymin = -Inf, ymax = Inf,
            alpha = 0.08, fill = b$color
        )
        # Add region label at top
        p <- p + annotate("text",
            x = (b$xmin + b$xmax) / 2, y = y_max * 0.94,
            label = b$label, size = 3.2,
            color = b$text_color, fontface = "bold"
        )
    }

    # Add Borzoi RF series label
    rf_series_in_plot <- rf_series[rf_series %in% name_order]
    if (length(rf_series_in_plot) > 0) {
        rf_positions <- which(name_order %in% rf_series_in_plot)
        if (length(rf_positions) > 0) {
            p <- p + annotate("text",
                x = mean(rf_positions), y = y_max * 0.90,
                label = "(Borzoi RF series)", size = 2.8,
                color = "gray50", fontface = "italic"
            )
        }
    }

    # Add vertical guide lines
    p <- p + geom_vline(
        xintercept = seq_along(name_order),
        linetype = "dashed", color = "gray80", alpha = 0.3
    )

    # Add line connections for series
    rf_series_in_plot <- rf_series[rf_series %in% name_order]
    if (length(rf_series_in_plot) > 0) {
        p <- p + geom_line(
            data = function(d) filter(d, col_name %in% rf_series_in_plot),
            linewidth = 1, alpha = 0.7
        )
    }
    minus_series_in_plot <- minus_series[minus_series %in% name_order]
    if (length(minus_series_in_plot) > 0) {
        p <- p + geom_line(
            data = function(d) filter(d, col_name %in% minus_series_in_plot),
            linewidth = 1, alpha = 0.7
        )
    }
    minus_mm10_in_plot <- minus_mm10_series[minus_mm10_series %in% name_order]
    if (length(minus_mm10_in_plot) > 0) {
        p <- p + geom_line(
            data = function(d) filter(d, col_name %in% minus_mm10_in_plot),
            linewidth = 1, alpha = 0.7
        )
    }

    # Add points and labels
    p <- p +
        geom_point(size = 3.5) +
        geom_text_repel(
            aes(label = overlap_label),
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
        scale_x_discrete(limits = name_order, labels = label_map) +
        scale_y_continuous(
            limits = c(0, y_max),
            breaks = seq(0, 1, 0.2),
            expand = c(0.02, 0)
        ) +
        labs(
            x = "Model Configuration",
            y = "Peak Overlap (proportion)",
            color = "Dataset",
            title = title,
            subtitle = "Overlap between top 2% predicted peaks and observed peaks"
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
# 7. Output Generation
# ==============================================================================

# Plot dimensions
PLOT_WIDTH <- 20   # Wider plot
PLOT_HEIGHT <- 10

#' Generate and save all plots
#'
#' @param results R^2 calculation results
#' @param track_tbl Track table
#' @param filters Filter definitions
#' @param output_dir Output directory
#' @param exclude_models Character vector of model_types to exclude from ALL plots
#' @param peak_overlap_results Peak overlap calculation results (optional)
generate_and_save_all <- function(results, track_tbl, filters, output_dir,
                                   exclude_models = DEFAULT_EXCLUDE_MODELS,
                                   peak_overlap_results = NULL) {
    plots <- list()

    # Generate R^2 plots for all filter/mark combinations
    for (fname in names(filters)) {
        for (mark in c("k27", "k4")) {
            is_k27 <- (mark == "k27")
            mark_label <- if (is_k27) "H3K27me3" else "H3K4me3"

            title <- paste0("Whole Genome ", mark_label, " Prediction Performance ",
                           filters[[fname]]$label)

            message(sprintf("Creating R^2 plot: %s %s", mark, fname))

            p <- create_rsqr_plot(
                cor_te = results[[fname]][[mark]]$test,
                cor_tr = results[[fname]][[mark]]$train,
                track_tbl = track_tbl,
                is_k27_filter = is_k27,
                title = title,
                exclude_models = exclude_models
            )

            # Save plot
            base_name <- paste0("gw_rsqr_", mark, filters[[fname]]$suffix, "_plot")
            ggsave(file.path(output_dir, paste0(base_name, ".png")), p,
                   width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)

            plots[[paste("rsqr", mark, fname, sep = "_")]] <- p
        }
    }

    # Save legacy filename for backwards compatibility
    if ("all" %in% names(filters)) {
        ggsave(file.path(output_dir, "gw_rsqr_plot.png"), plots[["rsqr_k27_all"]],
               width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)
    }

    # Generate peak overlap plots if results are provided (with filters like R^2)
    if (!is.null(peak_overlap_results)) {
        for (fname in names(filters)) {
            for (mark in c("k27", "k4")) {
                is_k27 <- (mark == "k27")
                mark_label <- if (is_k27) "H3K27me3" else "H3K4me3"
                overlap_df <- if (is_k27) peak_overlap_results[[fname]]$k27_df else peak_overlap_results[[fname]]$k4_df

                title <- paste0("Whole Genome ", mark_label, " Peak Overlap (Top 2%) ",
                               filters[[fname]]$label)

                message(sprintf("Creating peak overlap plot: %s %s", mark, fname))

                p <- create_peak_overlap_plot(
                    overlap_df = overlap_df,
                    track_tbl = track_tbl,
                    is_k27_filter = is_k27,
                    title = title,
                    exclude_models = exclude_models
                )

                # Save plot
                base_name <- paste0("gw_peak_overlap_", mark, filters[[fname]]$suffix, "_plot")
                ggsave(file.path(output_dir, paste0(base_name, ".png")), p,
                       width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)

                plots[[paste("peak_overlap", mark, fname, sep = "_")]] <- p
            }
        }
    }

    return(plots)
}

#' Print summary statistics
print_summary <- function(results, filters) {
    for (fname in names(filters)) {
        for (mark in c("k27", "k4")) {
            mark_label <- if (mark == "k27") "H3K27me3" else "H3K4me3"
            filter_label <- filters[[fname]]$label
            if (filter_label == "") filter_label <- "(All Regions)"

            message(sprintf("\n=== %s R^2 Summary %s ===", mark_label, filter_label))
            message("Test set (held out chromosomes and regions):")
            print(round(results[[fname]][[mark]]$test, 3))
            message("\nTrain set:")
            print(round(results[[fname]][[mark]]$train, 3))
        }
    }
}

# ==============================================================================
# 8. Main
# ==============================================================================

message("Starting pipeline...")

# Setup
link_dir <- setup_misha()

# Build track table
message("Building track table...")
track_tbl <- build_track_table()
track_names <- track_tbl$track_name
col_names <- track_tbl$col_name

# Validate tracks exist
stopifnot(all(gtrack.exists(track_names)))
stopifnot(length(track_names) == length(unique(col_names)))
message(sprintf("Track table built: %d tracks", nrow(track_tbl)))

# Load data
message("Loading data...")
data <- load_data(link_dir)
cg_trace_f <- data$cg_trace_f

# Setup train/test split
split_data <- setup_train_test_split(cg_trace_f)
cg_trace_f <- split_data$cg_trace_f
train_mod_gw <- split_data$train_mod_gw

# Extract prediction tracks
message("Extracting prediction tracks...")
gext2 <- gextract(track_names,
    iterator = cg_trace_f, intervals = cg_trace_f,
    colnames = col_names
) %>% arrange(intervalID) %fcache_rds% here("data/gext2.rds")
gext2[is.na(gext2)] <- 0

# Validate alignment
validate_alignment(gext2, cg_trace_f)

# Extract K4 response
message("Extracting K4 response...")
gvtrack.create("obs_k4", "jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3", "sum")
k4_resp <- gextract("obs_k4",
    intervals = cg_trace_f, iterator = cg_trace_f,
    colnames = "k4_1k"
)
k4_resp$lk4_1k <- log2(1 + k4_resp$k4_1k)

# Create filters
filters <- create_filters(cg_trace_f)

# Calculate R^2 for all combinations
message("Calculating R^2 values...")
results <- calc_all_rsqr(gext2, cg_trace_f, k4_resp, track_tbl, train_mod_gw, filters)

# Calculate peak overlap 
message("Calculating peak overlap...")
peak_overlap_results <- calc_all_peak_overlap(
    gext_data = gext2,
    obs_k27_vec = cg_trace_f$lk27_1k,
    obs_k4_vec = k4_resp$lk4_1k,
    track_tbl = track_tbl,
    filters = filters
)

# Generate and save plots
output_dir <- file.path(link_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Generating plots...")
plots <- generate_and_save_all(results, track_tbl, filters, output_dir,
                                peak_overlap_results = peak_overlap_results)

# Print summary
print_summary(results, filters)

# Print peak overlap summary
for (fname in names(filters)) {
    filter_label <- filters[[fname]]$label
    if (filter_label == "") filter_label <- "(All Regions)"
    message(sprintf("\n=== Peak Overlap Summary %s ===", filter_label))
    message("\nK27 Peak Overlap (Test set):")
    print(peak_overlap_results[[fname]]$k27_df %>% select(col_name, test_p) %>% arrange(desc(test_p)))
    message("\nK4 Peak Overlap (Test set):")
    print(peak_overlap_results[[fname]]$k4_df %>% select(col_name, test_p) %>% arrange(desc(test_p)))
}

message("\nDone!")

return(invisible(list(
    track_tbl = track_tbl,
    results = results,
    peak_overlap_results = peak_overlap_results,
    plots = plots
)))
