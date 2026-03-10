# ==============================================================================
# Borzoi Analysis Utilities
# ==============================================================================
#
# This file contains reusable utility functions for Borzoi model analysis,
# including peak overlap metrics, track configuration helpers, and data
# manipulation functions.
#
# Functions:
#   - compute_track_peak_overlap: Single track peak overlap metrics
#   - compute_all_track_peak_overlaps: Parallel computation for multiple tracks
#   - compute_peak_overlap_by_context: Peak overlap across multiple context sizes
#   - compute_track_auprc: Single track AUPRC metrics
#   - compute_all_track_auprcs: Parallel AUPRC computation for multiple tracks
#   - slugify_label: Convert labels to filesystem-safe names
#   - load_genome_variant_index: Load genome variant metadata
#   - add_nomenclature_columns: Add simplified nomenclature to track tables
#
# ==============================================================================

# ==============================================================================
# Peak Overlap Metrics (tst_symdif style)
# ==============================================================================

#' Compute train/test/val peak overlap metrics for a single track against observed data
#'
#' This function is similar to tst_symdif from borz_analysis.r but generalized to work
#' with a single track and optional genomic context.
#'
#' @param track_name Track name to compare against observed data
#' @param obs_track_name Observed data track name (default: K27me3 normalization track)
#' @param context_size Context size for gvtrack.iterator (default: 200)
#' @param T_q Quantile threshold for peak calling (default: 0.98)
#' @param test_chroms Test chromosome names (default: c("chr4", "chr10", "chr14", "chr15"))
#' @param val_chroms Validation chromosome names (default: c("chr8", "chr9"))
#' @param iterator_size Iterator size for gscreen (default: 100)
#'
#' @return Named vector with train_p, test_p, val_p and proportions
compute_track_peak_overlap <- function(track_name,
                                        obs_track_name = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm",
                                        context_size = 200,
                                        T_q = 0.98,
                                        test_chroms = c("chr4", "chr10", "chr14", "chr15"),
                                        val_chroms = c("chr8", "chr9"),
                                        exclude_chroms = c("chrX"),
                                        iterator_size = 100,
                                        intervals = gintervals.all()) {

    # Create virtual tracks with optional context
    gvtrack.create("obs_raw_tmp", obs_track_name, "avg")
    gvtrack.create("pred_raw_tmp", track_name, "avg")

    if (context_size > 0) {
        gvtrack.iterator("obs_raw_tmp", sshift = -context_size, eshift = context_size)
        gvtrack.iterator("pred_raw_tmp", sshift = -context_size, eshift = context_size)
    }

    # Build expressions that treat NA/NaN as 0 (using is.na which catches both NA and NaN)
    obs_expr <- "ifelse(is.na(obs_raw_tmp), 0, obs_raw_tmp)"
    pred_expr <- "ifelse(is.na(pred_raw_tmp), 0, pred_raw_tmp)"

    # Get intervals excluding chrX (and any other excluded chromosomes)
    valid_intervs <- intervals[!intervals$chrom %in% exclude_chroms, ]

    # Compute quantile thresholds genome-wide (excluding chrX)
    q_obs <- gquantiles(obs_expr, T_q, intervals = valid_intervs, iterator = iterator_size)
    q_pred <- gquantiles(pred_expr, T_q, intervals = valid_intervs, iterator = iterator_size)

    # Screen for peaks above threshold (excluding chrX)
    d_obs <- gscreen(paste0(obs_expr, " > ", q_obs), iterator = iterator_size, intervals = valid_intervs)
    d_pred <- gscreen(paste0(pred_expr, " > ", q_pred), iterator = iterator_size, intervals = valid_intervs)

    # Compute overlap
    d_both <- gintervals.intersect(d_obs, d_pred)

    # Clean up temporary tracks
    gvtrack.rm("obs_raw_tmp")
    gvtrack.rm("pred_raw_tmp")

    # Helper to compute overlap proportion for a chromosome set
    calc_overlap_p <- function(chroms) {
        f_pred <- d_pred$chrom %in% chroms
        f_both <- d_both$chrom %in% chroms
        pred_bp <- sum(d_pred$end[f_pred] - d_pred$start[f_pred])
        both_bp <- sum(d_both$end[f_both] - d_both$start[f_both])
        if (pred_bp > 0) both_bp / pred_bp else NA
    }

    # Define train chromosomes (everything except test, val, and excluded)
    all_chroms <- unique(d_pred$chrom)
    train_chroms <- setdiff(all_chroms, c(test_chroms, val_chroms, exclude_chroms))

    # Calculate metrics
    train_p <- calc_overlap_p(train_chroms)
    test_p <- calc_overlap_p(test_chroms)
    val_p <- calc_overlap_p(val_chroms)

    return(c(
        train_p = train_p,
        test_p = test_p,
        val_p = val_p
    ))
}

#' Compute peak overlap metrics for multiple tracks in parallel
#'
#' @param track_tbl Track table with track_name column
#' @param obs_track_name Observed data track name
#' @param context_size Context size for gvtrack.iterator (default: 200)
#' @param T_q Quantile threshold (default: 0.98)
#' @param test_chroms Test chromosome names
#' @param val_chroms Validation chromosome names
#' @param exclude_chroms Chromosomes to exclude from analysis (default: c("chrX"))
#' @param iterator_size Iterator size for gscreen (default: 100)
#' @param n_cores Number of cores for parallel processing (default: 50)
#'
#' @return Data frame with columns from track_tbl plus train_p, test_p, val_p
compute_all_track_peak_overlaps <- function(track_tbl,
                                             obs_track_name = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm",
                                             context_size = 200,
                                             T_q = 0.98,
                                             test_chroms = c("chr4", "chr10", "chr14", "chr15"),
                                             val_chroms = c("chr8", "chr9"),
                                             exclude_chroms = c("chrX"),
                                             iterator_size = 100,
                                             n_cores = 50) {

    message(sprintf("Computing peak overlap for %d tracks with context=%d, T_q=%g (excluding %s)...",
                   nrow(track_tbl), context_size, T_q, paste(exclude_chroms, collapse = ", ")))

    # Register parallel backend
    doMC::registerDoMC(cores = n_cores)

    # Process each track in parallel
    results <- plyr::adply(track_tbl, 1, function(row) {
        tryCatch({
            metrics <- compute_track_peak_overlap(
                track_name = row$track_name,
                obs_track_name = obs_track_name,
                context_size = context_size,
                T_q = T_q,
                test_chroms = test_chroms,
                val_chroms = val_chroms,
                exclude_chroms = exclude_chroms,
                iterator_size = iterator_size
            )

            data.frame(
                train_p = metrics["train_p"],
                test_p = metrics["test_p"],
                val_p = metrics["val_p"]
            )
        }, error = function(e) {
            message(sprintf("  ERROR processing %s: %s", row$col_name, e$message))
            data.frame(
                train_p = NA,
                test_p = NA,
                val_p = NA
            )
        })
    }, .parallel = TRUE, .id = NULL)

    return(results)
}

#' Compute peak overlap across multiple context sizes
#'
#' @param track_tbl Track table
#' @param context_sizes Vector of context sizes to test
#' @param obs_track_name Observed data track name
#' @param ... Additional parameters passed to compute_all_track_peak_overlaps
#'
#' @return List of data frames, one for each context size
compute_peak_overlap_by_context <- function(track_tbl,
                                             context_sizes = c(0, 500, 2000, 8000),
                                             obs_track_name = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm",
                                             ...) {

    results <- list()

    for (sc in context_sizes) {
        message(sprintf("\n=== Processing context size: %d ===", sc))
        results[[as.character(sc)]] <- compute_all_track_peak_overlaps(
            track_tbl = track_tbl,
            obs_track_name = obs_track_name,
            context_size = sc,
            ...
        )
    }

    return(results)
}

# ==============================================================================
# AUPRC Metrics
# ==============================================================================

#' Compute AUPRC for a single track against observed data
#'
#' AUPRC (Area Under Precision-Recall Curve) where:
#' - Positives are defined as bins where observed signal > T_q quantile (default: top 2%)
#' - Scores are the predicted signal values
#'
#' @param track_name Track name for predictions
#' @param obs_track_name Observed data track name
#' @param intervals Genomic intervals to evaluate on (required)
#' @param T_q Quantile threshold for defining positives (default: 0.98 = top 2%)
#' @param context_size Context size for gvtrack.iterator (default: 200)
#' @param test_chroms Test chromosome names
#' @param val_chroms Validation chromosome names
#' @param iterator_size Iterator size for gextract
#'
#' @return Named vector with auprc_train, auprc_test, auprc_val
compute_track_auprc <- function(track_name,
                                 obs_track_name = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm",
                                 intervals,
                                 T_q = 0.98,
                                 context_size = 200,
                                 context_size_obs = context_size,
                                 context_size_pred = context_size,
                                 test_chroms = c("chr4", "chr10", "chr14", "chr15"),
                                 val_chroms = c("chr8", "chr9"),
                                 exclude_chroms = c("chrX"),
                                 iterator_size = 100) {

    # Filter intervals to exclude chrX (and any other excluded chromosomes)
    intervals <- intervals[!intervals$chrom %in% exclude_chroms, ]

    # Create virtual tracks with context, treating NAs as zeros
    gvtrack.create("obs_auprc_raw", obs_track_name, "avg")
    gvtrack.create("pred_auprc_raw", track_name, "avg")

    if (context_size_obs > 0 || context_size_pred > 0) {
        gvtrack.iterator("obs_auprc_raw", sshift = -context_size_obs, eshift = context_size_obs)
        gvtrack.iterator("pred_auprc_raw", sshift = -context_size_pred, eshift = context_size_pred)
    }

    # Extract values at intervals using virtual tracks directly
    d <- gextract(c("obs_auprc_raw", "pred_auprc_raw"), intervals = intervals, iterator = iterator_size)

    # Clean up temporary tracks
    gvtrack.rm("obs_auprc_raw")
    gvtrack.rm("pred_auprc_raw")

    # Replace NaN with 0
    d$obs_auprc_raw[is.nan(d$obs_auprc_raw) | is.na(d$obs_auprc_raw)] <- 0
    d$pred_auprc_raw[is.nan(d$pred_auprc_raw) | is.na(d$pred_auprc_raw)] <- 0

    # Rename for clarity
    d$obs_val <- d$obs_auprc_raw
    d$pred_val <- d$pred_auprc_raw

    # Define chromosome splits (exclude_chroms already filtered out)
    all_chroms <- unique(d$chrom)
    train_chroms <- setdiff(all_chroms, c(test_chroms, val_chroms))

    # Compute global threshold for positives (top 2% of observed)
    q_obs <- quantile(d$obs_val, T_q, na.rm = TRUE)

    # Helper to compute AUPRC for a chromosome set
    calc_auprc <- function(chroms) {
        f <- d$chrom %in% chroms
        obs_vals <- d$obs_val[f]
        pred_vals <- d$pred_val[f]

        # Define labels: 1 if observed > threshold, 0 otherwise
        labels <- as.integer(obs_vals > q_obs)

        # Need both classes to compute AUPRC
        if (sum(labels == 1) == 0 || sum(labels == 0) == 0) {
            return(NA)
        }

        # Compute AUPRC using PRROC
        pr <- PRROC::pr.curve(
            scores.class0 = pred_vals[labels == 1],
            scores.class1 = pred_vals[labels == 0]
        )

        return(pr$auc.integral)
    }

    # Calculate AUPRC for each split
    auprc_train <- calc_auprc(train_chroms)
    auprc_test <- calc_auprc(test_chroms)
    auprc_val <- calc_auprc(val_chroms)

    return(c(
        auprc_train = auprc_train,
        auprc_test = auprc_test,
        auprc_val = auprc_val
    ))
}

#' Compute AUPRC metrics for multiple tracks in parallel
#'
#' @param track_tbl Track table with track_name column
#' @param obs_track_name Observed data track name
#' @param intervals Genomic intervals to evaluate on (required)
#' @param T_q Quantile threshold for defining positives (default: 0.98)
#' @param context_size Context size for gvtrack.iterator (default: 200)
#' @param test_chroms Test chromosome names
#' @param val_chroms Validation chromosome names
#' @param exclude_chroms Chromosomes to exclude from analysis (default: c("chrX"))
#' @param iterator_size Iterator size for gextract
#' @param n_cores Number of cores for parallel processing (default: 50)
#'
#' @return Data frame with columns from track_tbl plus auprc_train, auprc_test, auprc_val
compute_all_track_auprcs <- function(track_tbl,
                                      obs_track_name = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm",
                                      intervals,
                                      T_q = 0.98,
                                      context_size = 200,
                                      context_size_obs = context_size,
                                      context_size_pred = context_size,
                                      test_chroms = c("chr4", "chr10", "chr14", "chr15"),
                                      val_chroms = c("chr8", "chr9"),
                                      exclude_chroms = c("chrX"),
                                      iterator_size = 100,
                                      n_cores = 50) {

    message(sprintf("Computing AUPRC for %d tracks with context=%d, T_q=%g (excluding %s)...",
                   nrow(track_tbl), context_size, T_q, paste(exclude_chroms, collapse = ", ")))

    # Register parallel backend
    doMC::registerDoMC(cores = n_cores)

    if (!("obs_track_name" %in% colnames(track_tbl))) {
        track_tbl$obs_track_name <- obs_track_name
    }


    # Process each track in parallel
    results <- plyr::adply(track_tbl, 1, function(row) {
        tryCatch({
            metrics <- compute_track_auprc(
                track_name = row$track_name,
                obs_track_name = row$obs_track_name,
                intervals = intervals,
                T_q = T_q,
                context_size = context_size,
                context_size_obs = context_size_obs,
                context_size_pred = context_size_pred,
                test_chroms = test_chroms,
                val_chroms = val_chroms,
                exclude_chroms = exclude_chroms,
                iterator_size = iterator_size
            )

            data.frame(
                auprc_train = metrics["auprc_train"],
                auprc_test = metrics["auprc_test"],
                auprc_val = metrics["auprc_val"]
            )
        }, error = function(e) {
            message(sprintf("  ERROR processing %s: %s", row$col_name, e$message))
            data.frame(
                auprc_train = NA,
                auprc_test = NA,
                auprc_val = NA
            )
        })
    }, .parallel = TRUE, .id = NULL)

    return(results)
}

# ==============================================================================
# Track Configuration Utilities
# ==============================================================================

#' Slugify a label for filesystem-safe names
#'
#' @param x Character vector to slugify
#' @return Slugified character vector
slugify_label <- function(x) {
    x <- tolower(x)
    x <- gsub("%", "pct", x, fixed = TRUE)
    x <- gsub("[^a-z0-9]+", "-", x)
    x <- gsub("(^-|-$)", "", x)
    x
}

#' Load genome variant index mapping
#'
#' @param path Path to genome variant index file
#' @return Tibble with variant metadata
load_genome_variant_index <- function(path = here::here("analysis/genome_variant_index.tsv")) {
    if (!file.exists(path)) {
        message("Variant index not found at: ", path)
        return(tibble::tibble(
            variant_key = character(),
            genome_code = character(),
            description = character()
        ))
    }
    readr::read_tsv(path, show_col_types = FALSE)
}

#' Add simplified nomenclature columns to the track table
#'
#' @param track_tbl Track table
#' @param variant_index Genome variant index
#' @return Track table with nomenclature columns added
add_nomenclature_columns <- function(track_tbl, variant_index) {
    idx <- variant_index %>%
        dplyr::distinct(variant_key, genome_code, description)

    track_tbl <- track_tbl %>%
        dplyr::left_join(
            idx %>% dplyr::rename(
                train_variant_key = variant_key,
                train_code = genome_code,
                train_desc = description
            ),
            by = "train_variant_key"
        ) %>%
        dplyr::left_join(
            idx %>% dplyr::rename(
                infer_variant_key = variant_key,
                infer_code = genome_code,
                infer_desc = description
            ),
            by = "infer_variant_key"
        )

    track_tbl <- track_tbl %>%
        dplyr::mutate(
            train_code = dplyr::if_else(is.na(train_code), train_variant_key, train_code),
            infer_code = dplyr::if_else(is.na(infer_code), infer_variant_key, infer_code),
            features_slug = slugify_label(features),
            new_col_name = paste(model_type, train_code, infer_code, features_slug, sep = "__"),
            new_track_name = paste0("seq.IQ.pcg.borzoi.nomenclature_v2.", new_col_name)
        )

    return(track_tbl)
}

# ==============================================================================
# Data Validation Utilities
# ==============================================================================

#' Validate that extracted data aligns with intervals
#'
#' @param extracted_data Extracted genomic data with chrom, start, end columns
#' @param intervals Interval data frame
#' @param stop_on_error Whether to stop if validation fails (default TRUE)
#' @return TRUE if validation passes, FALSE otherwise
validate_alignment <- function(extracted_data, intervals, stop_on_error = TRUE) {
    checks <- list(
        chrom = all(extracted_data$chrom == intervals$chrom),
        start = all(extracted_data$start == intervals$start),
        end = all(extracted_data$end == intervals$end)
    )

    all_pass <- all(unlist(checks))

    if (!all_pass) {
        msg <- sprintf(
            "Alignment validation failed: chrom=%s, start=%s, end=%s",
            checks$chrom, checks$start, checks$end
        )
        if (stop_on_error) {
            stop(msg)
        } else {
            warning(msg)
        }
    }

    return(all_pass)
}

#' Load and filter cg_trace data
#'
#' This function loads only the cg_trace data (not the heavy feats, feats35, feats_iqdn),
#' and caches the filtered result for faster subsequent runs.
load_cg_trace_filtered <- function() {
    cache_file <- here("data/cg_trace_filtered.rds")

    if (file.exists(cache_file)) {
        message("Loading cached filtered cg_trace...")
        return(readRDS(cache_file))
    }

    message("Loading and filtering cg_trace (will cache result)...")
    cg_trace <- readRDS(here("data/cg_trace_mm10.rds"))

    cg_trace <- gextract.left_join("mapab.umap_k100", intervals = cg_trace, iterator = cg_trace)

    # Filter data
    f_norp <- cg_trace$d_ltr != 0 &
              cg_trace$d_line != 0 &            
              cg_trace$start > 3e6 &
              !is.na(cg_trace$mapab.umap_k100) &
              cg_trace$mapab.umap_k100 > 0.9 
    
    cg_trace_f <- cg_trace[f_norp, ]

    cg_trace_f <- cg_trace_f %>%
        select(chrom, start, end) %>%
        gintervals.neighbors("mapab.Encode_blacklist_v2") %>% 
        dplyr::filter(dist != 0) %>%
        dplyr::select(chrom, start, end)

    # Cache the result
    message(sprintf("Caching filtered cg_trace (%d rows) to %s", nrow(cg_trace_f), cache_file))
    saveRDS(cg_trace_f, cache_file)

    return(cg_trace_f)
}

# ==============================================================================
# Checkpoint Evaluation and Selection Utilities
# ==============================================================================

#' Get all metric checkpoints for a given genome directory
#'
#' Scans the metrics subdirectory structure to find all available checkpoint tracks
#' for a given genome variant.
#'
#' @param genome_dir Name of the genome directory (e.g., "silicusPlusCGD")
#' @param track_base_dir Base directory for flashzoi tracks
#' @param track_pattern Track name pattern to match (e.g., "rf524k_EB4_cnt")
#' @return Data frame with columns: metric_type, track_path, track_name, is_k4
get_metric_checkpoints <- function(genome_dir,
                                    track_base_dir = "/home/aviezerl/mm10/tracks/seq/IQ/pcg/flashzoi",
                                    track_pattern = NULL) {
    metrics_dir <- file.path(track_base_dir, genome_dir, "metrics")

    if (!dir.exists(metrics_dir)) {
        warning(sprintf("Metrics directory not found: %s", metrics_dir))
        return(tibble::tibble(
            metric_type = character(),
            track_path = character(),
            track_name = character(),
            is_k4 = logical()
        ))
    }

    # Find all metric subdirectories (iou, loss, pearson, etc.)
    metric_types <- list.dirs(metrics_dir, recursive = FALSE, full.names = FALSE)
    # Exclude 'legacy' directory
    metric_types <- setdiff(metric_types, "legacy")

    results <- list()

    for (metric in metric_types) {
        metric_path <- file.path(metrics_dir, metric)
        # Recursively find all .track directories under this metric type
        all_track_dirs <- list.dirs(metric_path, recursive = TRUE, full.names = TRUE)
        track_dirs <- all_track_dirs[grepl("\\.track$", all_track_dirs)]

        # Filter by pattern if provided
        if (!is.null(track_pattern)) {
            track_dirs <- track_dirs[grepl(track_pattern, basename(track_dirs))]
        }

        for (track_dir_full in track_dirs) {
            # Get the relative path from metric_path to track_dir
            rel_path <- gsub(paste0("^", gsub("([.()])", "\\\\\\1", metric_path), "/"), "", track_dir_full)
            track_dir_name <- basename(track_dir_full)
            track_dir_base <- gsub("\\.track$", "", track_dir_name)
            
            # Extract subdirectory path (if any) between metric and .track directory
            subdir_path <- dirname(rel_path)
            if (subdir_path == ".") {
                # No subdirectory, track is directly under metric type
                metric_suffix <- metric
            } else {
                # Include subdirectory in metric suffix (e.g., "iou.more_samples")
                metric_suffix <- paste(metric, subdir_path, sep = ".")
            }
            
            # Construct the misha track name
            track_name <- paste0("seq.IQ.pcg.flashzoi.", genome_dir, ".metrics.", metric_suffix, ".",
                                track_dir_base)

            results[[length(results) + 1]] <- tibble::tibble(
                metric_type = metric,
                track_path = track_dir_full,
                track_name = track_name,
                track_basename = track_dir_base,
                is_k4 = grepl("_k4", track_dir_base)
            )
        }
    }

    if (length(results) == 0) {
        return(tibble::tibble(
            metric_type = character(),
            track_path = character(),
            track_name = character(),
            track_basename = character(),
            is_k4 = logical()
        ))
    }

    dplyr::bind_rows(results)
}

#' Check for missing metric checkpoints across genomes
#'
#' Identifies which metric types (iou, loss, pearson) are missing for each
#' genome and track combination.
#'
#' @param genome_dirs Vector of genome directory names to check
#' @param track_pattern Track name pattern to check (e.g., "rf524k_EB4_cnt")
#' @param expected_metrics Vector of expected metric types (default: iou, loss, pearson)
#' @param track_base_dir Base directory for flashzoi tracks
#' @return List with: missing_summary (data frame), coverage_matrix (wide format)
check_missing_metrics <- function(genome_dirs,
                                   track_pattern = "rf524k_EB4_cnt",
                                   expected_metrics = c("iou", "loss", "pearson"),
                                   track_base_dir = "/home/aviezerl/mm10/tracks/seq/IQ/pcg/flashzoi") {

    all_checkpoints <- list()

    for (genome_dir in genome_dirs) {
        cp <- get_metric_checkpoints(genome_dir, track_base_dir, track_pattern)
        if (nrow(cp) > 0) {
            cp$genome_dir <- genome_dir
            all_checkpoints[[genome_dir]] <- cp
        }
    }

    if (length(all_checkpoints) == 0) {
        return(list(
            missing_summary = tibble::tibble(
                genome_dir = character(),
                track_basename = character(),
                missing_metrics = character()
            ),
            coverage_matrix = NULL
        ))
    }

    checkpoints_df <- dplyr::bind_rows(all_checkpoints)

    # Get all unique track basenames
    all_basenames <- unique(checkpoints_df$track_basename)

    # Build coverage matrix
    coverage_list <- list()
    missing_list <- list()

    for (genome in genome_dirs) {
        genome_cp <- checkpoints_df %>% dplyr::filter(genome_dir == genome)

        for (basename in all_basenames) {
            track_cp <- genome_cp %>% dplyr::filter(track_basename == basename)
            present_metrics <- unique(track_cp$metric_type)
            missing_metrics <- setdiff(expected_metrics, present_metrics)

            # Coverage entry
            for (metric in expected_metrics) {
                coverage_list[[length(coverage_list) + 1]] <- tibble::tibble(
                    genome_dir = genome,
                    track_basename = basename,
                    metric_type = metric,
                    present = metric %in% present_metrics
                )
            }

            # Missing entry (only if something is missing)
            if (length(missing_metrics) > 0 && length(present_metrics) > 0) {
                missing_list[[length(missing_list) + 1]] <- tibble::tibble(
                    genome_dir = genome,
                    track_basename = basename,
                    missing_metrics = paste(missing_metrics, collapse = ", "),
                    present_metrics = paste(present_metrics, collapse = ", "),
                    n_missing = length(missing_metrics)
                )
            }
        }
    }

    coverage_df <- dplyr::bind_rows(coverage_list)
    missing_df <- if (length(missing_list) > 0) dplyr::bind_rows(missing_list) else NULL

    # Create wide coverage matrix
    coverage_wide <- coverage_df %>%        
        tidyr::pivot_wider(
            names_from = metric_type,
            values_from = present,
            values_fill = FALSE
        )

    list(
        missing_summary = missing_df,
        coverage_matrix = coverage_wide,
        coverage_long = coverage_df
    )
}


#' Evaluate all metric checkpoints for a genome using overlap and R² metrics
#'
#' Computes peak overlap and R² for all available checkpoints across different
#' metric types (iou, loss, pearson).
#'
#' @param genome_dir Name of the genome directory (e.g., "silicusPlusCGD")
#' @param track_pattern Track name pattern to filter (e.g., "rf524k_EB4_cnt")
#' @param obs_k27_track Observed K27 track for overlap/R² calculation
#' @param obs_k4_track Observed K4 track for overlap/R² calculation
#' @param intervals Intervals for R² calculation (optional, will use cg_trace if NULL)
#' @param context_size Context size for overlap calculation
#' @param T_q Quantile threshold for peak calling
#' @param test_chroms Test chromosomes
#' @param val_chroms Validation chromosomes
#' @param iterator_size Iterator size for gscreen
#' @param track_base_dir Base directory for flashzoi tracks
#' @return Data frame with metrics for each checkpoint
evaluate_checkpoints <- function(genome_dir,
                                  track_pattern = "rf524k_EB4_cnt",
                                  obs_k27_track = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm",
                                  obs_k4_track = "jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3",
                                  intervals = NULL,
                                  context_size = 200,
                                  T_q = 0.98,
                                  test_chroms = c("chr4", "chr10", "chr14", "chr15"),
                                  val_chroms = c("chr8", "chr9"),
                                  iterator_size = 100,
                                  track_base_dir = "/home/aviezerl/mm10/tracks/seq/IQ/pcg/flashzoi") {

    # Get all checkpoints
    checkpoints <- get_metric_checkpoints(genome_dir, track_base_dir, track_pattern)

    if (nrow(checkpoints) == 0) {
        message(sprintf("No checkpoints found for %s with pattern '%s'", genome_dir, track_pattern))
        return(NULL)
    }

    message(sprintf("Evaluating %d checkpoints for %s...", nrow(checkpoints), genome_dir))

    # Load intervals if not provided
    if (is.null(intervals)) {
        cg_trace_f <- load_cg_trace_filtered()
        intervals <- cg_trace_f %>% dplyr::select(chrom, start, end)
    }

    # Split intervals
    train_chroms <- setdiff(unique(intervals$chrom), c(test_chroms, val_chroms))
    intervals_train <- intervals %>% dplyr::filter(chrom %in% train_chroms)
    intervals_test <- intervals %>% dplyr::filter(chrom %in% test_chroms)
    intervals_val <- intervals %>% dplyr::filter(chrom %in% val_chroms)

    # Use plyr::a_ply for parallel evaluation
    results <- plyr::alply(
        checkpoints, 
        1, 
        function(cp) {
            track_name <- cp$track_name
            is_k4 <- cp$is_k4
            obs_track <- if (is_k4) obs_k4_track else obs_k27_track
            obs_expr <- paste0("log2(1 + ", obs_track, ")")
            pred_expr <- paste0("log2(1 + ", track_name, ")")
            i <- which(checkpoints$track_name == track_name)
            message(sprintf("  [%d/%d] %s (%s)...", i, nrow(checkpoints),
                            cp$track_basename, cp$metric_type))

            # Check if track exists
            if (!gtrack.exists(track_name)) {
                message(sprintf("    Track not found: %s", track_name))
                return(NULL)
            }
            # Compute overlap metrics
            overlap <- compute_track_auprc(
                track_name = track_name,
                obs_track_name = obs_track,
                context_size = context_size,
                T_q = T_q,
                test_chroms = test_chroms,
                val_chroms = val_chroms,
                iterator_size = iterator_size,
                intervals = cg_trace_f 
            )
            # Compute R² metrics
            r2_train <- gcor(pred_expr, obs_expr, intervals = intervals_train,
                                iterator = iterator_size)^2
            r2_test <- gcor(pred_expr, obs_expr, intervals = intervals_test,
                            iterator = iterator_size)^2
            r2_val <- gcor(pred_expr, obs_expr, intervals = intervals_val,
                            iterator = iterator_size)^2

            tibble::tibble(
                genome_dir = genome_dir,
                metric_type = cp$metric_type,
                track_basename = cp$track_basename,
                track_name = track_name,
                is_k4 = is_k4,
                mark = if (is_k4) "K4" else "K27",
                auprc_train = overlap["auprc_train"],
                auprc_test = overlap["auprc_test"],
                auprc_val = overlap["auprc_val"],
                r2_train = r2_train,
                r2_test = r2_test,
                r2_val = r2_val
            )
        },
        .parallel = TRUE
    )

    # Gather non-empty results
    results <- Filter(Negate(is.null), results)

    if (length(results) == 0) {
        return(NULL)
    }

    dplyr::bind_rows(results)
}

#' Select best checkpoint based on specified metric
#'
#' @param eval_results Results from evaluate_checkpoints()
#' @param selection_metric Metric to use for selection: "auprc_test", "auprc_val",
#'        "r2_test", "r2_val"
#' @param mark Filter by mark ("K27", "K4", or NULL for both)
#' @return Data frame with best checkpoint per metric_type
select_best_checkpoint <- function(eval_results,
                                    selection_metric = "auprc_test",
                                    mark = NULL) {
    if (is.null(eval_results) || nrow(eval_results) == 0) {
        return(NULL)
    }

    df <- eval_results
    if (!is.null(mark)) {
        df <- df %>% dplyr::filter(mark == !!mark)
    }

    df %>%
        dplyr::group_by(metric_type, mark) %>%
        dplyr::slice_max(order_by = .data[[selection_metric]], n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(mark, dplyr::desc(.data[[selection_metric]]))
}

#' Update symlink to point to the best checkpoint
#'
#' @param genome_dir Name of the genome directory
#' @param track_basename Base name of the best track (e.g., "rf524k_EB4_cnt")
#' @param best_metric_type Metric type directory containing the best checkpoint
#' @param track_base_dir Base directory for flashzoi tracks
#' @param dry_run If TRUE, only print what would be done without making changes
#' @return TRUE if successful, FALSE otherwise
update_checkpoint_symlink <- function(genome_dir,
                                       track_basename,
                                       best_metric_type,
                                       track_base_dir = "/home/aviezerl/mm10/tracks/seq/IQ/pcg/flashzoi",
                                       dry_run = TRUE) {

    genome_path <- file.path(track_base_dir, genome_dir)
    symlink_name <- paste0(track_basename, ".track")
    symlink_path <- file.path(genome_path, symlink_name)
    target_rel <- file.path("metrics", best_metric_type, symlink_name)
    target_abs <- file.path(genome_path, target_rel)

    # Check target exists
    if (!dir.exists(target_abs)) {
        warning(sprintf("Target does not exist: %s", target_abs))
        return(FALSE)
    }

    # Check current symlink
    current_target <- NULL
    if (file.exists(symlink_path) || Sys.readlink(symlink_path) != "") {
        current_target <- Sys.readlink(symlink_path)
    }

    if (identical(current_target, target_rel)) {
        message(sprintf("  [%s] Symlink already points to %s", symlink_name, target_rel))
        return(TRUE)
    }

    if (dry_run) {
        message(sprintf("  [DRY RUN] Would update %s:", symlink_name))
        message(sprintf("    Current: %s", if (is.null(current_target)) "(none)" else current_target))
        message(sprintf("    New:     %s", target_rel))
        return(TRUE)
    }

    # Remove existing symlink if present
    if (!is.null(current_target)) {
        file.remove(symlink_path)
    }

    # Create new symlink (need to be in the directory for relative symlink)
    old_wd <- getwd()
    setwd(genome_path)
    result <- file.symlink(target_rel, symlink_name)
    setwd(old_wd)

    if (result) {
        message(sprintf("  [%s] Updated: %s -> %s", symlink_name, current_target, target_rel))
    } else {
        warning(sprintf("  [%s] Failed to create symlink", symlink_name))
    }

    return(result)
}


#' Evaluate and optionally update symlinks for multiple genomes
#'
#' @param genome_dirs Vector of genome directory names
#' @param track_pattern Track name pattern to filter
#' @param selection_metric Metric to use for selection
#' @param update_symlinks Whether to update symlinks (default FALSE)
#' @param dry_run If TRUE and update_symlinks is TRUE, only show what would change
#' @param ... Additional arguments passed to evaluate_checkpoints()
#' @return Combined evaluation results for all genomes
evaluate_and_select_checkpoints <- function(genome_dirs,
                                           track_pattern = "rf524k_EB4_cnt",
                                           selection_metric = "auprc_test",
                                           update_symlinks = FALSE,
                                           dry_run = TRUE,
                                           ...) {

    process_genome_dir <- function(genome_dir) {
        message(sprintf("\n=== Evaluating %s ===", genome_dir))

        eval_results <- evaluate_checkpoints(
            genome_dir = genome_dir,
            track_pattern = track_pattern,
            intervals = load_cg_trace_filtered(),
            ...
        )

        this_results <- NULL
        this_best <- list()

        if (!is.null(eval_results) && nrow(eval_results) > 0) {
            this_results <- eval_results

            # Select best for each mark
            for (m in c("K27", "K4")) {
                best <- select_best_checkpoint(eval_results, selection_metric, mark = m)
                if (!is.null(best) && nrow(best) > 0) {
                    best_row <- best[1, ]
                    this_best[[paste(genome_dir, m, sep = "_")]] <- best_row

                    message(sprintf("\n  Best %s checkpoint (%s = %.4f):",
                                   m, selection_metric, best_row[[selection_metric]]))
                    message(sprintf("    %s [%s]", best_row$track_basename, best_row$metric_type))

                    if (update_symlinks) {
                        update_checkpoint_symlink(
                            genome_dir = genome_dir,
                            track_basename = best_row$track_basename,
                            best_metric_type = best_row$metric_type,
                            dry_run = dry_run
                        )
                    }
                }
            }
        }

        # return as list to mimic the original structure
        list(
            eval_results = this_results,
            best = if (length(this_best) > 0) dplyr::bind_rows(this_best) else NULL
        )
    }

    res <- plyr::alply(genome_dirs, 1, process_genome_dir, .parallel = TRUE)

    # Combine across all genome_dirs
    all_results_list <- lapply(res, function(x) x$eval_results)
    all_best_list <- lapply(res, function(x) x$best)

    list(
        all_results = if (length(all_results_list) > 0 && !all(sapply(all_results_list, is.null))) 
            dplyr::bind_rows(all_results_list[!sapply(all_results_list, is.null)]) else NULL,
        best_checkpoints = if (length(all_best_list) > 0 && !all(sapply(all_best_list, is.null)))
            dplyr::bind_rows(all_best_list[!sapply(all_best_list, is.null)]) else NULL
    )
}

