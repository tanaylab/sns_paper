# ==============================================================================
# Borzoi Performance Plotting Functions
# ==============================================================================
#
# Plotting utilities for Borzoi model performance analysis.
# Source this file after loading CONFIG (from the notebook or borzoi_utils.R).
#
# Functions:
#   - plot_performance: Unified plot for overlap and R² metrics
#
# Constants:
#   - LEVELS: Predefined x-axis level orders for different series
#
# ==============================================================================

# ==============================================================================
# Predefined X-Axis Level Orders
# ==============================================================================

LEVELS <- list(
    mm10_minus = c(
        "mm10-CTCF", "mm10-TE", "mm10-Random", "mm10-CRE", "mm10-Exon",
        "mm10-CGD", "mm10-CGD-pad1k", "mm10-CGD-pad2k"
    ),
    silicus_plus = rev(c(
        "silicus", "silicus+CTCF", "silicus+TE", "silicus+Random",
        "silicus+Exon", "silicus+CRE", "silicus+CGD"
    )),
    silicus_telescope = c(
        "silicus", "silicus+CGD", "silicus+CGD+CRE", "silicus+CGD+CRE+CTCF",
        "silicus+CGD+CRE+CTCF+Exon", "silicus+CGD+CRE+CTCF+Exon+TE"
    ),
    fm_silicus = c(
        "silicus", "silicus+GC", "silicus+CGD+CTCF", "silicus+CGD+CTCF+CRE",
        "silicus+CGD+CTCF+CRE+Exon", "silicus+CGD+CTCF+CRE+Exon+TE",
        "silicus+CGD+CTCF+CRE+Exon+TE+Utr3"
    )
)

# ==============================================================================
# Unified Performance Plotting Function
# ==============================================================================

#' Unified performance plot for overlap, R², or AUPRC metrics
#'
#' @param data Data frame with metrics (train_p/test_p/val_p, r2_train/r2_test/r2_val, or auprc_train/auprc_test/auprc_val)
#' @param metric "overlap", "rsqr", or "auprc"
#' @param mark_label Label for histone mark (e.g., "H3K27me3")
#' @param series_filter Series to filter (NULL for no filtering)
#' @param title_suffix Suffix for the plot title
#' @param x_var Variable for x-axis ("genome_label" or "rf_numeric")
#' @param x_levels Factor levels for x-axis ordering
#' @param baseline_data Optional baseline data for reference lines (or "auto" to use series == "baseline")
#' @param show_val Whether to show validation data (default TRUE)
#' @param show_sns Whether to show SNS horizontal lines
#' @param sns_data SNS data frame for reference lines
#' @param x_transform Optional x-axis transformation ("log10", NULL)
#' @param style "standard" or "ppt" for presentation style
#' @param point_size Size of points
#' @param text_size Size of label text
#' @param line_width Width of lines
#' @param colors Named list with train, test, val colors (defaults to CONFIG$colors if available)
#'
#' @return ggplot object
#' @export
plot_performance <- function(data, metric = c("overlap", "rsqr", "auprc"),
                              mark_label, series_filter = NULL,
                              title_suffix = "",
                              x_var = "genome_label", x_levels = NULL,
                              baseline_data = "auto",
                              show_val = TRUE, show_sns = FALSE, sns_data = NULL,
                              x_transform = NULL,
                              style = c("standard", "ppt"),
                              point_size = 4, text_size = 3, line_width = 1.2,
                              colors = NULL, show_labels = NULL) {

    metric <- match.arg(metric)
    style <- match.arg(style)

    # Get colors from CONFIG if not provided
    if (is.null(colors)) {
        if (exists("CONFIG") && !is.null(CONFIG$colors)) {
            colors <- CONFIG$colors
        } else {
            colors <- list(train = "#2166AC", test = "#D6604D", val = "#7570B3")
        }
    }

    # Configure columns based on metric
    if (metric == "overlap") {
        train_col <- "train_p"
        test_col <- "test_p"
        val_col <- "val_p"
        y_label <- "Peak Overlap"
        y_limits <- c(0, 1)
    } else if (metric == "auprc") {
        train_col <- "auprc_train"
        test_col <- "auprc_test"
        val_col <- "auprc_val"
        y_label <- "AUPRC"
        y_limits <- c(0, 1)
    } else {
        train_col <- "r2_train"
        test_col <- "r2_test"
        val_col <- "r2_val"
        y_label <- expression(R^2)
        y_limits <- c(0, 1)
    }

    # Get baseline values for reference lines
    mm10_train <- mm10_test <- mm10_val <- NA
    if (identical(baseline_data, "auto") && "series" %in% names(data)) {
        baseline_df <- data %>% filter(series == "baseline")
        if (nrow(baseline_df) > 0) {
            mm10_train <- baseline_df[[train_col]][1]
            mm10_test <- baseline_df[[test_col]][1]
            mm10_val <- if (show_val && val_col %in% names(baseline_df)) baseline_df[[val_col]][1] else NA
        }
    } else if (is.data.frame(baseline_data) && nrow(baseline_data) > 0) {
        mm10_train <- baseline_data[[train_col]][1]
        mm10_test <- baseline_data[[test_col]][1]
        mm10_val <- if (show_val && val_col %in% names(baseline_data)) baseline_data[[val_col]][1] else NA
    }

    # Filter and prepare data
    plot_data <- data
    if (!is.null(series_filter) && "series" %in% names(data)) {
        plot_data <- data %>% filter(series == series_filter)
    }

    if ("order" %in% names(plot_data)) {
        plot_data <- plot_data %>% arrange(order)
    }

    if (x_var %in% names(plot_data)) {
        if (!is.null(x_levels)) {
            plot_data <- plot_data %>%
                mutate(!!sym(x_var) := factor(.data[[x_var]], levels = x_levels))
        } else if (is.character(plot_data[[x_var]])) {
            plot_data <- plot_data %>%
                mutate(!!sym(x_var) := factor(.data[[x_var]], levels = unique(.data[[x_var]])))
        }
    }

    # Pivot to long format
    cols_to_pivot <- c(train_col, test_col)
    name_map <- setNames(c("Train", "Test"), c(train_col, test_col))

    if (show_val && val_col %in% names(plot_data)) {
        cols_to_pivot <- c(cols_to_pivot, val_col)
        name_map <- c(name_map, setNames("Val", val_col))
    }

    plot_df <- plot_data %>%
        pivot_longer(cols = all_of(cols_to_pivot), names_to = "dataset", values_to = "value") %>%
        mutate(dataset = name_map[dataset])

    # Style adjustments
    if (style == "ppt") {
        base_size <- 16
        title_size <- 20
        axis_text_size <- 14
        legend_pos <- "top"
        if (is.null(show_labels)) {
            show_labels <- FALSE
        }
    } else {
        base_size <- 13
        title_size <- 16
        axis_text_size <- 11
        legend_pos <- "top"
        if (is.null(show_labels)) {
            show_labels <- TRUE
        }
    }

    # Build plot
    p <- ggplot(plot_df, aes(x = .data[[x_var]], y = value, color = dataset, group = dataset))

    # Add baseline reference lines if available
    if (!is.na(mm10_train)) {
        anchor_x <- if (is.factor(plot_df[[x_var]])) {
            levels(plot_df[[x_var]])[max(1, length(levels(plot_df[[x_var]])) - 1)]
        } else {
            max(plot_df[[x_var]], na.rm = TRUE)
        }

        p <- p +
            geom_hline(yintercept = mm10_train, linetype = "dashed",
                       color = colors$train, alpha = 0.4, linewidth = 0.8) +
            geom_hline(yintercept = mm10_test, linetype = "dashed",
                       color = colors$test, alpha = 0.4, linewidth = 0.8) +
            geom_text(
                data = data.frame(x = anchor_x, y = mm10_train),
                aes(x = x, y = y, label = paste0("mm10 Train: ", round(mm10_train, 3))),
                color = colors$train, vjust = -0.5, hjust = 0,
                fontface = "italic", size = text_size, inherit.aes = FALSE
            ) +
            geom_text(
                data = data.frame(x = anchor_x, y = mm10_test),
                aes(x = x, y = y, label = paste0("mm10 Test: ", round(mm10_test, 3))),
                color = colors$test, vjust = -0.5, hjust = 0,
                fontface = "italic", size = text_size, inherit.aes = FALSE
            )

        if (show_val && !is.na(mm10_val)) {
            p <- p +
                geom_hline(yintercept = mm10_val, linetype = "dashed",
                           color = colors$val, alpha = 0.4, linewidth = 0.8) +
                geom_text(
                    data = data.frame(x = anchor_x, y = mm10_val),
                    aes(x = x, y = y, label = paste0("mm10 Val: ", round(mm10_val, 3))),
                    color = colors$val, vjust = 1.5, hjust = 0,
                    fontface = "italic", size = text_size, inherit.aes = FALSE
                )
        }
    }

    # Add SNS reference lines for K27
    if (show_sns && !is.null(sns_data) && mark_label == "H3K27me3") {
        sns_colors <- c("sns_lm" = "#E7298A", "sns_brz2k" = "#66A61E")
        sns_labels <- c("sns_lm" = "SNS (linear)", "sns_brz2k" = "SNS (brz2k)")

        for (i in seq_len(nrow(sns_data))) {
            col <- sns_data$col_name[i]
            sns_test <- sns_data[[test_col]][i]
            col_color <- sns_colors[col]

            p <- p +
                geom_hline(yintercept = sns_test, linetype = "dotdash",
                           color = col_color, alpha = 0.6, linewidth = 0.8) +
                geom_text(
                    data = data.frame(x = 1, y = sns_test),
                    aes(x = x, y = y, label = paste0(sns_labels[col], " Test: ", round(sns_test, 3))),
                    color = col_color, vjust = -0.5, hjust = 0,
                    fontface = "italic", size = text_size - 0.5, inherit.aes = FALSE
                )

            if (show_val && val_col %in% names(sns_data)) {
                sns_val <- sns_data[[val_col]][i]
                p <- p +
                    geom_hline(yintercept = sns_val, linetype = "twodash",
                               color = col_color, alpha = 0.5, linewidth = 0.8)
            }
        }
    }

    # Add main plot elements
    p <- p +
        geom_line(linewidth = line_width, alpha = 0.8) +
        geom_point(size = point_size) 

    if (show_labels) {
        p <- p + 
            geom_text_repel(
                aes(label = round(value, 3)),
                size = text_size,
                max.overlaps = 20,
                box.padding = 0.3,
                show.legend = FALSE
            )
    }

    p <- p +
        scale_color_manual(values = c(
            "Train" = colors$train,
            "Test" = colors$test,
            "Val" = colors$val
        )) +
        scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.02, 0.1))) +
        labs(
            title = paste0(mark_label, " ", title_suffix),
            x = "",
            y = y_label,
            color = ""
        ) +
        theme_minimal(base_size = base_size) +
        theme(
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
            axis.title = element_text(size = base_size - 1, face = "bold"),
            legend.position = legend_pos,
            panel.grid.major.x = element_blank(),
            panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
        )

    # Add x-axis transformation if specified
    if (!is.null(x_transform) && x_transform == "log10" && is.numeric(plot_df[[x_var]])) {
        p <- p + scale_x_continuous(trans = "log10")
    }

    return(p)
}
