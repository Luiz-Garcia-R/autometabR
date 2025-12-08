#' Boxplot/Violin Plot for a Single Metabolite with Limma Significance
#'
#' @description
#' Generates a violin + boxplot visualization for a selected metabolite across
#' two experimental groups. The function integrates expression values,
#' metadata, and `limma` differential analysis results to annotate the plot
#' with significance levels (*, **, ***, ns).
#'
#' @details
#' This function expects `normalized_data` to follow the structure produced by
#' your metabolomics workflow, containing at minimum:
#'
#' - `expr_matrix`: a matrix or data frame with a column `"Sample"` and one
#'   column per metabolite.
#' - `metadata` (optional): a data frame with sample-level metadata, including
#'   a group column (default `"Group"`).
#' - `limma_results` (optional): a list containing a data frame `results`
#'   with columns `"Metabolite"` and `"P.Value"`.
#'
#' If metadata or limma results are not supplied directly, the function tries to
#' extract them from `normalized_data`.
#'
#' The function supports **exactly two groups** for plotting and significance
#' annotation.
#'
#' @param normalized_data A list containing at least `expr_matrix`, and optionally
#'   `metadata` and `limma_results`.
#' @param metadata Optional data frame with sample metadata. Must contain a `"Sample"`
#'   column and the grouping column specified in `group_col`.
#' @param metabolite Character string indicating the metabolite to plot.
#'   Must match a column name in `expr_matrix`.
#' @param group_col Character string indicating the grouping variable in metadata.
#'   Default is `"Group"`.
#' @param limma_results Optional list containing a `results` data frame with
#'   `Metabolite` and `P.Value` columns. If omitted, the function will try to
#'   use `normalized_data$limma_results`.
#'
#' @return
#' A `ggplot2` object containing a violin + boxplot with significance annotation.
#'
#' @examples
#' \dontrun{
#' metabol.boxplot(
#'   normalized_data = my_norm,
#'   metabolite = "Lactate",
#'   group_col = "Condition"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_point annotate
#'   scale_fill_manual theme_minimal labs element_text position_jitter
#' @importFrom dplyr left_join
#' @importFrom scales hue_pal
#'
#' @export


metabol.boxplot <- function(normalized_data,
                            metadata = NULL,
                            metabolite,
                            group_col = "Group",
                            limma_results = NULL) {

  # ------------------------------------------------------------
  # 0) Dependencies (CRAN-friendly)
  # ------------------------------------------------------------
  pkgs <- c("ggplot2", "dplyr", "scales")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))
  }

  # ------------------------------------------------------------
  # 1) Metadata fallback
  # ------------------------------------------------------------
  if (is.null(metadata)) {
    if (!is.null(normalized_data$metadata)) {
      metadata <- normalized_data$metadata
    } else {
      stop("metadata must be provided or available in normalized_data$metadata.")
    }
  }

  # ------------------------------------------------------------
  # 2) Extract expr_matrix
  # ------------------------------------------------------------
  if (is.null(normalized_data$expr_matrix)) {
    stop("normalized_data$expr_matrix not found.")
  }

  expr <- normalized_data$expr_matrix

  if (!"Sample" %in% colnames(expr)) {
    stop("normalized_data$expr_matrix must contain a column named 'Sample'.")
  }

  rownames(expr) <- expr[["Sample"]]
  expr <- expr[, setdiff(colnames(expr), "Sample"), drop = FALSE]

  if (!metabolite %in% colnames(expr)) {
    stop(sprintf("Metabolite '%s' not found in expr_matrix.", metabolite))
  }

  # ------------------------------------------------------------
  # 3) Merge with metadata
  # ------------------------------------------------------------
  df <- data.frame(
    Sample = rownames(expr),
    value = expr[[metabolite]],
    stringsAsFactors = FALSE
  )

  df <- dplyr::left_join(df, metadata, by = "Sample")

  if (!group_col %in% colnames(df)) {
    stop(sprintf("Group column '%s' not found in metadata.", group_col))
  }

  df$group <- df[[group_col]]

  groups <- unique(df$group)
  if (length(groups) != 2) {
    stop("Exactly two groups are required for the comparison.")
  }

  # ------------------------------------------------------------
  # 4) Get limma results
  # ------------------------------------------------------------
  if (is.null(limma_results)) {
    if (!is.null(normalized_data$limma_results)) {
      limma_results <- normalized_data$limma_results
    }
  }

  if (is.null(limma_results) || is.null(limma_results$results)) {
    stop("No limma results found. Run metabol.diff() or provide limma_results explicitly.")
  }

  lim <- limma_results$results

  lim_sub <- lim[lim$Metabolite == metabolite, , drop = FALSE]
  if (nrow(lim_sub) == 0) {
    stop(sprintf("Metabolite '%s' not found in limma results.", metabolite))
  }

  pval <- lim_sub$P.Value[1]

  signif_label <- if (pval < 0.001) {
    "***"
  } else if (pval < 0.01) {
    "**"
  } else if (pval < 0.05) {
    "*"
  }

  p_label <- signif(pval, 3)

  # ------------------------------------------------------------
  # 5) Aesthetic parameters
  # ------------------------------------------------------------
  y_max <- max(df$value, na.rm = TRUE)
  y_pos <- if (is.finite(y_max) && y_max > 0) y_max * 1.05 else 1

  colors <- scales::hue_pal()(length(groups))

  # ------------------------------------------------------------
  # 6) Final plot
  # ------------------------------------------------------------
  g <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.55, color = NA, adjust = 0.6) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA,
                          color = "gray20", linewidth = 0.4) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.1),
      alpha = 0.4, size = 1.8, color = "gray25"
    ) +
    ggplot2::annotate("text", x = 1.5, y = y_pos,
                      label = signif_label, size = 7) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = metabolite,
      subtitle = sprintf("Limma p-value: %s", p_label),
      x = NULL,
      y = "Normalized value"
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12)
    )

  g
}
