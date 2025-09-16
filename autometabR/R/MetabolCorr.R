#' Compute and visualize correlation between metabolites
#'
#' Computes pairwise correlations between metabolites and optionally plots
#' a heatmap with correlation coefficients. Can also summarize by group.
#'
#' @param normalized_data A numeric matrix or a \code{metabolNorm} object containing normalized metabolite data.
#'                        Rows are samples, columns are metabolites. May include a 'Sample' column.
#' @param metadata Optional metadata \code{data.frame} containing at least 'Sample' and 'Group' columns for aggregation.
#' @param group_level Logical, default TRUE. If TRUE and metadata is provided, averages data by group before computing correlations.
#' @param method Correlation method: "pearson", "spearman", "kendall", or "auto" (default). "auto" selects based on normality.
#' @param plot Logical, default TRUE. Whether to plot a heatmap (and scatterplot if 2 metabolites).
#'
#' @return A correlation matrix.
#'
#' @examples
#' \dontrun{
#' # Minimal example
#' set.seed(123)
#'
#' norm_data <- data.frame(
#'   Sample = paste0("S", 1:8),
#'   Met1   = rnorm(8, mean = 5, sd = 1),
#'   Met2   = rnorm(8, mean = 6, sd = 1.2),
#'   Met3   = rnorm(8, mean = 7, sd = 1.5),
#'   Met4   = rnorm(8, mean = 5.5, sd = 0.8),
#'   Met5   = rnorm(8, mean = 6.5, sd = 1.1)
#' )
#'
#' # Example metadata
#' meta <- data.frame(
#'   Sample = paste0("S", 1:8),
#'   Group  = rep(c("A", "B"), each = 4)
#' )
#'
#' correlation <- metabol.corr(norm_data, metadata = meta)
#' print(correlation)
#' }
#'
#' @references
#' Pearson, K. (1895). Note on regression and inheritance in the case of two parents.
#'   \emph{Proceedings of the Royal Society of London}, \strong{58}, 240–242.
#'   <doi:10.1098/rspl.1895.0041>
#'
#' Xia, J., & Wishart, D. S. (2016). Using MetaboAnalyst 3.0 for comprehensive metabolomics data analysis.
#'   \emph{Current Protocols in Bioinformatics}, \strong{55}(1), 14.10.1–14.10.91.
#'   <doi:10.1002/cpbi.11>
#'
#' @export

metabol.corr <- function(normalized_data, metadata = NULL, group_level = TRUE, method = "auto", plot = TRUE) {

  # --- Check required packages ---
  required_pkgs <- c("dplyr","ggplot2","pheatmap","tidyr")
  for(pkg in required_pkgs) {
    if(!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Please install '%s'", pkg))
  }

  # --- Extract numeric matrix if metabolNorm object ---
  if ("metabolNorm" %in% class(normalized_data)) {
    df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
    normalized_data <- normalized_data[[df_names]]
  }

  # --- Set rownames if 'Sample' column exists ---
  if ("Sample" %in% colnames(normalized_data)) {
    expr_df <- dplyr::select(normalized_data, -Sample)
    rownames(expr_df) <- normalized_data$Sample
  } else if ("sample" %in% colnames(normalized_data)) {
    expr_df <- dplyr::select(normalized_data, -sample)
    rownames(expr_df) <- normalized_data$sample
  } else {
    expr_df <- normalized_data
  }

  # --- Aggregate by group if requested ---
  if (!is.null(metadata) && group_level) {
    if (!all(c("Sample","Group") %in% colnames(metadata))) {
      stop("Metadata must contain 'Sample' and 'Group' columns")
    }

    df_long <- tidyr::pivot_longer(
      cbind(Sample = rownames(expr_df), expr_df),
      -Sample,
      names_to = "Metabolite",
      values_to = "Abundance"
    )

    df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

    df_grouped <- dplyr::summarise(
      dplyr::group_by(df_long, Metabolite, Group),
      Abundance = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

    df_grouped <- tidyr::pivot_wider(df_grouped, names_from = Group, values_from = Abundance)
    corr_df <- dplyr::select(df_grouped, -Metabolite)
  } else {
    corr_df <- expr_df
  }

  # --- Choose correlation method ---
  method_used <- if(method == "auto") {
    normality <- apply(corr_df, 2, function(x) tryCatch(stats::shapiro.test(x)$p.value > 0.05, error = function(e) FALSE))
    has_ties <- any(apply(corr_df, 2, function(x) anyDuplicated(x) > 0))
    if(all(normality) && !has_ties) "pearson" else "spearman"
  } else {
    method
  }

  # --- Correlation matrix ---
  corr_mat <- stats::cor(corr_df, method = method_used, use = "pairwise.complete.obs")

  # --- Create display text ---
  corr_text <- matrix(paste0("r = ", round(corr_mat, 2)),
                      nrow = nrow(corr_mat),
                      ncol = ncol(corr_mat),
                      dimnames = dimnames(corr_mat))

  # --- Plot heatmap ---
  if(plot) {
    pheatmap::pheatmap(
      corr_mat,
      display_numbers = corr_text,
      number_color = "black",
      main = paste("Correlation (", method_used, ")", sep = ""),
      cluster_rows = FALSE,
      cluster_cols = FALSE
    )

    # Scatterplot if 2 metabolites
    if(ncol(corr_df) == 2) {
      cols <- colnames(corr_df)
      test <- stats::cor.test(corr_df[[1]], corr_df[[2]], method = method_used)
      r_val <- round(test$estimate, 3)
      p_val <- if(test$p.value < 0.001) "<0.001" else signif(test$p.value, 3)

      strength <- cut(abs(r_val),
                      breaks = c(-Inf, 0.3, 0.5, 0.7, 0.9, Inf),
                      labels = c("very weak or none", "weak", "moderate", "strong", "very strong"),
                      right = FALSE)

      g <- ggplot2::ggplot(corr_df, ggplot2::aes(x = .data[[cols[1]]], y = .data[[cols[2]]])) +
        ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
        ggplot2::geom_smooth(method = ifelse(method_used=="pearson","lm","loess"),
                             se = FALSE, color = "red", linetype = "dashed") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Correlation scatterplot (", method_used, ")", sep=""),
          subtitle = paste0("r = ", r_val, " | p = ", p_val, " | ", strength),
          x = cols[1], y = cols[2]
        )
      print(g)
    }
  }

  return(corr_mat)
}




