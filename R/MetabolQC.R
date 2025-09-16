#' Metabolomics Quality Control (QC) Plots and Summary
#'
#' Generates QC plots and summaries for normalized metabolomics data,
#' including density plots, boxplots, correlation heatmaps, and mean ± SD per metabolite per group.
#'
#' @param normalized_data A \code{metabolNorm} object or list containing an \code{expr_matrix}.
#' @param metadata Data.frame containing at least 'Sample' column and grouping column.
#' @param group_col Character. Column name in \code{metadata} specifying the group/factor. Default "Group".
#' @param cluster_rows Logical. Whether to cluster rows in correlation heatmap. Default FALSE.
#' @param cluster_cols Logical. Whether to cluster columns in correlation heatmap. Default FALSE.
#' @param angle_col Numeric. Angle for x-axis labels. Default 45.
#' @param display_corr_values Logical. Whether to display correlation values in heatmap. Default FALSE.
#'
#' @return Invisibly returns a list containing:
#' \item{data_ready}{Normalized expression matrix used for QC.}
#' \item{presence_per_group}{Presence/absence matrix per metabolite per group.}
#' \item{corr_matrix_samples}{Pearson correlation matrix across samples.}
#' \item{plot_density}{Density plot of concentrations per group.}
#' \item{plot_boxplot}{Boxplot of concentrations per sample.}
#' \item{histograms_group}{List of mean ± SD barplots per metabolite per group.}
#' \item{summary}{QC summary statistics previously printed to console.}
#'
#' @examples
#' \dontrun{
#' # Minimal example
#' set.seed(123)
#' raw_data <- data.frame(
#'   Sample = paste0("S", 1:10),
#'   Met1 = rnorm(10),
#'   Met2 = rnorm(10)
#' )
#'
#' meta <- data.frame(
#'   Sample = paste0("S", 1:10),
#'   Group = rep(c("A","B"), each = 5)
#' )
#'
#' imp_data <- list(data = raw_data, metadata = meta)
#' class(imp_data) <- "metabolR"
#'
#' norm_obj <- metabol.normalize(imp_data)
#'
#' metabol.qc(norm_obj, metadata = meta)
#' }
#'
#' @references
#' Dunn WB, Broadhurst D, Begley P, et al. (2011). Procedures for large-scale metabolic profiling of serum and plasma using gas chromatography and liquid chromatography coupled to mass spectrometry. \emph{Nature Protocols}, 6(7): 1060–1083. <doi:10.1038/nprot.2011.335>
#' Broadhurst D, Goodacre R, Reinke SN, et al. (2018). Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. \emph{Metabolomics}, 14: 72. <doi:10.1007/s11306-018-1367-3>
#'
#' @export

metabol.qc <- function(normalized_data, metadata,
                       group_col = "Group",
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       angle_col = 45,
                       display_corr_values = FALSE) {

  # --- Required packages ---
  req_pkgs <- c("ggplot2", "pheatmap", "tidyr")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # --- Extract numeric matrix ---
  if (!is.list(normalized_data) || is.null(normalized_data$expr_matrix)) {
    stop("Input must be a list returned by metabol.normalize() (containing expr_matrix).")
  }
  data_mat <- normalized_data$expr_matrix
  if (!is.matrix(data_mat) && !is.data.frame(data_mat)) data_mat <- as.matrix(data_mat)

  # --- Align metadata ---
  if (!"Sample" %in% colnames(metadata)) stop("`metadata` must contain a column named 'Sample'.")
  if (!group_col %in% colnames(metadata)) stop(paste0("`metadata` must contain a column named '", group_col, "'."))
  common_samples <- intersect(rownames(data_mat), metadata$Sample)
  if (length(common_samples) == 0) stop("No sample names match between metadata$Sample and expression matrix rownames.")
  metadata_sub <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]
  data_mat <- data_mat[common_samples, , drop = FALSE]

  # --- Correlation matrix across samples ---
  corr_mat_samples <- stats::cor(t(as.matrix(data_mat)), use = "pairwise.complete.obs", method = "pearson")

  # --- Optional display of correlation values ---
  corr_text <- if(display_corr_values) {
    matrix(paste0("r = ", formatC(round(corr_mat_samples, 2), format = "f", digits = 2)),
           nrow = nrow(corr_mat_samples),
           ncol = ncol(corr_mat_samples),
           dimnames = dimnames(corr_mat_samples))
  } else FALSE

  # --- Long format for ggplot ---
  df_long <- data.frame(Sample = rownames(data_mat), data_mat, check.names = FALSE, stringsAsFactors = FALSE)
  df_long <- tidyr::pivot_longer(df_long, cols = -Sample, names_to = "Metabolite", values_to = "Concentration")
  df_long <- merge(df_long, metadata_sub, by = "Sample", all.x = TRUE, sort = FALSE)

  # --- Density plot per group ---
  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(x = Concentration, color = .data[[group_col]], fill = .data[[group_col]])) +
    ggplot2::geom_density(alpha = 0.3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Distribution per group", x = "Log2 Concentration", y = "Density", color = group_col, fill = group_col) +
    ggplot2::theme(legend.position = "bottom")
  print(p_density)

  # --- Boxplot per sample ---
  p_boxplot <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = Concentration, fill = Sample)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Boxplot per sample", x = "", y = "Log2 Concentration") +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = angle_col, hjust = 1))
  print(p_boxplot)

  # --- Correlation heatmap ---
  pheatmap::pheatmap(
    corr_mat_samples,
    main = "Correlation heatmap across all samples",
    silent = FALSE,
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    clustering_method = "complete",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    fontsize = 10,
    legend = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = ifelse(cluster_rows, 30, 0),
    treeheight_col = ifelse(cluster_cols, 30, 0),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    display_numbers = corr_text,
    number_color = "black",
    angle_col = angle_col
  )

  # --- Presence/absence per group ---
  presence_per_group <- t(sapply(colnames(data_mat), function(metab) {
    tapply(data_mat[, metab], metadata_sub[[group_col]], function(x) as.integer(any(x > 0)))
  }))

  # --- Mean +/- SD histogram per group ---
  histo_list <- list()
  for (grp in unique(metadata_sub[[group_col]])) {
    samples_in_group <- metadata_sub$Sample[metadata_sub[[group_col]] == grp]
    df_group <- data_mat[samples_in_group, , drop = FALSE]
    df_plot <- data.frame(
      Metabolite = colnames(df_group),
      Media = apply(df_group, 2, mean, na.rm = TRUE),
      SD = apply(df_group, 2, sd, na.rm = TRUE)
    )
    p_hist <- ggplot2::ggplot(df_plot, ggplot2::aes(x = reorder(Metabolite, -Media), y = Media)) +
      ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Media - SD, ymax = Media + SD), width = 0.4, color = "black") +
      ggplot2::labs(title = paste("Mean +/- SD per metabolite -", group_col, grp), x = "", y = "Log2 Concentration") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    histo_list[[grp]] <- p_hist
    print(p_hist)
  }

  # --- QC summary ---
  n_metabs <- ncol(data_mat)
  n_samples <- nrow(data_mat)
  mean_global <- mean(as.matrix(data_mat), na.rm = TRUE)
  mean_corr_all <- mean(corr_mat_samples[lower.tri(corr_mat_samples)], na.rm = TRUE)

  summary_info <- list(
    n_metabolites = n_metabs,
    n_samples = n_samples,
    mean_global = mean_global,
    mean_corr_all = mean_corr_all,
    mean_corr_by_group = setNames(
      lapply(unique(metadata_sub[[group_col]]), function(grp) {
        samples_in_group <- metadata_sub$Sample[metadata_sub[[group_col]] == grp]
        if(length(samples_in_group) > 1) {
          corr_grp <- stats::cor(t(data_mat[samples_in_group, , drop = FALSE]), use = "pairwise.complete.obs")
          mean(corr_grp[lower.tri(corr_grp)], na.rm = TRUE)
        } else NA
      }),
      unique(metadata_sub[[group_col]])
    )
  )

  if (mean_corr_all < 0.7) message("Warning: Low average correlation may indicate variability or normalization issues.")

  # --- Return invisibly ---
  invisible(list(
    data_ready = data_mat,
    presence_per_group = presence_per_group,
    corr_matrix_samples = corr_mat_samples,
    plot_density = p_density,
    plot_boxplot = p_boxplot,
    histograms_group = histo_list,
    summary = summary_info
  ))
}
