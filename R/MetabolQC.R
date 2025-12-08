#' Metabolomics Quality Control (QC) Plots and Summary
#'
#' Generates QC plots and summaries for normalized metabolomics data,
#' including density plots, boxplots, correlation heatmaps, and mean ± SD per metabolite per group.
#'
#' @param normalized_data A \code{metabolNorm} object or list containing an \code{expr_matrix}.
#' @param metadata Data.frame containing at least 'Sample' column and grouping column.
#' @param group_col Character. Column name in \code{metadata} specifying the group/factor. Default "Group".
#' @param corr_type Character. Select betwenn heatmap, correlogram or both.
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

metabol.qc <- function(normalized_data,
                       metadata = NULL,
                       group_col = "Group",
                       corr_type = c("heatmap", "correlogram", "both"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       angle_col = 45,
                       display_corr_values = FALSE) {

  corr_type <- match.arg(corr_type)

  # ============================================================
  # Required packages
  # ============================================================
  req_pkgs <- c("ggplot2", "pheatmap", "tidyr", "reshape2")
  if (corr_type %in% c("correlogram", "both"))
    req_pkgs <- c(req_pkgs, "GGally")

  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs))
    stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # ============================================================
  # Validate input object
  # ============================================================
  if (!is.list(normalized_data) ||
      is.null(normalized_data$expr_matrix)) {
    stop("Input must be the list returned by metabol.normalize() (missing $expr_matrix).")
  }

  data_mat <- normalized_data$expr_matrix

  # ============================================================
  # Metadata auto-retrieval
  # ============================================================
  if (is.null(metadata)) {
    if (!is.null(normalized_data$metadata)) {
      metadata <- normalized_data$metadata
      message("Using metadata stored inside normalized_data.")
    } else {
      stop("Metadata not provided and none found inside normalized_data.")
    }
  }

  if (!"Sample" %in% colnames(metadata))
    stop("metadata must contain a column called 'Sample'.")
  if (!group_col %in% colnames(metadata))
    stop(paste0("metadata must contain a group column '", group_col, "'."))

  # ============================================================
  # Expression matrix cleanup
  # ============================================================
  if ("Sample" %in% colnames(data_mat)) {
    rownames(data_mat) <- data_mat$Sample
    data_mat <- data_mat[, setdiff(colnames(data_mat), "Sample"), drop = FALSE]
  }

  data_mat <- as.matrix(data_mat)
  storage.mode(data_mat) <- "numeric"

  # ============================================================
  # Align metadata and expression matrix
  # ============================================================
  common_samples <- intersect(rownames(data_mat), metadata$Sample)
  if (length(common_samples) == 0)
    stop("No overlapping sample names between metadata and expr_matrix.")

  data_mat <- data_mat[common_samples, , drop = FALSE]
  metadata_sub <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]

  # ============================================================
  # SAMPLE–SAMPLE CORRELATION MATRIX
  # ============================================================
  corr_mat_samples <- stats::cor(
    t(data_mat),
    use = "pairwise.complete.obs",
    method = "pearson"
  )

  corr_labels <- if (display_corr_values) {
    matrix(
      paste0("r=", sprintf("%.2f", corr_mat_samples)),
      nrow = nrow(corr_mat_samples),
      dimnames = dimnames(corr_mat_samples)
    )
  } else {
    FALSE
  }

  # ============================================================
  # Long format for ggplot
  # ============================================================
  df_long <- data.frame(
    Sample = rownames(data_mat),
    data_mat,
    check.names = FALSE
  )

  df_long <- tidyr::pivot_longer(
    df_long,
    cols = -Sample,
    names_to = "Metabolite",
    values_to = "Concentration"
  )

  df_long <- merge(df_long, metadata_sub, by = "Sample", all.x = TRUE, sort = FALSE)

  # ============================================================
  # Density plot
  # ============================================================
  p_density <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Concentration, color = .data[[group_col]], fill = .data[[group_col]])
  ) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Distribution per group",
      x = "Normalized Concentration",
      y = "Density"
    ) +
    ggplot2::theme(legend.position = "bottom")

  print(p_density)

  # ============================================================
  # Boxplot per sample
  # ============================================================
  p_boxplot <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Sample, y = Concentration, fill = Sample)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Boxplot per sample", y = "Normalized Concentration", x = "") +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = angle_col, hjust = 1)
    )

  print(p_boxplot)

  # ============================================================
  # Mean +/- SD barplots per group
  # ============================================================
  histo_list <- list()

  groups_unique <- unique(metadata_sub[[group_col]])

  for (grp in groups_unique) {

    samples_grp <- metadata_sub$Sample[metadata_sub[[group_col]] == grp]
    df_group <- data_mat[samples_grp, , drop = FALSE]

    df_stats <- data.frame(
      Metabolite = colnames(df_group),
      Mean  = apply(df_group, 2, mean),
      SD    = apply(df_group, 2, sd)
    )

    df_stats <- df_stats[order(-df_stats$Mean), ]
    topn <- min(30, nrow(df_stats))
    df_stats2 <- df_stats[seq_len(topn), ]

    p_hist <- ggplot2::ggplot(
      df_stats2,
      ggplot2::aes(x = stats::reorder(Metabolite, -Mean), y = Mean)
    ) +
      ggplot2::geom_col(fill = "skyblue") +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD),
        width = 0.3
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0("Mean +/- SD (Top ", topn, ") - ", group_col, ": ", grp),
        x = "", y = "Normalized concentration"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    histo_list[[grp]] <- p_hist
    print(p_hist)
  }

  # ============================================================
  # HEATMAP
  # ============================================================
  if (corr_type %in% c("heatmap", "both")) {

    cm <- corr_mat_samples

    # Cluster rows/cols
    if (cluster_rows) {
      cm <- cm[stats::hclust(stats::dist(cm))$order, , drop = FALSE]
    }
    if (cluster_cols) {
      cm <- cm[, stats::hclust(stats::dist(t(cm)))$order, drop = FALSE]
    }

    var_names <- colnames(cm)

    corr_long <- reshape2::melt(cm)
    colnames(corr_long) <- c("Var1", "Var2", "value")

    corr_long$Var1 <- factor(corr_long$Var1, levels = var_names)
    corr_long$Var2 <- factor(corr_long$Var2, levels = var_names)

    # keep lower triangle
    corr_long <- corr_long[
      as.numeric(corr_long$Var1) >= as.numeric(corr_long$Var2),
    ]

    p_heatmap <- ggplot2::ggplot(
      corr_long,
      ggplot2::aes(x = Var2, y = Var1, fill = value)
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-1, 1), name = "Correlation"
      ) +
      ggplot2::labs(
        title = "Samples Correlation",
        x = "", y = ""
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::coord_fixed() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text = ggplot2::element_text(size = 9),
        plot.title = ggplot2::element_text(size = 14, hjust = 0.5)
      )

    print(p_heatmap)
  }

  # ============================================================
  # Correlogram (GGally)
  # ============================================================
  correlogram_obj <- NULL

  if (corr_type %in% c("correlogram", "both")) {

    df_corr <- as.data.frame(corr_mat_samples)

    correlogram_obj <- GGally::ggpairs(
      df_corr,
      upper = list(continuous = GGally::wrap("cor", size = 3)),
      diag  = list(continuous = GGally::wrap("densityDiag")),
      lower = list(continuous = GGally::wrap("points", alpha = 0.4, size = 1.5))
    ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = ggplot2::element_text(size = 8),
        panel.grid = ggplot2::element_blank()
      )

    print(correlogram_obj)
  }

  # ============================================================
  # Summary
  # ============================================================
  summary_info <- list(
    n_metabolites = ncol(data_mat),
    n_samples = nrow(data_mat),
    mean_corr_all_samples = mean(corr_mat_samples[base::lower.tri(corr_mat_samples)]),
    mean_corr_by_group = tapply(
      rownames(data_mat),
      metadata_sub[[group_col]],
      function(smp) {
        if (length(smp) < 2) return(NA_real_)
        m <- stats::cor(t(data_mat[smp, , drop = FALSE]))
        mean(m[base::lower.tri(m)])
      }
    )
  )

  invisible(list(
    data_matrix = data_mat,
    corr_matrix = corr_mat_samples,
    plot_density = p_density,
    plot_boxplot = p_boxplot,
    barplots = histo_list,
    correlogram = correlogram_obj,
    summary = summary_info
  ))
}


utils::globalVariables(c("Mean", "Var1", "Var2"))
