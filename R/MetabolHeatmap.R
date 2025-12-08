#' Heatmap of Top Variable Metabolites
#'
#' Generates a heatmap of the top N most variable metabolites from a normalized metabolomics dataset.
#' Useful for visualizing patterns of variation across samples or groups.
#'
#' @param normalized_data A \code{metabolNorm} object, a list containing a normalized data.frame,
#'                        or a numeric data.frame/matrix (samples x metabolites).
#' @param metadata A data.frame with at least columns 'Sample' and a group column.
#' @param group_colname Character, name of the group column in metadata. Default "Group".
#' @param top_n Integer, number of top variable metabolites to include in heatmap. Default 100.
#' @param cluster_rows Logical, whether to cluster rows. Default TRUE.
#' @param cluster_cols Logical, whether to cluster columns. Default TRUE.
#' @param show_rownames Logical, show metabolite names. Default FALSE.
#' @param show_colnames Logical, show sample names. Default TRUE.
#' @param diff_results optional, use DEGs identified by \code{metabol.diff}.
#' @param palette Color palette for heatmap. Default blue-white-red.
#' @param results Logical, whether to return a list with z-score matrix and top metabolites. Default FALSE.
#' @param assign_result Logical, whether to assign the result object to the environment. Default FALSE.
#' @param assign_name Character, name for the assigned object. Default "heatmap_data".
#' @param envir Environment in which to assign result object. Default parent.frame().
#'
#' @return Invisibly returns a list (if \code{results = TRUE}) containing:
#' \describe{
#'   \item{top_metabolites}{Character vector of top variable metabolites.}
#'   \item{expr_z_score}{Z-score matrix of top metabolites (rows = metabolites, cols = samples).}
#'   \item{annotation_col}{Data frame with sample group annotations.}
#'   \item{var_table}{Data frame of variance values for selected metabolites.}
#' }
#'
#' @examples
#' # Minimal example
#' set.seed(123)
#'
#' expr_mat <- matrix(rnorm(48), nrow = 6)
#' rownames(expr_mat) <- paste0("S", 1:6)
#' colnames(expr_mat) <- paste0("M", 1:8)
#'
#' meta <- data.frame(
#'   Sample = rownames(expr_mat),
#'   Group  = rep(c("A","B"), each = 3)
#' )
#'
#' minimal_norm <- list(
#'   expr_matrix = data.frame(
#'     Sample = rownames(expr_mat),
#'     expr_mat,
#'     check.names = FALSE
#'   ),
#'   metadata = meta,
#'   method = "none"
#' )
#' class(minimal_norm) <- "metabolNorm"
#'
#' metabol.heatmap(
#'   normalized_data = minimal_norm,
#'   top_n = 3,
#'   results = FALSE
#' )
#'
#'
#' @references
#' Eisen MB, Spellman PT, Brown PO, Botstein D (1998). Cluster analysis and display of genome-wide expression patterns. \emph{Proceedings of the National Academy of Sciences}, 95(25):14863–14868. <doi:10.1073/pnas.95.25.14863>
#' Xia J, Wishart DS (2016). Using MetaboAnalyst 3.0 for Comprehensive Metabolomics Data Analysis. \emph{Current Protocols in Bioinformatics}, 55:14.10.1–14.10.91. <doi:10.1002/cpbi.11>
#'
#' @export

metabol.heatmap <- function(normalized_data,
                            metadata = NULL,
                            group_colname = "Group",
                            top_n = 100,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            show_rownames = FALSE,
                            show_colnames = TRUE,
                            palette = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
                            diff_results = NULL,
                            results = FALSE,
                            assign_result = FALSE,
                            assign_name = "heatmap_data",
                            envir = parent.frame()) {

  # ============================================================
  # 0) Dependencies
  # ============================================================
  if (!requireNamespace("pheatmap", quietly = TRUE))
    stop("Please install the 'pheatmap' package.")

  # ============================================================
  # 1) Extract metadata if not provided
  # ============================================================
  if (is.null(metadata)) {
    if (!is.null(normalized_data$metadata)) {
      metadata <- normalized_data$metadata
      message("Using metadata stored inside normalized_data.")
    } else {
      stop("No metadata provided and none found inside normalized_data.")
    }
  }

  # ============================================================
  # 2) Extract expression matrix
  # ============================================================
  if (!inherits(normalized_data, "metabolNorm"))
    stop("normalized_data must be a metabolNorm object.")

  if (!"expr_matrix" %in% names(normalized_data))
    stop("normalized_data is metabolNorm but contains no expr_matrix.")

  expr_mat <- normalized_data$expr_matrix

  if (!"Sample" %in% colnames(expr_mat))
    stop("expr_matrix must contain a 'Sample' column.")

  # Standardize
  rownames(expr_mat) <- expr_mat$Sample
  expr_mat$Sample <- NULL

  # numeric-only
  expr_mat <- expr_mat[, vapply(expr_mat, is.numeric, logical(1)), drop = FALSE]
  expr_mat <- as.matrix(expr_mat)

  # ============================================================
  # 3) Ranking using diff_results or fallback to variance
  # ============================================================
  top_metabs <- NULL

  if (is.null(diff_results)) {
    if (!is.null(normalized_data$limma_results) &&
        !is.null(normalized_data$limma_results$results)) {

      cand <- normalized_data$limma_results$results

      if (all(c("Metabolite", "logFC") %in% colnames(cand))) {
        diff_results <- cand
        message("Using log2FC from normalized_data$limma_results to rank metabolites.")
      }
    }
  }

  if (!is.null(diff_results)) {

    if (all(c("Metabolite", "logFC") %in% colnames(diff_results))) {

      ranked <- diff_results[order(abs(diff_results$logFC), decreasing = TRUE), ]
      top_metabs <- utils::head(ranked$Metabolite, top_n)

    } else {
      warning("diff_results provided but missing Metabolite/logFC columns. Ignoring.")
      diff_results <- NULL
    }
  }

  # Fallback = variance-based ranking
  if (is.null(diff_results)) {

    message("Tip: you can provide diff_results or run metabol.diff() to rank metabolites by log2FC.")

    var_metabs <- apply(expr_mat, 2, stats::var, na.rm = TRUE)
    var_metabs <- sort(var_metabs, decreasing = TRUE)

    top_metabs <- names(var_metabs)[seq_len(min(top_n, length(var_metabs)))]
  }

  # ============================================================
  # 4) Align metadata and samples
  # ============================================================
  if (!"Sample" %in% colnames(metadata))
    stop("metadata must contain a 'Sample' column.")

  if (!group_colname %in% colnames(metadata))
    stop(sprintf("metadata must contain the column '%s'.", group_colname))

  common_samples <- intersect(rownames(expr_mat), as.character(metadata$Sample))

  if (length(common_samples) == 0)
    stop("No overlapping sample names between metadata and expression matrix.")

  metadata_sub <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]
  expr_mat <- expr_mat[common_samples, , drop = FALSE]

  # ============================================================
  # 5) Z-score transform and transpose
  # ============================================================
  z_mat <- scale(expr_mat[, top_metabs, drop = FALSE], center = TRUE, scale = TRUE)
  expr_z_top <- t(z_mat)

  # ============================================================
  # 6) Column annotation
  # ============================================================
  sample_groups <- metadata_sub[[group_colname]]
  names(sample_groups) <- metadata_sub$Sample

  ann_col <- data.frame(Group = as.factor(sample_groups[colnames(expr_z_top)]))
  rownames(ann_col) <- colnames(expr_z_top)

  # ============================================================
  # 7) Heatmap
  # ============================================================
  pheatmap::pheatmap(expr_z_top,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     annotation_col = ann_col,
                     color = palette,
                     border_color = NA,
                     main = sprintf("Top %d Metabolites", length(top_metabs)),
                     fontsize = 10)

  # ============================================================
  # 8) Return object
  # ============================================================
  result <- list(
    top_metabolites = top_metabs,
    expr_z_score = expr_z_top,
    annotation_col = ann_col
  )

  if (isTRUE(assign_result)) {
    assign(assign_name, result, envir = envir)
  }

  if (isTRUE(results))
    return(result)

  invisible(NULL)
}
