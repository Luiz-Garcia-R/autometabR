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
#' expr_mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' rownames(expr_mat) <- paste0("S", 1:5)
#' colnames(expr_mat) <- paste0("Met", 1:4)
#' meta <- data.frame(
#'   Sample = paste0("S", 1:5),
#'   Group = c("A","A","B","B","B")
#' )
#'
#' metabol.heatmap(expr_mat, meta, top_n = 3, results = TRUE)
#'
#' @references
#' Eisen MB, Spellman PT, Brown PO, Botstein D (1998). Cluster analysis and display of genome-wide expression patterns. \emph{Proceedings of the National Academy of Sciences}, 95(25):14863–14868. <doi:10.1073/pnas.95.25.14863>
#' Xia J, Wishart DS (2016). Using MetaboAnalyst 3.0 for Comprehensive Metabolomics Data Analysis. \emph{Current Protocols in Bioinformatics}, 55:14.10.1–14.10.91. <doi:10.1002/cpbi.11>
#'
#' @export

metabol.heatmap <- function(normalized_data, metadata,
                            group_colname = "Group",
                            top_n = 100,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            show_rownames = FALSE,
                            show_colnames = TRUE,
                            palette = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
                            results = FALSE,
                            assign_result = FALSE,
                            assign_name = "heatmap_data",
                            envir = parent.frame()) {

  # --- Required package ---
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install the 'pheatmap' package.")

  # --- Extract expression matrix ---
  expr_mat <- NULL
  if (is.list(normalized_data)) {
    if ("expr_matrix" %in% names(normalized_data)) {
      expr_mat <- normalized_data$expr_matrix
    } else {
      df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
      if (length(df_names) >= 1) {
        tmp <- normalized_data[[df_names[1]]]
        if (is.data.frame(tmp) && "Sample" %in% colnames(tmp)) {
          expr_mat <- tmp
          rownames(expr_mat) <- expr_mat$Sample
          expr_mat$Sample <- NULL
        }
      }
      if (is.null(expr_mat)) {
        df_candidates <- Filter(function(x) is.data.frame(x) && "Sample" %in% colnames(x), normalized_data)
        if (length(df_candidates) >= 1) {
          tmp <- df_candidates[[1]]
          expr_mat <- tmp
          rownames(expr_mat) <- expr_mat$Sample
          expr_mat$Sample <- NULL
        }
      }
    }
  } else if (is.data.frame(normalized_data) || is.matrix(normalized_data)) {
    expr_mat <- as.data.frame(normalized_data, stringsAsFactors = FALSE, check.names = FALSE)
    if ("Sample" %in% colnames(expr_mat)) {
      rownames(expr_mat) <- expr_mat$Sample
      expr_mat$Sample <- NULL
    }
  }

  if (is.null(expr_mat)) stop("Invalid input: cannot extract expression matrix.")

  # --- Ensure numeric matrix ---
  numeric_cols <- vapply(expr_mat, is.numeric, logical(1))
  expr_mat <- expr_mat[, numeric_cols, drop = FALSE]
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  if (is.null(rownames(expr_mat))) stop("Expression matrix must have sample names as rownames.")

  # --- Metadata checks ---
  if (!"Sample" %in% colnames(metadata)) stop("`metadata` must contain a column named 'Sample'.")
  if (!group_colname %in% colnames(metadata)) stop(sprintf("`metadata` must contain the column '%s'.", group_colname))

  # --- Align samples ---
  common_samples <- intersect(rownames(expr_mat), metadata$Sample)
  if (length(common_samples) == 0) stop("No sample names match between metadata$Sample and expression matrix rownames.")
  metadata_sub <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]
  expr_mat <- expr_mat[common_samples, , drop = FALSE]

  # --- Select top N metabolites ---
  var_metabs <- apply(expr_mat, 2, stats::var, na.rm = TRUE)
  n_sel <- min(top_n, length(var_metabs))
  top_metabs <- names(sort(var_metabs, decreasing = TRUE))[1:n_sel]

  # --- Z-score and transpose for pheatmap ---
  z_mat <- scale(expr_mat[, top_metabs, drop = FALSE], center = TRUE, scale = TRUE)
  expr_z_top <- t(z_mat)

  # --- Column annotation ---
  sample_groups <- metadata_sub[[group_colname]]
  names(sample_groups) <- metadata_sub$Sample
  ann_col <- data.frame(Group = as.factor(sample_groups[colnames(expr_z_top)]))
  rownames(ann_col) <- colnames(expr_z_top)

  # --- Generate heatmap ---
  pheatmap::pheatmap(expr_z_top,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     annotation_col = ann_col,
                     color = palette,
                     border_color = NA,
                     main = sprintf("Top %d Most Variable Metabolites", n_sel),
                     fontsize = 10)

  # --- Return object ---
  result <- list(
    top_metabolites = top_metabs,
    expr_z_score = expr_z_top,
    annotation_col = ann_col,
    var_table = data.frame(Metabolite = names(var_metabs), Variance = var_metabs)[match(top_metabs, names(var_metabs)), , drop = FALSE]
  )

  if (assign_result) {
    assign(assign_name, result, envir = envir)
  }

  if (results) return(result) else invisible(NULL)
}
