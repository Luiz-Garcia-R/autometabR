#' Dimensionality Reduction for Metabolomics Data
#'
#' Performs PCA, UMAP, and PLS-DA on normalized metabolomics data.
#' Returns plots and a list with results including VIP scores for PLS-DA.
#'
#' @param normalized_data A \code{metabolNorm} object or a numeric data.frame/matrix.
#'                        Rows are samples, columns are metabolites.
#' @param metadata A data.frame containing at least a 'Sample' column and a group/condition column.
#' @param n_neighbors Integer, number of neighbors for UMAP. Defaults to half of the number of samples if NULL.
#' @param ncomp_plsda Integer, number of PLS-DA components to compute. Default 2.
#' @param cross_val Integer, number of folds for PLS-DA cross-validation. Default 5.
#' @param verbose Logical, whether to show messages and plots. Default TRUE.
#'
#' @return A list with the following elements (invisible):
#' \describe{
#'   \item{expr_matrix}{The numeric matrix used for dimensionality reduction.}
#'   \item{metadata}{Metadata used in the analysis.}
#'   \item{res_dimred}{A list containing \code{plsda_res}, \code{pca_res}, \code{umap_res}, \code{vip_scores}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Minimal example
#' set.seed(123)
#' expr_mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' rownames(expr_mat) <- paste0("S", 1:20)
#' colnames(expr_mat) <- paste0("Met", 1:5)
#'
#' meta$Sample <- factor(meta$Sample)
#' meta <- data.frame(
#'  Sample = paste0("S", 1:20),
#'  Group  = rep(c("A","B"), each = 10)
#' )
#'
#' res <- metabol.dimred(expr_mat, metadata = meta, cross_val = 3)
#' }
#' @references
#' Jolliffe IT (2002). Principal Component Analysis. 2nd edition. Springer Series in Statistics.
#' Barker M, Rayens W (2003). Partial least squares for discrimination. \emph{Journal of Chemometrics}, 17:166–173. <doi:10.1002/cem.785>
#' McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. \emph{arXiv:1802.03426}.
#'
#' @export

metabol.dimred <- function(normalized_data, metadata,
                           n_neighbors = NULL, ncomp_plsda = 2,
                           cross_val = 5, verbose = TRUE) {

  # --- Check packages ---
  pkgs <- c("ggplot2", "umap", "dplyr", "mixOmics", "caret")
  for(p in pkgs){
    if(!requireNamespace(p, quietly = TRUE)){
      stop(sprintf("Please install '%s'", p))
    }
  }

  # --- Extract expression matrix ---
  if (inherits(normalized_data, "metabolNorm")) {
    expr_mat <- normalized_data$expr_matrix
  } else if (is.data.frame(normalized_data) || is.matrix(normalized_data)) {
    expr_mat <- as.matrix(normalized_data)
  } else {
    stop("Input must be a 'metabolNorm' object or a numeric matrix/data.frame.")
  }

  # --- Detect metadata columns ---
  col_sample <- grep("^Sample$", names(metadata), ignore.case = TRUE, value = TRUE)
  col_group  <- grep("^(Group|Condition|Cluster|Class)$", names(metadata), ignore.case = TRUE, value = TRUE)
  if (length(col_sample) != 1) stop("Unable to detect 'Sample' column in metadata.")
  if (length(col_group) != 1) stop("Unable to detect 'Group/Condition' column in metadata.")

  samples <- metadata[[col_sample]]
  groups  <- as.factor(metadata[[col_group]])

  # --- Ensure sample order matches metadata ---
  if (!all(rownames(expr_mat) == samples)) {
    expr_mat <- expr_mat[samples, , drop = FALSE]
  }

  n_samples <- nrow(expr_mat)

  # --- PCA ---
  pca_res <- stats::prcomp(expr_mat, center = TRUE, scale. = TRUE)
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  df_pca <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                       Group = groups, Sample = samples)
  p_pca <- ggplot2::ggplot(df_pca, ggplot2::aes(PC1, PC2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PCA - Dimensionality Reduction",
                  x = sprintf("PC1 (%.1f%%)", 100*var_explained[1]),
                  y = sprintf("PC2 (%.1f%%)", 100*var_explained[2]))

  # --- UMAP ---
  if (is.null(n_neighbors) || n_neighbors >= n_samples) {
    n_neighbors <- max(2, floor(n_samples / 2))
    if (verbose) message(sprintf("n_neighbors for UMAP set to %d", n_neighbors))
  }
  umap_res <- umap::umap(expr_mat, n_neighbors = n_neighbors, init = "random")
  df_umap <- data.frame(UMAP1 = umap_res$layout[,1], UMAP2 = umap_res$layout[,2],
                        Group = groups, Sample = samples)
  p_umap <- ggplot2::ggplot(df_umap, ggplot2::aes(UMAP1, UMAP2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "UMAP - Dimensionality Reduction")

  # --- PLS-DA ---
  plsda_res <- mixOmics::plsda(X = expr_mat, Y = groups, ncomp = ncomp_plsda)
  df_plsda <- data.frame(plsda_res$variates$X[,1:2], Group = groups, Sample = samples)
  colnames(df_plsda)[1:2] <- c("Comp1", "Comp2")
  vip_scores <- mixOmics::vip(plsda_res)

  # --- Cross-validation ---
  acc <- NA
  try({
    suppressWarnings({
      set.seed(123)
      cv_res <- mixOmics::perf(plsda_res, validation = "Mfold", folds = cross_val, progressBar = FALSE)
    })
    acc <- round(mean(cv_res$error.rate$overall), 3)
  }, silent = TRUE)

  subtitle_plsda <- paste0("Accuracy: ", acc)

  p_plsda <- ggplot2::ggplot(df_plsda, ggplot2::aes(Comp1, Comp2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::stat_ellipse(level = 0.95, linetype = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PLS-DA Score Plot",
                  subtitle = subtitle_plsda,
                  x = "Component 1", y = "Component 2")

  # --- Print plots if verbose ---
  if(verbose){
    print(p_pca)
    print(p_umap)
    print(p_plsda)
  }

  # --- Return results invisibly ---
  invisible(list(
    expr_matrix = expr_mat,
    metadata = metadata,
    res_dimred = list(
      plsda_res = plsda_res,
      pca_res = pca_res,
      umap_res = umap_res,
      vip_scores = vip_scores
    )
  ))
}
