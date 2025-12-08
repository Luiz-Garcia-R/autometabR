#' OPLS-DA Analysis for Metabolomics Data
#'
#' Performs Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#' on normalized metabolomics data and optionally plots the score plot.
#'
#' @param normalized_data A \code{metabolNorm} object or numeric matrix/data.frame of metabolite expression (rows = samples, cols = metabolites).
#' @param metadata Data.frame containing at least a 'Sample' column and a grouping column.
#' @param group_col Character. Column name in \code{metadata} specifying the group/factor. Default "Group".
#' @param predI Integer. Number of predictive components. Default 1.
#' @param orthoI Integer or NA. Number of orthogonal components. Default NA (automatic selection).
#' @param crossvalI Integer. Number of cross-validation folds. Default 5.
#' @param permI Integer. Number of permutations for validation. Default 100.
#' @param plot Logical. Whether to plot the OPLS-DA scores. Default TRUE.
#'
#' @return Invisibly returns a list containing:
#' \item{oplsda_res}{The ropls::opls model object.}
#' \item{scores_df}{Data.frame of predictive and orthogonal scores with group labels.}
#' \item{plot}{ggplot object (if \code{plot = TRUE}), else NULL.}
#'
#' @examples
#' \dontrun{
#' # Minimal example
#' if (requireNamespace("ropls", quietly = TRUE)) {
#'   set.seed(123)
#'   raw_data <- data.frame(
#'     Sample = paste0("S", 1:20),
#'     Met1 = c(rep(2,10), rep(-2,10)),
#'     Met2 = c(rep(1,10), rep(-1,10))
#'   )
#'
#'   meta <- data.frame(
#'     Sample = paste0("S", 1:20),
#'     Group  = rep(c("A","B"), each = 10)
#'   )
#'
#'   imp_data <- list(data = raw_data, metadata = meta)
#'   class(imp_data) <- "metabolR"
#'
#'   norm_obj <- metabol.normalize(imp_data, assign_results = FALSE)
#'
#'   metabol.oplsda(norm_obj, meta)
#'   }
#' }
#' @references
#' Trygg J, Wold S. (2002). Orthogonal projections to latent structures (O-PLS). \emph{Journal of Chemometrics}, 16(3): 119–128. <doi:10.1002/cem.695>
#' Worley B, Powers R. (2013). Multivariate Analysis in Metabolomics. \emph{Current Metabolomics}, 1(1): 92–107. <doi:10.2174/2213235X11301010092>
#'
#' @export

metabol.oplsda <- function(normalized_data, metadata = NULL,
                           group_col = "Group",
                           predI = 1, orthoI = NA,
                           crossvalI = 5, permI = 100,
                           plot = TRUE) {

  # --- Packages ---
  pkgs <- c("ggplot2", "dplyr", "ropls")
  for(p in pkgs){
    if(!requireNamespace(p, quietly = TRUE)){
      stop(sprintf("Please install '%s'", p))
    }
  }

  # --- Expression matrix ---
  if (inherits(normalized_data, "metabolNorm")) {
    expr_mat <- normalized_data$expr_matrix
  } else if (is.matrix(normalized_data) || is.data.frame(normalized_data)) {
    expr_mat <- normalized_data
  } else {
    stop("Input must be a metabolNorm object or matrix/data.frame")
  }

  expr_mat <- as.matrix(expr_mat)

  # --- Remove Sample column ---
  if ("Sample" %in% colnames(expr_mat)) {
    expr_mat <- expr_mat[, setdiff(colnames(expr_mat), "Sample"), drop = FALSE]
  }

  # --- Detect metadata automatically ---
  if (is.null(metadata)) {
    if (!is.null(normalized_data$metadata)) {
      metadata <- normalized_data$metadata
      message("Using metadata stored inside normalized_data.")
    } else {
      stop("Metadata not provided and not found inside normalized_data.")
    }
  }

  # --- Checks ---
  if (!"Sample" %in% colnames(metadata))
    stop("metadata must contain a 'Sample' column.")

  if (!group_col %in% colnames(metadata))
    stop(sprintf("metadata must contain '%s'.", group_col))

  # --- Align samples ---
  common_samples <- intersect(rownames(expr_mat), metadata$Sample)
  if (length(common_samples) == 0) {
    stop("No shared sample names between expression matrix and metadata.")
  }

  expr_mat <- expr_mat[common_samples, , drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]

  groups <- as.factor(metadata[[group_col]])

  # --- Ensure numeric expression matrix ---
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"

  # --- Remove zero-variance variables ---
  nzv <- apply(expr_mat, 2, sd, na.rm = TRUE) > 0
  expr_mat <- expr_mat[, nzv, drop = FALSE]

  # --- Run OPLS-DA ---
  set.seed(123)
  oplsda_res <- ropls::opls(
    x = expr_mat,
    y = groups,
    predI = predI,
    orthoI = orthoI,
    crossvalI = crossvalI,
    permI = permI
  )

  # --- Scores dataframe ---
  scores_df <- as.data.frame(oplsda_res@scoreMN)

  if (ncol(scores_df) == 1) {
    colnames(scores_df) <- "comp1"
    scores_df$ortho1 <- 0
  } else {
    colnames(scores_df) <- c("comp1", "ortho1")
  }

  scores_df$Group <- groups

  # --- Plot ---
  p <- NULL
  if (plot) {
    p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = comp1, y = ortho1, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(
        title = "OPLS-DA Scores",
        x = "Predictive Component",
        y = "Orthogonal Component",
        color = "Group"
      )
    print(p)
  }

  invisible(list(
    oplsda_res = oplsda_res,
    scores_df = scores_df,
    plot = p
  ))
}
