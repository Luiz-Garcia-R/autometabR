#' ROC Curve for One or More Metabolites
#'
#' Computes a combined ROC curve and AUC for one or more metabolites,
#' using logistic regression when multiple metabolites are specified.
#'
#' @param normalized_data A list returned by \code{metabol.normalize()}, containing \code{$log2_normalized}.
#' @param metabolites Character vector. Name(s) of the metabolite(s) to evaluate.
#' @param metadata Data.frame containing at least 'Sample' and the grouping column.
#' @param group Character. Column in \code{metadata} specifying class labels. Default "Group".
#' @param plot Logical. Whether to plot the ROC curve. Default TRUE.
#'
#' @return Invisibly returns a list with:
#' \item{roc_object}{Object of class \code{roc} from the pROC package.}
#' \item{auc}{Area under the curve (AUC) value.}
#' \item{metabolites}{Metabolite(s) used.}
#' \item{summary}{Summary information about samples, group sizes, and direction.}
#'
#' @examples
#' # Minimal example
#' set.seed(123)
#' raw_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Creatine = c(rnorm(10, 5, 1), rnorm(10, 6, 1)),
#'   Glucose  = c(rnorm(10, 10, 2), rnorm(10, 10.5, 2))
#' )
#'
#' meta <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Group = rep(c("A", "B"), each = 10)
#' )
#'
#' imp_data <- list(expr_matrix = raw_data)
#'
#' res_roc <- metabol.roc(
#'   normalized_data = imp_data,
#'   metabolites = c("Creatine", "Glucose"),
#'   metadata = meta
#' )
#'
#' res_roc$auc
#'
#' @references
#' Hanley JA, McNeil BJ. (1982). The meaning and use of the area under a receiver operating characteristic (ROC) curve. \emph{Radiology}, 143(1): 29–36. <doi:10.1148/radiology.143.1.7063747>
#' Zou KH, O'Malley AJ, Mauri L. (2007). Receiver-operating characteristic analysis for evaluating diagnostic tests and predictive models. \emph{Circulation}, 115: 654–657. <doi:10.1161/CIRCULATIONAHA.105.594929>
#'
#' @export

metabol.roc <- function(normalized_data,
                        metadata = NULL,
                        metabolites,
                        group = "Group",
                        plot = TRUE) {

  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' required.")

  # ------------------------------------------------------------
  # 1) validate metabolites
  # ------------------------------------------------------------
  metabolites <- as.character(metabolites)

  # ------------------------------------------------------------
  # 2) Detect expr_matrix
  # ------------------------------------------------------------
  if (!is.list(normalized_data) || !"expr_matrix" %in% names(normalized_data))
    stop("normalized_data must be the list returned by metabol.normalize() and contain 'expr_matrix'.")

  data_mat <- normalized_data$expr_matrix

  # Remover coluna Sample se existir
  if ("Sample" %in% colnames(data_mat)) {
    rownames(data_mat) <- data_mat$Sample
    data_mat <- data_mat[, setdiff(colnames(data_mat), "Sample"), drop = FALSE]
  }

  # ------------------------------------------------------------
  # 3) Detect metadata
  # ------------------------------------------------------------
  if (is.null(metadata)) {
    if (!is.null(normalized_data$metadata)) {
      metadata <- normalized_data$metadata
      message("Using metadata stored inside normalized_data.")
    } else {
      stop("metadata not provided and none found inside normalized_data.")
    }
  }

  if (!"Sample" %in% colnames(metadata))
    stop("metadata must contain a column called 'Sample'.")

  if (!group %in% colnames(metadata))
    stop(paste0("metadata must contain column '", group, "'."))

  # ------------------------------------------------------------
  # 4) Align metadata <-> expr_matrix
  # ------------------------------------------------------------
  common_samples <- intersect(rownames(data_mat), metadata$Sample)
  if (length(common_samples) == 0)
    stop("No overlapping sample names between expr_matrix and metadata$Sample.")

  metadata <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]
  data_mat <- data_mat[common_samples, , drop = FALSE]

  # ------------------------------------------------------------
  # 5) Prepare model
  # ------------------------------------------------------------
  classe <- as.factor(metadata[[group]])
  group_levels <- levels(classe)

  if (length(metabolites) > 3) {
    message("Warning: Using more than 3 metabolites may inflate AUC due to overfitting.")
  }

  # Check metabolites
  missing_met <- setdiff(metabolites, colnames(data_mat))
  if (length(missing_met))
    stop("These metabolites were not found in the expression matrix: ",
         paste(missing_met, collapse = ", "))

  combined_data <- data.frame(classe, data_mat[, metabolites, drop = FALSE])

  # ------------------------------------------------------------
  # 6) Adjusts logistic model
  # ------------------------------------------------------------
  formula_str <- paste("classe ~", paste(metabolites, collapse = " + "))
  logit_model <- stats::glm(stats::as.formula(formula_str),
                            data = combined_data,
                            family = stats::binomial)

  pred <- stats::predict(logit_model, type = "response")

  # ------------------------------------------------------------
  # 7) ROC e AUC
  # ------------------------------------------------------------
  roc_comb <- pROC::roc(classe, pred, levels = group_levels, direction = "<")
  auc_comb <- pROC::auc(roc_comb)

  # ------------------------------------------------------------
  # 8) Plot
  # ------------------------------------------------------------
  if (plot) {
    plot(1 - roc_comb$specificities, roc_comb$sensitivities,
         type = "l", col = "#1f77b4", lwd = 3,
         main = "Combined ROC",
         xlab = "False Positive Rate",
         ylab = "True Positive Rate",
         xlim = c(0,1), ylim = c(0,1),
         bty = "n", las = 1)

    mtext(side = 3, line = 0.5,
          at = mean(par("usr")[1:2]),
          text = paste("Metabolites:", paste(metabolites, collapse = " + ")),
          cex = 1)

    abline(a = 0, b = 1, col = "#d3d3d3", lty = 2, lwd = 2)
    text(0.8, 0.05, paste("AUC =", round(auc_comb, 3)),
         col = "#1f77b4", cex = 1.2, font = 2)
    grid(col = "lightgray", lty = "dotted")
  }

  # ------------------------------------------------------------
  # 9) Return
  # ------------------------------------------------------------
  invisible(list(
    roc_object = roc_comb,
    auc = auc_comb,
    metabolites = metabolites,
    summary = list(
      metabolites = metabolites,
      n_samples = length(classe),
      group_sizes = table(classe),
      auc = auc_comb,
      direction = paste(group_levels[1], "<", group_levels[2])
    )
  ))
}
