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
#' imp_data <- list(log2_normalized = raw_data)
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

metabol.roc <- function(normalized_data, metabolites, metadata, group = "Group", plot = TRUE) {

  # --- Check required package ---
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' required. Install with install.packages('pROC').")

  metabolites <- as.character(metabolites)

  # --- Extract data and metadata ---
  if (!"log2_normalized" %in% names(normalized_data)) stop("normalized_data must contain 'log2_normalized'.")
  data_mat <- normalized_data$log2_normalized

  if (!group %in% colnames(metadata)) stop(paste("The column", group, "was not found in metadata"))
  classe <- as.factor(metadata[[group]])
  group_levels <- levels(classe)

  # --- Warn for too many metabolites ---
  if (length(metabolites) > 3) {
    message("Warning: Using more than 3 metabolites may lead to overfitting and artificially high AUC values.")
  }

  # --- Subset metabolites ---
  combined_data <- data_mat[, metabolites, drop = FALSE]

  # --- Fit logistic regression ---
  logit_model <- stats::glm(classe ~ ., data = combined_data, family = stats::binomial)
  pred <- stats::predict(logit_model, type = "response")

  # --- Compute ROC ---
  roc_comb <- pROC::roc(classe, pred, levels = group_levels, direction = "<")
  auc_comb <- pROC::auc(roc_comb)

  # --- Plot ROC if requested ---
  if (plot) {
    plot(1 - roc_comb$specificities, roc_comb$sensitivities,
         type = "l", col = "#1f77b4", lwd = 3,
         main = "Combined ROC",
         xlab = "False Positive Rate",
         ylab = "True Positive Rate",
         xlim = c(0,1), ylim = c(0,1),
         bty = "n", las = 1)
    mtext(side = 3, line = 0.5, at = mean(par("usr")[1:2]),
          text = paste("Metabolites:", paste(metabolites, collapse = " + ")), cex = 1)
    abline(a = 0, b = 1, col = "#d3d3d3", lty = 2, lwd = 2)
    text(0.8, 0.05, paste("AUC =", round(auc_comb, 3)), col = "#1f77b4", cex = 1.2, font = 2)
    grid(col = "lightgray", lty = "dotted")
  }

  # --- Prepare summary info instead of printing ---
  summary_info <- list(
    metabolites = metabolites,
    n_samples = length(classe),
    group_sizes = table(classe),
    auc = auc_comb,
    direction = paste(group_levels[1], "<", group_levels[2])
  )

  # --- Return invisibly ---
  invisible(list(
    roc_object = roc_comb,
    auc = auc_comb,
    metabolites = metabolites,
    summary = summary_info
  ))
}
