#' Normalize Metabolomics Data
#'
#' Normalizes metabolomics data, applies filtering, outlier removal, log transformation,
#' global imputation, and optionally collapses highly correlated metabolites.
#'
#' @param imp_data A \code{metabolR} object containing raw data (\code{data}) and metadata (\code{metadata}).
#' @param method Normalization method: "log2" (default) or "log10".
#' @param remove_zeros Logical, whether to remove metabolites that are all zeros. Default FALSE.
#' @param remove_outliers Logical, whether to replace outliers with NA. Default TRUE.
#' @param impute Logical, whether to perform global imputation. Default TRUE.
#' @param impute_shift Shift applied to median for imputation. Default 1.8.
#' @param impute_scale Scale applied to SD for imputation. Default 0.3.
#' @param max_missing_prop_per_group Maximum allowed missing proportion per group. Default 0.5.
#' @param pseudo_count Numeric vector or NULL, pseudo-count per metabolite. Default NULL (computed automatically).
#' @param clip_negative_log Logical, whether to clip negative values to zero before log transformation. Default TRUE.
#' @param collapse_correlated Logical, whether to collapse highly correlated metabolites. Default TRUE.
#' @param cor_threshold Correlation threshold for collapsing metabolites. Default 0.95.
#' @param assign_result Logical, whether to assign the resulting object to the environment. Default TRUE.
#' @param assign_name Character, name of the object to assign if \code{assign_result = TRUE}. Default "normalized_data".
#' @param envir Environment where to assign the result if \code{assign_result = TRUE}. Default parent.frame().
#'
#' @return Invisibly returns a \code{metabolNorm} object with the normalized data, presence info, checks, and summary.
#'
#' @examples
#' # Minimal example
#' raw_data <- data.frame(Sample = paste0("S", 1:5), Met1 = rnorm(5), Met2 = rnorm(5))
#' meta <- data.frame(Sample = paste0("S", 1:5), Group = c("A","A","B","B","B"))
#' imp_data <- list(data = raw_data, metadata = meta)
#'
#' class(imp_data) <- "metabolR"
#'
#' norm_obj <- metabol.normalize(imp_data)
#' print(norm_obj)
#'
#' @references
#' Xia J, Psychogios N, Young N, Wishart DS. (2009). MetaboAnalyst: a web server for metabolomic data analysis and interpretation. \emph{Nucleic Acids Research}, 37(Web Server issue): W652–W660. <doi:10.1093/nar/gkp356>
#' van den Berg RA, Hoefsloot HCJ, Westerhuis JA, Smilde AK, van der Werf MJ. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. \emph{BMC Genomics}, 7:142. <doi:10.1186/1471-2164-7-142>
#'
#' @export

metabol.normalize <- function(imp_data, method = "log2",
                              remove_zeros = FALSE, remove_outliers = TRUE,
                              impute = TRUE, impute_shift = 1.8, impute_scale = 0.3,
                              max_missing_prop_per_group = 0.5, pseudo_count = NULL,
                              clip_negative_log = TRUE, collapse_correlated = TRUE,
                              cor_threshold = 0.95, assign_result = TRUE, assign_name = "normalized_data",
                              envir = parent.frame()) {

  checks <- list()

  # --- Input validation ---
  if (!inherits(imp_data, "metabolR")) stop("Input must be a 'metabolR' object.")
  raw_data <- imp_data$data
  metadata <- imp_data$metadata

  col_sample <- which(names(metadata) %in% c("Sample","sample","SAMPLE"))
  col_group <- which(names(metadata) %in% c("Group","group","GROUP","Condition","condition","COND"))
  if (length(col_sample) != 1) stop("Sample column not detected in metadata.")
  if (length(col_group) != 1) stop("Group column not detected in metadata.")

  samples <- metadata[[col_sample]]
  groups <- metadata[[col_group]]

  sample_idx <- match(samples, raw_data$Sample)
  if (any(is.na(sample_idx))) stop("Some samples in metadata are not present in raw_data$Sample")
  expr_data <- as.data.frame(raw_data[sample_idx, -1, drop = FALSE])
  rownames(expr_data) <- raw_data$Sample[sample_idx]

  # --- Presence per group ---
  presence_per_group <- t(apply(expr_data, 2, function(x)
    tapply(x, groups, function(vals) as.integer(any(vals > 0)))
  ))
  presence_per_group <- as.data.frame(presence_per_group)
  presence_per_group$Metabolite <- colnames(expr_data)
  presence_per_group <- presence_per_group[, c("Metabolite", setdiff(names(presence_per_group), "Metabolite"))]

  # --- Remove zero metabolites ---
  if (remove_zeros) {
    nonzero_cols <- colSums(expr_data != 0, na.rm = TRUE) > 0
    expr_data <- expr_data[, nonzero_cols, drop = FALSE]
    checks$removed_all_zeros <- sum(!nonzero_cols)
  }

  # --- Remove outliers ---
  if (remove_outliers) {
    remove_outlier_col <- function(x) {
      qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- qnt[2] - qnt[1]
      x[x < qnt[1] - 3*iqr | x > qnt[2] + 3*iqr] <- NA
      x
    }
    expr_data <- apply(expr_data, 2, remove_outlier_col)
    checks$outliers <- TRUE
  }

  # --- Filter metabolites with too many missing ---
  remove_idx <- sapply(1:ncol(expr_data), function(i) {
    by_group <- tapply(expr_data[, i], groups, function(vals) mean(is.na(vals) | vals == 0))
    all(by_group > max_missing_prop_per_group)
  })
  expr_data <- expr_data[, !remove_idx, drop = FALSE]
  presence_per_group <- presence_per_group[!remove_idx, , drop = FALSE]
  checks$missing_filter <- sum(remove_idx)

  # --- Clip negatives ---
  if (clip_negative_log && ncol(expr_data) > 0) {
    n_neg <- sum(expr_data < 0, na.rm = TRUE)
    if (n_neg > 0) expr_data[expr_data < 0] <- 0
    checks$clip_neg <- n_neg
  }

  # --- Pseudo-count ---
  if (is.null(pseudo_count)) {
    pseudo_count <- apply(expr_data, 2, function(x) {
      vals <- x[x > 0 & is.finite(x)]
      if (length(vals) == 0) return(1e-6)
      min(vals)/2
    })
  }

  # --- Normalization ---
  norm_matrix <- expr_data
  if (method %in% c("log2", "log10")) {
    for (i in 1:ncol(norm_matrix)) {
      norm_matrix[, i] <- switch(method,
                                 log2 = log2(norm_matrix[, i] + pseudo_count[i]),
                                 log10 = log10(norm_matrix[, i] + pseudo_count[i]))
    }
  }

  # --- Imputation ---
  if (impute) {
    observed <- norm_matrix[norm_matrix > 0 & is.finite(norm_matrix)]
    imp_mean <- median(observed, na.rm = TRUE) - impute_shift * sd(observed, na.rm = TRUE)
    imp_sd <- sd(observed, na.rm = TRUE) * impute_scale
    miss_idx <- which(norm_matrix <= 0 | is.na(norm_matrix))
    if (length(miss_idx) > 0) norm_matrix[miss_idx] <- rnorm(length(miss_idx), mean = imp_mean, sd = imp_sd)
    checks$imputation <- c(mean = imp_mean, sd = imp_sd)
  }

  # --- Collapse correlated metabolites ---
  if (collapse_correlated && ncol(norm_matrix) > 1) {
    cor_mat <- cor(norm_matrix, use = "pairwise.complete.obs")
    high_cor <- which(abs(cor_mat) >= cor_threshold & row(cor_mat) < col(cor_mat), arr.ind = TRUE)
    collapsed <- 0
    if (nrow(high_cor) > 0) {
      for (idx in 1:nrow(high_cor)) {
        c1 <- high_cor[idx, 1]
        c2 <- high_cor[idx, 2]
        norm_matrix[, c1] <- rowMeans(norm_matrix[, c(c1, c2)], na.rm = TRUE)
        norm_matrix <- norm_matrix[, -c2, drop = FALSE]
        collapsed <- collapsed + 1
      }
    }
    checks$collapsed <- collapsed
  }

  # --- Final normalized object ---
  norm_name <- paste0(method, "_normalized")
  normalized_df <- cbind(Sample = rownames(norm_matrix), as.data.frame(norm_matrix))
  expr_matrix <- as.matrix(norm_matrix)
  rownames(expr_matrix) <- rownames(norm_matrix)

  obj <- list()
  obj[[norm_name]] <- normalized_df
  obj$raw_preserved <- raw_data
  obj$expr_matrix <- expr_matrix
  obj$presence_per_group <- presence_per_group
  obj$method <- method
  obj$checks <- checks
  class(obj) <- "metabolNorm"

  if (assign_result) assign(assign_name, obj, envir = envir)
  invisible(obj)
}

#' @export
print.metabolNorm <- function(x, ...) {
  summary_list <- list(
    method = x$method,
    n_metabolites = ncol(x[[grep("_normalized$", names(x))]]) - 1,
    n_samples = nrow(x[[grep("_normalized$", names(x))]]),
    checks = x$checks
  )
  invisible(summary_list)
}
