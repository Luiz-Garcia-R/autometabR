#' Normalize Metabolomics Data
#'
#' Normalizes metabolomics data, applies filtering, outlier removal, log transformation,
#' global imputation, and optionally collapses highly correlated metabolites.
#'
#' @param imp_data A \code{metabolR} object containing raw data (\code{data}) and metadata (\code{metadata}).
#' @param method Normalization method: "log2" (default) or "log10".
#' @param remove_zeros Logical, whether to remove metabolites that are all zeros. Default FALSE.
#' @param remove_outliers Logical, whether to replace outliers with NA. Default TRUE.
#' @param help Logical, if TRUE shows detailed information about the normalization steps instead of running them. Default FALSE.
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
#' set.seed(123)
#'
#' # expression data
#' expr_data <- data.frame(
#'   Sample = paste0('S', 1:5),
#'   Met1 = c(5.2, 3.1, 4.8, 5.5, 4.2),
#'   Met2 = c(2.3, 2.7, 2.5, 2.6, 2.4),
#'   Met3 = c(8.1, 7.8, 8.5, 8.0, 8.3)
#' )
#'
#' # metadata
#' metadata <- data.frame(
#'   Sample = paste0('S', 1:5),
#'   Group = c('A','A','B','B','A')
#' )
#'
#' # Create metabolR object
#' imp_data <- list(data = expr_data, metadata = metadata)
#' class(imp_data) <- "metabolR"
#'
#' # Test normalization
#' norm_res <- metabol.normalize(imp_data)
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
                              clip_negative_log = TRUE,
                              collapse_correlated = TRUE, cor_threshold = 0.95,
                              assign_result = TRUE, assign_name = "normalized_data",
                              envir = parent.frame(),
                              help = FALSE) {

  if (help || missing(imp_data)) {
    message("Function metabol.normalize(): normalize metabolites in a 'metabolR' object.\nSee help for full description.")
    return(invisible(NULL))
  }

  if (!inherits(imp_data, "metabolR")) stop("Input must be a 'metabolR' object.")

  # Ensure data.frame
  raw_data <- as.data.frame(imp_data$data)
  metadata <- as.data.frame(imp_data$metadata)

  # Identify Sample and Group columns
  col_sample <- which(names(metadata) %in% c("Sample","sample","SAMPLE"))
  col_group  <- which(names(metadata) %in% c("Group","group","GROUP","Condition","condition","COND"))
  if (length(col_sample) != 1) stop("Sample column not detected in metadata.")
  if (length(col_group) != 1) stop("Group column not detected in metadata.")

  samples <- metadata[[col_sample]]
  groups  <- metadata[[col_group]]

  # --- Ensure sample matching ---
  sample_idx <- match(samples, raw_data[[names(raw_data)[col_sample]]])
  if (any(is.na(sample_idx))) stop("Some samples in metadata are not present in raw_data$Sample")

  # --- Extract numeric columns only ---
  numeric_cols <- setdiff(names(raw_data), names(raw_data)[col_sample])
  expr_data <- raw_data[sample_idx, numeric_cols, drop = FALSE]
  expr_data <- as.data.frame(expr_data) # ensure compatibility

  # --- Remove zero metabolites ---
  if (remove_zeros && ncol(expr_data) > 0) {
    nonzero_cols <- colSums(expr_data != 0, na.rm = TRUE) > 0
    expr_data <- expr_data[, nonzero_cols, drop = FALSE]
  }

  # --- Remove outliers ---
  if (remove_outliers && ncol(expr_data) > 0) {
    remove_outlier_col <- function(x) {
      qnt <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- qnt[2] - qnt[1]
      x[x < qnt[1] - 3*iqr | x > qnt[2] + 3*iqr] <- NA
      x
    }
    expr_data <- as.data.frame(lapply(expr_data, remove_outlier_col))
  }

  # --- Filter metabolites with too many NA or zeros ---
  if (ncol(expr_data) > 0) {
    remove_idx <- sapply(seq_len(ncol(expr_data)), function(i) {
      by_group <- tapply(expr_data[[i]], groups, function(vals) mean(is.na(vals) | vals == 0))
      all(by_group > max_missing_prop_per_group)
    })
    expr_data <- expr_data[, !remove_idx, drop = FALSE]
  }

  # --- Initialize normalized matrix ---
  if (ncol(expr_data) == 0) {
    norm_matrix <- data.frame(matrix(nrow = nrow(expr_data), ncol = 0))
  } else {

    # Pseudo-count
    if (is.null(pseudo_count)) {
      pseudo_count <- sapply(expr_data, function(x) {
        vals <- x[x > 0 & is.finite(x)]
        if (length(vals) == 0) return(1e-6)
        min(vals)/2
      })
    }

    # Logarithmic normalization
    norm_matrix <- expr_data
    if (method %in% c("log2","log10")) {
      for (i in seq_len(ncol(norm_matrix))) {
        norm_matrix[[i]] <- switch(method,
                                   log2 = log2(norm_matrix[[i]] + pseudo_count[i]),
                                   log10 = log10(norm_matrix[[i]] + pseudo_count[i]))
      }
    }

    # Global imputation
    if (impute) {
      observed <- unlist(norm_matrix)
      observed <- observed[observed > 0 & is.finite(observed)]
      imp_mean <- median(observed, na.rm = TRUE) - impute_shift * sd(observed, na.rm = TRUE)
      imp_sd   <- sd(observed, na.rm = TRUE) * impute_scale
      miss_idx <- which(norm_matrix <= 0 | is.na(norm_matrix), arr.ind = TRUE)
      if (nrow(miss_idx) > 0) norm_matrix[miss_idx] <- stats::rnorm(nrow(miss_idx), mean = imp_mean, sd = imp_sd)
    }

    # Clip negative values
    if (clip_negative_log) {
      norm_matrix[norm_matrix < 0] <- 0
    }

    # Collapse highly correlated metabolites
    if (collapse_correlated && ncol(norm_matrix) > 1) {
      cor_mat <- stats::cor(norm_matrix, use = "pairwise.complete.obs")
      high_cor <- which(abs(cor_mat) >= cor_threshold & row(cor_mat) < col(cor_mat), arr.ind = TRUE)
      if (nrow(high_cor) > 0) {
        for (idx in seq_len(nrow(high_cor))) {
          c1 <- high_cor[idx, 1]
          c2 <- high_cor[idx, 2]
          norm_matrix[[c1]] <- rowMeans(norm_matrix[, c(c1, c2)], na.rm = TRUE)
          norm_matrix <- norm_matrix[, -c2, drop = FALSE]
        }
      }
    }
  }

  # --- Build final object ---
  normalized_df <- cbind(Sample = raw_data[[names(raw_data)[col_sample]]][sample_idx], norm_matrix)

  obj <- list(
    expr_matrix = normalized_df,
    raw_preserved = raw_data,
    method = method
  )
  class(obj) <- "metabolNorm"

  if (assign_result) {
    assign(assign_name, obj, envir = envir)
    print(obj)
    invisible(obj)
  } else {
    return(obj)
  }
}

# --- Safe print method ---
print.metabolNorm <- function(x, ...) {
  message("Object of class 'metabolNorm'")
  message("------------------------------")
  message("Samples: ", nrow(x$expr_matrix))
  message("Metabolites: ", ncol(x$expr_matrix) - 1)
  message("Normalization method: ", x$method)
  message("------------------------------")
  invisible(x)
}
