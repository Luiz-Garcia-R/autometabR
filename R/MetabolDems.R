#' Identify Differentially Expressed Metabolites (DEMs)
#'
#' Computes log2 fold change and t-test p-values between two groups for each metabolite.
#' Flags metabolites as significant based on p-value and fold-change cutoff.
#'
#' @param normalized_data A \code{metabolNorm} object or a list containing \code{log2_normalized} data.frame/matrix.
#'                        Rows are samples, columns are metabolites, must include 'Sample' column.
#' @param metadata A data.frame containing at least a column with sample names and a group assignment.
#' @param group_col Character, default "Group". Column name in metadata representing group labels.
#' @param p_value Numeric, default 0.05. Threshold for significance.
#' @param fc_cutoff Numeric, default 0.5. Minimum absolute log2 fold-change for significance.
#' @param results Logical, default FALSE. If TRUE, returns full results table; if FALSE, returns invisible NULL.
#'
#' @return If \code{results = TRUE}, a data.frame with columns: \code{Metabolite}, \code{log2FC}, \code{p_value}, \code{Significant}.
#'         Otherwise, prints a brief summary and returns \code{invisible(NULL)}.
#'
#' @examples
#' # Minimal example for metabol.dems
#' expr_matrix <- data.frame(
#'   Sample = c("S1","S2","S3","S4"),
#'   Met1   = c(5,6,1,2),
#'   Met2   = c(2,3,2,3)
#' )
#'
#' norm_data <- list(expr_matrix = expr_matrix)
#' class(norm_data) <- "metabolNorm"
#'
#' meta <- data.frame(
#'   Sample = c("S1","S2","S3","S4"),
#'   Group  = c("A","A","B","B")
#' )
#'
#' metabol.dems(norm_data, metadata = meta)
#'
#' @references
#' Xia J, Wishart DS (2016). Using MetaboAnalyst 3.0 for Comprehensive
#' Metabolomics Data Analysis. \emph{Current Protocols in Bioinformatics},
#' 55: 14.10.1–14.10.91. <doi:10.1002/cpbi.11>
#'
#' Storey JD, Tibshirani R (2003). Statistical significance for genomewide
#' studies. \emph{Proceedings of the National Academy of Sciences},
#' 100(16):9440–9445. <doi:10.1073/pnas.1530509100>
#'
#' @export

metabol.dems <- function(normalized_data, metadata, group_col = "Group",
                         p_value = 0.05, fc_cutoff = 0.5, results = FALSE) {

  # --- Check packages ---
  if(!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")

  # --- Extract normalized data ---
  if ("metabolNorm" %in% class(normalized_data)) {
    data <- normalized_data$expr_matrix
  } else {
    stop("normalized_data must be an object of class 'metabolNorm'")
  }

  rownames(data) <- data$Sample
  data <- data[, setdiff(colnames(data), "Sample")]

  # --- Prepare groups ---
  df_long <- data.frame(Sample = rownames(data))
  df_long <- dplyr::left_join(df_long, metadata, by = "Sample")
  groups <- unique(df_long[[group_col]])
  if(length(groups) != 2) stop("Exactly 2 groups are required.")

  # --- Compute DEMs ---
  results_tbl <- data.frame(
    Metabolite = colnames(data),
    log2FC = NA,
    p_value = NA,
    Significant = NA
  )

  for(i in seq_along(colnames(data))) {
    g1 <- as.numeric(data[df_long[[group_col]] == groups[1], i])
    g2 <- as.numeric(data[df_long[[group_col]] == groups[2], i])

    g1 <- g1[is.finite(g1)]
    g2 <- g2[is.finite(g2)]

    log2fc <- log2(mean(g2) / mean(g1))
    t_res <- stats::t.test(g1, g2)

    results_tbl$log2FC[i] <- log2fc
    results_tbl$p_value[i] <- t_res$p.value
    results_tbl$Significant[i] <- (t_res$p.value < p_value) & (abs(log2fc) >= fc_cutoff)
  }

  # --- Output ---
  if (!results) {
    sig_count <- sum(results_tbl$Significant, na.rm = TRUE)
    msg <- paste0(
      "\n==============================\n",
      " Differential Expression Report\n",
      "==============================\n",
      "Groups compared: ", groups[1], " vs ", groups[2], "\n",
      "Significant metabolites: ", sig_count, "\n",
      "Cutoffs: p-value: ", p_value,
      " | log2FC: ", fc_cutoff, "\n",
      "==============================\n",
      "Tip: To access detailed results, use 'results = TRUE'.\n"
    )
    message(msg)
    invisible(NULL)
  } else {
    return(results_tbl)
  }
}
