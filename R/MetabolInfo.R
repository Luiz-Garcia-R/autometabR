#' Metabolomics Dataset Summary
#'
#' Provides a summary of a metabolomics dataset including sample numbers, metabolite counts,
#' missing values, concentration statistics, group sizes, and detection rates per metabolite.
#'
#' @param normalized_data A list returned by \code{metabol.normalize()} containing \code{expr_matrix}.
#' @param metadata A data.frame with at least columns \code{Sample} and a grouping column.
#' @param group_col Character, name of the grouping column in \code{metadata}. Default "Group".
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{n_samples}{Number of samples.}
#'   \item{n_metabolites}{Number of metabolites.}
#'   \item{groups}{Unique groups found in metadata.}
#'   \item{missing_total}{Total number of missing values.}
#'   \item{missing_perc}{Percentage of missing values.}
#'   \item{mean_concentration}{Mean concentration of all metabolites.}
#'   \item{sd_concentration}{Standard deviation of metabolite concentrations.}
#'   \item{range_concentration}{Minimum and maximum concentrations.}
#'   \item{group_summary}{Data frame with number of samples per group.}
#'   \item{metabolite_summary}{Data frame with metabolite detection counts and rates.}
#' }
#'
#' @examples
#' # Minimal example
#' expr_mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' rownames(expr_mat) <- paste0("S", 1:5)
#' colnames(expr_mat) <- paste0("Met", 1:4)
#' norm_data <- list(expr_matrix = expr_mat)
#' meta <- data.frame(Sample = paste0("S", 1:5), Group = c("A","A","B","B","B"))
#'
#' info <- metabol.info(norm_data, meta)
#'
#' print(info)
#'
#' @references
#' Want EJ, Masson P, Michopoulos F, Wilson ID, Theodoridis G, Plumb RS, Shockcor J, Loftus N, Holmes E, Nicholson JK. (2013). Global metabolic profiling of animal and human tissues via UPLC-MS. \emph{Nature Protocols}, 8(1): 17–32. <doi:10.1038/nprot.2012.135>
#' Xia J, Psychogios N, Young N, Wishart DS. (2009). MetaboAnalyst: a web server for metabolomic data analysis and interpretation. \emph{Nucleic Acids Research}, 37(Web Server issue): W652–W660. <doi:10.1093/nar/gkp356>
#'
#' @export

metabol.info <- function(normalized_data, metadata, group_col = "Group") {
  # ---- check input ----
  if (!is.list(normalized_data) || is.null(normalized_data$expr_matrix)) {
    stop("Input must be a list returned by metabol.normalize() containing 'expr_matrix'.")
  }

  data_mat <- normalized_data$expr_matrix

  # --- Remove Sample column se existir ---
  if ("Sample" %in% colnames(data_mat)) {
    rownames(data_mat) <- data_mat$Sample
    data_mat <- data_mat[, setdiff(colnames(data_mat), "Sample"), drop = FALSE]
  }

  # --- Ensure numeric ---
  data_mat <- as.matrix(data_mat)
  storage.mode(data_mat) <- "numeric"

  if (!"Sample" %in% colnames(metadata)) stop("metadata must contain the column 'Sample'.")
  if (!group_col %in% colnames(metadata)) stop(sprintf("metadata must contain the column '%s'.", group_col))

  # ---- align samples ----
  common_samples <- intersect(rownames(data_mat), metadata$Sample)
  data_mat <- data_mat[common_samples, , drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]

  # ---- basic info ----
  n_samples <- nrow(data_mat)
  n_metabolites <- ncol(data_mat)
  groups <- unique(metadata[[group_col]])

  # ---- missing values ----
  missing_total <- sum(is.na(data_mat))
  missing_perc <- 100 * missing_total / (n_samples * n_metabolites)

  # ---- concentration summaries ----
  mean_conc <- mean(data_mat, na.rm = TRUE)
  sd_conc <- sd(data_mat, na.rm = TRUE)
  range_conc <- range(data_mat, na.rm = TRUE)

  # ---- group-level summaries ----
  group_summary <- aggregate(metadata$Sample, list(metadata[[group_col]]), length)
  colnames(group_summary) <- c(group_col, "N_samples")

  # ---- metabolite-level summaries ----
  metab_nonzero <- colSums(data_mat > 0, na.rm = TRUE)
  metab_detect_rate <- metab_nonzero / n_samples
  metab_summary <- data.frame(
    Metabolite = colnames(data_mat),
    Detected_in_samples = metab_nonzero,
    Detection_rate = round(100 * metab_detect_rate, 1)
  )

  # ---- create metabolInfo object ----
  out <- list(
    n_samples = n_samples,
    n_metabolites = n_metabolites,
    groups = groups,
    missing_total = missing_total,
    missing_perc = missing_perc,
    mean_concentration = mean_conc,
    sd_concentration = sd_conc,
    range_concentration = range_conc,
    group_summary = group_summary,
    metabolite_summary = metab_summary
  )
  class(out) <- "metabolInfo"
  return(out)
}

# --- custom print method ---
#' @method print metabolInfo
#' @export
print.metabolInfo <- function(x, ...) {
  message("===============================")
  message("Object of class 'metabolInfo'")
  message("===============================")
  message("Samples:     ", x$n_samples)
  message("Metabolites: ", x$n_metabolites)
  message("Groups:      ", paste(x$groups, collapse = ", "))
  message("Missing values: ", x$missing_total, " (", round(x$missing_perc,1), "%)")
  message("Mean concentration: ", round(x$mean_concentration, 3))
  message("SD concentration:   ", round(x$sd_concentration, 3))
  message("Range:             [", round(x$range_concentration[1],3), ", ", round(x$range_concentration[2],3), "]")
  message("===============================")
  message("Tip: assign to a variable and access details with $group_summary or $metabolite_summary")
  invisible(x)
}


