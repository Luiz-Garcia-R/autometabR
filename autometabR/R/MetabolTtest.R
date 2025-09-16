#' T-test or Mann-Whitney Test for a Metabolite
#'
#' Performs a t-test (or Mann-Whitney test) comparing a single metabolite between two groups.
#'
#' @param normalized_data List returned by \code{metabol.normalize()}, containing \code{$log2_normalized}.
#' @param metadata Data.frame containing at least columns 'Sample' and the grouping column.
#' @param metabolite Character. Name of the metabolite to test.
#' @param group_col Character. Column in \code{metadata} defining groups (default = "Group").
#' @param utest Logical. If TRUE, performs Mann-Whitney instead of t-test (default FALSE).
#' @param return_type Character. "htest" (returns only test object) or "all" (returns test, data, and plot). Default "htest".
#'
#' @return Invisibly returns a list or htest object depending on \code{return_type}.
#'
#' @examples
#' raw_data <- data.frame(Sample = paste0("S", 1:6),
#'                        Lactate = rnorm(6))
#' meta <- data.frame(Sample = paste0("S", 1:6), Group = c("A","A","B","B","A","B"))
#' imp_data <- list(log2_normalized = raw_data)
#'
#' metabol.ttest(imp_data, metadata = meta, metabolite = "Lactate", return_type = "all")
#'
#' @references
#' Student. (1908). The probable error of a mean. \emph{Biometrika}, 6(1): 1–25. <doi:10.2307/2331554>
#' Mann, H.B., & Whitney, D.R. (1947). On a test of whether one of two random variables is stochastically larger than the other. \emph{Annals of Mathematical Statistics}, 18: 50–60. <doi:10.1214/aoms/1177730491>
#'
#' @export

metabol.ttest <- function(normalized_data, metadata, metabolite,
                          group_col = "Group", utest = FALSE,
                          return_type = c("htest", "all")) {

  # --- Check packages ---
  pkgs <- c("ggplot2", "dplyr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  return_type <- match.arg(return_type)

  # --- Extract metabolite data ---
  data_mat <- normalized_data$log2_normalized
  if (!"Sample" %in% colnames(data_mat)) stop("normalized_data must have a 'Sample' column")

  rownames(data_mat) <- data_mat$Sample
  data_mat <- data_mat[, setdiff(colnames(data_mat), "Sample"), drop = FALSE]

  if (!metabolite %in% colnames(data_mat)) stop(paste0("Metabolite '", metabolite, "' not found"))

  # --- Merge with metadata ---
  df_long <- data.frame(Sample = rownames(data_mat), Value = data_mat[[metabolite]])
  df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

  # --- Check exactly 2 groups ---
  groups <- unique(df_long[[group_col]])
  if (length(groups) != 2) stop("Exactly 2 groups are required")

  # --- Extract numeric vectors ---
  g1 <- as.numeric(df_long$Value[df_long[[group_col]] == groups[1]])
  g2 <- as.numeric(df_long$Value[df_long[[group_col]] == groups[2]])
  g1 <- g1[is.finite(g1)]
  g2 <- g2[is.finite(g2)]

  if (length(g1) < 2 || length(g2) < 2) stop("Each group must have at least 2 non-NA values")

  # --- Check for many zeros ---
  if (mean(g1 == 0) > 0.5 || mean(g2 == 0) > 0.5) stop("More than 50% zeros in one group; test not recommended")

  # --- Statistical test ---
  if (utest) {
    test_name <- "Mann-Whitney"
    warning("Running Mann-Whitney on normalized data is not recommended")
    res <- suppressWarnings(stats::wilcox.test(g1, g2, exact = FALSE))
  } else {
    test_name <- "t-test"
    res <- stats::t.test(g1, g2)
  }

  # --- P-values & significance labels ---
  pval <- res$p.value
  signif_label <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"

  # --- Plot ---
  y_max <- max(df_long$Value, na.rm = TRUE)
  y_pos <- y_max * 1.02
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[group_col]], y = Value, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.65, color = "black") +
    ggplot2::annotate("text", x = 1.5, y = y_pos, label = signif_label, size = 6) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(title = metabolite, subtitle = paste0(test_name, " p-value: ", signif(pval, 3)),
                  x = "", y = "Normalized value")
  print(p)

  # --- Return invisibly ---
  if (return_type == "all") {
    return(invisible(list(test = res, data = df_long, plot = p)))
  } else {
    return(invisible(res))
  }
}
