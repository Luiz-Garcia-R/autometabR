#' Comparative Statistical Test for a Metabolite (t-test, ANOVA, Wilcoxon, or Kruskal-Wallis)
#'
#' Performs an automatic comparison of a single metabolite across two or more groups,
#' selecting the appropriate statistical test based on data normality.
#' For two groups, it performs a t-test or Wilcoxon test.
#' For more than two groups, it performs ANOVA or Kruskal-Wallis followed by post-hoc tests (Tukey or Dunn).
#' Significant differences are represented by grouping letters in the plot.
#'
#' @param normalized_data List returned by \code{metabol.normalize()}, containing \code{$expr_matrix}.
#' @param metadata Data.frame containing at least columns 'Sample' and the grouping column.
#' @param metabolite Character. Name of the metabolite to test.
#' @param group_col Character. Column in \code{metadata} defining groups (default = "Group").
#' @param adjust Character. P-value adjustment method for non-parametric post-hoc tests (default = "BH").
#' @param verbose Logical. Whether to print group means, SDs, and test summaries (default = TRUE).
#' @param theme_style ggplot2 theme applied to the plot (default = \code{ggplot2::theme_minimal()}).
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{model}: test or ANOVA model object.
#'   \item \code{letters}: data.frame with group letters indicating significant differences.
#'   \item \code{summary}: data.frame with group means and standard deviations.
#'   \item \code{plot}: ggplot object showing the boxplot with significance letters.
#' }
#'
#' @examples
#' # Example dataset
#' set.seed(1)
#' expr <- data.frame(
#'   Glucose = c(rnorm(3, 5, 1), rnorm(3, 6, 1), rnorm(3, 8, 1))
#' )
#' rownames(expr) <- paste0("S", 1:9)  # força os nomes de amostra
#' meta <- data.frame(
#'   Sample = paste0("S", 1:9),
#'   Group = rep(c("A", "B", "C"), each = 3)
#' )
#'
#' norm_data <- list(expr_matrix = expr)
#' class(norm_data) <- "normalized_data"
#'
#' metabol.anova(norm_data, metadata = meta, metabolite = "Glucose")
#'
#' @references
#' Student. (1908). The probable error of a mean. \emph{Biometrika}, 6(1): 1–25. <doi:10.2307/2331554>
#' Kruskal, W.H., & Wallis, W.A. (1952). Use of ranks in one-criterion variance analysis. \emph{Journal of the American Statistical Association}, 47(260): 583–621. <doi:10.1080/01621459.1952.10483441>
#' Dunn, O.J. (1964). Multiple comparisons using rank sums. \emph{Technometrics}, 6(3): 241–252. <doi:10.1080/00401706.1964.10490181>
#'
#' @export
metabol.anova <- function(normalized_data, metadata, metabolite,
                            group_col = "Group",
                            adjust = "BH",
                            verbose = TRUE,
                            theme_style = ggplot2::theme_minimal()) {

  # --- Check classes ---
  if (!inherits(normalized_data, "metabolNorm") && !inherits(normalized_data, "normalized_data"))
    stop("Input object must be of class 'metabolNorm' or 'normalized_data'.")

  if (!is.data.frame(metadata))
    stop("'metadata' must be a data.frame.")

  data <- normalized_data$expr_matrix

  # --- Check grouping column ---
  if (!group_col %in% colnames(metadata))
    stop("The grouping column '", group_col, "' was not found in 'metadata'.")

  # --- Check metabolite ---
  if (!metabolite %in% colnames(data))
    stop("Specified metabolite ('", metabolite, "') not found in normalized data.")

  # --- Prepare data frame ---
  df <- data.frame(
    Sample = rownames(data),
    value = as.numeric(data[, metabolite])
  )
  df <- merge(df, metadata[, c("Sample", group_col)], by = "Sample", all.x = TRUE)
  colnames(df)[which(colnames(df) == group_col)] <- "group"
  df$group <- factor(df$group)

  # --- Normality check ---
  shapiro_p <- tapply(df$value, df$group, function(x) {
    if(length(x) < 3) return(NA)
    tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
  })

  is_normal <- all(shapiro_p > 0.05, na.rm = TRUE)

  groups <- unique(df$group)

  # --- Main test and post-hoc ---
  if (length(groups) == 2) {
    if (is_normal) {
      model <- stats::t.test(value ~ group, data = df)
      p_global <- model$p.value
      p_label <- if (p_global < 0.001) "t-test: p < 0.001" else paste0("t-test: p = ", signif(p_global, 3))
      letters_df <- data.frame(group = groups, letter = c("a", "b"))
    } else {
      model <- stats::wilcox.test(value ~ group, data = df)
      p_global <- model$p.value
      p_label <- if (p_global < 0.001) "Wilcoxon: p < 0.001" else paste0("Wilcoxon: p = ", signif(p_global, 3))
      letters_df <- data.frame(group = groups, letter = c("a", "b"))
    }

  } else {
    if (is_normal) {
      model <- stats::aov(value ~ group, data = df)
      p_global <- summary(model)[[1]][["Pr(>F)"]][1]
      p_label <- if (p_global < 0.001) "ANOVA: p < 0.001" else paste0("ANOVA: p = ", signif(p_global, 3))
      tukey_res <- stats::TukeyHSD(model)
      letters <- multcompView::multcompLetters4(model, tukey_res)
      letters_df <- data.frame(group = names(letters$group$Letters),
                               letter = letters$group$Letters)
    } else {
      model <- rstatix::kruskal_test(df, value ~ group)
      p_global <- model$p
      p_label <- if (p_global < 0.001) "Kruskal-Wallis: p < 0.001" else paste0("Kruskal-Wallis: p = ", signif(p_global, 3))
      posthoc <- rstatix::dunn_test(df, value ~ group, p.adjust.method = adjust)
      pvals_posthoc <- setNames(posthoc$p.adj, paste(posthoc$group1, posthoc$group2, sep = "-"))
      letters_raw <- multcompView::multcompLetters(pvals_posthoc)$Letters
      letters_df <- data.frame(group = names(letters_raw), letter = letters_raw)
    }
  }

  # --- Mean and SD summary ---
  mean_sd <- stats::aggregate(value ~ group, data = df,
                              function(x) c(mean = mean(x), sd = sd(x)))
  mean_sd <- do.call(data.frame, mean_sd)
  colnames(mean_sd)[2:3] <- c("mean", "sd")

  if (verbose) {
    message("Group means and SDs:")
    print(mean_sd)
  }

  # --- Letter placement ---
  max_vals <- stats::aggregate(value ~ group, data = df, max)
  letters_df <- merge(max_vals, letters_df, by = "group")
  letters_df$value <- letters_df$value + 0.05 * diff(range(df$value))

  # --- Plot ---
  g <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.6, color = "black") +
    ggplot2::geom_text(data = letters_df,
                       ggplot2::aes(x = group, y = value, label = letter),
                       size = 4.5, vjust = 0) +
    theme_style +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = metabolite,
      subtitle = p_label,
      x = "Group",
      y = "Normalized value"
    )

  print(g)

  invisible(list(
    model = model,
    letters = letters_df,
    summary = mean_sd,
    plot = g
  ))
}
