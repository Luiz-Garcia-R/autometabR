#' VIP Scores and Heatmap for Dimensionality Reduction Results
#'
#' Computes Variable Importance in Projection (VIP) scores from a PLS-DA
#' object contained in a dimensionality reduction result, plots top VIPs and
#' a heatmap across groups, and optionally calculates AUC for each variable.
#'
#' @param dimred_obj List returned by \code{metabol.dimred()} containing \code{expr_matrix}, \code{metadata}, and \code{res_dimred$plsda_res}.
#' @param n_top Number of top VIP variables to show (default = 20).
#' @param component Component to extract VIP scores from (default = 1).
#' @param plot Logical, whether to plot VIP and heatmap (default = TRUE).
#' @param calc_auc Logical, whether to calculate AUC for each VIP variable (default = FALSE).
#' @param group_col Name of column in metadata defining groups (default = "Group").
#'
#' @return Invisibly returns a list containing VIP scores, top VIPs, AUC (if requested), and plot objects.
#'
#' @examples
#' \dontrun{
#' # Minimal example
#' set.seed(123)
#'
#' expr_mat <- data.frame(
#'   Met1 = rnorm(6),
#'   Met2 = rnorm(6),
#'   Met3 = rnorm(6)
#' )
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:6),
#'   Group  = rep(c("A","B"), each = 3)
#' )
#'
#' plsda_res <- mixOmics::plsda(expr_mat, metadata$Group, ncomp = 2)
#'
#' dimred_obj <- list(
#'   expr_matrix = expr_mat,
#'   metadata    = metadata,
#'   res_dimred  = list(plsda_res = plsda_res)
#' )
#'
#' metabol.vip(dimred_obj, n_top = 3, component = 1, plot = TRUE)
#' }
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics. \emph{Chemometrics and Intelligent Laboratory Systems}, 58(2), 109–130. <doi:10.1016/S0169-7439(01)00155-1>
#' Chong, I.G., & Jun, C.H. (2005). Performance of some variable selection methods when multicollinearity is present. \emph{Chemometrics and Intelligent Laboratory Systems}, 78(1–2), 103–112. <doi:10.1016/j.chemolab.2005.04.003>
#'
#' @export

metabol.vip <- function(dimred_obj, n_top = 20, component = 1, plot = TRUE,
                        calc_auc = FALSE, group_col = "Group") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'")
  if (calc_auc && !requireNamespace("pROC", quietly = TRUE)) stop("Please install 'pROC'")
  if (plot && !requireNamespace("patchwork", quietly = TRUE)) stop("Please install 'patchwork'")

  # --- Check input ---
  if (!is.list(dimred_obj) || !all(c("expr_matrix","metadata","res_dimred") %in% names(dimred_obj))) {
    stop("Input must be an object returned by metabol.dimred()")
  }

  expr_mat <- dimred_obj$expr_matrix
  metadata <- dimred_obj$metadata
  plsda_res <- dimred_obj$res_dimred$plsda_res

  # --- VIP scores ---
  vip_scores <- mixOmics::vip(plsda_res)
  vip_comp <- vip_scores[, component]
  top_vip <- sort(vip_comp, decreasing = TRUE)[1:n_top]

  df_vip <- data.frame(Variable = names(top_vip), VIP = top_vip)

  # --- Heatmap preparation ---
  cluster_groups <- as.factor(metadata[[group_col]])
  medias_por_cluster <- aggregate(
    expr_mat[, df_vip$Variable, drop = FALSE],
    by = list(Cluster = cluster_groups),
    FUN = mean
  )

  heatmap_long <- tidyr::pivot_longer(
    medias_por_cluster,
    cols = -Cluster,
    names_to = "Variable",
    values_to = "ValorMedio"
  )

  heatmap_long <- heatmap_long %>%
    dplyr::group_by(Variable) %>%
    dplyr::mutate(ValorZ = scale(ValorMedio)[,1]) %>%
    dplyr::ungroup()

  df_vip$Variable <- factor(df_vip$Variable, levels = rev(df_vip$Variable))
  heatmap_long$Variable <- factor(heatmap_long$Variable, levels = levels(df_vip$Variable))

  # --- Plots ---
  if (plot) {
    plot_vip <- ggplot2::ggplot(df_vip, ggplot2::aes(x = Variable, y = VIP)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(title = paste("VIP Scores - Component", component), y = "VIP Score", x = NULL) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

    plot_heatmap <- ggplot2::ggplot(heatmap_long, ggplot2::aes(x = as.factor(Cluster), y = Variable, fill = ValorZ)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Z-score") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(x = "", y = NULL, title = "") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )

    print(plot_vip + plot_heatmap + patchwork::plot_layout(widths = c(2.5, 1)))
  }

  # --- Optional AUC calculation ---
  resultados_auc <- NULL
  if (calc_auc) {
    resultados_auc <- data.frame(Metabolito = character(), AUC = numeric())
    classe <- as.factor(metadata[[group_col]])

    for (metab in df_vip$Variable) {
      roc_obj <- pROC::roc(classe, expr_mat[[metab]], levels = levels(classe), direction = "<")
      auc_val <- as.numeric(pROC::auc(roc_obj))
      resultados_auc <- rbind(resultados_auc, data.frame(Metabolito = metab, AUC = auc_val))
    }
    resultados_auc <- resultados_auc[order(-resultados_auc$AUC), ]
  }

  invisible(list(
    vip_scores = vip_scores,
    top_vip = df_vip,
    auc = resultados_auc,
    plot_vip = if(plot) plot_vip else NULL,
    plot_heatmap = if(plot) plot_heatmap else NULL
  ))
}
