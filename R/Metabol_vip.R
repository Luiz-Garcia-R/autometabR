#' VIP Scores and Heatmap for Dimensionality Reduction Results
#'
#' Computes Variable Importance in Projection (VIP) scores from a PLS-DA
#' object contained in a dimensionality reduction result, plots top VIPs and
#' a heatmap across groups, and optionally calculates AUC for each variable.
#'
#' @param normalized_data Object from metabol.dimred().
#' @param top_n Integer, number of top metabolites to display. Default = 20.
#' @param component Integer, which component to extract VIP scores from. Default = 1.
#' @param plot Logical, whether to generate a bar plot. Default = TRUE.
#' @param calc_auc Logical, whether to compute AUC for classification performance. Default = FALSE.
#' @param group_col Character, metadata column indicating sample groups. Default = "Group".
#' @param assign_result Logical, whether to assign the result to the environment. Default = TRUE.
#' @param assign_name Character, name of the object if assigned. Default = "vip_data".
#' @param envir Environment where the object will be assigned. Default = parent.frame().
#'
#' @return A data frame of VIP scores, optionally plotted or assigned to the environment.
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
#' dimred_data <- list(
#'   expr_matrix = expr_mat,
#'   metadata    = metadata,
#'   res_dimred  = list(plsda_res = plsda_res)
#' )
#'
#' metabol.vip(dimred_data, n_top = 3, component = 1, plot = TRUE)
#' }
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics. \emph{Chemometrics and Intelligent Laboratory Systems}, 58(2), 109–130. <doi:10.1016/S0169-7439(01)00155-1>
#' Chong, I.G., & Jun, C.H. (2005). Performance of some variable selection methods when multicollinearity is present. \emph{Chemometrics and Intelligent Laboratory Systems}, 78(1–2), 103–112. <doi:10.1016/j.chemolab.2005.04.003>
#'
#' @export

metabol.vip <- function(
    normalized_data,
    top_n = 20,
    component = 1,
    plot = TRUE,
    calc_auc = FALSE,
    group_col = "Group",
    assign_result = TRUE,
    assign_name = "vip_data",
    envir = parent.frame()
) {

  # ------------------------------
  # 1) Required packages
  # ------------------------------
  pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "pROC")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p != "pROC" || calc_auc) {
        stop(sprintf("Please install '%s'", p))
      }
    }
  }

  # ------------------------------
  # 2) Detect object structure
  # ------------------------------
  obj <- if ("dimred_obj" %in% names(normalized_data)) {
    normalized_data$dimred_obj
  } else {
    normalized_data
  }

  if (is.null(obj$res_dimred$vip_scores)) {
    stop("This object has no VIP scores. Run metabol.dimred() first.")
  }

  if (!is.list(obj) ||
      !all(c("expr_matrix", "metadata", "res_dimred") %in% names(obj))) {
    stop("Input must be an object returned by metabol.dimred().")
  }

  expr_mat <- obj$expr_matrix
  metadata <- obj$metadata

  # ------------------------------
  # 3) Extract VIP scores
  # ------------------------------
  vip_scores <- obj$res_dimred$vip_scores
  if (is.null(vip_scores)) {
    plsda_res <- obj$res_dimred$plsda_res
    vip_scores <- mixOmics::vip(plsda_res)
  }

  if (is.null(vip_scores)) {
    stop("No VIP scores found. Run metabol.dimred() first.")
  }

  if (component > ncol(vip_scores)) {
    stop("Requested component exceeds number of components in VIP matrix.")
  }

  vip_comp <- vip_scores[, component]

  n_top <- min(top_n, length(vip_comp))
  top_vip <- sort(vip_comp, decreasing = TRUE)[seq_len(n_top)]

  df_vip <- stats::setNames(
    data.frame(Metabolite = names(top_vip), VIP = as.numeric(top_vip)),
    c("Metabolite", "VIP")
  )

  # ------------------------------
  # 4) Heatmap preparation
  # ------------------------------
  if (!group_col %in% colnames(metadata)) {
    stop(sprintf("Column '%s' not found in metadata.", group_col))
  }

  cluster_groups <- as.factor(metadata[[group_col]])

  # Ensure metabolites exist in expression matrix
  valid_metabs <- df_vip$Metabolite[df_vip$Metabolite %in% colnames(expr_mat)]
  if (length(valid_metabs) == 0) {
    stop("None of the top VIP metabolites were found in expr_matrix.")
  }
  df_vip <- df_vip[df_vip$Metabolite %in% valid_metabs, , drop = FALSE]

  # Mean by cluster
  medias_por_cluster <- stats::aggregate(
    expr_mat[, df_vip$Metabolite, drop = FALSE],
    by = list(Cluster = cluster_groups),
    FUN = base::mean
  )

  heatmap_long <- tidyr::pivot_longer(
    medias_por_cluster,
    cols = -Cluster,
    names_to = "Metabolite",
    values_to = "ValorMedio"
  )

  heatmap_long <- dplyr::group_by(heatmap_long, .data$Metabolite)
  heatmap_long <- dplyr::mutate(
    heatmap_long,
    ValorZ = as.numeric(scale(.data$ValorMedio)[, 1])
  )
  heatmap_long <- dplyr::ungroup(heatmap_long)

  df_vip$Metabolite <- factor(df_vip$Metabolite, levels = rev(df_vip$Metabolite))
  heatmap_long$Metabolite <- factor(
    heatmap_long$Metabolite,
    levels = levels(df_vip$Metabolite)
  )

  # ------------------------------
  # 5) Plots
  # ------------------------------
  plot_vip <- plot_heatmap <- NULL

  if (plot) {
    plot_vip <- ggplot2::ggplot(df_vip, ggplot2::aes(x = .data$Metabolite, y = .data$VIP)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(
        title = paste("VIP Scores - Component", component),
        y = "VIP Score",
        x = NULL
      ) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

    plot_heatmap <- ggplot2::ggplot(
      heatmap_long,
      ggplot2::aes(
        x = as.factor(.data$Cluster),
        y = .data$Metabolite,
        fill = .data$ValorZ
      )
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        low = "#4575b4",
        mid = "white",
        high = "#d73027",
        midpoint = 0,
        name = "Z-score"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(x = "", y = NULL, title = "") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )

    # explicit print() is fine inside interactive plotting
    print(plot_vip + plot_heatmap + patchwork::plot_layout(widths = c(2.5, 1)))
  }

  # ------------------------------
  # 6) Optional AUC calculation
  # ------------------------------
  resultados_auc <- NULL

  if (calc_auc) {
    resultados_auc <- data.frame(Metabolite = character(), AUC = numeric())
    classe <- as.factor(metadata[[group_col]])

    for (metab in df_vip$Metabolite) {
      roc_obj <- pROC::roc(
        response = classe,
        predictor = expr_mat[, as.character(metab)],
        levels = levels(classe),
        direction = "<"
      )

      auc_val <- as.numeric(pROC::auc(roc_obj))
      resultados_auc <- rbind(
        resultados_auc,
        data.frame(Metabolite = metab, AUC = auc_val)
      )
    }

    resultados_auc <- resultados_auc[order(-resultados_auc$AUC), , drop = FALSE]
  }

  # ------------------------------
  # 7) Assign + Return
  # ------------------------------
  result <- list(
    vip_scores = vip_scores,
    top_vip = df_vip,
    auc = resultados_auc,
    plot_vip = plot_vip,
    plot_heatmap = plot_heatmap
  )

  if (assign_result) {
    base::assign(assign_name, result, envir = envir)
  }

  invisible(result)
}
