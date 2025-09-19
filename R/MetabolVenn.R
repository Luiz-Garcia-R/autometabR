#' Venn Diagram for Metabolite Presence Across Groups
#'
#' Generates a Venn diagram showing shared and unique metabolites across groups.
#'
#' @param normalized_data Object returned by \code{metabol.qc()}, containing \code{$presence_per_group}.
#' @param cutoff Numeric. Minimum value to consider a metabolite as present (default = 0).
#' @param results Logical. If TRUE, returns detailed lists of metabolites; default FALSE.
#'
#' @return Invisibly returns a list with metabolites by group, unique metabolites, and shared metabolites (if results = TRUE). Also prints a Venn diagram.
#'
#' @examples
#' \dontrun{
#' presence_per_group <- data.frame(
#'   MetaboliteID = paste0("Metab", 1:7),
#'   GroupA = c(1,1,0,2,2,0,1),
#'   GroupB = c(0,1,3,1,0,1,3)
#' )
#'
#' qc_res <- list(presence_per_group = presence_per_group)
#'
#' metabol.venn(qc_res, cutoff = 0, results = TRUE)
#' }
#'
#' @references
#' Venn, J. (1880). On the diagrammatic and mechanical representation of propositions and reasoning. \emph{The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science}, 10: 1–18. <doi:10.1080/14786448008626994>
#' Chen, H., & Boutros, P.C. (2011). VennDiagram: a package for the generation of highly-customizable Venn and Euler diagrams in R. \emph{BMC Bioinformatics}, 12: 35. <doi:10.1186/1471-2105-12-35>
#'
#' @export

metabol.venn <- function(normalized_data, cutoff = 0, results = FALSE) {

  # --- Check packages ---
  req_pkgs <- c("ggvenn", "RColorBrewer", "ggplot2")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # --- Extract presence matrix ---
  if (is.list(normalized_data) && !is.null(normalized_data$presence_per_group)) {
    presence_per_group <- normalized_data$presence_per_group
  } else if (!is.null(attr(normalized_data, "presence_per_group"))) {
    presence_per_group <- attr(normalized_data, "presence_per_group")
  } else {
    stop("Object does not contain 'presence_per_group'. Use metabol.qc(..., remove_zeros = FALSE).")
  }

  presence_per_group <- as.data.frame(presence_per_group)

  # --- Groups and metabolites ---
  groups <- colnames(presence_per_group)[-1]  # first column is MetaboliteID
  metabolite_ids <- presence_per_group[[1]]

  # --- Metabolites per group ---
  metabolites_by_group <- lapply(groups, function(g) metabolite_ids[presence_per_group[[g]] > cutoff])
  names(metabolites_by_group) <- groups

  # --- Unique and shared metabolites ---
  metabolites_unique <- lapply(groups, function(g) setdiff(metabolites_by_group[[g]], unlist(metabolites_by_group[groups != g])))
  names(metabolites_unique) <- groups
  metabolites_shared <- Reduce(intersect, metabolites_by_group)

  # --- Venn diagram ---
  n_colors <- max(3, length(groups))
  brewer_colors <- RColorBrewer::brewer.pal(n_colors, "Set1")[1:length(groups)]

  p <- ggvenn::ggvenn(
    metabolites_by_group,
    fill_color = brewer_colors,
    stroke_size = 0,
    set_name_size = 4,
    text_size = 4
  ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "right"
    )

  print(p)

  # --- Optional message ---
  if (!results) {
    message("Tip: To access detailed lists of metabolites, call with 'results = TRUE'.")
  }

  # --- Return results invisibly ---
  if (results) {
    return(invisible(list(
      metabolites_by_group = metabolites_by_group,
      metabolites_unique = metabolites_unique,
      metabolites_shared = metabolites_shared
    )))
  } else {
    invisible(NULL)
  }
}
