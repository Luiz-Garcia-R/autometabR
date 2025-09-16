#' Print Method for metabNorm Objects
#'
#' Custom print method for objects of class \code{metabolNorm}.
#'
#' @param x An object of class \code{metabolNorm}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.metabolNorm <- function(x, ...) {
  norm_name <- names(x)[grep("_normalized$", names(x))]
  norm_df <- x[[norm_name]]

  cat("Object of class 'metabolNorm'\n")
  cat("------------------------------\n")
  cat("Normalization method:", x$method, "\n")
  cat("Metabolites:", ncol(norm_df) - 1, "\n")
  cat("Samples:", nrow(norm_df), "\n\n")

  if(length(x$checks) > 0){
    cat("Checks/Transformations performed:\n")
    for(name in names(x$checks)){
      val <- x$checks[[name]]
      if(is.numeric(val) && length(val) > 1){
        cat(" -", name, ":", paste(round(val, 2), collapse=", "), "\n")
      } else cat(" -", name, ":", val, "\n")
    }
  }

  cat("------------------------------\n")
  cat("Use functions like metabol.dimred(), metabol.corr(), metabol.heatmap() on this object.\n")
  invisible(x)
}
