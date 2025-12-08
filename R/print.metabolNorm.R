#' Print Method for metabolNorm Objects
#'
#' Custom print method for objects of class \code{metabolNorm}.
#'
#' @param x An object of class \code{metabolNorm}.
#' @param ... Additional arguments (ignored).
#'
#' @export

print.metabolNorm <- function(x, ...) {
  message("==============================")
  message("Object of class 'metabolNorm'")
  message("==============================")
  message("Normalization method: ", x$method)
  message("Metabolites: ", ncol(x$expr_matrix) - 1)
  message("Samples: ", nrow(x$expr_matrix))

  if(!is.null(x$checks) && length(x$checks) > 0){
    message("\nChecks/Transformations performed:")
    for(name in names(x$checks)){
      val <- x$checks[[name]]
      if(is.numeric(val) && length(val) > 1){
        message(" - ", name, ": ", paste(round(val, 2), collapse=", "))
      } else {
        message(" - ", name, ": ", val)
      }
    }
  }
  message("==============================")
  message("Use functions like metabol.dimred(), metabol.corr(), metabol.heatmap() on this object.")
  invisible(x)
}

