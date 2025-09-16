.onLoad <- function(libname, pkgname) {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    invisible(NULL)
  }
}
