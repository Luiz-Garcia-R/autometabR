.onLoad <- function(libname, pkgname) {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    invisible(NULL)
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    crayon::green("autometabR "), "loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined metabolomics data analysis.\n",
    "Use ", crayon::green("?autometabR"), " for general help.\n",
    "--------------------------------------------------\n"
  )
}
