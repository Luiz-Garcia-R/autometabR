.onAttach <- function(libname, pkgname) {

  version_text <- paste0(
    "autometabR v.",
    utils::packageVersion("autometabR")
  )

  if (requireNamespace("crayon", quietly = TRUE)) {
    version_text <- crayon::green(version_text)
  }

  packageStartupMessage(
    "\n",
    version_text, " loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined metabolomics data analysis.\n",
    "GitHub: https://github.com/Luiz-Garcia-R/autometabR\n"
  )
}
