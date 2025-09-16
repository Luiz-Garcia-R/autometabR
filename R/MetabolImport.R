#' Import metabolomics data and metadata
#'
#' Prepares raw metabolomics data and optional metadata into a standardized \code{metabolR} object.
#' Ensures consistent formatting, validation, and easy downstream processing.
#'
#' @param raw_data A data.frame or tibble containing a 'Sample' column and metabolite quantifications in the remaining columns.
#' @param metadata Optional data.frame or tibble containing at least a 'Sample' column. A 'Group' column is recommended.
#' @param strict Logical; if TRUE, stops execution on data validation errors. Default is TRUE.
#' @param assign_result Logical; if TRUE, assigns the object to the environment. Default is FALSE.
#' @param assign_name Character; name to assign the object if assign_result = TRUE.
#' @param envir Environment in which to assign the object if assign_result = TRUE. Default is parent.frame().
#'
#' @return An object of class \code{metabolR} containing data, metadata, warnings, and summary counts.
#'
#' @examples
#' raw_data <- data.frame(
#'   Sample = c("S1", "S2", "S3"),
#'   Met1 = c(10, 15, 12),
#'   Met2 = c(5, 3, 6)
#' )
#' metadata <- data.frame(
#'   Sample = c("S1", "S2", "S3"),
#'   Group  = c("Control", "Treatment", "Control")
#' )
#'
#' imp <- metabol.import(raw_data, metadata)
#' print(imp)
#'
#' @references
#' Xia J, Sinelnikov IV, Han B, Wishart DS (2015). MetaboAnalyst 3.0—making metabolomics more meaningful. \emph{Nucleic Acids Research}, 43(W1): W251–W257. <doi:10.1093/nar/gkv380>
#' Pang Z, Chong J, Zhou G, de Lima Morais DA, Chang L, Barrette M, et al. (2021). MetaboAnalyst 5.0: narrowing the gap between raw spectra and functional insights. \emph{Nucleic Acids Research}, 49(W1): W388–W396. <doi:10.1093/nar/gkab382>
#'
#' @export

metabol.import <- function(raw_data, metadata = NULL, strict = TRUE,
                           assign_result = FALSE, assign_name = "imp_data",
                           envir = parent.frame()) {

  errors <- character()
  warnings_list <- character()

  add_error   <- function(txt) errors <<- c(errors, txt)
  add_warning <- function(txt) warnings_list <<- c(warnings_list, txt)

  ## 1. Type check
  if (!is.data.frame(raw_data)) add_error("'raw_data' must be a data.frame or tibble.")

  ## 2. Sample column
  if (!"Sample" %in% colnames(raw_data)) add_error("'raw_data' must contain a 'Sample' column as the first column.")

  ## 3. Clean metabolite names
  n_fixed <- 0
  if ("Sample" %in% colnames(raw_data)) {
    metabolite_cols <- setdiff(colnames(raw_data), "Sample")
    cleaned_names <- gsub("[^A-Za-z0-9]", "_", metabolite_cols)
    cleaned_names <- make.unique(cleaned_names, sep = "_")

    if (!identical(metabolite_cols, cleaned_names)) {
      colnames(raw_data)[-1] <- cleaned_names
      n_fixed <- sum(metabolite_cols != cleaned_names)
      add_warning(paste0(n_fixed, " metabolite names adjusted."))
    }
  }

  ## 4. Check numeric columns
  if ("Sample" %in% colnames(raw_data)) {
    non_numeric <- names(raw_data[-1])[!sapply(raw_data[-1], is.numeric)]
    if (length(non_numeric) > 0) add_error(paste("The following columns are not numeric:", paste(non_numeric, collapse = ", ")))
  }

  ## 5. Metadata checks
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) add_error("'metadata' must be a data.frame or tibble.") else {
      if (!"Sample" %in% colnames(metadata)) add_error("Metadata must contain a column named 'Sample'.") else {
        missing_samples <- setdiff(metadata$Sample, raw_data$Sample)
        if (length(missing_samples) > 0) add_warning("Some samples in metadata are not present in 'raw_data'.")
      }
      if (!"Group" %in% colnames(metadata)) add_warning("Metadata does not contain a 'Group' column.")
    }
  }

  ## 6. Stop if errors
  if (length(errors) > 0) {
    if (strict) stop(paste("Data validation error(s):\n-", paste(errors, collapse = "\n- ")), call. = FALSE)
    warnings_list <- c(warnings_list, errors)
  }

  ## 7. Build object
  obj <- list(
    data = raw_data,
    metadata = metadata,
    warnings = warnings_list,
    n_metabolites = ncol(raw_data) - 1,
    n_samples = nrow(raw_data),
    n_fixed = n_fixed
  )
  class(obj) <- "metabolR"

  ## 8. Assign if requested
  if(assign_result) assign(assign_name, obj, envir = envir)

  return(obj)
}

#' @export
print.metabolR <- function(x, ...) {
  msg <- paste0(
    "Object of class 'metabolR'\n",
    "------------------------------\n",
    "Samples:     ", x$n_samples, "\n",
    "Metabolites: ", x$n_metabolites, "\n"
  )
  if (!is.null(x$metadata)) msg <- paste0(msg, "Metadata:    ", nrow(x$metadata), " rows\n")
  if (x$n_fixed > 0) msg <- paste0(msg, "Note:        ", x$n_fixed, " metabolite names were adjusted.\n")
  message(msg)
  invisible(x)
}
