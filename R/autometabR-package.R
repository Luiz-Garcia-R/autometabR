#' autometabR: Quick Reference and Main Functions for Metabolomics Analysis
#'
#' autometabR provides a comprehensive set of functions for metabolomics data analysis,
#' from raw data import to normalization, quality control, exploratory and differential
#' analysis. This help topic serves as a quick reference and guide for users.
#'
#'
#' ## Main Workflow
#'
#' The recommended workflow guides users from raw data import to quality control and initial metabolite evaluation:
#'  - metabol.import() – Import and validate raw metabolomics data (raw_data + metadata).
#'  - metabol.normalize() – Normalize imported data (imp_data), filter features, impute missing values, and remove outliers.
#'  - metabol.qc() – Generate QC plots including boxplots, PCA, and density distributions.
#'  - metabol.info() – Summarize general metabolite characteristics and group-level statistics.
#'
#'
#'
#' Main Functions Overview
#' | Function | Description |
#' |----------|-------------|
#' | `metabol.import()`   | Import and validate metabolomics data |
#' | `metabol.normalize()`| Normalize, filter, and impute missing values |
#' | `metabol.qc()`       | QC plots: boxplots, PCA, density distributions |
#' | `metabol.info()`     | Summarize metabolites and group statistics |
#' | `metabol.corr()`     | Sample/group correlation matrices |
#' | `metabol.diff()`     | Identify differentially expressed metabolites |
#' | `metabol.dimred()`   | PCA, UMAP and PLS-DA for dimensionality reduction |
#' | `metabol.enrich()`   | Metabolite set enrichment analysis |
#' | `metabol.heatmap()`  | Heatmap of top variable metabolites |
#' | `metabol.oplsda()`   | Group separation via OPLS-DA |
#' | `metabol.roc()`      | ROC curves for discriminant metabolites |
#' | `metabol.boxplot()`  | Boxplot based on log2FC difference |
#' | `metabol.vip()`      | Highlight top contributing metabolites |
#'
#' Contact and Contributions
#' For suggestions, bug reports, or contributions, see the
#' [GitHub repository](https://github.com/Luiz-Garcia-R/autometabR).
#'
#' Example Workflow
#'
#' A small synthetic metabolomics dataset with 50 metabolites and 5 samples per group.
#' This dataset is for demonstration purposes only.
#'
#' @examples
#' \dontrun{
#'
#' # Small synthetic metabolomics dataset
#' metabolites <- paste0("M", sprintf("%03d", 1:50))
#' raw_data <- data.frame(
#'   Sample = paste0("S", 1:10),
#'   matrix(
#'     rnorm(50*10, mean = 100, sd = 20),
#'     nrow = 10,
#'     ncol = 50,
#'     dimnames = list(NULL, metabolites)
#'   )
#' )
#'
#' metadata <- data.frame(
#'   Sample = raw_data$Sample,
#'   Group  = rep(c("Control","Treatment"), each = 5)
#' )
#'
#' # Convert to metabolR object
#' imp_data <- metabol.import(raw_data, metadata)
#'
#' # Normalize
#' normalized <- metabol.normalize(imp_data)
#'
#' # QC
#' metabol.qc(normalized)
#'
#' # Exploratory and differential analysis
#' metabol.corr(normalized)
#' metabol.dimred(normalized)
#'
#' # Print summary of normalization
#' print(normalized)
#'}
#'
#' @name autometabR
#'
"_PACKAGE"
