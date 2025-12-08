# autometabR

<!-- badges: start --> <!-- badges: end -->

**autometabR** is an R package designed to streamline metabolomics data analysis, from raw data import to normalization, quality control, statistical testing, and visualization.
It provides functions for filtering, imputation, outlier detection, dimensionality reduction, correlation, enrichment, heatmaps, volcano plots, and statistical evaluation.

**Note:** The package is primarily optimized for pairwise analyses (e.g., Control vs. Treatment).
Workflows involving multiple groups can be explored, but most functions are tuned for two-group comparisons.

## Installation

You can install the development version directly from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install autometabR from GitHub
devtools::install_github("Luiz-Garcia-R/autometabR")
```

## Input Data Format (Very Important)

To use autometabR, you need two input data frames:

  1. raw_data – metabolite quantifications (from COLMAR, Bayesil, etc.)
    - Must contain a first column called Sample
    - Other columns must have metabolites with numeric quantification
  
Example of raw_data:
    
| Sample       | Metabolite\_1 | Metabolite\_2 | Metabolite\_3 | Metabolite\_4 | Metabolite\_5 | Metabolite\_6 |
| ------------ | ----------    | ----------    | ----------    | ------------  | ------------  | ------------  |
| Sample1      | 1.2           | 1.5           | 1.1           | 2.0           | 2.1           | 0.0           |
| Sample2      | 0.0           | 1.0           | 3.6           | 3.5           | 3.7           | 0.0           |
| Sample3      | 2.1           | 2.3           | 2.2           | 1.9           | 0.0           | 0.0           |
| Sample4      | 2.2           | 0.0           | 0.0           | 2.5           | 0.0           | 0.0           |


  2. metadata – describes your samples and experimental groups
    Must contain two columns:
      - Sample: matches exactly the first column of raw_data
      - Group: experimental group/condition

Example of metadata:

| Sample       | Group     |
| ------------ | --------- |
| Control\_1   | Control   |
| Control\_2   | Control   |
| Control\_3   | Control   |
| Treatment\_1 | Treatment |
| Treatment\_2 | Treatment |
| Treatment\_3 | Treatment |

## Main Workflow

The recommended workflow guides users from raw data import to quality control and initial metabolite evaluation:
  - metabol.import() – Import and validate raw metabolomics data (raw_data + metadata).
  - metabol.normalize() – Normalize data, filter features, impute missing values, and remove outliers.
  - metabol.qc() – Generate QC plots including boxplots, PCA, and density distributions.
  - metabol.info() – Summarize general metabolite characteristics and group-level statistics.


# Exploratory and differential analysis
For deeper insights and group comparisons:
  - metabol.corr() – Compute correlations among samples or experimental groups.
  - metabol.diff() – Identify differentially expressed metabolites between conditions.
  - metabol.dimred() – Perform dimensionality reduction using PCA and UMAP.
  - metabol.heatmap() – Visualize top variable metabolites with heatmaps.
  - metabol.oplsda() – Evaluate group differences via OPLS-DA.
  - metabol.roc() – Assess the discriminatory power of specific metabolites.
  - metabol.vip() – Highlight the most important metabolites contributing to group separation.
  - metabol.volcano() – Generate volcano plots for differential metabolite analysis.


## Example (Minimal)

# Load example raw data and metadata

```r
# Example raw data
raw_data <- data.frame(
  MetaboliteID = c("M001","M002","M003"),
  Control1 = c(1.0, 2.0, 1.5),
  Control2 = c(1.2, 2.1, 1.7),
  Treatment1 = c(2.3, 1.8, 2.0),
  Treatment2 = c(2.1, 1.9, 2.2)
)

metadata <- data.frame(
  Sample = c("Control1","Control2","Treatment1","Treatment2"),
  Group = c("Control","Control","Treatment","Treatment")
)

# Import data
obj <- metabol.import(raw_data, metadata) # returns a metabR object

# Normalize
normalized_data <- metabol.normalize(obj)

# QC plots
metabol.qc(normalized_data)

# Correlation
corr_mat <- metabol.corr(normalized_data)
```

## Contact

For questions, suggestions, or contributions, open an issue or pull request on GitHub.

Thank you for using autometabR!
