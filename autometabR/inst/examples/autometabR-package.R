## ---- autometabR-package.R ----

# Small synthetic metabolomics dataset
metabolites <- paste0("M", sprintf("%03d", 1:50))
raw_data <- data.frame(
  Sample = paste0("S", 1:10),
  matrix(
    rnorm(50*10, mean = 100, sd = 20),
    nrow = 10,
    ncol = 50,
    dimnames = list(NULL, metabolites)
  )
)

metadata <- data.frame(
  Sample = raw_data$Sample,
  Group  = rep(c("Control","Treatment"), each = 5)
)

# Convert to metabolR object
imp_data <- metabol.import(raw_data, metadata)

# Normalize
normalized <- metabol.normalize(imp_data)

# QC
metabol.qc(normalized, metadata)

# Exploratory and differential analysis
metabol.corr(normalized, metadata)
metabol.dems(normalized, metadata)

# Print summary of normalization
print(normalized)

