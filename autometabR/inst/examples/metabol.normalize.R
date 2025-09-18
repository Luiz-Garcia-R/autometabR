## ---- metabol.normalize.R ----

# Minimal example
set.seed(123)

# expression data
expr_data <- data.frame(
  Sample = paste0('S', 1:5),
  Met1 = c(5.2, 3.1, 4.8, 5.5, 4.2),
  Met2 = c(2.3, 2.7, 2.5, 2.6, 2.4),
  Met3 = c(8.1, 7.8, 8.5, 8.0, 8.3)
)

# metadata
metadata <- data.frame(
  Sample = paste0('S', 1:5),
  Group = c('A','A','B','B','A')
)

# Create metabolR object
imp_data <- list(data = expr_data, metadata = metadata)
class(imp_data) <- "metabolR"

# Test normalization
norm_res <- metabol.normalize(imp_data)
