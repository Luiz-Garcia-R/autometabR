## ---- metabol.corr.R ----

# Minimal example
set.seed(123)
norm_data <- data.frame(
  Sample = paste0("S", 1:8),
  Met1   = rnorm(8, mean = 5, sd = 1),
  Met2   = rnorm(8, mean = 6, sd = 1.2),
  Met3   = rnorm(8, mean = 7, sd = 1.5),
  Met4   = rnorm(8, mean = 5.5, sd = 0.8),
  Met5   = rnorm(8, mean = 6.5, sd = 1.1)
)

meta <- data.frame(
  Sample = paste0("S", 1:8),
  Group  = rep(c("A", "B"), each = 4)
)

correlation <- metabol.corr(norm_data, metadata = meta)
print(correlation)
