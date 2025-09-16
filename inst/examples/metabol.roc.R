## ---- metabol.roc.R ----

# Minimal example
set.seed(123)
raw_data <- data.frame(
  Sample = paste0("S", 1:20),
  Creatine = c(rnorm(10, 5, 1), rnorm(10, 6, 1)),
  Glucose  = c(rnorm(10, 10, 2), rnorm(10, 10.5, 2))
)

meta <- data.frame(
  Sample = paste0("S", 1:20),
  Group = rep(c("A", "B"), each = 10)
)

imp_data <- list(log2_normalized = raw_data)

res_roc <- metabol.roc(
  normalized_data = imp_data,
  metabolites = c("Creatine", "Glucose"),
  metadata = meta
)
res_roc$auc
