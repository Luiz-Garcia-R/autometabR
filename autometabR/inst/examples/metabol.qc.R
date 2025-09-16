## ---- metabol.qc.R ----

# Minimal example
set.seed(123)
raw_data <- data.frame(
  Sample = paste0("S", 1:10),
  Met1 = rnorm(10),
  Met2 = rnorm(10)
)

meta <- data.frame(
  Sample = paste0("S", 1:10),
  Group = rep(c("A","B"), each = 5)
)

imp_data <- list(data = raw_data, metadata = meta)
class(imp_data) <- "metabolR"

norm_obj <- metabol.normalize(imp_data)

metabol.qc(norm_obj, metadata = meta)
