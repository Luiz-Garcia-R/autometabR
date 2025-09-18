## ---- metabol.oplsda.R ----

# Minimal example
if (requireNamespace("ropls", quietly = TRUE)) {
  set.seed(123)
  raw_data <- data.frame(
    Sample = paste0("S", 1:20),
    Met1 = c(rep(2,10), rep(-2,10)),
    Met2 = c(rep(1,10), rep(-1,10))
  )

  meta <- data.frame(
    Sample = paste0("S", 1:20),
    Group  = rep(c("A","B"), each = 10)
  )

  imp_data <- list(data = raw_data, metadata = meta)
  class(imp_data) <- "metabolR"

  norm_obj <- metabol.normalize(imp_data)
  metabol.oplsda(norm_obj, meta)
}
