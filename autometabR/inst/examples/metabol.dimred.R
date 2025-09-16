## ---- metabol.dimred.R ----

# Minimal example
set.seed(123)
expr_mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
rownames(expr_mat) <- paste0("S", 1:20)
colnames(expr_mat) <- paste0("Met", 1:5)

meta$Sample <- factor(meta$Sample)
meta <- data.frame(
  Sample = paste0("S", 1:20),
  Group  = rep(c("A","B"), each = 10)
)

res <- metabol.dimred(expr_mat, metadata = meta, cross_val = 3)
