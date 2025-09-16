## ---- metabol.heatmap.R ----

# Minimal example
expr_mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
rownames(expr_mat) <- paste0("S", 1:5)
colnames(expr_mat) <- paste0("Met", 1:4)
meta <- data.frame(
  Sample = paste0("S", 1:5),
  Group = c("A","A","B","B","B")
)
metabol.heatmap(expr_mat, meta, top_n = 3, results = TRUE)
