## ---- metabol.info.R ----

# Minimal example
expr_mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
rownames(expr_mat) <- paste0("S", 1:5)
colnames(expr_mat) <- paste0("Met", 1:4)

norm_data <- list(expr_matrix = expr_mat)
meta <- data.frame(Sample = paste0("S", 1:5), Group = c("A","A","B","B","B"))

info <- metabol.info(norm_data, meta)

print(info)
