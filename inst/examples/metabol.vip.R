## ---- metabol.vip.R ----

# Minimal example
set.seed(123)

expr_mat <- data.frame(
  Met1 = rnorm(6),
  Met2 = rnorm(6),
  Met3 = rnorm(6)
)
metadata <- data.frame(
  Sample = paste0("S", 1:6),
  Group  = rep(c("A","B"), each = 3)
)

plsda_res <- mixOmics::plsda(expr_mat, metadata$Group, ncomp = 2)

dimred_data <- list(
  expr_matrix = expr_mat,
  metadata    = metadata,
  res_dimred  = list(plsda_res = plsda_res)
)

metabol.vip(dimred_data, n_top = 3, component = 1, plot = TRUE)
