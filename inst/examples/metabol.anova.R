## ---- metabol.anova.R ----

# Minimal example
set.seed(1)
expr <- data.frame(
  Glucose = c(rnorm(3, 5, 1), rnorm(3, 6, 1), rnorm(3, 8, 1))
)
rownames(expr) <- paste0("S", 1:9)  # força os nomes de amostra
meta <- data.frame(
  Sample = paste0("S", 1:9),
  Group = rep(c("A", "B", "C"), each = 3)
)

norm_data <- list(expr_matrix = expr)
class(norm_data) <- "normalized_data"

metabol.anova(norm_data, metadata = meta, metabolite = "Glucose")
