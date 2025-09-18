## ---- metabol.dems.R ----

# Minimal example for metabol.dems
expr_matrix <- data.frame(
  Sample = c("S1","S2","S3","S4"),
  Met1   = c(5,6,1,2),
  Met2   = c(2,3,2,3)
)

norm_data <- list(expr_matrix = expr_matrix)
class(norm_data) <- "metabolNorm"

meta <- data.frame(
  Sample = c("S1","S2","S3","S4"),
  Group  = c("A","A","B","B")
)

metabol.dems(norm_data, metadata = meta)
