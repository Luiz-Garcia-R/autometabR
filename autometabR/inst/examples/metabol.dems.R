## ---- metabol.dems.R ----

# Minimal example
norm_data <- list(
  log2_normalized = data.frame(
    Sample = c("S1","S2","S3","S4"),
    Met1 = c(1,2,3,4),
    Met2 = c(2,3,1,5)
  )
)
meta <- data.frame(
  Sample = c("S1","S2","S3","S4"),
  Group  = c("A","A","B","B")
)
metabol.dems(norm_data, metadata = meta)
