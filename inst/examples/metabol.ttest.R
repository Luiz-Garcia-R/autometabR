## ---- metabol.ttest.R ----

# Minimal example
raw_data <- data.frame(Sample = paste0("S", 1:6),
                       Lactate = rnorm(6))
meta <- data.frame(Sample = paste0("S", 1:6), Group = c("A","A","B","B","A","B"))
imp_data <- list(expr_matrix = raw_data)

metabol.ttest(imp_data, metadata = meta, metabolite = "Lactate", return_type = "all")
