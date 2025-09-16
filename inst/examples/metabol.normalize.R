## ---- metabol.normalize.R ----

# Minimal example
raw_data <- data.frame(Sample = paste0("S", 1:5), Met1 = rnorm(5), Met2 = rnorm(5))
meta <- data.frame(Sample = paste0("S", 1:5), Group = c("A","A","B","B","B"))
imp_data <- list(data = raw_data, metadata = meta)

class(imp_data) <- "metabolR"

norm_obj <- metabol.normalize(imp_data)

print(norm_obj)
