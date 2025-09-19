## ---- metabol.import.R ----

raw_data <- data.frame(
  Sample = c("S1", "S2", "S3"),
  Met1 = c(10, 15, 12),
  Met2 = c(5, 3, 6)
)
metadata <- data.frame(
  Sample = c("S1", "S2", "S3"),
  Group  = c("Control", "Treatment", "Control")
)

imp <- metabol.import(raw_data, metadata)

print(imp)
