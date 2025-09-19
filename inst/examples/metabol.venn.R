## ---- metabol.venn.R ----

# Minimal example
presence_per_group <- data.frame(
  MetaboliteID = paste0("Metab", 1:7),
  GroupA = c(1,1,0,2,2,0,1),
  GroupB = c(0,1,3,1,0,1,3)
)

qc_res <- list(presence_per_group = presence_per_group)

metabol.venn(qc_res, cutoff = 0, results = TRUE)
