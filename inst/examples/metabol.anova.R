## ---- metabol.anova.R ----

# Minimal example
set.seed(1)
expr <- data.frame(
        Sample = paste0("S", 1:9),
        Glucose = c(rnorm(3, 5, 1), rnorm(3, 6, 1), rnorm(3, 8, 1))
)

meta <- data.frame(
        Sample = paste0("S", 1:9),
        Group = rep(c("A", "B", "C"), each = 3)
)

norm_data <- list(expr_matrix = expr)
class(norm_data) <- "normalized_data"

metabol.anova(norm_data, metadata = meta, metabolite = "Glucose")
