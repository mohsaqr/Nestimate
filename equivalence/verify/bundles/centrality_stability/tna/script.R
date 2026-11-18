library(tna)
sequences <- read.csv("../data/sequences.csv")
model <- tna(sequences)
set.seed(1)
stability <- estimate_cs(model, measures = c("InStrength", "OutStrength", "Betweenness"), iter = 200)
stability
result <- do.call(rbind, lapply(names(stability), function(measure) {
  correlations <- stability[[measure]]$correlations
  data.frame(drop = paste0(measure, ", dropping ", attr(stability, "drop_prop")),
             mean_cor = colMeans(correlations), sd_cor = apply(correlations, 2, sd),
             prop_above = colMeans(correlations >= 0.7))
}))
result
