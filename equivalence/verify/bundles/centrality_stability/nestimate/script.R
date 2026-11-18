library(Nestimate)
sequences <- read.csv("../data/sequences.csv")
network <- build_network(sequences, method = "relative")
stability <- centrality_stability(network, measures = c("InStrength", "OutStrength", "Betweenness"),
                                  iter = 200, seed = 1)
stability
result <- subset(transform(summary(stability), drop = paste0(measure, ", dropping ", drop_prop)),
                 select = c(drop, mean_cor, sd_cor, prop_above))
result
