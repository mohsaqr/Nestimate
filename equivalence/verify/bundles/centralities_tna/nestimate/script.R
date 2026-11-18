library(Nestimate)
sequences <- read.csv("../data/sequences.csv")
network <- build_network(sequences, method = "relative")
result <- net_centrality(network, measures = c("BetweennessRSP", "Diffusion", "Clustering"),
                         normalize_diffusion = FALSE)
result
