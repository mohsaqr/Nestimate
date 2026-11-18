library(Nestimate)
sequences <- read.csv("../data/sequences.csv")
network <- build_network(sequences, method = "relative")
result <- net_centrality(network, measures = c("OutStrength", "InStrength", "ClosenessIn",
                                               "ClosenessOut", "Closeness", "Betweenness"))
result
