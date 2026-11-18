library(igraph)
sequences <- read.csv("../data/sequences.csv")
states <- sort(unique(na.omit(unlist(sequences))))
from <- factor(unlist(sequences[-ncol(sequences)]), levels = states)
to <- factor(unlist(sequences[-1]), levels = states)
weights <- unclass(prop.table(table(from, to), margin = 1))
diag(weights) <- 0
graph <- graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE)
lengths <- 1 / E(graph)$weight
result <- data.frame(state = V(graph)$name,
                     OutStrength = strength(graph, mode = "out"),
                     InStrength = strength(graph, mode = "in"),
                     ClosenessIn = closeness(graph, mode = "in", weights = lengths),
                     ClosenessOut = closeness(graph, mode = "out", weights = lengths),
                     Closeness = closeness(graph, mode = "all", weights = lengths),
                     Betweenness = betweenness(graph, weights = lengths))
result
