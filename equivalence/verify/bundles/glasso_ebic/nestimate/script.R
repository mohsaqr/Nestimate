library(Nestimate)
strategies <- read.csv("../data/srl_strategies.csv")
network <- build_network(strategies, method = "glasso")
network
result <- transform(as.data.frame(network), quantity = paste(from, "--", to),
                    kind = "partial correlation", value = weight)
result
