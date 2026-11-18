library(Nestimate)
strategies <- read.csv("../data/srl_strategies.csv")
correlation <- build_network(strategies, method = "cor")
partial <- build_network(strategies, method = "pcor")
result <- transform(rbind(as.data.frame(correlation), as.data.frame(partial)),
                    quantity = paste0(network, ": ", from, " -- ", to), kind = network,
                    value = weight)
result
