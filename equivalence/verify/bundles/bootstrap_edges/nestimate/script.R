library(Nestimate)
sequences <- read.csv("../data/sequences.csv")
network <- build_network(sequences, method = "relative")
boot <- bootstrap_network(network, iter = 200, seed = 1)
boot
result <- subset(transform(summary(boot), edge = paste(from, "->", to)),
                 select = c(edge, weight, mean, sd, p_value, sig, ci_lower, ci_upper, cr_lower, cr_upper))
result
