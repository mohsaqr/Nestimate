library(qgraph)
strategies <- read.csv("../data/srl_strategies.csv")
fit <- EBICglasso(cor(strategies), n = nrow(strategies), gamma = 0.5, returnAllResults = TRUE)
fit$lambda[which.min(fit$ebic)]
network <- fit$optnet
pairs <- which(upper.tri(network) & network != 0, arr.ind = TRUE)
result <- data.frame(quantity = paste(colnames(network)[pairs[, "row"]], "--",
                                      colnames(network)[pairs[, "col"]]),
                     kind = "partial correlation", value = network[pairs])
result
