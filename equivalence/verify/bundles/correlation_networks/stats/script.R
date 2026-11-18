library(stats)
strategies <- read.csv("../data/srl_strategies.csv")
correlation <- cor(strategies)
partial <- -cov2cor(solve(correlation))
pairs <- which(upper.tri(correlation), arr.ind = TRUE)
labels <- paste(colnames(correlation)[pairs[, "row"]], "--", colnames(correlation)[pairs[, "col"]])
result <- data.frame(quantity = c(paste0("cor: ", labels), paste0("pcor: ", labels)),
                     kind = rep(c("cor", "pcor"), each = nrow(pairs)),
                     value = c(correlation[pairs], partial[pairs]))
result
