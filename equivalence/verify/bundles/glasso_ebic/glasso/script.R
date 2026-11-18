library(glasso)
strategies <- read.csv("../data/srl_strategies.csv")
S <- cor(strategies)
n <- nrow(strategies)
p <- ncol(S)
largest <- max(abs(S[upper.tri(S)]))
penalties <- exp(seq(log(0.01 * largest), log(largest), length.out = 100))
path <- glassopath(S, rholist = penalties, thr = 1e-10, maxit = 1e5, penalize.diagonal = FALSE,
                   trace = 0)
ebic <- apply(path$wi, 3, function(K) {
  edges <- sum(K[lower.tri(K)] != 0)
  loglik <- n / 2 * (log(det(K)) - sum(diag(K %*% S)))
  -2 * loglik + edges * log(n) + 4 * edges * 0.5 * log(p)
})
penalties[which.min(ebic)]
network <- -cov2cor(path$wi[, , which.min(ebic)])
dimnames(network) <- dimnames(S)
pairs <- which(upper.tri(network) & network != 0, arr.ind = TRUE)
result <- data.frame(quantity = paste(colnames(network)[pairs[, "row"]], "--",
                                      colnames(network)[pairs[, "col"]]),
                     kind = "partial correlation", value = network[pairs])
result
