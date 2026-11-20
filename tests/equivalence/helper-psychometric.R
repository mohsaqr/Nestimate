# Inputs and reference computations shared by the psychometric-network equivalence tests (correlation, partial
# correlation, graphical lasso, Ising, mixed graphical models, network comparison, predictability), written in
# base R. Nothing here calls Nestimate.

#' Continuous data with a random correlation structure
#'
#' @param n,p Rows and variables.
#' @param rho Strength of the random off-diagonal structure.
#' @param seed Seed of the draw.
#' @return A data frame with columns V1..Vp.
psy_continuous <- function(n, p, rho, seed) {
  set.seed(seed)
  sigma <- diag(p) + rho * matrix(stats::rnorm(p * p, sd = 0.3), p, p)
  sigma <- (sigma + t(sigma)) / 2
  diag(sigma) <- 1
  sigma <- sigma + max(0, 0.05 - min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values)) * diag(p)
  draws <- matrix(stats::rnorm(n * p), n, p) %*% chol(sigma)
  colnames(draws) <- paste0("V", seq_len(p))
  as.data.frame(draws)
}

#' Binary data from thresholded correlated normals (a chain of correlations rho)
#'
#' @param n,p Rows and variables.
#' @param rho Correlation of neighbouring variables.
#' @param seed Seed of the draw.
#' @return A data frame of 0/1 integers, columns V1..Vp.
psy_binary <- function(n, p, rho, seed) {
  set.seed(seed)
  sigma <- diag(p)
  sigma[cbind(seq_len(p - 1L), seq_len(p - 1L) + 1L)] <- rho
  sigma[cbind(seq_len(p - 1L) + 1L, seq_len(p - 1L))] <- rho
  draws <- matrix(stats::rnorm(n * p), n, p) %*% chol(sigma)
  binary <- as.data.frame(ifelse(draws > 0, 1L, 0L))
  colnames(binary) <- paste0("V", seq_len(p))
  binary
}

#' A symmetric matrix as one row per unordered pair (i < j in the order of its names)
#'
#' @param ... Named symmetric matrices with the same dimnames; each becomes a column.
#' @return A data frame with `from`, `to` and one column per matrix.
psy_pairs <- function(...) {
  matrices <- list(...)
  nodes <- rownames(matrices[[1L]])
  index <- which(upper.tri(matrices[[1L]]), arr.ind = TRUE)
  index <- index[order(index[, "col"], index[, "row"]), , drop = FALSE]
  pairs <- data.frame(from = nodes[index[, "row"]], to = nodes[index[, "col"]])
  pairs[names(matrices)] <- lapply(matrices, \(m) m[index])
  pairs
}

#' Partial correlations of a precision matrix: -K_ij / sqrt(K_ii K_jj), zero diagonal
#'
#' @param precision A precision matrix.
#' @return The matrix of partial correlations.
psy_partial <- function(precision) {
  scale <- sqrt(diag(precision))
  partial <- -precision / outer(scale, scale)
  diag(partial) <- 0
  partial
}

#' The extended BIC of the graphical lasso along a lambda path (Foygel and Drton 2010, as qgraph::EBICglasso)
#'
#' EBIC = -2 loglik + E log(n) + 4 E gamma log(p), loglik = n / 2 (log det K - tr(S K)), E the number of edges.
#' Each K is glasso::glasso() run to thr = 1e-12 without a penalty on the diagonal. At lambda >= max|S_ij| the
#' solution is known in closed form, the diagonal K = diag(1 / S_ii) (Friedman, Hastie and Tibshirani 2008), and
#' that is used there: the coordinate descent stops one step short at that boundary and leaves one edge.
#'
#' @param cor_matrix A correlation matrix.
#' @param n Number of observations.
#' @param gamma EBIC hyperparameter.
#' @param lambdas The penalties.
#' @return A numeric vector, one EBIC per lambda.
psy_ebic_path <- function(cor_matrix, n, gamma, lambdas) {
  p <- ncol(cor_matrix)
  boundary <- max(abs(cor_matrix[upper.tri(cor_matrix)]))
  vapply(lambdas, \(lambda) {
    # a path built as exp(seq(log(.))) can return lambda_max a rounding error below max|S_ij|
    precision <- if (lambda >= boundary * (1 - 1e-12)) diag(1 / diag(cor_matrix)) else
      glasso::glasso(cor_matrix, rho = lambda, penalize.diagonal = FALSE, thr = 1e-12)$wi
    edges <- sum(precision[upper.tri(precision)] != 0)
    loglik <- n / 2 * (as.numeric(determinant(precision)$modulus) - sum(diag(cor_matrix %*% precision)))
    -2 * loglik + edges * log(n) + 4 * edges * gamma * log(p)
  }, numeric(1))
}

#' A simulated panel: a stationary VAR(1) per person around a person-specific mean, occasions numbered by day and beep
#'
#' @param n_id,days,beeps People, days per person and occasions per day.
#' @param p Number of variables (V1..Vp).
#' @param seed Seed of the draw.
#' @return A data frame `id`, `day`, `beep`, V1..Vp.
psy_panel <- function(n_id, days, beeps, p, seed) {
  set.seed(seed)
  transition <- diag(0.3, p) + matrix(stats::runif(p * p, -0.1, 0.1), p, p)
  people <- lapply(seq_len(n_id), \(person) {
    means <- stats::rnorm(p, 0, 0.6)
    shocks <- matrix(stats::rnorm(days * beeps * p), ncol = p)
    states <- Reduce(\(previous, shock) as.vector(transition %*% previous) + shock,
                     asplit(shocks[-1L, , drop = FALSE], 1L), accumulate = TRUE, init = shocks[1L, ])
    series <- sweep(do.call(rbind, states), 2L, means, "+")
    colnames(series) <- paste0("V", seq_len(p))
    data.frame(id = person, day = rep(seq_len(days), each = beeps), beep = rep(seq_len(beeps), days), series)
  })
  do.call(rbind, people)
}

#' mlVAR's fixed temporal effects as one row per outcome and predictor: estimate, standard error, p-value
#'
#' @param fit An mlVAR::mlVAR() result with one lag.
#' @return A data frame `outcome`, `predictor`, `beta`, `se`, `p`.
psy_mlvar_coefs <- function(fit) {
  beta <- fit$results$Beta
  vars <- rownames(beta$mean)
  cells <- expand.grid(predictor = vars, outcome = vars, stringsAsFactors = FALSE)
  index <- cbind(match(cells$outcome, vars), match(cells$predictor, vars), 1L)
  data.frame(outcome = cells$outcome, predictor = cells$predictor,
             beta = beta$mean[index], se = beta$SE[index], p = beta$P[index])
}
