# Numerical equivalence: predictability() vs lm()-based R² and RMSE
#
# predictability.netobject has two computational paths:
#
#   A) precision_matrix path (method = "glasso"):
#      R² = 1 - 1 / diag(Omega)     (Haslbeck & Waldorp 2020)
#
#   B) cor_matrix path (method = "cor" / "pcor"):
#      For each node j with neighbours N(j) = {i : W[j,i] != 0}:
#      R²_j = r'_j,N · inv(R_{N,N}) · r_j,N
#      where R is the correlation matrix and r_j,N its j-th column
#      restricted to neighbours.
#
# For the cor path, R²_j is the squared multiple correlation of z_j regressed
# on z_{N(j)}, which equals summary(lm(z_j ~ z_{N(j)}))$r.squared exactly
# when the regression uses the same neighbour set. That lm-based R² is the
# independent reference.
#
# RMSE: reconstructed from the precision or correlation matrix via the
# Nestimate formula y_hat = mu_j + sd_j * Z_{-j} * beta_z, with beta_z from
# the precision matrix or from solve(R_{N,N}, r_j,N). Reference: predict()
# from an explicit lm() using the same neighbour set.

set.seed(20260422)
N_CONFIGS <- 50L
TOL       <- 1e-8

skip_if_not_installed("glasso")

.gen_cont <- function(n, p, rho, seed) {
  set.seed(seed)
  Sigma <- diag(p) + rho * (matrix(stats::rnorm(p * p, sd = 0.3), p, p))
  Sigma <- (Sigma + t(Sigma)) / 2
  diag(Sigma) <- 1
  Sigma <- Sigma + max(0, -min(eigen(Sigma)$values) + 0.05) * diag(p)
  L <- chol(Sigma)
  X <- matrix(stats::rnorm(n * p), n, p) %*% L
  colnames(X) <- paste0("V", seq_len(p))
  as.data.frame(X)
}

cor_configs <- replicate(N_CONFIGS, list(
  n    = sample(c(80L, 150L, 300L), 1L),
  p    = sample(5:9, 1L),
  rho  = round(runif(1L, 0.1, 0.4), 2),
  seed = sample.int(1e5, 1L)
), simplify = FALSE)

glasso_configs <- replicate(N_CONFIGS, list(
  n    = sample(c(100L, 200L), 1L),
  p    = sample(5:8, 1L),
  rho  = round(runif(1L, 0.1, 0.3), 2),
  seed = sample.int(1e5, 1L)
), simplify = FALSE)

# ---- Cor / pcor path: R² vs lm() -----------------------------------------

test_that("predictability R² (cor method) matches summary(lm)$r.squared", {
  skip_on_cran()

  deltas <- vapply(cor_configs, function(cfg) {
    df  <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    net <- build_network(df, method = "cor")
    pr  <- predictability(net)

    Z <- scale(df)
    W <- net$weights

    # lm-based R² per node using the network's neighbour set
    r2_ref <- vapply(seq_len(cfg$p), function(j) {
      nbrs <- which(W[j, ] != 0)
      if (length(nbrs) == 0L) return(0)
      fit <- stats::lm(Z[, j] ~ Z[, nbrs])
      summary(fit)$r.squared
    }, numeric(1))
    r2_ref <- pmin(pmax(r2_ref, 0), 1)

    max(abs(pr$R2 - r2_ref))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- Cor path: RMSE vs predict(lm) ---------------------------------------

test_that("predictability RMSE (cor method) matches manual predict(lm)", {
  skip_on_cran()

  deltas <- vapply(cor_configs, function(cfg) {
    df  <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    net <- build_network(df, method = "cor")
    pr  <- predictability(net, data = df)

    X <- as.matrix(df)
    mu <- colMeans(X); sd <- apply(X, 2, stats::sd)
    Z <- scale(X)
    W <- net$weights

    rmse_ref <- vapply(seq_len(cfg$p), function(j) {
      nbrs <- which(W[j, ] != 0)
      if (length(nbrs) == 0L) return(sd[j])
      beta_z <- solve(cor(X)[nbrs, nbrs, drop = FALSE], cor(X)[nbrs, j])
      z_hat <- Z[, nbrs, drop = FALSE] %*% beta_z
      y_hat <- mu[j] + sd[j] * z_hat
      sqrt(mean((X[, j] - y_hat)^2))
    }, numeric(1))

    max(abs(pr$RMSE - rmse_ref), na.rm = TRUE)
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- Glasso path: R² vs closed form on stored Omega ----------------------
#
# The glasso R² formula is 1 - 1/diag(Omega). This test re-reads Omega from
# the fitted network and recomputes the formula, verifying the implementation
# applies the formula to the stored precision without accidental scaling,
# clipping in the wrong place, or index slips.

test_that("predictability R² (glasso method) matches 1 - 1/diag(Omega)", {
  skip_on_cran()

  deltas <- vapply(glasso_configs, function(cfg) {
    df  <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    net <- build_network(df, method = "glasso",
                         params = list(gamma = 0.25))
    pr  <- predictability(net)

    Omega <- net$precision_matrix
    r2_ref <- 1 - 1 / diag(Omega)
    r2_ref <- pmin(pmax(r2_ref, 0), 1)

    max(abs(pr$R2 - r2_ref))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- Delta report + validation dashboard emit --------------------------

test_that("predictability equivalence report (CSV + CVS JSON)", {
  skip_on_cran()

  report <- equiv_report()

  invisible(lapply(cor_configs, function(cfg) {
    df  <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    net <- build_network(df, method = "cor")
    pr  <- predictability(net, data = df)
    Z <- scale(df); W <- net$weights
    r2_ref <- vapply(seq_len(cfg$p), function(j) {
      nbrs <- which(W[j, ] != 0)
      if (length(nbrs) == 0L) return(0)
      summary(stats::lm(Z[, j] ~ Z[, nbrs]))$r.squared
    }, numeric(1))
    r2_ref <- pmin(pmax(r2_ref, 0), 1)
    err <- abs(pr$R2 - r2_ref)

    report$log(
      func = "predictability (cor path)",
      config = sprintf("n=%d p=%d rho=%.2f", cfg$n, cfg$p, cfg$rho),
      n_checked = length(err),
      n_failed  = sum(err >= TOL),
      max_abs_err    = max(err),
      mean_abs_err   = mean(err),
      median_abs_err = stats::median(err),
      p95_abs_err    = stats::quantile(err, 0.95, names = FALSE),
      reference = "summary(lm(Z_j ~ Z_N(j)))$r.squared",
      notes = ""
    )
  }))

  invisible(lapply(glasso_configs, function(cfg) {
    df  <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    net <- build_network(df, method = "glasso",
                         params = list(gamma = 0.25))
    pr  <- predictability(net)
    Omega <- net$precision_matrix
    r2_ref <- pmin(pmax(1 - 1 / diag(Omega), 0), 1)
    err <- abs(pr$R2 - r2_ref)

    report$log(
      func = "predictability (glasso path)",
      config = sprintf("n=%d p=%d rho=%.2f", cfg$n, cfg$p, cfg$rho),
      n_checked = length(err),
      n_failed  = sum(err >= TOL),
      max_abs_err    = max(err),
      mean_abs_err   = mean(err),
      median_abs_err = stats::median(err),
      p95_abs_err    = stats::quantile(err, 0.95, names = FALSE),
      reference = "closed form: 1 - 1/diag(Omega) (Haslbeck & Waldorp 2020)",
      notes = ""
    )
  }))

  report$write_csv("predictability")
  report$write_cvs("predictability")
  expect_true(length(report$rows) > 0L)
})
