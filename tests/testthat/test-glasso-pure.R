# Pure-R graphical lasso: correctness is certified against the estimand itself
# (the stationarity / KKT conditions of the convex objective), so these tests
# need NO reference solver -- glasso does not have to be installed. Numerical
# equivalence to glasso::glasso is checked separately in
# local_testing_and_equivalence/ (gated, glasso installed locally).

# --- synthetic correlation matrices ------------------------------------------
ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
compound <- function(p, rho) {
  m <- matrix(rho, p, p); diag(m) <- 1; m
}

test_that(".glasso_fit returns the unique optimum (KKT < 1e-8)", {
  cases <- list(
    ar1(6, 0.5), ar1(8, 0.6), compound(5, 0.4), compound(7, -0.2)
  )
  rhos <- c(0.05, 0.1, 0.2, 0.35)
  for (S in cases) {
    for (rho in rhos) {
      fit <- .glasso_fit(S, rho)
      v <- .glasso_kkt_violation(fit$wi, S, rho)
      expect_true(v < 1e-8,
                  info = sprintf("KKT violation %.2e at rho=%g (p=%d)",
                                 v, rho, ncol(S)))
    }
  }
})

test_that(".glasso_kkt_violation flags a non-optimal precision matrix", {
  S <- ar1(5, 0.5)
  fit <- .glasso_fit(S, 0.1)
  bad <- fit$wi
  bad[1, 2] <- bad[2, 1] <- bad[1, 2] + 0.3
  expect_gt(.glasso_kkt_violation(bad, S, 0.1),
            .glasso_kkt_violation(fit$wi, S, 0.1))
})

test_that(".glasso_fit warm start reaches the same optimum as cold", {
  S <- ar1(7, 0.55)
  cold <- .glasso_fit(S, 0.15)
  warm <- .glasso_fit(S, 0.15, w_init = S, beta_init = matrix(0, 7, 7))
  expect_equal(cold$wi, warm$wi, tolerance = 1e-8)
})

test_that(".glasso_fit zero constraint forces the named edges to exactly zero", {
  S <- compound(6, 0.3)
  zero <- rbind(c(1L, 4L), c(2L, 5L), c(3L, 6L))   # (i, j) pairs, i < j
  fit <- .glasso_fit(S, rho = 0, zero = zero)
  expect_identical(fit$wi[1, 4], 0)
  expect_identical(fit$wi[4, 1], 0)
  expect_identical(fit$wi[2, 5], 0)
  expect_identical(fit$wi[3, 6], 0)
  # an unconstrained off-diagonal is generally non-zero
  expect_true(abs(fit$wi[1, 2]) > 1e-8)
})

test_that(".glassopath_fit returns a wi array, each slice near-optimal", {
  S <- ar1(6, 0.5)
  rholist <- exp(seq(log(0.4), log(0.02), length.out = 12))
  # Default tolerance is glasso's own thr = 1e-4 (loose-by-design for the
  # bootstrap path): slices are valid glasso fits to ~1e-4, not machine
  # precision. Tightening the tolerance drives every slice to the optimum.
  gp <- .glassopath_fit(S, rholist)
  expect_equal(dim(gp$wi), c(6L, 6L, 12L))
  worst_loose <- max(vapply(seq_along(rholist), function(k)
    .glasso_kkt_violation(gp$wi[, , k], S, rholist[k]), numeric(1)))
  expect_true(worst_loose < 1e-3,
              info = sprintf("worst loose-path KKT %.2e", worst_loose))

  gp_tight <- .glassopath_fit(S, rholist, tol_outer = 1e-8, tol_inner = 1e-10)
  worst_tight <- max(vapply(seq_along(rholist), function(k)
    .glasso_kkt_violation(gp_tight$wi[, , k], S, rholist[k]), numeric(1)))
  expect_true(worst_tight < 1e-8,
              info = sprintf("worst tight-path KKT %.2e", worst_tight))
})

# --- end-to-end through build_network (no glasso dependency) ------------------
test_that("build_network(method='glasso') returns a KKT-optimal network", {
  set.seed(1)
  p <- 8
  Sigma <- ar1(p, 0.5)
  L <- chol(Sigma)
  X <- matrix(stats::rnorm(400 * p), 400, p) %*% L
  colnames(X) <- paste0("V", seq_len(p))
  df <- as.data.frame(X)

  fit <- build_network(df, method = "glasso", gamma = 0.5)
  expect_s3_class(fit, "netobject")

  est <- .estimator_glasso(df, gamma = 0.5)
  v <- .glasso_kkt_violation(est$precision_matrix, est$cor_matrix,
                             est$lambda_selected)
  expect_true(v < 1e-7,
              info = sprintf("end-to-end KKT violation %.2e", v))
})
