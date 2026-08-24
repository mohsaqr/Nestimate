# Delegated EBIC-glasso: since 0.8.6 the solver lives in psychnets
# (Imports), which certifies its own kernel against the KKT conditions.
# What Nestimate still owns -- and what these tests pin -- is the public
# surface: the estimator's return contract, the branch mapping
# (refit -> "unregularized"), and that the returned precision matrix is
# KKT-optimal for the selected lambda, checked with psychnets::glasso_kkt.
# (The kernel-level tests that lived here moved to psychnets with the math.)

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

.glasso_test_df <- function(p = 8L, n = 400L, seed = 1L) {
  set.seed(seed)
  L <- chol(ar1(p, 0.5))
  X <- matrix(stats::rnorm(n * p), n, p) %*% L
  colnames(X) <- paste0("V", seq_len(p))
  as.data.frame(X)
}

test_that("build_network(method='glasso') returns a KKT-optimal network", {
  df <- .glasso_test_df()

  fit <- build_network(df, method = "glasso", gamma = 0.5)
  expect_s3_class(fit, "netobject")

  est <- .estimator_glasso(df, gamma = 0.5)
  v <- psychnets::glasso_kkt(est$precision_matrix, est$cor_matrix,
                             est$lambda_selected)
  expect_true(v < 1e-7,
              info = sprintf("end-to-end KKT violation %.2e", v))
})

test_that(".estimator_glasso keeps the full pre-delegation field contract", {
  est <- .estimator_glasso(.glasso_test_df(), gamma = 0.5)
  expect_named(est, c("matrix", "nodes", "directed", "cleaned_data",
                      "precision_matrix", "cor_matrix", "lambda_selected",
                      "ebic_selected", "lambda_path", "ebic_path",
                      "gamma", "n", "p"))
  expect_false(est$directed)
  expect_identical(dim(est$matrix), dim(est$precision_matrix))
  expect_length(est$ebic_path, length(est$lambda_path))
  # the selected lambda is on the path and matches the minimal EBIC
  expect_true(any(abs(est$lambda_path - est$lambda_selected) < 1e-12))
  expect_equal(est$ebic_selected, min(est$ebic_path))
})

test_that("refit = TRUE maps to the unregularized-refit branch and changes weights", {
  df <- .glasso_test_df()
  base  <- .estimator_glasso(df, gamma = 0.5)
  refit <- .estimator_glasso(df, gamma = 0.5, refit = TRUE)
  # identical support, different (unpenalized) magnitudes
  expect_identical(base$matrix != 0, refit$matrix != 0)
  expect_gt(max(abs(base$matrix - refit$matrix)), 1e-6)
})

test_that("cor and pcor delegation keeps the return contract and symmetry", {
  df <- .glasso_test_df(p = 6L, seed = 2L)
  ec <- .estimator_cor(df)
  expect_named(ec, c("matrix", "nodes", "directed", "cleaned_data",
                     "cor_matrix", "n", "p"))
  expect_identical(unname(diag(ec$matrix)), rep(0, 6L))
  expect_identical(ec$matrix, t(ec$matrix))

  ep <- .estimator_pcor(df)
  expect_named(ep, c("matrix", "nodes", "directed", "cleaned_data",
                     "precision_matrix", "cor_matrix", "n", "p"))
  expect_identical(ep$matrix, t(ep$matrix))
  # precision inverts the correlation matrix (delegated solve)
  expect_equal(ep$precision_matrix %*% ep$cor_matrix, diag(6L),
               tolerance = 1e-8, ignore_attr = TRUE)
})
