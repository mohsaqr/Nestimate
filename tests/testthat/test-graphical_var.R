# Tests for graphical_var.R — own implementation (no graphicalVAR package)

# ============================================
# Helper
# ============================================

.make_gvar_data <- function(n = 200, d = 3, seed = 42) {
  set.seed(seed)
  data <- data.frame(matrix(rnorm(n * d), n, d))
  names(data) <- paste0("V", seq_len(d))
  # Add temporal structure: V1 -> V2
  for (t in 2:n) {
    data[t, "V2"] <- 0.4 * data[t - 1, "V1"] + rnorm(1, sd = 0.5)
  }
  data
}


# ============================================
# Basic functionality
# ============================================

test_that("graphical_var produces correct output structure", {
  data <- .make_gvar_data()
  res <- graphical_var(data, vars = c("V1", "V2", "V3"), n_lambda = 5)

  expect_s3_class(res, "gvar_result")
  expect_equal(dim(res$beta), c(3, 3))
  expect_equal(dim(res$kappa), c(3, 3))
  expect_equal(dim(res$PCC), c(3, 3))
  expect_equal(dim(res$PDC), c(3, 3))
  expect_equal(res$labels, c("V1", "V2", "V3"))
  expect_true(is.numeric(res$n_obs))
  expect_true(is.numeric(res$EBIC))
  expect_true(is.numeric(res$lambda_beta))
  expect_true(is.numeric(res$lambda_kappa))
  expect_equal(res$gamma, 0.5)
  # Aliases
  expect_identical(res$temporal, res$beta)
  expect_identical(res$contemporaneous, res$PCC)
})

test_that("graphical_var recovers known temporal structure", {
  set.seed(99)
  n <- 300
  data <- data.frame(V1 = rnorm(n), V2 = rnorm(n), V3 = rnorm(n))
  for (t in 2:n) {
    data$V2[t] <- 0.5 * data$V1[t - 1] + rnorm(1, sd = 0.5)
  }
  res <- graphical_var(data, vars = c("V1", "V2", "V3"), n_lambda = 8)
  # V1 -> V2 should be the strongest temporal edge (beta[outcome, predictor])
  expect_true(abs(res$beta["V2", "V1"]) > 0.2)
})

test_that("graphical_var beta has named rows/cols", {
  set.seed(1)
  d <- data.frame(X = rnorm(150), Y = rnorm(150))
  res <- graphical_var(d, vars = c("X", "Y"), n_lambda = 5)

  expect_equal(rownames(res$beta), c("X", "Y"))
  expect_equal(colnames(res$beta), c("X", "Y"))
  expect_equal(rownames(res$PCC), c("X", "Y"))
})

test_that("graphical_var kappa is symmetric", {
  data <- .make_gvar_data()
  res <- graphical_var(data, vars = c("V1", "V2", "V3"), n_lambda = 5)
  expect_true(isSymmetric(res$kappa))
})

test_that("graphical_var PCC has zero diagonal", {
  data <- .make_gvar_data()
  res <- graphical_var(data, vars = c("V1", "V2", "V3"), n_lambda = 5)
  expect_equal(unname(diag(res$PCC)), c(0, 0, 0))
})


# ============================================
# Panel data (id/day/beep)
# ============================================

test_that("graphical_var handles panel data with id", {
  set.seed(42)
  data_panel <- do.call(rbind, lapply(1:3, function(id) {
    data.frame(subj = id, V1 = rnorm(60), V2 = rnorm(60))
  }))
  res <- graphical_var(data_panel, vars = c("V1", "V2"),
                       id = "subj", n_lambda = 5)
  expect_s3_class(res, "gvar_result")
})

test_that("graphical_var handles panel data with id/day/beep", {
  set.seed(42)
  data_panel <- do.call(rbind, lapply(1:3, function(id) {
    data.frame(
      subj = id,
      day  = rep(1:4, each = 10),
      beep = rep(1:10, 4),
      V1 = rnorm(40), V2 = rnorm(40)
    )
  }))
  res <- graphical_var(data_panel, vars = c("V1", "V2"),
                       id = "subj", day = "day", beep = "beep",
                       n_lambda = 5)
  expect_s3_class(res, "gvar_result")
})


# ============================================
# Input validation
# ============================================

test_that("graphical_var errors on too few variables", {
  d <- data.frame(V1 = rnorm(50))
  expect_error(graphical_var(d, vars = "V1"), "length\\(vars\\) >= 2")
})

test_that("graphical_var errors on too few observations", {
  d <- data.frame(V1 = rnorm(3), V2 = rnorm(3))
  expect_error(graphical_var(d, vars = c("V1", "V2"), n_lambda = 5),
               "Fewer than 3|Too few")
})

test_that("graphical_var accepts matrix input", {
  set.seed(1)
  mat <- matrix(rnorm(200), 100, 2)
  colnames(mat) <- c("A", "B")
  res <- graphical_var(as.data.frame(mat), vars = c("A", "B"), n_lambda = 5)
  expect_s3_class(res, "gvar_result")
})


# ============================================
# S3 methods
# ============================================

test_that("print.gvar_result works", {
  data <- .make_gvar_data(n = 150)
  res <- graphical_var(data, vars = c("V1", "V2", "V3"), n_lambda = 5)
  expect_output(print(res), "Graphical VAR Result")
  expect_output(print(res), "Variables:")
  expect_output(print(res), "EBIC:")
})

test_that("summary.gvar_result works", {
  data <- .make_gvar_data(n = 150)
  res <- graphical_var(data, vars = c("V1", "V2", "V3"), n_lambda = 5)
  expect_output(summary(res), "Temporal Network")
  expect_output(summary(res), "Contemporaneous Network")
  expect_output(summary(res), "Partial Directed Correlations")
})


# ============================================
# Internal functions
# ============================================

test_that(".gvar_compute_pcc returns symmetric matrix with zero diagonal", {
  set.seed(1)
  A <- matrix(runif(9), 3, 3)
  kappa <- t(A) %*% A + diag(3)
  result <- Nestimate:::.gvar_compute_pcc(kappa)
  expect_true(isSymmetric(result))
  expect_true(all(diag(result) == 0))
  expect_equal(dim(result), c(3, 3))
})

test_that(".gvar_compute_pdc returns matrix with correct dimensions", {
  set.seed(2)
  p <- 3L
  A <- matrix(runif(p * p), p, p)
  kappa <- t(A) %*% A + diag(p)
  beta <- matrix(rnorm(p * p, sd = 0.3), p, p)
  result <- Nestimate:::.gvar_compute_pdc(beta, kappa)
  expect_equal(dim(result), c(p, p))
  expect_true(is.numeric(result))
})

test_that(".gvar_compute_pdc falls back to identity when kappa is singular", {
  kappa <- matrix(c(1, 1, 1, 1), 2, 2)
  beta <- matrix(c(0.3, 0.1, 0.2, 0.4), 2, 2)
  result <- Nestimate:::.gvar_compute_pdc(beta, kappa)
  expect_equal(dim(result), c(2, 2))
  expect_true(is.numeric(result))
})

test_that(".gvar_build_lag_pairs errors on too few rows", {
  d <- data.frame(V1 = 1, V2 = 2)
  expect_error(Nestimate:::.gvar_build_lag_pairs(d, c("V1", "V2"), NULL, NULL, NULL),
               "Fewer than 3")
})

test_that(".gvar_multivariate_lasso returns sparse beta for large lambda", {
  set.seed(1)
  X <- matrix(rnorm(100 * 3), 100, 3)
  Y <- matrix(rnorm(100 * 3), 100, 3)
  beta <- Nestimate:::.gvar_multivariate_lasso(X, Y, lambda = 10)
  expect_equal(dim(beta), c(3, 3))
  expect_true(all(beta == 0))  # large lambda -> all zeros
})

test_that(".gvar_multivariate_lasso returns non-sparse beta for small lambda", {
  set.seed(1)
  X <- matrix(rnorm(200 * 3), 200, 3)
  Y <- X %*% matrix(c(0.5, 0.1, 0, 0, 0.3, 0.2, 0, 0, 0.4), 3, 3) +
    matrix(rnorm(200 * 3, sd = 0.1), 200, 3)
  beta <- Nestimate:::.gvar_multivariate_lasso(X, Y, lambda = 0.001)
  expect_true(sum(beta != 0) > 0)
})
