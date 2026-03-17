# Tests for graphical_var.R: graphical_var, ml_graphical_var

# ============================================
# Single-subject graphical VAR
# ============================================

test_that("graphical_var produces correct output structure", {
  skip_if_not_installed("graphicalVAR")

  set.seed(42)
  n <- 200
  d <- data.frame(V1 = rnorm(n), V2 = rnorm(n), V3 = rnorm(n))
  # Add autoregressive signal
  for (i in 2:n) {
    d$V1[i] <- 0.3 * d$V1[i - 1] + rnorm(1, sd = 0.5)
    d$V2[i] <- 0.2 * d$V2[i - 1] + 0.1 * d$V1[i - 1] + rnorm(1, sd = 0.5)
  }

  res <- graphical_var(d, vars = c("V1", "V2", "V3"), n_lambda = 5)

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
})

test_that("graphical_var beta has named rows/cols", {
  skip_if_not_installed("graphicalVAR")

  set.seed(1)
  d <- data.frame(X = rnorm(150), Y = rnorm(150))
  res <- graphical_var(d, vars = c("X", "Y"), n_lambda = 5)

  expect_equal(rownames(res$beta), c("X", "Y"))
  expect_equal(colnames(res$beta), c("X", "Y"))
  expect_equal(rownames(res$PCC), c("X", "Y"))
})

test_that("graphical_var print method works", {
  skip_if_not_installed("graphicalVAR")

  set.seed(1)
  d <- data.frame(V1 = rnorm(150), V2 = rnorm(150))
  res <- graphical_var(d, vars = c("V1", "V2"), n_lambda = 5)
  expect_output(print(res), "Graphical VAR")
})

test_that("graphical_var summary method works", {
  skip_if_not_installed("graphicalVAR")

  set.seed(1)
  d <- data.frame(V1 = rnorm(150), V2 = rnorm(150))
  res <- graphical_var(d, vars = c("V1", "V2"), n_lambda = 5)
  expect_output(summary(res), "Temporal Network")
})

# ============================================
# Multilevel graphical VAR
# ============================================

test_that("ml_graphical_var produces correct output structure", {
  skip_if_not_installed("graphicalVAR")

  set.seed(42)
  data_panel <- do.call(rbind, lapply(1:5, function(id) {
    data.frame(id = id, V1 = rnorm(40), V2 = rnorm(40), V3 = rnorm(40))
  }))

  res <- ml_graphical_var(data_panel, vars = c("V1", "V2", "V3"),
                          id = "id")

  expect_s3_class(res, "ml_graphical_var_result")
  expect_equal(dim(res$temporal), c(3, 3))
  expect_equal(dim(res$PCC), c(3, 3))
  expect_equal(dim(res$between), c(3, 3))
  expect_equal(res$labels, c("V1", "V2", "V3"))
  expect_equal(res$n_subjects, 5)
  expect_equal(length(res$ids), 5)
  # Subject networks should be lists
  expect_true(is.list(res$subject_PCC) || is.null(res$subject_PCC))
})

test_that("ml_graphical_var print method works", {
  skip_if_not_installed("graphicalVAR")

  set.seed(1)
  data_panel <- do.call(rbind, lapply(1:3, function(id) {
    data.frame(id = id, V1 = rnorm(40), V2 = rnorm(40))
  }))

  res <- ml_graphical_var(data_panel, vars = c("V1", "V2"), id = "id")
  expect_output(print(res), "Panel Graphical VAR")
})

test_that("ml_graphical_var summary method works", {
  skip_if_not_installed("graphicalVAR")

  set.seed(1)
  data_panel <- do.call(rbind, lapply(1:3, function(id) {
    data.frame(id = id, V1 = rnorm(40), V2 = rnorm(40))
  }))

  res <- ml_graphical_var(data_panel, vars = c("V1", "V2"), id = "id")
  expect_output(summary(res), "Temporal Network")
})
