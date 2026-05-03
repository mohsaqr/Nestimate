# Equivalence test_that() blocks extracted from
# tests/testthat/test-estimator_ising.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ---- Tests for Ising Model Estimator ----

# ---- Synthetic data generators ----

#' Generate chain-structure binary data with conditional dependencies
#'
#' X1 -> X2 -> X3 -> ... -> Xp (chain)
#' Each X[j] depends on X[j-1] via logistic model with coupling strength.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Number of variables.
#' @param seed Integer. Random seed.
#' @param coupling Numeric. Logistic regression coupling strength between
#'   adjacent nodes.
#' @return Data frame of 0/1 variables.
#' @noRd
.make_ising_data <- function(n = 300, p = 5, seed = 42, coupling = 1.5) {
  set.seed(seed)
  mat <- matrix(0L, nrow = n, ncol = p)
  colnames(mat) <- paste0("V", seq_len(p))

  # First node: independent with ~50% probability

  mat[, 1] <- rbinom(n, 1, 0.5)

  # Chain: each subsequent node depends on previous
  for (j in seq(2, p)) {
    eta <- coupling * mat[, j - 1] - coupling / 2
    prob <- 1 / (1 + exp(-eta))
    mat[, j] <- rbinom(n, 1, prob)
  }

  as.data.frame(mat)
}


#' Generate independent binary data (null model)
#'
#' All columns are independent Bernoulli(0.5).
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Number of variables.
#' @param seed Integer. Random seed.
#' @return Data frame of 0/1 variables.
#' @noRd
.make_ising_null <- function(n = 300, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rbinom(n * p, 1, 0.5), nrow = n, ncol = p)
  colnames(mat) <- paste0("V", seq_len(p))
  as.data.frame(mat)
}


# ---- Input validation tests ----

test_that("Ising: exact match with IsingFit::IsingFit (AND rule)", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("IsingFit")

  df <- .make_ising_data(300, 5, seed = 42, coupling = 1.5)
  our <- .estimator_ising(df, gamma = 0.25, rule = "AND")
  ref <- IsingFit::IsingFit(df, gamma = 0.25, AND = TRUE,
                             progressbar = FALSE, plot = FALSE)

  expect_equal(our$matrix, ref$weiadj, tolerance = 1e-10)
  expect_equal(unname(our$thresholds), unname(ref$thresholds), tolerance = 1e-10)
})

test_that("Ising: exact match with IsingFit::IsingFit (OR rule)", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("IsingFit")

  df <- .make_ising_data(300, 5, seed = 77, coupling = 1.5)
  our <- .estimator_ising(df, gamma = 0.25, rule = "OR")
  ref <- IsingFit::IsingFit(df, gamma = 0.25, AND = FALSE,
                             progressbar = FALSE, plot = FALSE)

  expect_equal(our$matrix, ref$weiadj, tolerance = 1e-10)
  expect_equal(unname(our$thresholds), unname(ref$thresholds), tolerance = 1e-10)
})

test_that("Ising: exact match with IsingFit across 20 random configs", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("IsingFit")

  max_diffs <- vapply(seq_len(20), function(seed) {
    set.seed(seed)
    n <- sample(100:400, 1)
    p <- sample(3:6, 1)
    gamma <- sample(c(0, 0.25, 0.5), 1)
    rule <- sample(c("AND", "OR"), 1)

    mat <- matrix(0L, nrow = n, ncol = p)
    colnames(mat) <- paste0("V", seq_len(p))
    mat[, 1] <- rbinom(n, 1, 0.5)
    for (j in 2:p) {
      eta <- 1.5 * mat[, j - 1] - 0.75
      mat[, j] <- rbinom(n, 1, 1 / (1 + exp(-eta)))
    }
    df <- as.data.frame(mat)

    our <- .estimator_ising(df, gamma = gamma, rule = rule)
    ref <- IsingFit::IsingFit(df, gamma = gamma, AND = (rule == "AND"),
                               progressbar = FALSE, plot = FALSE)

    max(abs(our$matrix - ref$weiadj))
  }, numeric(1))

  expect_true(all(max_diffs == 0))
})


