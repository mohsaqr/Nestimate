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

test_that("Ising: non-binary data raises error", {
  skip_if_not_installed("glmnet")
  df <- data.frame(X1 = c(0, 1, 2, 0, 1), X2 = c(1, 0, 1, 0, 1))
  expect_error(.prepare_ising_input(df), "binary.*0/1")
})

test_that("Ising: non-numeric columns are dropped", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  df$Name <- letters[seq_len(50)]
  prepared <- .prepare_ising_input(df)
  expect_equal(prepared$p, 3L)
  expect_false("Name" %in% prepared$nodes)
})

test_that("Ising: id_col is excluded", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  df$id <- seq_len(50)
  prepared <- .prepare_ising_input(df, id_col = "id")
  expect_false("id" %in% prepared$nodes)
  expect_equal(prepared$p, 3L)
})

test_that("Ising: zero-variance columns are dropped", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  df$Constant <- 1L
  expect_message(
    prepared <- .prepare_ising_input(df),
    "zero-variance"
  )
  expect_false("Constant" %in% prepared$nodes)
})

test_that("Ising: all-NA columns are dropped", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  df$Empty <- NA_real_
  expect_message(
    prepared <- .prepare_ising_input(df),
    "all-NA"
  )
  expect_false("Empty" %in% prepared$nodes)
})

test_that("Ising: rows with NAs are dropped", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  df$V1[1:3] <- NA
  expect_message(
    prepared <- .prepare_ising_input(df),
    "3 rows with NA"
  )
  expect_equal(prepared$n, 47L)
})

test_that("Ising: fewer than 2 columns after cleaning errors", {
  skip_if_not_installed("glmnet")
  df <- data.frame(X1 = c(0, 1, 0, 1, 0))
  expect_error(.prepare_ising_input(df), "At least 2")
})

test_that("Ising: non-syntactic column names are dropped", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  df[["100%"]] <- rbinom(50, 1, 0.5)
  expect_message(
    prepared <- .prepare_ising_input(df),
    "non-syntactic"
  )
  expect_equal(prepared$p, 3L)
})

test_that("Ising: continuous data raises error", {
  skip_if_not_installed("glmnet")
  df <- data.frame(X1 = rnorm(20), X2 = rnorm(20))
  expect_error(.prepare_ising_input(df), "binary.*0/1")
})


# ---- .log1pexp tests ----

test_that("log1pexp is numerically stable", {
  # Large positive: log(1 + exp(100)) ~ 100
  expect_equal(.log1pexp(100), 100, tolerance = 1e-10)
  # Large negative: log(1 + exp(-100)) ~ exp(-100)
  expect_equal(.log1pexp(-100), exp(-100), tolerance = 1e-30)
  # Zero: log(1 + exp(0)) = log(2)
  expect_equal(.log1pexp(0), log(2), tolerance = 1e-12)
  # Vectorized
  result <- .log1pexp(c(-100, 0, 100))
  expect_length(result, 3)
})


# ---- Basic structure tests ----

test_that("Ising: output has correct structure", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)

  expect_type(result, "list")
  expect_true(all(c("matrix", "nodes", "directed", "cleaned_data",
                     "thresholds", "asymm_weights", "rule", "gamma",
                     "n", "p", "lambda_selected") %in% names(result)))
})

test_that("Ising: output matrix is symmetric", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_true(isSymmetric(result$matrix, tol = 1e-10))
})

test_that("Ising: output matrix has zero diagonal", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_true(all(diag(result$matrix) == 0))
})

test_that("Ising: directed is FALSE", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_false(result$directed)
})

test_that("Ising: nodes match column names", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_equal(result$nodes, paste0("V", 1:4))
})

test_that("Ising: thresholds and lambda have correct length", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_length(result$thresholds, 4)
  expect_length(result$lambda_selected, 4)
  expect_true(all(result$lambda_selected > 0))
})

test_that("Ising: asymm_weights is p x p", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_equal(dim(result$asymm_weights), c(4, 4))
})

test_that("Ising: n and p are correct", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  result <- .estimator_ising(df)
  expect_equal(result$n, 100)
  expect_equal(result$p, 4)
})


# ---- Algorithm correctness tests ----

test_that("Ising: chain data detects adjacent edges", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(500, 4, seed = 123, coupling = 2.0)
  result <- .estimator_ising(df, gamma = 0.25)

  # Adjacent pairs should have nonzero edges
  expect_true(result$matrix[1, 2] != 0)  # V1-V2
  expect_true(result$matrix[2, 3] != 0)  # V2-V3
  expect_true(result$matrix[3, 4] != 0)  # V3-V4
})

test_that("Ising: null data produces sparse/empty network", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_null(200, 5, seed = 99)
  result <- .estimator_ising(df, gamma = 0.25)

  # Should be very sparse (most or all edges zero)
  n_nonzero <- sum(result$matrix[upper.tri(result$matrix)] != 0)
  total_possible <- choose(5, 2)
  expect_true(n_nonzero <= total_possible / 2)
})

test_that("Ising: chain non-adjacent edges weaker than adjacent", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(500, 4, seed = 123, coupling = 2.0)
  result <- .estimator_ising(df, gamma = 0.25)

  # Mean adjacent edge weight should exceed mean non-adjacent
  adj_weights <- c(
    abs(result$matrix[1, 2]),
    abs(result$matrix[2, 3]),
    abs(result$matrix[3, 4])
  )
  nonadj_weight <- abs(result$matrix[1, 4])  # V1-V4 (furthest apart)
  expect_true(mean(adj_weights) > nonadj_weight)
})


# ---- Parameter tests ----

test_that("Ising: higher gamma produces sparser network", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(200, 5, seed = 77)

  result_low <- .estimator_ising(df, gamma = 0)
  result_high <- .estimator_ising(df, gamma = 1.0)

  n_edges_low <- sum(result_low$matrix[upper.tri(result_low$matrix)] != 0)
  n_edges_high <- sum(result_high$matrix[upper.tri(result_high$matrix)] != 0)

  expect_true(n_edges_high <= n_edges_low)
})

test_that("Ising: AND rule produces <= OR rule edges", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(200, 5, seed = 55)

  result_and <- .estimator_ising(df, rule = "AND")
  result_or <- .estimator_ising(df, rule = "OR")

  n_and <- sum(result_and$matrix[upper.tri(result_and$matrix)] != 0)
  n_or <- sum(result_or$matrix[upper.tri(result_or$matrix)] != 0)

  expect_true(n_and <= n_or)
})

test_that("Ising: rule is stored in result", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 3)

  result_and <- .estimator_ising(df, rule = "AND")
  expect_equal(result_and$rule, "AND")

  result_or <- .estimator_ising(df, rule = "OR")
  expect_equal(result_or$rule, "OR")
})

test_that("Ising: gamma is stored in result", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 3)
  result <- .estimator_ising(df, gamma = 0.5)
  expect_equal(result$gamma, 0.5)
})

test_that("Ising: invalid rule errors", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  expect_error(.estimator_ising(df, rule = "MEAN"))
})

test_that("Ising: invalid gamma errors", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(50, 3)
  expect_error(.estimator_ising(df, gamma = -1))
})

test_that("Ising: nlambda is respected", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 3)
  # Just verify it runs without error with different nlambda

  result_20 <- .estimator_ising(df, nlambda = 20L)
  result_50 <- .estimator_ising(df, nlambda = 50L)
  expect_true(is.matrix(result_20$matrix))
  expect_true(is.matrix(result_50$matrix))
})


# ---- Symmetrization unit tests ----

test_that("symmetrize AND: edge only if both directions nonzero", {
  coef <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  coef[1, 2] <- 0.5   # a -> b
  coef[2, 1] <- 0.3   # b -> a (both nonzero: keep)
  coef[1, 3] <- 0.4   # a -> c
  coef[3, 1] <- 0.0   # c -> a = 0 (one zero: drop with AND)

  sym <- .symmetrize_ising(coef, rule = "AND")

  # a-b: both nonzero -> average
  expect_equal(sym[1, 2], (0.5 + 0.3) / 2)
  expect_equal(sym[2, 1], (0.5 + 0.3) / 2)
  # a-c: one zero -> dropped

  expect_equal(sym[1, 3], 0)
  expect_equal(sym[3, 1], 0)
})

test_that("symmetrize OR: simple average of both directions", {
  coef <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  coef[1, 2] <- 0.5
  coef[2, 1] <- 0.3
  coef[1, 3] <- 0.4
  coef[3, 1] <- 0.0

  sym <- .symmetrize_ising(coef, rule = "OR")

  # a-b: (0.5 + 0.3) / 2
  expect_equal(sym[1, 2], (0.5 + 0.3) / 2)
  # a-c: (0.4 + 0.0) / 2 = 0.2 (average including zeros, matching IsingFit)
  expect_equal(sym[1, 3], 0.2)
  expect_equal(sym[3, 1], 0.2)
})


# ---- Integration tests: build_network ----

test_that("Ising: build_network(method='ising') works", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  net <- build_network(df, method = "ising")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "ising")
  expect_false(net$directed)
  expect_equal(net$n_nodes, 4)
  expect_true(isSymmetric(net$matrix, tol = 1e-10))
})

test_that("Ising: alias 'isingfit' resolves correctly", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 3)
  net <- build_network(df, method = "isingfit")

  expect_equal(net$method, "ising")
})

test_that("Ising: params are passed through build_network", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 3)
  net <- build_network(df, method = "ising",
                       params = list(gamma = 0.5, rule = "OR"))

  expect_equal(net$gamma, 0.5)
  expect_equal(net$rule, "OR")
})

test_that("Ising: print.netobject works for ising", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  net <- build_network(df, method = "ising")

  output <- capture.output(print(net))
  expect_true(any(grepl("Ising", output)))
  expect_true(any(grepl("Gamma", output)))
  expect_true(any(grepl("Rule", output)))
  expect_true(any(grepl("Threshold", output)))
})

test_that("Ising: list_estimators includes ising", {
  skip_if_not_installed("glmnet")
  # Force registry reload
  .register_builtin_estimators()
  est_list <- list_estimators()
  expect_true("ising" %in% est_list$name)
})

test_that("Ising: scaling works with ising method", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  net <- build_network(df, method = "ising", scaling = "max")

  max_abs <- max(abs(net$matrix))
  if (max_abs > 0) {
    expect_equal(max_abs, 1, tolerance = 1e-10)
  }
})

test_that("Ising: threshold works with ising method", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 4)
  net <- build_network(df, method = "ising", threshold = 0.1)

  nonzero <- net$matrix[net$matrix != 0]
  if (length(nonzero) > 0) {
    expect_true(all(abs(nonzero) >= 0.1))
  }
})


# ---- Bootstrap and permutation integration ----

test_that("Ising: bootstrap_network runs without error", {
  skip_if_not_installed("glmnet")
  df <- .make_ising_data(100, 3)

  boot <- bootstrap_network(build_network(df, method = "ising"), iter = 10)
  expect_s3_class(boot, "net_bootstrap")
})

# ---- Equivalence with IsingFit ----

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


test_that("Ising: permutation_test runs without error", {
  skip_if_not_installed("glmnet")
  df1 <- .make_ising_data(80, 3, seed = 1)
  df2 <- .make_ising_data(80, 3, seed = 2)

  net1 <- build_network(df1, method = "ising")
  net2 <- build_network(df2, method = "ising")

  perm <- permutation_test(net1, net2, iter = 10)
  expect_s3_class(perm, "net_permutation")
})
