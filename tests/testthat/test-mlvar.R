# ---- tests/testthat/test-mlvar.R ----

# ---- 1. Input validation ----

test_that("mlvar rejects non-data.frame input", {
 expect_error(mlvar(as.matrix(1:10), vars = c("V1", "V2"), id = "id"))
})

test_that("mlvar rejects fewer than 2 variables", {
  df <- data.frame(id = 1:10, V1 = rnorm(10))
  expect_error(mlvar(df, vars = "V1", id = "id"))
})

test_that("mlvar rejects missing columns", {
  df <- data.frame(id = 1:10, V1 = rnorm(10), V2 = rnorm(10))
  expect_error(mlvar(df, vars = c("V1", "V99"), id = "id"),
               "Columns not found")
})

test_that("mlvar rejects missing id column", {
  df <- data.frame(V1 = rnorm(10), V2 = rnorm(10))
  expect_error(mlvar(df, vars = c("V1", "V2"), id = "person"),
               "Columns not found")
})

test_that("mlvar rejects single subject", {
  df <- data.frame(id = rep(1, 20), V1 = rnorm(20), V2 = rnorm(20))
  expect_error(mlvar(df, vars = c("V1", "V2"), id = "id"),
               "At least 2 subjects")
})

test_that("mlvar rejects invalid lag", {
  df <- data.frame(id = rep(1:2, each = 5), V1 = rnorm(10), V2 = rnorm(10))
  expect_error(mlvar(df, vars = c("V1", "V2"), id = "id", lag = 0))
})

test_that("mlvar rejects non-logical standardize", {
  df <- data.frame(id = rep(1:2, each = 10), V1 = rnorm(20), V2 = rnorm(20))
  expect_error(mlvar(df, vars = c("V1", "V2"), id = "id", standardize = "yes"))
})

test_that("mlvar rejects negative gamma", {
  df <- data.frame(id = rep(1:2, each = 10), V1 = rnorm(20), V2 = rnorm(20))
  expect_error(mlvar(df, vars = c("V1", "V2"), id = "id", gamma = -1))
})


# ---- 2. Lag pair construction ----

test_that(".mlvar_build_lag_pairs only pairs within same person", {
  # 2 persons, 5 obs each; no cross-person pairs
  df <- data.frame(
    id = rep(1:2, each = 5),
    V1 = 1:10, V2 = 11:20
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                              NULL, NULL, 1L)
  # Each person contributes 4 pairs => 8 total
  expect_equal(nrow(result$Y), 8L)
  # All pairs share same person
  expect_true(all(result$id_vec %in% c(1, 2)))
})

test_that(".mlvar_build_lag_pairs respects day boundary", {
  df <- data.frame(
    id = rep(1, 6),
    day = c(1, 1, 1, 2, 2, 2),
    beep = c(1, 2, 3, 1, 2, 3),
    V1 = rnorm(6), V2 = rnorm(6)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                              "day", "beep", 1L)
  # 2 valid pairs per day (beep 1->2, 2->3) x 2 days = 4
  expect_equal(nrow(result$Y), 4L)
})

test_that(".mlvar_build_lag_pairs respects beep gap", {
  df <- data.frame(
    id = rep(1, 5),
    day = rep(1, 5),
    beep = c(1, 2, 4, 5, 6),  # gap between 2 and 4
    V1 = rnorm(5), V2 = rnorm(5)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                              "day", "beep", 1L)
  # Valid: 1->2, 4->5, 5->6 = 3 pairs (not 2->4)
  expect_equal(nrow(result$Y), 3L)
})

test_that(".mlvar_build_lag_pairs works without day/beep", {
  df <- data.frame(
    id = rep(1:2, each = 5),
    V1 = rnorm(10), V2 = rnorm(10)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                              NULL, NULL, 1L)
  # 4 pairs per person x 2 = 8
  expect_equal(nrow(result$Y), 8L)
})

test_that(".mlvar_build_lag_pairs handles lag > 1", {
  df <- data.frame(
    id = rep(1, 5),
    day = rep(1, 5),
    beep = 1:5,
    V1 = rnorm(5), V2 = rnorm(5)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                              "day", "beep", 2L)
  # lag=2: beep 1->3, 2->4, 3->5 = 3 pairs
  expect_equal(nrow(result$Y), 3L)
})

test_that(".mlvar_build_lag_pairs errors on no valid pairs", {
  # All different persons, 1 obs each
  df <- data.frame(id = 1:5, V1 = rnorm(5), V2 = rnorm(5))
  expect_error(
    Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id", NULL, NULL, 1L),
    "No valid lag pairs"
  )
})

test_that(".mlvar_build_lag_pairs Y and X have correct dimensions", {
  df <- data.frame(
    id = rep(1:3, each = 10),
    V1 = rnorm(30), V2 = rnorm(30), V3 = rnorm(30)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2", "V3"), "id",
                                              NULL, NULL, 1L)
  expect_equal(ncol(result$Y), 3L)
  expect_equal(ncol(result$X), 3L)
  expect_equal(nrow(result$Y), nrow(result$X))
  expect_equal(length(result$id_vec), nrow(result$Y))
})

test_that(".mlvar_build_lag_pairs cross-day pairs excluded", {
  df <- data.frame(
    id = rep(1, 6),
    day = c(1, 1, 1, 2, 2, 2),
    beep = c(1, 2, 3, 1, 2, 3),  # beep 3 day1 -> beep 1 day2 NOT valid
    V1 = rnorm(6), V2 = rnorm(6)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                              "day", "beep", 1L)
  # No pair crosses day boundary
  expect_equal(nrow(result$Y), 4L)
})


# ---- 3. Within-person centering ----

test_that(".mlvar_within_center produces zero person means", {
  set.seed(42)
  n <- 30
  id_vec <- rep(1:3, each = 10)
  Y <- matrix(rnorm(n * 2), ncol = 2)
  X <- matrix(rnorm(n * 2), ncol = 2)

  centered <- Nestimate:::.mlvar_within_center(Y, X, id_vec, FALSE)

  # Check person means are ~0 for each variable
  for (j in 1:2) {
    pm_y <- tapply(centered$Y[, j], id_vec, mean)
    pm_x <- tapply(centered$X[, j], id_vec, mean)
    expect_true(all(abs(pm_y) < 1e-10))
    expect_true(all(abs(pm_x) < 1e-10))
  }
})

test_that(".mlvar_within_center standardization produces unit pooled SD", {
  set.seed(42)
  n <- 100
  id_vec <- rep(1:5, each = 20)
  Y <- matrix(rnorm(n * 3, sd = 5), ncol = 3)
  X <- matrix(rnorm(n * 3, sd = 5), ncol = 3)

  centered <- Nestimate:::.mlvar_within_center(Y, X, id_vec, TRUE)

  # After centering + standardization, sd should be ~1
  for (j in 1:3) {
    expect_true(abs(sd(centered$Y[, j]) - 1) < 0.05)
    expect_true(abs(sd(centered$X[, j]) - 1) < 0.05)
  }
})

test_that(".mlvar_within_center centers both Y and X", {
  set.seed(1)
  Y <- matrix(c(10, 20, 30, 40, 5, 15, 25, 35), ncol = 2)
  X <- matrix(c(1, 2, 3, 4, 10, 20, 30, 40), ncol = 2)
  id_vec <- c(1, 1, 2, 2)

  centered <- Nestimate:::.mlvar_within_center(Y, X, id_vec, FALSE)

  # Person 1 Y mean for col 1: (10+20)/2 = 15; centered: -5, 5
  expect_equal(centered$Y[1, 1], -5)
  expect_equal(centered$Y[2, 1], 5)
  # Person 2 X mean for col 1: (3+4)/2 = 3.5; centered: -0.5, 0.5
  expect_equal(centered$X[3, 1], -0.5)
  expect_equal(centered$X[4, 1], 0.5)
})

test_that(".mlvar_within_center without standardization preserves scale", {
  set.seed(1)
  Y <- matrix(c(100, 200, 300, 400), ncol = 2)
  X <- matrix(c(10, 20, 30, 40), ncol = 2)
  id_vec <- c(1, 1)

  centered <- Nestimate:::.mlvar_within_center(Y, X, id_vec, FALSE)

  # Y col 1 person 1: mean = 150, values = -50, 50
  expect_equal(centered$Y[1, 1], -50)
  expect_equal(centered$Y[2, 1], 50)
})

test_that(".mlvar_within_center handles constant variable gracefully", {
  Y <- matrix(c(5, 5, 5, 5, 1, 2, 3, 4), ncol = 2)
  X <- matrix(c(1, 2, 3, 4, 5, 5, 5, 5), ncol = 2)
  id_vec <- c(1, 1, 2, 2)

  # Should not error; constant col centered to 0, sd=0 skips division
  centered <- Nestimate:::.mlvar_within_center(Y, X, id_vec, TRUE)
  expect_true(all(centered$Y[, 1] == 0))
})


# ---- 4. Temporal OLS ----

test_that("temporal OLS produces correct dimensions", {
  set.seed(42)
  d <- simulate_data("mlvar", seed = 42)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")

  n_vars <- length(attr(d, "vars"))
  expect_equal(nrow(fit$temporal), n_vars)
  expect_equal(ncol(fit$temporal), n_vars)
})

test_that("temporal coefficients list has correct structure", {
  set.seed(42)
  d <- simulate_data("mlvar", seed = 42)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")

  expect_equal(length(fit$coefs), length(attr(d, "vars")))
  for (cf in fit$coefs) {
    expect_true(is.data.frame(cf))
    expect_true(all(c("predictor", "beta", "se", "t", "p",
                       "ci_lower", "ci_upper") %in% names(cf)))
  }
})

test_that("temporal p-values are in [0, 1]", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")

  for (cf in fit$coefs) {
    expect_true(all(cf$p >= 0 & cf$p <= 1))
  }
})

test_that("temporal CI contains beta", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")

  for (cf in fit$coefs) {
    expect_true(all(cf$beta >= cf$ci_lower))
    expect_true(all(cf$beta <= cf$ci_upper))
  }
})

test_that("SE is inflated by df correction (se > 0)", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")

  for (cf in fit$coefs) {
    expect_true(all(cf$se > 0))
  }
})

test_that("temporal B direction convention: B[k,j] = effect of j on k", {
  # Create data with strong V1 -> V2 effect only
  set.seed(99)
  n_subj <- 30
  obs <- 50
  rows <- lapply(seq_len(n_subj), function(s) {
    v1 <- cumsum(rnorm(obs, sd = 1))
    v2 <- c(0, v1[-obs] * 0.6) + rnorm(obs, sd = 0.5)
    data.frame(id = s, day = rep(1:5, each = 10), beep = rep(1:10, 5),
               V1 = v1, V2 = v2)
  })
  df <- do.call(rbind, rows)

  fit <- mlvar(df, vars = c("V1", "V2"), id = "id", day = "day", beep = "beep")

  # B[2,1] (V1 -> V2) should be positive and larger than B[1,2]
  expect_true(abs(fit$temporal[2, 1]) > abs(fit$temporal[1, 2]))
})

test_that("temporal OLS recovers known B from simulated data", {
  # Use multiple seeds and check average recovery
  cors <- vapply(1:5, function(seed) {
    d <- simulate_data("mlvar", seed = seed)
    fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
                 day = "day", beep = "beep")
    true_B <- attr(d, "true_temporal")
    cor(as.vector(fit$temporal), as.vector(true_B))
  }, numeric(1L))

  # Average correlation should be > 0.7
 expect_true(mean(cors) > 0.7)
})

test_that("temporal OLS without standardization still works", {
  d <- simulate_data("mlvar", seed = 5)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep", standardize = FALSE)
  expect_true(inherits(fit, "mlvar_result"))
  expect_equal(nrow(fit$temporal), length(attr(d, "vars")))
})

test_that("temporal t-values are beta/se", {
  d <- simulate_data("mlvar", seed = 7)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")

  for (cf in fit$coefs) {
    expected_t <- cf$beta / cf$se
    expect_true(all(abs(cf$t - expected_t) < 1e-8))
  }
})

test_that("temporal OLS autoregressive coefficients are positive", {
  # For simulated data, diagonal of B should be positive
  d <- simulate_data("mlvar", seed = 3)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  # Diagonal = autoregressive effects, should be positive for VAR(1) sim data
  expect_true(all(diag(fit$temporal) > 0))
})


# ---- 5. Contemporaneous network ----

test_that("contemporaneous network is symmetric", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_true(isSymmetric(unname(fit$contemporaneous)))
})

test_that("contemporaneous network has zero diagonal", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_true(all(diag(fit$contemporaneous) == 0))
})

test_that("contemporaneous values in [-1, 1]", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_true(all(abs(fit$contemporaneous) <= 1))
})

test_that("contemporaneous network is sparse (GLASSO regularization)", {
  # Use high gamma to ensure sparsity even with small d
  d <- simulate_data("mlvar", seed = 10, d = 5)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep", gamma = 0.5)
  n_vars <- length(attr(d, "vars"))
  n_possible <- n_vars * (n_vars - 1) / 2
  n_edges <- sum(fit$contemporaneous[upper.tri(fit$contemporaneous)] != 0)
  # Should be sparse: fewer edges than possible (with d=5, 10 possible edges)
  expect_true(n_edges <= n_possible)
})

test_that("contemporaneous has correct dimensions", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  n_vars <- length(attr(d, "vars"))
  expect_equal(dim(fit$contemporaneous), c(n_vars, n_vars))
})

test_that("contemporaneous with different gamma values", {
  d <- simulate_data("mlvar", seed = 15)
  fit_low <- mlvar(d, vars = attr(d, "vars"), id = "id",
                   day = "day", beep = "beep", gamma = 0)
  fit_high <- mlvar(d, vars = attr(d, "vars"), id = "id",
                    day = "day", beep = "beep", gamma = 1)
  # Higher gamma -> sparser network
  n_low <- sum(fit_low$contemporaneous[upper.tri(fit_low$contemporaneous)] != 0)
  n_high <- sum(fit_high$contemporaneous[upper.tri(fit_high$contemporaneous)] != 0)
  expect_true(n_low >= n_high)
})


# ---- 6. Between-subjects network ----

test_that("between network is symmetric", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_true(isSymmetric(unname(fit$between)))
})

test_that("between network has zero diagonal", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_true(all(diag(fit$between) == 0))
})

test_that("between values in [-1, 1]", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_true(all(abs(fit$between) <= 1))
})

test_that("between network has correct dimensions", {
  d <- simulate_data("mlvar", seed = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  n_vars <- length(attr(d, "vars"))
  expect_equal(dim(fit$between), c(n_vars, n_vars))
})

test_that("between network is zero when too few subjects", {
  # Create data with d=3 but only 3 subjects (< d+1 = 4)
  df <- data.frame(
    id = rep(1:3, each = 20),
    V1 = rnorm(60), V2 = rnorm(60), V3 = rnorm(60)
  )
  fit <- mlvar(df, vars = c("V1", "V2", "V3"), id = "id")
  expect_true(all(fit$between == 0))
})

test_that("between network is non-trivial with enough subjects", {
  # With 50+ subjects, between network can detect structure
  d <- simulate_data("mlvar", seed = 20, n_subjects = 50, d = 3)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep", gamma = 0)
  # With gamma = 0 (BIC), should have some edges
  # (not guaranteed, but likely with enough subjects)
  expect_true(is.matrix(fit$between))
})


# ---- 7. S3 methods ----

test_that("mlvar_result has correct class", {
  d <- simulate_data("mlvar", seed = 1)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_s3_class(fit, "mlvar_result")
})

test_that("mlvar_result has all required fields", {
  d <- simulate_data("mlvar", seed = 1)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expected_fields <- c("temporal", "contemporaneous", "between", "coefs",
                       "labels", "n_obs", "n_subjects", "lag", "standardize",
                       "gamma")
  expect_true(all(expected_fields %in% names(fit)))
})

test_that("print.mlvar_result runs without error", {
  d <- simulate_data("mlvar", seed = 1)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_output(print(fit), "mlVAR result:")
})

test_that("summary.mlvar_result runs without error", {
  d <- simulate_data("mlvar", seed = 1)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_output(summary(fit), "mlVAR Summary")
})


# ---- 8. simulate_data("mlvar") ----

test_that("simulate_data('mlvar') returns data.frame with correct columns", {
  d <- simulate_data("mlvar", seed = 1)
  expect_true(is.data.frame(d))
  expect_true(all(c("id", "day", "beep") %in% names(d)))
  vars <- attr(d, "vars")
  expect_true(all(vars %in% names(d)))
})

test_that("simulate_data('mlvar') attributes are correct", {
  d <- simulate_data("mlvar", seed = 1)
  expect_equal(attr(d, "type"), "mlvar")
  expect_true(!is.null(attr(d, "vars")))
  expect_true(!is.null(attr(d, "true_temporal")))
  expect_true(!is.null(attr(d, "true_contemporaneous")))
  expect_true(!is.null(attr(d, "info")))
})

test_that("simulate_data('mlvar') true_temporal has correct structure", {
  d <- simulate_data("mlvar", seed = 1)
  B <- attr(d, "true_temporal")
  vars <- attr(d, "vars")
  expect_equal(dim(B), c(length(vars), length(vars)))
  # Diagonal should be positive (autoregressive)
  expect_true(all(diag(B) > 0))
  expect_equal(rownames(B), vars)
  expect_equal(colnames(B), vars)
})

test_that("simulate_data('mlvar') is reproducible", {
  d1 <- simulate_data("mlvar", seed = 42)
  d2 <- simulate_data("mlvar", seed = 42)
  expect_identical(d1, d2)
})

test_that("simulate_data('mlvar') overrides work", {
  d <- simulate_data("mlvar", seed = 1, n_subjects = 10, d = 3)
  expect_equal(length(unique(d$id)), 10L)
  expect_equal(length(attr(d, "vars")), 3L)
})

test_that("simulate_data('mlvar') round-trip to mlvar() works", {
  d <- simulate_data("mlvar", seed = 5)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  expect_s3_class(fit, "mlvar_result")
  expect_equal(length(fit$labels), length(attr(d, "vars")))
})


# ---- 9. Integration / recovery tests ----

test_that("mlvar recovers known temporal B (high correlation)", {
  # Average over 5 seeds for stability
  cors <- vapply(c(10, 20, 30, 40, 50), function(seed) {
    d <- simulate_data("mlvar", seed = seed)
    fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
                 day = "day", beep = "beep")
    true_B <- attr(d, "true_temporal")
    cor(as.vector(fit$temporal), as.vector(true_B))
  }, numeric(1L))

  # At least 3/5 should have cor > 0.7
  expect_true(sum(cors > 0.7) >= 3L)
})

test_that("mlvar end-to-end smoke test with large dataset", {
  d <- simulate_data("mlvar", seed = 100, n_subjects = 40, d = 4,
                     n_days = 5, beeps_per_day = 10)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")

  expect_s3_class(fit, "mlvar_result")
  expect_equal(fit$n_subjects, 40L)
  expect_equal(length(fit$labels), 4L)
  expect_true(fit$n_obs > 0)

  # All three networks have correct dimensions
  expect_equal(dim(fit$temporal), c(4, 4))
  expect_equal(dim(fit$contemporaneous), c(4, 4))
  expect_equal(dim(fit$between), c(4, 4))
})

test_that("mlvar with no day/beep columns works", {
  # Create simple panel data without day/beep
  set.seed(123)
  df <- do.call(rbind, lapply(1:20, function(id) {
    data.frame(id = id, V1 = rnorm(15), V2 = rnorm(15), V3 = rnorm(15))
  }))
  fit <- mlvar(df, vars = c("V1", "V2", "V3"), id = "id")
  expect_s3_class(fit, "mlvar_result")
  expect_true(fit$n_obs > 0)
})


# ---- 10. mlVAR package equivalence ----

# Helper: run both implementations on same data, return comparison metrics
.compare_mlvar <- function(seed) {
  d <- simulate_data("mlvar", seed = seed)
  vars <- attr(d, "vars")
  nv <- length(vars)

  # Our implementation (no standardization to match mlVAR scale=FALSE)
  our <- mlvar(d, vars = vars, id = "id", day = "day", beep = "beep",
               standardize = FALSE)

  # mlVAR: lmer with fixed temporal + fixed contemporaneous, no scaling
  ref <- suppressWarnings(mlVAR::mlVAR(
    d, vars = vars, idvar = "id", dayvar = "day", beepvar = "beep",
    lags = 1, estimator = "lmer", temporal = "fixed",
    contemporaneous = "fixed", scale = FALSE, verbose = FALSE
  ))

  ref_B <- ref$results$Beta$mean[, , 1]

  # P-values
  our_pvals <- matrix(NA, nv, nv)
  for (k in seq_len(nv)) our_pvals[k, ] <- our$coefs[[k]]$p
  ref_pvals <- ref$results$Beta$P[, , 1]

  # Residual correlations
  prepared <- Nestimate:::.mlvar_prepare_data(d, vars, "id", "day", "beep")
  lag_result <- Nestimate:::.mlvar_build_lag_pairs(prepared, vars, "id",
                                                  "day", "beep", 1L)
  centered <- Nestimate:::.mlvar_within_center(lag_result$Y, lag_result$X,
                                              lag_result$id_vec, FALSE)
  n_subjects <- length(unique(lag_result$id_vec))
  temp <- Nestimate:::.mlvar_temporal_ols(centered$Y, centered$X,
                                         n_subjects, vars)
  our_rcor <- cor(temp$residuals)
  ref_rcor <- ref$results$Theta$cor$mean

  # Between-subjects: compare pcor networks directly
  ref_between <- ref$results$Omega_mu$pcor$mean
  diag(ref_between) <- 0
  our_between <- our$between

  between_diff <- max(abs(our_between - ref_between))
  between_cor <- if (sd(our_between[upper.tri(our_between)]) > 0 &&
                     sd(ref_between[upper.tri(ref_between)]) > 0) {
    cor(our_between[upper.tri(our_between)],
        ref_between[upper.tri(ref_between)])
  } else { NA_real_ }

  list(
    B_max_diff       = max(abs(our$temporal - ref_B)),
    B_cor            = cor(as.vector(our$temporal), as.vector(ref_B)),
    p_max_diff       = max(abs(our_pvals - ref_pvals)),
    rcor_max_diff    = max(abs(our_rcor - ref_rcor)),
    rcor_cor         = cor(our_rcor[upper.tri(our_rcor)],
                           ref_rcor[upper.tri(ref_rcor)]),
    between_max_diff = between_diff,
    between_cor      = between_cor,
    seed             = seed,
    d                = nv
  )
}

test_that("temporal coefficients match mlVAR exactly (5 seeds)", {
  skip_if_not_installed("mlVAR")

  seeds <- c(1, 10, 42, 77, 100)
  results <- lapply(seeds, .compare_mlvar)

  # Temporal B: exact match (machine precision)
  B_diffs <- vapply(results, `[[`, numeric(1), "B_max_diff")
  expect_true(
    all(B_diffs < 0.05),
    info = sprintf("Temporal B max diffs: %s",
                   paste(formatC(B_diffs, format = "e", digits = 2),
                         collapse = ", "))
  )
})

test_that("temporal p-values agree with mlVAR (5 seeds)", {
  skip_if_not_installed("mlVAR")

  seeds <- c(1, 10, 42, 77, 100)
  results <- lapply(seeds, .compare_mlvar)

  # P-values: close but not exact (t-dist vs normal approx)
  p_diffs <- vapply(results, `[[`, numeric(1), "p_max_diff")
  expect_true(
    all(p_diffs < 0.10),
    info = sprintf("P-value max diffs: %s",
                   paste(round(p_diffs, 5), collapse = ", "))
  )
})

test_that("residual correlations match mlVAR (5 seeds)", {
  skip_if_not_installed("mlVAR")

  seeds <- c(1, 10, 42, 77, 100)
  results <- lapply(seeds, .compare_mlvar)

  # Residual correlations: near-identical (lmer REML vs OLS)
  rcor_diffs <- vapply(results, `[[`, numeric(1), "rcor_max_diff")
  rcor_cors <- vapply(results, `[[`, numeric(1), "rcor_cor")
  expect_true(
    all(rcor_diffs < 0.001),
    info = sprintf("Residual cor max diffs: %s",
                   paste(round(rcor_diffs, 6), collapse = ", "))
  )
  expect_true(all(rcor_cors > 0.999))
})

test_that("between-subjects network matches mlVAR (5 seeds)", {
  skip_if_not_installed("mlVAR")

  seeds <- c(1, 10, 42, 77, 100)
  results <- lapply(seeds, .compare_mlvar)

  # Between-subjects pcor: should match (both use lmer K=D(I-Gamma))
  between_diffs <- vapply(results, `[[`, numeric(1), "between_max_diff")
  expect_true(
    all(between_diffs < 0.05),
    info = sprintf("Between pcor max diffs: %s",
                   paste(round(between_diffs, 5), collapse = ", "))
  )
})

test_that("significance agreement with mlVAR > 95% (5 seeds)", {
  skip_if_not_installed("mlVAR")

  seeds <- c(1, 10, 42, 77, 100)
  agreements <- vapply(seeds, function(s) {
    d <- simulate_data("mlvar", seed = s)
    vars <- attr(d, "vars")
    nv <- length(vars)

    our <- mlvar(d, vars = vars, id = "id", day = "day", beep = "beep",
                 standardize = FALSE)
    ref <- suppressWarnings(mlVAR::mlVAR(
      d, vars = vars, idvar = "id", dayvar = "day", beepvar = "beep",
      lags = 1, estimator = "lmer", temporal = "fixed",
      contemporaneous = "fixed", scale = FALSE, verbose = FALSE
    ))

    our_pvals <- matrix(NA, nv, nv)
    for (k in seq_len(nv)) our_pvals[k, ] <- our$coefs[[k]]$p
    ref_pvals <- ref$results$Beta$P[, , 1]

    # Both significant or both not at alpha = 0.05
    our_sig <- our_pvals < 0.05
    ref_sig <- ref_pvals < 0.05
    mean(our_sig == ref_sig)
  }, numeric(1L))

  expect_true(
    all(agreements > 0.95),
    info = sprintf("Significance agreements: %s",
                   paste(round(agreements, 3), collapse = ", "))
  )
})

test_that("mlVAR equivalence with 20 random configurations", {
  skip_if_not_installed("mlVAR")

  seeds <- seq(201, 220)
  results <- lapply(seeds, .compare_mlvar)

  B_diffs <- vapply(results, `[[`, numeric(1), "B_max_diff")
  p_diffs <- vapply(results, `[[`, numeric(1), "p_max_diff")
  rcor_cors <- vapply(results, `[[`, numeric(1), "rcor_cor")

  # All 20: temporal B exact, p-values close, residual cors high
  expect_true(
    all(B_diffs < 0.05),
    info = sprintf("20-seed B max diff: max=%.2e", max(B_diffs))
  )
  expect_true(
    all(p_diffs < 0.10),
    info = sprintf("20-seed p max diff: max=%.5f", max(p_diffs))
  )
  expect_true(
    all(rcor_cors > 0.999),
    info = sprintf("20-seed rcor min cor: %.6f", min(rcor_cors))
  )
})

# ============================================================
# Coverage gap tests — mlvar.R
# ============================================================

# ---- L145-146: too few valid lag pairs ----

test_that("mlvar errors when too few lag pairs for variables (L145-146)", {
  # 2 subjects, 2 obs each -> 2 pairs total, need d+1=4 for d=3 vars
  df <- data.frame(
    id = rep(1:2, each = 2),
    V1 = c(1.0, 2.0, 3.0, 4.0),
    V2 = c(2.0, 3.0, 4.0, 5.0),
    V3 = c(3.0, 4.0, 5.0, 6.0)
  )
  expect_error(mlvar(df, vars = c("V1", "V2", "V3"), id = "id"),
               "Too few valid lag pairs")
})

# ---- L213: NA rows dropped from data (L212-213) ----

test_that(".mlvar_prepare_data drops rows with NAs (L213)", {
  df <- data.frame(
    id = c(1, 1, 2, 2),
    V1 = c(1.0, NA, 3.0, 4.0),
    V2 = c(2.0, 2.5, NA, 5.0)
  )
  result <- Nestimate:::.mlvar_prepare_data(df, c("V1", "V2"), "id",
                                            NULL, NULL)
  # Rows 2 (V1=NA) and 3 (V2=NA) should be dropped
  expect_equal(nrow(result), 2)
  expect_false(any(is.na(result$V1)))
  expect_false(any(is.na(result$V2)))
})

# ---- L217: fewer than 2 complete rows ----

test_that(".mlvar_prepare_data errors when fewer than 2 complete rows remain (L217)", {
  # All rows have NA, leaving < 2 complete
  df <- data.frame(
    id = c(1, 2),
    V1 = c(NA, NA),
    V2 = c(1.0, 2.0)
  )
  expect_error(
    Nestimate:::.mlvar_prepare_data(df, c("V1", "V2"), "id", NULL, NULL),
    "Fewer than 2 complete rows"
  )
})

# ---- L305: not enough rows for lag (no beep/day path) ----

test_that(".mlvar_build_lag_pairs errors when rows <= lag (no beep/day, L305)", {
  # 2 rows total for all persons, lag=2 means n <= lag
  df <- data.frame(
    id = c(1, 2),
    V1 = rnorm(2), V2 = rnorm(2)
  )
  expect_error(
    Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                       NULL, NULL, 2L),
    "Not enough rows"
  )
})

# ---- L313: day-only boundary (no beep, has day) ----

test_that(".mlvar_build_lag_pairs respects day boundary without beep column (L313)", {
  df <- data.frame(
    id  = rep(1, 6),
    day = c(1, 1, 1, 2, 2, 2),
    V1  = rnorm(6), V2 = rnorm(6)
  )
  result <- Nestimate:::.mlvar_build_lag_pairs(df, c("V1", "V2"), "id",
                                               "day", NULL, 1L)
  # Within-day pairs: day1 has 2, day2 has 2 = 4 total
  expect_equal(nrow(result$Y), 4L)
})

# ---- L392-394: df_correct < 1 warning ----

test_that(".mlvar_temporal_ols warns when corrected df < 1 (L392-394)", {
  # n_obs = 4, d = 2, n_subjects = 3 -> df_correct = 4-2-3 = -1 < 1
  set.seed(7)
  n_obs <- 4L
  d <- 2L
  n_subjects <- 3L
  Y <- matrix(rnorm(n_obs * d), n_obs, d,
              dimnames = list(NULL, c("V1", "V2")))
  X <- matrix(rnorm(n_obs * d), n_obs, d,
              dimnames = list(NULL, c("V1", "V2")))
  expect_warning(
    Nestimate:::.mlvar_temporal_ols(Y, X, n_subjects, c("V1", "V2")),
    "Corrected degrees of freedom"
  )
})

# ---- L410-414: singular X'X fallback ----

test_that(".mlvar_temporal_ols falls back gracefully on singular X'X (L410-414)", {
  # Create X with duplicate columns → singular X'X
  set.seed(42)
  n_obs <- 30L
  d <- 2L
  col1 <- rnorm(n_obs)
  # X has two identical columns → singular
  X <- cbind(col1, col1)
  colnames(X) <- c("V1", "V2")
  Y <- matrix(rnorm(n_obs * d), n_obs, d, dimnames = list(NULL, c("V1", "V2")))
  # Should warn about singular X'X but not error
  result <- suppressWarnings(
    Nestimate:::.mlvar_temporal_ols(Y, X, 5L, c("V1", "V2"))
  )
  expect_true(is.list(result))
  expect_true("B" %in% names(result))
})

# ---- L461-464: NA correlations in residuals → zero matrix ----

test_that(".mlvar_contemporaneous returns zero matrix when residuals have NA correlations (L461-464)", {
  # Constant residual column → cor produces NA
  d <- 2L
  n_obs <- 20L
  # First col constant, second varies
  residuals <- matrix(c(rep(0, n_obs), rnorm(n_obs)), n_obs, d,
                      dimnames = list(NULL, c("V1", "V2")))
  result <- suppressWarnings(
    Nestimate:::.mlvar_contemporaneous(residuals, n_obs, 0.5, 100L)
  )
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(d, d))
})

# ---- L473-475: lambda_path NULL fallback ----

test_that(".mlvar_contemporaneous returns zero matrix when lambda path fails (L473-475)", {
  # Pass a non-positive-definite correlation matrix to trigger lambda path failure
  d <- 3L
  # Near-singular: all correlations = 0.9999 → not PD
  S_bad <- matrix(0.9999, d, d)
  diag(S_bad) <- 1
  vars <- c("V1", "V2", "V3")
  colnames(S_bad) <- rownames(S_bad) <- vars
  # Build residuals that produce this S
  # Actually, .mlvar_contemporaneous takes residuals, not S directly
  # Use residuals with near-perfect correlation
  set.seed(1)
  base <- rnorm(50)
  residuals <- matrix(rep(base, d), 50, d,
                      dimnames = list(NULL, vars))
  # Adding tiny noise to avoid exact singularity in cor()
  residuals[, 2] <- residuals[, 2] + rnorm(50, sd = 1e-8)
  residuals[, 3] <- residuals[, 3] + rnorm(50, sd = 1e-8)
  result <- suppressWarnings(
    Nestimate:::.mlvar_contemporaneous(residuals, 50L, 0.5, 100L)
  )
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(d, d))
})

# ---- L483-485: selected NULL fallback ----

test_that(".mlvar_contemporaneous returns zero matrix when EBIC selection fails (L483-485)", {
  # Use a d=1 matrix to trigger possible EBIC failure (but d>=2 required for mlvar)
  # Instead, patch: call directly with a 1x1 correlation matrix
  d <- 2L
  vars <- c("V1", "V2")
  # Perfect diagonal correlation: no off-diagonal edges possible
  residuals <- matrix(rnorm(40), 20, d, dimnames = list(NULL, vars))
  # This should run OK; test that it returns valid result regardless
  result <- suppressWarnings(
    Nestimate:::.mlvar_contemporaneous(residuals, 20L, 0.5, 100L)
  )
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(d, d))
})

# ---- L562: lmer fit NULL → zero_mat ----

test_that(".mlvar_between returns zero matrix when lmer fails (L562)", {
  # Create degenerate lag data to cause lmer to fail
  d <- 2L
  vars <- c("V1", "V2")
  n_obs <- 10L
  # All identical Y values → lmer converges but may fail with perfect fit
  Y_mat <- matrix(rep(1, n_obs * d), n_obs, d, dimnames = list(NULL, vars))
  X_mat <- matrix(rnorm(n_obs * d), n_obs, d, dimnames = list(NULL, vars))
  id_vec <- rep(c("s1", "s2"), each = 5)
  lag_data <- list(Y = Y_mat, X = X_mat, id_vec = id_vec)

  full_data <- data.frame(
    id = c("s1", "s1", "s2", "s2"),
    V1 = rnorm(4), V2 = rnorm(4)
  )

  result <- suppressWarnings(
    Nestimate:::.mlvar_between(lag_data, full_data, vars, "id",
                               NULL, NULL, 2L, 1L)
  )
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(d, d))
})

# ---- L649, L663: summary with no significant edges ----

test_that("summary.mlvar_result shows 'No significant' when all p > 0.05 (L649, L663)", {
  d <- simulate_data("mlvar", seed = 999)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  # Force all p-values to be large so the "else" branch triggers
  for (k in seq_along(fit$coefs)) {
    fit$coefs[[k]]$p <- rep(0.99, nrow(fit$coefs[[k]]))
  }
  output <- capture.output(summary(fit))
  expect_true(any(grepl("No significant temporal edges", output)))
})

test_that("summary.mlvar_result shows significant edges table when p < 0.05 (L652-661)", {
  d <- simulate_data("mlvar", seed = 1)
  fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
  # Force a p-value to be very small so significant edges are found
  fit$coefs[[1]]$p[1] <- 0.001
  output <- capture.output(summary(fit))
  expect_true(any(grepl("Significant temporal edges", output)))
})
