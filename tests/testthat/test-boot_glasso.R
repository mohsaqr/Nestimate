# ---- Tests for boot_glasso() ----

# ---- Helper: generate structured data with known edges ----
.make_test_data <- function(n = 200, p = 5, seed = 42) {
  set.seed(seed)
  # Create correlated data via Cholesky decomposition
  Sigma <- diag(p)
  # Add some correlations between adjacent variables
  for (i in seq_len(p - 1)) {
    Sigma[i, i + 1] <- 0.4
    Sigma[i + 1, i] <- 0.4
  }
  L <- chol(Sigma)
  mat <- matrix(rnorm(n * p), n, p) %*% L
  df <- as.data.frame(mat)
  names(df) <- paste0("V", seq_len(p))
  df
}

# Small run settings for speed
SMALL_ITER <- 20L
SMALL_CS_ITER <- 10L
SMALL_CS_DROP <- c(0.25, 0.5, 0.75)
FAST_CENT <- c("strength", "expected_influence")


# ========================================================================
# Input validation
# ========================================================================

test_that("boot_glasso rejects non-numeric input", {
  df <- data.frame(A = letters[1:10], B = LETTERS[1:10])
  expect_error(boot_glasso(df, iter = 10), "2 numeric")
})

test_that("boot_glasso rejects too few rows", {
  df <- data.frame(V1 = 1:2, V2 = 3:4)
  expect_error(boot_glasso(df, iter = 10), "3 complete rows|3 observations")
})

test_that("boot_glasso rejects too few columns", {
  df <- data.frame(V1 = rnorm(20))
  expect_error(boot_glasso(df, iter = 10), "2 numeric|2 variable")
})

test_that("boot_glasso rejects invalid iter", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 1))
  expect_error(boot_glasso(df, iter = "abc"))
})

test_that("boot_glasso rejects invalid alpha", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, alpha = 0))
  expect_error(boot_glasso(df, iter = 10, alpha = 1))
  expect_error(boot_glasso(df, iter = 10, alpha = -0.5))
})

test_that("boot_glasso rejects invalid gamma", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, gamma = -1))
})

test_that("boot_glasso rejects invalid nlambda", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, nlambda = 1))
})

test_that("boot_glasso rejects invalid cor_method", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, cor_method = "invalid"))
})

test_that("boot_glasso rejects invalid centrality measures", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, centrality = "invalid"))
})

test_that("boot_glasso rejects invalid cs_drop values", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, cs_drop = c(0, 0.5)))
  expect_error(boot_glasso(df, iter = 10, cs_drop = c(0.5, 1)))
})

test_that("boot_glasso rejects invalid ncores", {
  df <- .make_test_data(50, 3)
  expect_error(boot_glasso(df, iter = 10, ncores = 0))
})

test_that("boot_glasso rejects non-glasso netobject", {
  df <- .make_test_data(50, 3)
  net <- build_network(df, method = "cor")
  expect_error(boot_glasso(net, iter = 10), "method.*glasso")
})

test_that("boot_glasso rejects invalid input types", {
  expect_error(boot_glasso("not a data frame", iter = 10), "data frame")
  expect_error(boot_glasso(list(a = 1), iter = 10), "data frame")
})


# ========================================================================
# Structure
# ========================================================================

test_that("boot_glasso returns correct class", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_s3_class(result, "boot_glasso")
})

test_that("boot_glasso result has all expected fields", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)

  expected_fields <- c(
    "original_pcor", "original_precision", "original_centrality",
    "original_predictability", "edge_ci", "edge_inclusion",
    "thresholded_pcor", "centrality_ci", "cs_coefficient", "cs_data",
    "edge_diff_p", "centrality_diff_p", "predictability_ci",
    "boot_edges", "boot_centrality", "boot_predictability",
    "nodes", "n", "p", "iter", "cs_iter", "cs_drop", "alpha",
    "gamma", "nlambda", "centrality_measures", "cor_method",
    "lambda_path", "lambda_selected", "timing"
  )
  for (f in expected_fields) {
    expect_true(f %in% names(result), info = sprintf("Missing field: %s", f))
  }
})

test_that("boot_glasso matrices have correct dimensions", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- 4
  n_upper <- p * (p - 1) / 2

  expect_equal(nrow(result$original_pcor), p)
  expect_equal(ncol(result$original_pcor), p)
  expect_equal(nrow(result$original_precision), p)
  expect_equal(ncol(result$original_precision), p)
  expect_equal(nrow(result$thresholded_pcor), p)
  expect_equal(ncol(result$thresholded_pcor), p)
  expect_equal(nrow(result$boot_edges), SMALL_ITER)
  expect_equal(ncol(result$boot_edges), n_upper)
  expect_equal(nrow(result$boot_predictability), SMALL_ITER)
  expect_equal(ncol(result$boot_predictability), p)
})

test_that("boot_glasso centrality dimensions are correct", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- 4

  for (m in FAST_CENT) {
    expect_equal(length(result$original_centrality[[m]]), p)
    expect_equal(nrow(result$boot_centrality[[m]]), SMALL_ITER)
    expect_equal(ncol(result$boot_centrality[[m]]), p)
  }
})

test_that("boot_glasso edge_ci has correct structure", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  ci <- result$edge_ci
  expect_true(is.data.frame(ci))
  expect_true(all(c("edge", "weight", "ci_lower", "ci_upper", "inclusion")
                    %in% names(ci)))
  expect_equal(nrow(ci), 4 * 3 / 2)
})

test_that("boot_glasso stores correct metadata", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         alpha = 0.10, gamma = 0.25, seed = 1)
  expect_equal(result$iter, SMALL_ITER)
  expect_equal(result$cs_iter, SMALL_CS_ITER)
  expect_equal(result$alpha, 0.10)
  expect_equal(result$gamma, 0.25)
  expect_equal(result$p, 4)
  expect_equal(result$n, 100)
  expect_equal(result$cor_method, "pearson")
  expect_equal(result$centrality_measures, FAST_CENT)
})

test_that("boot_glasso timing vector has correct names", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_true(all(c("total", "phase1", "bootstrap", "case_drop", "statistics")
                    %in% names(result$timing)))
  expect_true(all(result$timing >= 0))
})

test_that("boot_glasso lambda_path is decreasing", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  lp <- result$lambda_path
  expect_true(all(diff(lp) < 0))
  expect_true(length(lp) == result$nlambda)
})


# ========================================================================
# Correctness
# ========================================================================

test_that("boot_glasso original_pcor is symmetric with zero diagonal", {
  df <- .make_test_data(150, 5)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  pcor <- result$original_pcor
  expect_true(isSymmetric(pcor, tol = 1e-10))
  expect_equal(unname(diag(pcor)), rep(0, 5))
})

test_that("boot_glasso CIs bracket original edge weights (mostly)", {
  df <- .make_test_data(200, 5)
  result <- boot_glasso(df, iter = 100, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 42)
  ci <- result$edge_ci
  # At alpha = 0.05, we expect ~95% of original weights within CIs
  within_ci <- ci$weight >= ci$ci_lower & ci$weight <= ci$ci_upper
  expect_true(mean(within_ci) >= 0.7)  # Conservative check
})

test_that("boot_glasso inclusion probabilities are in [0, 1]", {
  df <- .make_test_data(150, 5)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_true(all(result$edge_inclusion >= 0))
  expect_true(all(result$edge_inclusion <= 1))
})

test_that("boot_glasso CS-coefficient is in [0, max(cs_drop)]", {
  df <- .make_test_data(200, 5)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  for (m in FAST_CENT) {
    cs <- result$cs_coefficient[[m]]
    expect_true(cs >= 0)
    expect_true(cs <= max(SMALL_CS_DROP))
  }
})

test_that("boot_glasso edge_diff_p is symmetric with zero diagonal", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  if (!is.null(result$edge_diff_p)) {
    p_mat <- result$edge_diff_p
    expect_true(isSymmetric(p_mat, tol = 1e-10))
    expect_equal(unname(diag(p_mat)), rep(0, ncol(p_mat)))
    expect_true(all(p_mat >= 0))
    expect_true(all(p_mat <= 1))
  }
})

test_that("boot_glasso centrality_diff_p is symmetric with zero diagonal", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  for (m in FAST_CENT) {
    p_mat <- result$centrality_diff_p[[m]]
    expect_true(isSymmetric(p_mat, tol = 1e-10))
    expect_equal(unname(diag(p_mat)), rep(0, ncol(p_mat)))
    expect_true(all(p_mat >= 0))
    expect_true(all(p_mat <= 1))
  }
})

test_that("boot_glasso thresholded network has fewer or equal edges", {
  df <- .make_test_data(200, 5)
  result <- boot_glasso(df, iter = 50, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 42)
  orig_edges <- sum(result$original_pcor[upper.tri(result$original_pcor)] != 0)
  thresh_edges <- sum(result$thresholded_pcor[upper.tri(
    result$thresholded_pcor)] != 0)
  expect_true(thresh_edges <= orig_edges)
})

test_that("boot_glasso predictability is in [0, 1]", {
  df <- .make_test_data(150, 5)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_true(all(result$original_predictability >= 0))
  expect_true(all(result$original_predictability <= 1))
  # Bootstrap predictability also in [0,1]
  valid_pred <- result$boot_predictability[!is.na(result$boot_predictability)]
  expect_true(all(valid_pred >= 0))
  expect_true(all(valid_pred <= 1))
})

test_that("boot_glasso strength equals rowSums(abs(pcor))", {
  df <- .make_test_data(150, 5)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP,
                         centrality = c("strength"), seed = 1)
  expected_strength <- rowSums(abs(result$original_pcor))
  expect_equal(unname(result$original_centrality$strength),
               unname(expected_strength), tolerance = 1e-10)
})

test_that("boot_glasso expected_influence equals rowSums(pcor)", {
  df <- .make_test_data(150, 5)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP,
                         centrality = c("expected_influence"), seed = 1)
  expected_ei <- rowSums(result$original_pcor)
  expect_equal(unname(result$original_centrality$expected_influence),
               unname(expected_ei), tolerance = 1e-10)
})

test_that("boot_glasso cs_data has correct structure", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  cs <- result$cs_data
  expect_true(is.data.frame(cs))
  expect_true(all(c("drop_prop", "measure", "mean_cor") %in% names(cs)))
  expect_equal(nrow(cs), length(SMALL_CS_DROP) * length(FAST_CENT))
  expect_true(all(cs$drop_prop %in% SMALL_CS_DROP))
})

test_that("boot_glasso predictability_ci has correct structure", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  pred <- result$predictability_ci
  expect_true(is.data.frame(pred))
  expect_true(all(c("node", "r2", "ci_lower", "ci_upper") %in% names(pred)))
  expect_equal(nrow(pred), 4)
  expect_true(all(pred$ci_lower <= pred$r2 | abs(pred$ci_lower - pred$r2) < 0.01))
  expect_true(all(pred$ci_upper >= pred$r2 | abs(pred$ci_upper - pred$r2) < 0.01))
})


# ========================================================================
# Reproducibility
# ========================================================================

test_that("boot_glasso with same seed gives identical results", {
  df <- .make_test_data(100, 4)
  r1 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = 123)
  r2 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = 123)
  expect_equal(r1$boot_edges, r2$boot_edges)
  expect_equal(r1$cs_coefficient, r2$cs_coefficient)
  expect_equal(r1$original_pcor, r2$original_pcor)
})

test_that("boot_glasso with different seeds gives different results", {
  df <- .make_test_data(100, 4)
  r1 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = 1)
  r2 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = 2)
  expect_false(identical(r1$boot_edges, r2$boot_edges))
})

test_that("boot_glasso with NULL seed uses random state", {
  df <- .make_test_data(100, 4)
  set.seed(999)
  r1 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = NULL)
  set.seed(888)
  r2 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = NULL)
  # Different seeds should give different results
  expect_false(identical(r1$boot_edges, r2$boot_edges))
})


# ========================================================================
# netobject input
# ========================================================================

test_that("boot_glasso accepts glasso netobject", {
  df <- .make_test_data(100, 4)
  net <- build_network(df, method = "glasso",
                        params = list(gamma = 0.5, nlambda = 100L))
  result <- boot_glasso(net, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$gamma, 0.5)
  expect_equal(result$p, 4)
})

test_that("boot_glasso extracts params from netobject", {
  df <- .make_test_data(100, 4)
  net <- build_network(df, method = "glasso",
                        params = list(gamma = 0.25, nlambda = 50L,
                                       cor_method = "spearman"))
  result <- boot_glasso(net, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_equal(result$gamma, 0.25)
  expect_equal(result$nlambda, 50L)
  expect_equal(result$cor_method, "spearman")
})

test_that("boot_glasso from netobject matches dataframe input", {
  df <- .make_test_data(100, 4)
  net <- build_network(df, method = "glasso",
                        params = list(gamma = 0.5, nlambda = 100L))

  r1 <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = 42, gamma = 0.5, nlambda = 100L)
  r2 <- boot_glasso(net, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                     cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                     seed = 42)

  # Original pcors should match
  expect_equal(r1$original_pcor, r2$original_pcor, tolerance = 1e-8)
  # Boot edges should match since same seed
  expect_equal(r1$boot_edges, r2$boot_edges, tolerance = 1e-8)
})


# ========================================================================
# S3 methods
# ========================================================================

test_that("print.boot_glasso runs without error", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_output(print(result), "GLASSO Bootstrap")
  expect_output(print(result), "Centrality Stability")
  expect_output(print(result), "Timing")
})

test_that("print.boot_glasso shows correct iteration count", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  output <- capture.output(print(result))
  expect_true(any(grepl(as.character(SMALL_ITER), output)))
})

test_that("print.boot_glasso shows CS labels", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  output <- capture.output(print(result))
  # Should show one of the labels
  expect_true(any(grepl("Stable|Acceptable|Unstable", output)))
})

test_that("summary.boot_glasso returns edges by default", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  s <- summary(result)
  expect_true(is.data.frame(s))
  expect_true("edge" %in% names(s))
})

test_that("summary.boot_glasso type='centrality' returns list", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  s <- summary(result, type = "centrality")
  expect_true(is.list(s))
  expect_equal(names(s), FAST_CENT)
  for (m in FAST_CENT) {
    expect_true(is.data.frame(s[[m]]))
    expect_true(all(c("node", "value", "ci_lower", "ci_upper")
                      %in% names(s[[m]])))
  }
})

test_that("summary.boot_glasso type='cs' returns cs_data", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  s <- summary(result, type = "cs")
  expect_true(is.data.frame(s))
  expect_true("drop_prop" %in% names(s))
})

test_that("summary.boot_glasso type='predictability' works", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  s <- summary(result, type = "predictability")
  expect_true(is.data.frame(s))
  expect_true(all(c("node", "r2", "ci_lower", "ci_upper") %in% names(s)))
})

test_that("summary.boot_glasso type='all' returns list of all", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  s <- summary(result, type = "all")
  expect_true(is.list(s))
  expect_true(all(c("edges", "centrality", "cs", "predictability")
                    %in% names(s)))
})

test_that("plot.boot_glasso type='edges' produces ggplot", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "edges")
  expect_s3_class(p, "ggplot")
})

test_that("plot.boot_glasso type='stability' produces ggplot", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "stability")
  expect_s3_class(p, "ggplot")
})

test_that("plot.boot_glasso type='edge_diff' produces ggplot", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "edge_diff")
  expect_s3_class(p, "ggplot")
})

test_that("plot.boot_glasso type='centrality_diff' produces ggplot", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "centrality_diff")
  expect_s3_class(p, "ggplot")
})

test_that("plot.boot_glasso edge_diff order='sample' is default", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p1 <- plot(result, type = "edge_diff")
  p2 <- plot(result, type = "edge_diff", order = "sample")
  # Both should be identical ggplot objects (same default)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  # Same x-axis levels
  expect_equal(levels(p1$data$edge1), levels(p2$data$edge1))
})

test_that("plot.boot_glasso edge_diff order='id' sorts alphabetically", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "edge_diff", order = "id")
  expect_s3_class(p, "ggplot")
  # Levels should be alphabetically sorted
  edge_levels <- levels(p$data$edge1)
  expect_equal(edge_levels, sort(edge_levels))
})

test_that("plot.boot_glasso edge_diff uses discrete fill", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "edge_diff")
  # Data should have fill column with discrete values
  expect_true("fill" %in% names(p$data))
  expect_true(all(p$data$fill %in% c("significant", "non-significant", "diagonal", "blank")))
})

test_that("plot.boot_glasso edge_diff full matrix has n^2 tiles", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "edge_diff")
  n_edges <- ncol(result$edge_diff_p)
  expect_equal(nrow(p$data), n_edges * n_edges)
})

test_that("plot.boot_glasso centrality_diff order='sample' is default", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p1 <- plot(result, type = "centrality_diff")
  p2 <- plot(result, type = "centrality_diff", order = "sample")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_equal(levels(p1$data$node1), levels(p2$data$node1))
})

test_that("plot.boot_glasso centrality_diff order='id' sorts alphabetically", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "centrality_diff", order = "id")
  expect_s3_class(p, "ggplot")
  node_levels <- levels(p$data$node1)
  expect_equal(node_levels, sort(node_levels))
})

test_that("plot.boot_glasso centrality_diff uses discrete fill", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "centrality_diff")
  expect_true("fill" %in% names(p$data))
  expect_true(all(p$data$fill %in% c("significant", "non-significant", "diagonal", "blank")))
})

test_that("plot.boot_glasso centrality_diff full matrix has p^2 tiles", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "centrality_diff")
  n_nodes <- result$p
  expect_equal(nrow(p$data), n_nodes * n_nodes)
})

test_that("plot.boot_glasso type='inclusion' produces ggplot", {
  skip_if_not_installed("ggplot2")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  p <- plot(result, type = "inclusion")
  expect_s3_class(p, "ggplot")
})

test_that("plot.boot_glasso type='network' runs without error", {
  skip_if_not_installed("cograph")
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_no_error(plot(result, type = "network"))
})


# ========================================================================
# Edge cases
# ========================================================================

test_that("boot_glasso works with p=2", {
  set.seed(42)
  df <- data.frame(A = rnorm(100), B = 0.5 * rnorm(100) + rnorm(100))
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$p, 2)
  expect_equal(ncol(result$boot_edges), 1)  # Only 1 edge for p=2
})

test_that("boot_glasso works with small n", {
  set.seed(42)
  df <- data.frame(V1 = rnorm(20), V2 = rnorm(20), V3 = rnorm(20))
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = c(0.25, 0.5), centrality = FAST_CENT,
                         seed = 1)
  expect_s3_class(result, "boot_glasso")
})

test_that("boot_glasso works with iter=2 (minimum)", {
  df <- .make_test_data(50, 3)
  result <- boot_glasso(df, iter = 2, cs_iter = 2,
                         cs_drop = c(0.5), centrality = FAST_CENT, seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(nrow(result$boot_edges), 2)
})

test_that("boot_glasso works with only strength centrality", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP,
                         centrality = "strength", seed = 1)
  expect_equal(result$centrality_measures, "strength")
  expect_equal(length(result$cs_coefficient), 1)
  expect_true("strength" %in% names(result$cs_coefficient))
})

test_that("boot_glasso works with single cs_drop value", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = 0.5, centrality = FAST_CENT, seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$cs_drop, 0.5)
  expect_equal(nrow(result$cs_data), length(FAST_CENT))
})

test_that("boot_glasso handles data with NAs (dropped rows)", {
  set.seed(42)
  df <- .make_test_data(100, 4)
  df$V1[1:5] <- NA
  expect_message(
    result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                           cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                           seed = 1),
    "Dropping.*rows"
  )
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$n, 95)
})

test_that("boot_glasso handles matrix input", {
  set.seed(42)
  mat <- matrix(rnorm(200 * 5), 200, 5)
  colnames(mat) <- paste0("V", 1:5)
  result <- boot_glasso(mat, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$p, 5)
})


# ========================================================================
# Parallel execution
# ========================================================================

test_that("boot_glasso with ncores=2 runs without error", {
  skip_on_os("windows")  # mclapply is fork-based, not available on Windows
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         ncores = 2, seed = 1)
  expect_s3_class(result, "boot_glasso")
})

test_that("boot_glasso parallel gives similar structure to sequential", {
  skip_on_os("windows")
  df <- .make_test_data(100, 4)
  r_seq <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                        cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                        ncores = 1, seed = 1)
  r_par <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                        cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                        ncores = 2, seed = 1)
  # Same original network (seed applied before bootstrap)
  expect_equal(r_seq$original_pcor, r_par$original_pcor)
  # Same dimensions
  expect_equal(dim(r_seq$boot_edges), dim(r_par$boot_edges))
})


# ========================================================================
# Correlation methods
# ========================================================================

test_that("boot_glasso works with spearman correlation", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         cor_method = "spearman", seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$cor_method, "spearman")
})

test_that("boot_glasso works with kendall correlation", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         cor_method = "kendall", seed = 1)
  expect_s3_class(result, "boot_glasso")
  expect_equal(result$cor_method, "kendall")
})


# ========================================================================
# All four centrality measures
# ========================================================================

test_that("boot_glasso computes all four centrality measures", {
  df <- .make_test_data(100, 4)
  all_cent <- c("strength", "expected_influence", "closeness", "betweenness")
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = all_cent,
                         seed = 1)
  expect_equal(result$centrality_measures, all_cent)
  expect_equal(length(result$original_centrality), 4)
  expect_equal(length(result$cs_coefficient), 4)
  expect_equal(length(result$centrality_ci), 4)
  expect_equal(length(result$centrality_diff_p), 4)
})

test_that("boot_glasso closeness and betweenness are non-negative", {
  df <- .make_test_data(100, 4)
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP,
                         centrality = c("closeness", "betweenness"),
                         seed = 1)
  # Some values may be NaN/Inf for disconnected networks, check finite ones
  cl <- result$original_centrality$closeness
  bt <- result$original_centrality$betweenness
  expect_true(all(cl[is.finite(cl)] >= 0))
  expect_true(all(bt[is.finite(bt)] >= 0))
})


# ========================================================================
# Integration: matches build_network(method="glasso")
# ========================================================================

test_that("boot_glasso original matches build_network glasso", {
  df <- .make_test_data(150, 5)
  net <- build_network(df, method = "glasso",
                        params = list(gamma = 0.5, nlambda = 100L))
  result <- boot_glasso(df, iter = SMALL_ITER, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         gamma = 0.5, nlambda = 100L, seed = 1)

  # The partial correlation matrices should match closely
  # (both use the same lambda path and EBIC selection)
  expect_equal(result$original_pcor, net$matrix, tolerance = 1e-6)
})

test_that("boot_glasso CS mean correlations are high for structured data", {
  # Strong edges should give stable centrality (high mean cors)
  set.seed(42)
  n <- 300
  x1 <- rnorm(n)
  x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5)
  x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  result <- boot_glasso(df, iter = 100, cs_iter = 50,
                         cs_drop = seq(0.1, 0.9, by = 0.1),
                         centrality = c("strength"), seed = 42)

  # With strong structure, mean correlation at 50% drop should be > 0.7
  cs_50 <- result$cs_data[result$cs_data$drop_prop == 0.5, "mean_cor"]
  expect_true(cs_50 > 0.7)

  # CS coefficient should be numeric and within valid range
  cs_val <- result$cs_coefficient[["strength"]]
  expect_true(is.numeric(cs_val))
  expect_true(cs_val >= 0)
  expect_true(cs_val <= 0.9)
})

test_that("boot_glasso detects known strong edges", {
  set.seed(42)
  n <- 300
  x1 <- rnorm(n)
  x2 <- 0.8 * x1 + rnorm(n, sd = 0.3)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4)

  result <- boot_glasso(df, iter = 100, cs_iter = SMALL_CS_ITER,
                         cs_drop = SMALL_CS_DROP, centrality = FAST_CENT,
                         seed = 42)

  # V1--V2 should have high inclusion probability
  v1v2_name <- "V1 -- V2"
  v1v2_incl <- result$edge_inclusion[v1v2_name]
  expect_true(v1v2_incl > 0.9)

  # V1--V2 should be in thresholded network
  expect_true(abs(result$thresholded_pcor["V1", "V2"]) > 0)
})


# ========================================================================
# Internal helper tests
# ========================================================================

test_that(".bg_upper_tri_indices returns correct indices", {
  ut <- .bg_upper_tri_indices(4)
  expect_equal(length(ut$row_idx), 6)  # 4*3/2
  expect_equal(length(ut$col_idx), 6)
  # All row < col
  expect_true(all(ut$row_idx < ut$col_idx))
})

test_that(".bg_build_edge_names creates correct labels", {
  nodes <- c("A", "B", "C")
  ut <- .bg_upper_tri_indices(3)
  names <- .bg_build_edge_names(nodes, ut$row_idx, ut$col_idx)
  expect_equal(names, c("A -- B", "A -- C", "B -- C"))
})

test_that(".bg_cs_label returns correct labels", {
  expect_equal(.bg_cs_label(0.7), "Stable")
  expect_equal(.bg_cs_label(0.5), "Stable")
  expect_equal(.bg_cs_label(0.4), "Acceptable")
  expect_equal(.bg_cs_label(0.25), "Acceptable")
  expect_equal(.bg_cs_label(0.2), "Unstable")
  expect_equal(.bg_cs_label(0), "Unstable")
  expect_equal(.bg_cs_label(NA), "Unknown")
})

test_that(".bg_cs_coefficient returns 0 when all below threshold", {
  # cors_per_prop is a list of vectors (one per drop proportion)
  cors <- list(rep(0.3, 10), rep(0.3, 10), rep(0.3, 10))  # All below 0.7
  cs_drop <- c(0.25, 0.5, 0.75)
  result <- .bg_cs_coefficient(cors, cs_drop, 0.7, 0.95)
  expect_equal(result, 0)
})

test_that(".bg_cs_coefficient returns max drop_prop when all above", {
  cors <- list(rep(0.9, 10), rep(0.9, 10), rep(0.9, 10))  # All above 0.7
  cs_drop <- c(0.25, 0.5, 0.75)
  result <- .bg_cs_coefficient(cors, cs_drop, 0.7, 0.95)
  expect_equal(result, 0.75)
})

test_that(".bg_threshold_network zeros edges with CI spanning zero", {
  pcor <- matrix(c(0, 0.3, -0.2, 0.3, 0, 0.1, -0.2, 0.1, 0), 3, 3)
  colnames(pcor) <- rownames(pcor) <- c("A", "B", "C")
  ut <- .bg_upper_tri_indices(3)

  # CI for A-B: [0.1, 0.5] -> significant
  # CI for A-C: [-0.4, 0.1] -> spans zero
  # CI for B-C: [-0.1, 0.3] -> spans zero
  ci_lower <- c(0.1, -0.4, -0.1)
  ci_upper <- c(0.5, 0.1, 0.3)

  result <- .bg_threshold_network(pcor, ci_lower, ci_upper, ut)
  expect_true(result["A", "B"] != 0)
  expect_equal(result["A", "C"], 0)
  expect_equal(result["B", "C"], 0)
})

test_that(".bg_edge_diff_test returns correct dimensions", {
  set.seed(42)
  boot_edges <- matrix(rnorm(50 * 6), 50, 6)
  colnames(boot_edges) <- paste0("e", 1:6)
  result <- .bg_edge_diff_test(boot_edges)
  expect_equal(nrow(result), 6)
  expect_equal(ncol(result), 6)
  expect_true(isSymmetric(result))
  expect_equal(unname(diag(result)), rep(0, 6))
})

test_that(".bg_centrality_diff_test returns correct dimensions", {
  set.seed(42)
  boot_cent <- matrix(rnorm(50 * 4), 50, 4)
  colnames(boot_cent) <- paste0("V", 1:4)
  result <- .bg_centrality_diff_test(boot_cent)
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
  expect_true(isSymmetric(result))
  expect_equal(unname(diag(result)), rep(0, 4))
})


# ========================================================================
# bootnet equivalence
# ========================================================================

test_that("boot_glasso original network matches bootnet EBICglasso", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  # Run across 5 seeds
  for (s in c(100, 200, 300, 400, 500)) {
    set.seed(s)
    n <- 200
    x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
    x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
    x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
    df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

    bn_net <- suppressWarnings(suppressMessages(
      bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                                corMethod = "cor", missing = "listwise",
                                verbose = FALSE)
    ))

    bg <- boot_glasso(df, iter = 10, cs_iter = 10, cs_drop = c(0.5),
                       centrality = "strength", gamma = 0.5, seed = 1)

    # Original networks should be very close (same EBIC/glasso)
    max_diff <- max(abs(bn_net$graph - bg$original_pcor))
    expect_true(max_diff < 0.01,
                info = sprintf("seed=%d, max_diff=%.6f", s, max_diff))
  }
})

test_that("boot_glasso edge CIs correlate highly with bootnet", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  bn_net <- suppressWarnings(suppressMessages(
    bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                              corMethod = "cor", missing = "listwise",
                              verbose = FALSE)
  ))

  set.seed(1)
  bn_boot <- suppressWarnings(suppressMessages(
    bootnet::bootnet(bn_net, nBoots = 200, type = "nonparametric",
                      nCores = 1, verbose = FALSE,
                      statistics = c("edge", "strength"))
  ))

  bg <- boot_glasso(df, iter = 200, cs_iter = 10, cs_drop = c(0.5),
                     centrality = "strength", gamma = 0.5, seed = 1)

  # Extract bootnet edge CIs
  bn_edges <- bn_boot$bootTable[bn_boot$bootTable$type == "edge", ]
  edge_ids <- unique(bn_edges$id)
  bn_ci <- do.call(rbind, lapply(edge_ids, function(eid) {
    vals <- bn_edges$value[bn_edges$id == eid]
    data.frame(edge = eid, ci_lower = quantile(vals, 0.025),
               ci_upper = quantile(vals, 0.975), stringsAsFactors = FALSE)
  }))

  bg_ci <- bg$edge_ci
  bg_ci$bn_id <- gsub(" -- ", "--", bg_ci$edge)
  merged <- merge(bn_ci, bg_ci, by.x = "edge", by.y = "bn_id")

  # Edge CIs should correlate highly (r > 0.99)
  edge_ci_r <- cor(merged$ci_lower.x, merged$ci_lower.y)
  expect_true(edge_ci_r > 0.99,
              info = sprintf("Edge CI r = %.4f", edge_ci_r))

  # Inclusion probabilities should correlate highly
  bn_incl <- vapply(edge_ids, function(eid) {
    mean(bn_edges$value[bn_edges$id == eid] != 0)
  }, numeric(1))
  bg_incl <- bg$edge_inclusion[gsub("--", " -- ", names(bn_incl))]
  incl_r <- cor(bn_incl, bg_incl)
  expect_true(incl_r > 0.99,
              info = sprintf("Inclusion r = %.4f", incl_r))
})

test_that("boot_glasso strength CIs correlate highly with bootnet", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  bn_net <- suppressWarnings(suppressMessages(
    bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                              corMethod = "cor", missing = "listwise",
                              verbose = FALSE)
  ))

  set.seed(1)
  bn_boot <- suppressWarnings(suppressMessages(
    bootnet::bootnet(bn_net, nBoots = 200, type = "nonparametric",
                      nCores = 1, verbose = FALSE,
                      statistics = c("edge", "strength"))
  ))

  bg <- boot_glasso(df, iter = 200, cs_iter = 10, cs_drop = c(0.5),
                     centrality = "strength", gamma = 0.5, seed = 1)

  bn_str <- bn_boot$bootTable[bn_boot$bootTable$type == "strength", ]
  bn_str_ci <- do.call(rbind, lapply(unique(bn_str$node1), function(nd) {
    vals <- bn_str$value[bn_str$node1 == nd]
    data.frame(node = nd, ci_lower = quantile(vals, 0.025),
               ci_upper = quantile(vals, 0.975), stringsAsFactors = FALSE)
  }))

  bg_str_ci <- bg$centrality_ci$strength
  str_merged <- merge(bn_str_ci, bg_str_ci, by = "node")
  str_ci_r <- cor(str_merged$ci_lower.x, str_merged$ci_lower.y)
  expect_true(str_ci_r > 0.99,
              info = sprintf("Strength CI r = %.4f", str_ci_r))
})

test_that("boot_glasso CS-coefficient matches bootnet corStability", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  # Run 5 seeds, require 80% agreement within 0.15
  cs_diffs <- numeric(5)

  for (i in seq_len(5)) {
    s <- i * 100
    set.seed(s)
    n <- 200
    x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
    x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
    x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
    df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

    bn_net <- suppressWarnings(suppressMessages(
      bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                                corMethod = "cor", missing = "listwise",
                                verbose = FALSE)
    ))

    set.seed(s + 1000)
    bn_case <- suppressWarnings(suppressMessages(
      bootnet::bootnet(bn_net, nBoots = 500, type = "case", nCores = 1,
                        verbose = FALSE, statistics = c("strength"),
                        caseMin = 0.05, caseMax = 0.75, caseN = 10)
    ))
    bn_cs <- suppressWarnings(suppressMessages(
      bootnet::corStability(bn_case, verbose = FALSE)
    ))

    bg <- boot_glasso(df, iter = 50, cs_iter = 500,
                       cs_drop = seq(0.05, 0.75, length.out = 10),
                       centrality = "strength", gamma = 0.5, seed = s)

    cs_diffs[i] <- abs(as.numeric(bn_cs["strength"]) -
                        as.numeric(bg$cs_coefficient["strength"]))
  }

  # At least 80% should match within 0.15
  pct_match <- mean(cs_diffs <= 0.15)
  expect_true(pct_match >= 0.8,
              info = sprintf("CS match: %.0f%%, diffs: %s",
                             pct_match * 100,
                             paste(round(cs_diffs, 2), collapse = ", ")))

  # Mean diff should be small
  expect_true(mean(cs_diffs) < 0.1,
              info = sprintf("Mean CS diff = %.3f", mean(cs_diffs)))
})

test_that("boot_glasso strength matches bootnet strength definition", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  bn_net <- suppressWarnings(suppressMessages(
    bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                              corMethod = "cor", missing = "listwise",
                              verbose = FALSE)
  ))

  bg <- boot_glasso(df, iter = 10, cs_iter = 10, cs_drop = c(0.5),
                     centrality = "strength", gamma = 0.5, seed = 1)

  # bootnet uses qgraph::centrality()$InDegree for strength
  bn_str <- qgraph::centrality(bn_net$graph)$InDegree
  bg_str <- bg$original_centrality$strength

  # Should match within tolerance (tiny diff from pcor rounding)
  expect_equal(unname(bg_str), unname(bn_str), tolerance = 0.01)
})
