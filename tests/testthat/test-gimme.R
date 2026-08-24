# ===========================================================================
# Tests for build_gimme() — GIMME Network Analysis
#
# Since 0.9.0 the search is delegated to idiographic::fit_gimme(), whose
# algorithm is upstream-gimme-exact (>= 10.0) and carries its own kernel
# test suite. What Nestimate owns — and what this file pins — is the
# public surface: the argument contract, the error surface, the net_gimme
# field contract, seed reproducibility, and that print/summary/plot
# dispatch (to idiographic's methods) without error.
# ===========================================================================

# gimme fits a lavaan SEM per subject and is inherently slow. Under CRAN's
# 10-min check budget that's half the budget on one file. Skip on CRAN; the
# full suite runs locally and in CI.
testthat::skip_on_cran()

# --- Helper: generate test data ---
.make_gimme_data <- function(n_subjects = 10, n_time = 80, n_vars = 3,
                              seed = 42) {
  set.seed(seed)
  vars <- paste0("V", seq_len(n_vars))
  data_list <- lapply(seq_len(n_subjects), function(i) {
    # Simple AR(1) + some cross-effects; the loop is the natural recursion
    mat <- matrix(0, n_time, n_vars)
    mat[1, ] <- stats::rnorm(n_vars)
    for (t in 2:n_time) {
      mat[t, ] <- 0.3 * mat[t - 1, ] + stats::rnorm(n_vars, sd = 0.7)
      if (i <= n_subjects / 2) {
        mat[t, 2] <- mat[t, 2] + 0.2 * mat[t - 1, 1]
      }
    }
    df <- as.data.frame(mat)
    colnames(df) <- vars
    df$id <- i
    df$time <- seq_len(n_time)
    df
  })
  long_data <- do.call(rbind, data_list)
  list(data = long_data, vars = vars)
}


# ===========================================================================
# Section 1: Input validation (Nestimate's own error surface)
# ===========================================================================
test_that("gimme rejects non-data.frame input", {
  expect_error(build_gimme(matrix(1:10, 2, 5), vars = "V1", id = "id"),
               "data.frame")
})

test_that("gimme rejects missing variables", {
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = c("V1", "NONEXISTENT"), id = "id"),
               "not found")
})

test_that("gimme rejects single variable", {
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = "V1", id = "id"))
})

test_that("gimme rejects missing id column", {
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "nope"))
})

test_that("gimme rejects single subject", {
  sim <- .make_gimme_data(n_subjects = 1, n_time = 50)
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id"),
               "at least 2 individuals")
})

test_that("build_gimme rejects invalid groupcutoff", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id",
                            groupcutoff = 0), "groupcutoff")
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id",
                            groupcutoff = 1.5), "groupcutoff")
})

test_that("build_gimme rejects exogenous names outside vars", {
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id",
                           exogenous = "NOPE"), "exogenous")
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id",
                           exogenous = sim$vars), "endogenous")
})


# ===========================================================================
# Section 2: Return contract
# ===========================================================================
test_that("gimme returns net_gimme class", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_s3_class(res, "net_gimme")
})

test_that("gimme keeps the full pre-delegation field contract", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  legacy_fields <- c("temporal", "temporal_avg", "contemporaneous",
                     "contemporaneous_avg", "coefs", "psi", "fit",
                     "path_counts", "paths", "group_paths",
                     "individual_paths", "syntax", "labels",
                     "n_subjects", "n_obs", "config")
  expect_true(all(legacy_fields %in% names(res)))
})

test_that("gimme temporal and contemporaneous matrices have correct dimensions", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  p <- length(sim$vars)
  expect_identical(dim(res$temporal), c(p, p))
  expect_identical(dim(res$contemporaneous), c(p, p))
  expect_identical(rownames(res$temporal), sim$vars)
})

test_that("gimme per-person coef matrices have correct dimensions", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  p <- length(sim$vars)
  expect_length(res$coefs, res$n_subjects)
  for (m in res$coefs) expect_identical(dim(m), c(p, 2L * p))
})


# ===========================================================================
# Section 3: AR paths
# ===========================================================================
test_that("gimme with ar=TRUE has autoregressive paths for all subjects", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     ar = TRUE, seed = 1)
  # AR entries of path_counts (diagonal of the lagged block) equal n_subjects
  ar_counts <- diag(res$path_counts[, seq_along(sim$vars), drop = FALSE])
  expect_true(all(ar_counts == res$n_subjects))
})

test_that("gimme with ar=FALSE does not force autoregressive paths", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     ar = FALSE, seed = 1)
  expect_s3_class(res, "net_gimme")
})


# ===========================================================================
# Section 4: Path structures
# ===========================================================================
test_that("gimme path counts are non-negative integers bounded by n_subjects", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_true(all(res$path_counts >= 0))
  expect_true(all(res$path_counts <= res$n_subjects))
  expect_true(all(res$path_counts == round(res$path_counts)))
})

test_that("gimme group_paths is character and individual_paths a named list", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_type(res$group_paths, "character")
  expect_type(res$individual_paths, "list")
  expect_length(res$individual_paths, res$n_subjects)
})


# ===========================================================================
# Section 5: Fit indices
# ===========================================================================
test_that("gimme fit data.frame has correct structure and valid ranges", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_s3_class(res$fit, "data.frame")
  expect_identical(nrow(res$fit), as.integer(res$n_subjects))
  expect_true(all(c("rmsea", "srmr", "cfi") %in% names(res$fit)))
  ok <- stats::complete.cases(res$fit[, c("rmsea", "cfi")])
  expect_true(all(res$fit$rmsea[ok] >= 0))
  expect_true(all(res$fit$cfi[ok] <= 1 + 1e-8))
})


# ===========================================================================
# Section 6: Reproducibility
# ===========================================================================
test_that("gimme produces identical results with same seed", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 5, n_time = 60, seed = 3)
  r1 <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                    seed = 7)
  r2 <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                    seed = 7)
  expect_identical(r1$path_counts, r2$path_counts)
  expect_identical(r1$group_paths, r2$group_paths)
  expect_equal(r1$coefs, r2$coefs, tolerance = 0)
})


# ===========================================================================
# Section 7: S3 dispatch (idiographic's methods since 0.9.0)
# ===========================================================================
test_that("print, summary, and plot dispatch without error", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_output(print(res))
  # idiographic's summary returns a tidy one-row-per-network data.frame
  expect_s3_class(summary(res), "data.frame")
  expect_no_error(plot(res))
  expect_no_error(plot(res, layer = "temporal"))
})


# ===========================================================================
# Section 8: Argument pass-through
# ===========================================================================
test_that("gimme works with 4 variables", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 70, n_vars = 4, seed = 2)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 2)
  expect_identical(dim(res$temporal), c(4L, 4L))
})

test_that("gimme stores config correctly", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     ar = FALSE, groupcutoff = 0.6, seed = 1)
  expect_false(res$config$ar)
  expect_identical(res$config$groupcutoff, 0.6)
})

test_that("build_gimme hybrid=TRUE runs and returns the contract", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 5, n_time = 60, seed = 4)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     hybrid = TRUE, seed = 4)
  expect_s3_class(res, "net_gimme")
})

test_that("build_gimme standardize=TRUE runs without error", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 5, n_time = 60, seed = 5)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     standardize = TRUE, seed = 5)
  expect_s3_class(res, "net_gimme")
})

test_that("build_gimme runs with exogenous specification", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  expect_no_error(
    res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                       exogenous = "V3", ar = TRUE, seed = 1)
  )
  expect_s3_class(res, "net_gimme")
})

test_that("build_gimme works without time argument", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", seed = 1)
  expect_s3_class(res, "net_gimme")
})


# ===========================================================================
# Section 9: Degenerate subjects
# ===========================================================================
test_that("build_gimme drops subjects with too few time points", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 5, n_time = 60, seed = 1)
  bad_subj <- data.frame(V1 = c(1, 2), V2 = c(1, 2), V3 = c(1, 2),
                          id = 999, time = 1:2)
  combined <- rbind(sim$data, bad_subj)
  res <- suppressWarnings(
    build_gimme(combined, vars = sim$vars, id = "id", time = "time",
                seed = 1)
  )
  expect_equal(res$n_subjects, 5)
  expect_false("999" %in% names(res$coefs))
})

test_that("build_gimme errors when every subject is too short", {
  skip_if_not_installed("lavaan")
  short <- data.frame(V1 = rep(1:2, 2), V2 = rep(1:2, 2),
                      id = rep(1:2, each = 2), time = rep(1:2, 2))
  expect_error(
    suppressWarnings(
      build_gimme(short, vars = c("V1", "V2"), id = "id", time = "time")
    )
  )
})
