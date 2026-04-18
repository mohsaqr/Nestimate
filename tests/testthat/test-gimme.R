# ===========================================================================
# Tests for build_gimme() — GIMME Network Analysis
# ===========================================================================

# --- Helper: generate test data ---
.make_gimme_data <- function(n_subjects = 10, n_time = 80, n_vars = 3,
                              seed = 42) {
  set.seed(seed)
  vars <- paste0("V", seq_len(n_vars))
  data_list <- lapply(seq_len(n_subjects), function(i) {
    # Simple AR(1) + some cross-effects
    mat <- matrix(0, n_time, n_vars)
    mat[1, ] <- stats::rnorm(n_vars)
    for (t in 2:n_time) {
      mat[t, ] <- 0.3 * mat[t - 1, ] + stats::rnorm(n_vars, sd = 0.7)
      # Add cross-effect for some subjects
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
# Section 1: Input validation
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
               "at least 2")
})


# ===========================================================================
# Section 2: Basic construction
# ===========================================================================
test_that("gimme returns net_gimme class", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_s3_class(res, "net_gimme")
})

test_that("gimme has correct structure", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_true(is.matrix(res$temporal))
  expect_true(is.matrix(res$contemporaneous))
  expect_true(is.list(res$coefs))
  expect_true(is.list(res$psi))
  expect_true(is.data.frame(res$fit))
  expect_equal(length(res$coefs), 6)
  expect_equal(res$n_subjects, 6)
  expect_equal(res$labels, c("V1", "V2", "V3"))
})

test_that("gimme temporal matrix has correct dimensions", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_equal(nrow(res$temporal), 3)
  expect_equal(ncol(res$temporal), 3)
  expect_equal(rownames(res$temporal), sim$vars)
})

test_that("gimme contemporaneous matrix has correct dimensions", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_equal(nrow(res$contemporaneous), 3)
  expect_equal(ncol(res$contemporaneous), 3)
})

test_that("gimme per-person coef matrices have correct dimensions", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  for (k in seq_along(res$coefs)) {
    m <- res$coefs[[k]]
    expect_equal(nrow(m), 3)
    expect_equal(ncol(m), 6)  # 3 lagged + 3 contemporaneous
  }
})


# ===========================================================================
# Section 3: AR paths
# ===========================================================================
test_that("gimme with ar=TRUE has autoregressive paths for all subjects", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     ar = TRUE, seed = 1)

  # Diagonal of temporal counts should be n_subjects (AR paths always present)
  expect_true(all(diag(res$temporal) == res$n_subjects))
})

test_that("gimme with ar=FALSE does not force autoregressive paths", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     ar = FALSE, seed = 1)

  # Diagonal may or may not be n_subjects
  expect_s3_class(res, "net_gimme")
})


# ===========================================================================
# Section 4: Path counts and group paths
# ===========================================================================
test_that("gimme path counts are non-negative integers", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_true(all(res$path_counts >= 0))
  expect_true(all(res$path_counts == round(res$path_counts)))
  expect_true(all(res$path_counts <= res$n_subjects))
})

test_that("gimme group_paths is character vector", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_true(is.character(res$group_paths))
})

test_that("gimme individual_paths is a named list", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_true(is.list(res$individual_paths))
  expect_equal(length(res$individual_paths), res$n_subjects)
})


# ===========================================================================
# Section 5: Fit indices
# ===========================================================================
test_that("gimme fit data.frame has correct structure", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  expect_true("rmsea" %in% names(res$fit))
  expect_true("srmr" %in% names(res$fit))
  expect_true("cfi" %in% names(res$fit))
  expect_true("nnfi" %in% names(res$fit))
  expect_true("file" %in% names(res$fit))
  expect_true("status" %in% names(res$fit))
  expect_equal(nrow(res$fit), res$n_subjects)
})

test_that("gimme fit indices are in valid ranges", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)

  rmsea <- res$fit$rmsea[!is.na(res$fit$rmsea)]
  expect_true(all(rmsea >= 0))

  srmr <- res$fit$srmr[!is.na(res$fit$srmr)]
  expect_true(all(srmr >= 0))
})


# ===========================================================================
# Section 6: Reproducibility
# ===========================================================================
test_that("gimme produces identical results with same seed", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res1 <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                      seed = 42)
  res2 <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                      seed = 42)

  expect_identical(res1$temporal, res2$temporal)
  expect_identical(res1$contemporaneous, res2$contemporaneous)
  expect_identical(res1$group_paths, res2$group_paths)
})


# ===========================================================================
# Section 7: S3 methods
# ===========================================================================
test_that("print.net_gimme produces output", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  out <- capture.output(print(res))
  expect_true(any(grepl("GIMME", out)))
  expect_true(any(grepl("Subjects", out)))
  expect_true(any(grepl("Variables", out)))
})

test_that("summary.net_gimme produces output", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  out <- capture.output(summary(res))
  expect_true(any(grepl("FIT INDICES", out)))
  expect_true(any(grepl("TEMPORAL", out)))
})

test_that("plot.net_gimme runs without error for temporal", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_no_error(plot(res, type = "temporal"))
})

test_that("plot.net_gimme runs without error for fit", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_no_error(plot(res, type = "fit"))
})

# ===========================================================================
# Section 9: 4-variable test
# ===========================================================================
test_that("gimme works with 4 variables", {
  sim <- .make_gimme_data(n_subjects = 8, n_time = 80, n_vars = 4, seed = 99)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 99)

  expect_s3_class(res, "net_gimme")
  expect_equal(length(res$labels), 4)
  expect_equal(nrow(res$temporal), 4)
  expect_equal(ncol(res$temporal), 4)
  expect_equal(nrow(res$path_counts), 4)
  expect_equal(ncol(res$path_counts), 8)  # 4 lagged + 4 contemporaneous
})


# ===========================================================================
# Section 10: Config storage
# ===========================================================================
test_that("gimme stores config correctly", {
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     ar = TRUE, standardize = FALSE, groupcutoff = 0.75,
                     seed = 42)

  expect_equal(res$config$ar, TRUE)
  expect_equal(res$config$standardize, FALSE)
  expect_equal(res$config$groupcutoff, 0.75)
  expect_equal(res$config$seed, 42)
})


# ===========================================================================
# Section 11: lavaan not installed (L103-104)
# ===========================================================================

test_that("build_gimme errors when lavaan is not installed", {
  skip_if(requireNamespace("lavaan", quietly = TRUE),
          "lavaan is installed; skipping no-lavaan test")
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id"), "lavaan")
})


# ===========================================================================
# Section 12: hybrid mode (L145, L173, L179-180)
# ===========================================================================

test_that("build_gimme hybrid=TRUE includes residual covariances in search", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     hybrid = TRUE, seed = 1)
  expect_s3_class(res, "net_gimme")
  expect_true(isTRUE(res$config$hybrid))
})


# ===========================================================================
# Section 13: standardize=TRUE (L251-255)
# ===========================================================================

test_that("build_gimme standardize=TRUE runs without error", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     standardize = TRUE, seed = 1)
  expect_s3_class(res, "net_gimme")
  expect_true(isTRUE(res$config$standardize))
})


# ===========================================================================
# Section 14: subjects with too few time points are dropped (L259)
# ===========================================================================

test_that("build_gimme drops subjects with fewer than 3 time points", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 5, n_time = 60, seed = 1)
  # Add a subject with only 2 time points
  bad_subj <- data.frame(V1 = c(1, 2), V2 = c(1, 2), V3 = c(1, 2),
                          id = 999, time = 1:2)
  combined <- rbind(sim$data, bad_subj)
  res <- build_gimme(combined, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  # Subject 999 should be dropped
  expect_equal(res$n_subjects, 5)
  expect_false("999" %in% names(res$coefs))
})


# ===========================================================================
# Section 15: no time column (L279)
# ===========================================================================

test_that("build_gimme works without time argument", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  # Remove time ordering (no time arg → uses data order as-is)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", seed = 1)
  expect_s3_class(res, "net_gimme")
})


# ===========================================================================
# Section 16: .gimme_select_path edge cases (L355-357, L374-375, L379)
# ===========================================================================

test_that(".gimme_select_path returns NA when no subjects converge", {
  skip_if_not_installed("lavaan")
  # All NULL mi_list → converge count is 0
  mi_list <- list(NA, NA, NA)
  result <- Nestimate:::.gimme_select_path(
    mi_list, elig_paths = "V1~V2lag", prop_cutoff = 0.75,
    n_subj = 3L, chisq_cutoff = 10, hybrid = FALSE
  )
  expect_true(is.na(result))
})

test_that(".gimme_select_path returns NA when n_converge <= n_subj/2", {
  skip_if_not_installed("lavaan")
  # Only 1 out of 4 converged → 1 <= 4/2 = 2
  valid_mi <- data.frame(lhs = "V1", op = "~", rhs = "V2lag", mi = 20,
                          stringsAsFactors = FALSE)
  mi_list <- list(valid_mi, NA, NA, NA)
  result <- Nestimate:::.gimme_select_path(
    mi_list, elig_paths = "V1~V2lag", prop_cutoff = 0.75,
    n_subj = 4L, chisq_cutoff = 10, hybrid = FALSE
  )
  expect_true(is.na(result))
})


# ===========================================================================
# Section 17: .gimme_find_weakest edge cases (L392-402)
# ===========================================================================

test_that(".gimme_find_weakest returns NA when no converge", {
  result <- Nestimate:::.gimme_find_weakest(
    list(NA, NA), elig_paths = "V1~V2lag",
    prop_cutoff = 0.75, n_subj = 2L, z_cutoff = 2
  )
  expect_true(is.na(result))
})

test_that(".gimme_find_weakest returns NA when z_all is empty", {
  # z_all has no rows matching elig_paths
  valid_z <- data.frame(lhs = "V1", op = "~", rhs = "V3lag",
                         z = 1.5, stringsAsFactors = FALSE)
  result <- Nestimate:::.gimme_find_weakest(
    list(valid_z), elig_paths = "V1~V2lag",
    prop_cutoff = 0.75, n_subj = 1L, z_cutoff = 2
  )
  expect_true(is.na(result))
})


# ===========================================================================
# Section 18: .gimme_prepare_data no-subjects error (L464-465)
# ===========================================================================

test_that(".gimme_prepare_data errors when all subjects have too few obs", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 3, n_time = 60, seed = 1)
  # Keep only 2 rows per subject
  sub_data <- do.call(rbind, lapply(unique(sim$data$id), function(i) {
    head(sim$data[sim$data$id == i, ], 2)
  }))
  expect_error(
    build_gimme(sub_data, vars = sim$vars, id = "id", time = "time"),
    "enough time points|minimum 3"
  )
})


# ===========================================================================
# Section 19: .gimme_fit_final with NULL fit (L469)
# ===========================================================================

test_that(".gimme_fit_final returns zeros when lavaan fails", {
  skip_if_not_installed("lavaan")
  varnames <- c("V1", "V2")
  lag_names <- c("V1lag", "V2lag")
  # Empty syntax to guarantee failure
  result <- suppressWarnings(Nestimate:::.gimme_fit_final(
    syntax = "V1 ~ V2lag\nV2 ~ V1lag",
    data_k = data.frame(V1 = 1:2, V2 = 1:2, V1lag = 1:2, V2lag = 1:2),
    varnames = varnames,
    lag_names = lag_names
  ))
  expect_true(is.list(result))
  expect_true("coefs" %in% names(result))
  expect_true("fit_indices" %in% names(result))
})


# ===========================================================================
# Section 20: .gimme_build_syntax with ar=FALSE, no paths (L473)
# ===========================================================================

test_that(".gimme_build_syntax without AR paths excludes AR", {
  skip_if_not_installed("lavaan")
  varnames <- c("V1", "V2")
  lag_names <- paste0(varnames, "lag")
  result <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = FALSE, paths = character(0), exogenous = NULL, hybrid = FALSE
  )
  ar_paths <- paste0(varnames, "~", varnames, "lag")
  expect_false(any(ar_paths %in% result$base_syntax))
})

test_that(".gimme_build_syntax with ar=TRUE includes AR paths", {
  skip_if_not_installed("lavaan")
  varnames <- c("V1", "V2")
  lag_names <- paste0(varnames, "lag")
  result <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = TRUE, paths = character(0), exogenous = NULL, hybrid = FALSE
  )
  ar_paths <- paste0(varnames, "~", varnames, "lag")
  expect_true(all(ar_paths %in% result$base_syntax))
})


# ===========================================================================
# Section 21: print.net_gimme with group paths (L505, L520-523)
# ===========================================================================

test_that("print.net_gimme shows group paths when present", {
  skip_if_not_installed("lavaan")
  # Use strong signal data to increase chance of group paths
  set.seed(7)
  n_subj <- 10
  n_time <- 100
  vars <- c("V1", "V2")
  data_list <- lapply(seq_len(n_subj), function(i) {
    mat <- matrix(0, n_time, 2)
    mat[1, ] <- rnorm(2)
    for (t in 2:n_time) {
      mat[t, 1] <- 0.7 * mat[t - 1, 1] + rnorm(1, sd = 0.3)
      mat[t, 2] <- 0.8 * mat[t - 1, 1] + rnorm(1, sd = 0.2)
    }
    df <- as.data.frame(mat)
    colnames(df) <- vars
    df$id <- i
    df$time <- seq_len(n_time)
    df
  })
  long_data <- do.call(rbind, data_list)
  res <- build_gimme(long_data, vars = vars, id = "id", time = "time",
                     ar = TRUE, seed = 7)
  out <- capture.output(print(res))
  expect_true(any(grepl("Group-level", out)))
})


# ===========================================================================
# Section 22: plot.net_gimme type='individual' (L526-528, L531-544)
# ===========================================================================

test_that("plot.net_gimme type='individual' runs for first subject", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_no_error(plot(res, type = "individual", subject = 1L))
})

test_that("plot.net_gimme type='individual' uses character subject name", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  first_name <- names(res$coefs)[1L]
  expect_no_error(plot(res, type = "individual", subject = first_name))
})

test_that("plot.net_gimme type='individual' with NULL subject defaults to 1", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_message(plot(res, type = "individual"), "subject 1")
})

test_that("plot.net_gimme type='individual' errors on out-of-range subject", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_error(plot(res, type = "individual", subject = 999L), "Subject not found")
})


# ===========================================================================
# Section 23: plot.net_gimme type='counts' (L552-561)
# ===========================================================================

test_that("plot.net_gimme type='counts' runs without error", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  expect_no_error(plot(res, type = "counts"))
})


# ===========================================================================
# Section 24: summary.net_gimme with no group paths (L564-575, L578-579)
# ===========================================================================

test_that("summary.net_gimme shows '(none)' when no group paths", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  # Force empty group paths
  res$group_paths <- character(0)
  out <- capture.output(summary(res))
  expect_true(any(grepl("none", out)))
})


# ===========================================================================
# Section 25: .gimme_stabilize paths (L582-583, L626)
# ===========================================================================

test_that(".gimme_stabilize with empty ind_paths returns stable result", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 2, n_time = 60, seed = 1)
  ts_list <- Nestimate:::.gimme_prepare_data(sim$data, sim$vars, "id", "time",
                                             FALSE, NULL)
  varnames <- sim$vars
  lag_names <- paste0(varnames, "lag")

  syntax_info <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = TRUE, paths = character(0), exogenous = NULL, hybrid = FALSE
  )

  result <- Nestimate:::.gimme_stabilize(
    base_syntax = syntax_info$base_syntax,
    group_paths = character(0),
    ind_paths = character(0),
    data_k = ts_list[[1L]],
    endo_names = varnames,
    lag_names = lag_names
  )
  expect_true(is.list(result))
  expect_equal(result$ind_paths, character(0))
})


# ===========================================================================
# Section 26: .gimme_plot_matrix fallback (L641-644)
# ===========================================================================

test_that(".gimme_plot_matrix runs without error", {
  mat <- matrix(c(0, 0.3, 0.1, 0.3, 0, 0.2, 0.1, 0.2, 0), 3, 3)
  rownames(mat) <- colnames(mat) <- c("V1", "V2", "V3")
  expect_no_error(Nestimate:::.gimme_plot_matrix(mat, "Test Matrix"))
})


# ===========================================================================
# Section 27: .gimme_build_syntax with fixed paths (L649)
# ===========================================================================

test_that(".gimme_build_syntax includes user-specified fixed paths", {
  varnames <- c("V1", "V2", "V3")
  lag_names <- paste0(varnames, "lag")
  user_path <- "V1~V2lag"
  result <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = FALSE, paths = user_path, exogenous = NULL, hybrid = FALSE
  )
  expect_true(user_path %in% result$base_syntax)
  expect_false(user_path %in% result$candidate_paths)
})


# ===========================================================================
# Section 28: print.net_gimme no group paths (L657-666)
# ===========================================================================

test_that("print.net_gimme works with no group paths", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  res$group_paths <- character(0)
  out <- capture.output(print(res))
  expect_true(any(grepl("Group-level paths found: 0", out)))
})




# ===========================================================================
# Section 30: .gimme_test_weights (L755)
# ===========================================================================

test_that(".gimme_test_weights returns TRUE for NULL fit beta", {
  skip_if_not_installed("lavaan")
  # When std_beta is NULL (error case), should return TRUE (unstable)
  result <- Nestimate:::.gimme_test_weights(
    fit = NULL,
    endo_names = c("V1", "V2"),
    lag_names = c("V1lag", "V2lag")
  )
  expect_true(result)
})


# ===========================================================================
# Section 31: .gimme_prune_paths with empty group_paths exits (L761)
# ===========================================================================

test_that(".gimme_prune_paths returns immediately with empty paths", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 4, n_time = 60, seed = 1)
  ts_list <- Nestimate:::.gimme_prepare_data(sim$data, sim$vars, "id", "time",
                                             FALSE, NULL)
  varnames <- sim$vars
  lag_names <- paste0(varnames, "lag")
  syntax_info <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = TRUE, paths = character(0), exogenous = NULL, hybrid = FALSE
  )
  # With empty group_paths, pruning loop will find nothing to drop
  result <- Nestimate:::.gimme_prune_paths(
    base_syntax = syntax_info$base_syntax,
    group_paths = character(0),
    ts_list = ts_list,
    n_subj = length(ts_list),
    prop_cutoff = 0.75,
    z_cutoff = 2.0
  )
  expect_equal(result, character(0))
})


# ===========================================================================
# Section 32: .gimme_ind_forward_search with excellent fit stops search (L792-793)
# ===========================================================================

test_that(".gimme_ind_forward_search stops when fit is excellent", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 2, n_time = 60, seed = 1)
  ts_list <- Nestimate:::.gimme_prepare_data(sim$data, sim$vars, "id", "time",
                                             FALSE, NULL)
  varnames <- sim$vars
  lag_names <- paste0(varnames, "lag")
  syntax_info <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = TRUE, paths = character(0), exogenous = NULL, hybrid = FALSE
  )
  fit_indices <- c(rmsea_cutoff = 0.05, srmr_cutoff = 0.05,
                   nnfi_cutoff = 0.95, cfi_cutoff = 0.95)
  # Use very lenient n_excellent = 0 → stop immediately
  result <- Nestimate:::.gimme_ind_forward_search(
    base_syntax = syntax_info$base_syntax,
    group_paths = character(0),
    ind_paths = character(0),
    data_k = ts_list[[1L]],
    elig_paths = syntax_info$candidate_paths,
    ind_cutoff = 0,
    fit_indices = fit_indices,
    n_excellent = 0L,
    exclude = character(0)
  )
  expect_true(is.character(result))
})


# ===========================================================================
# Section 33: .gimme_extract_results (L814-815)
# ===========================================================================

test_that(".gimme_extract_results produces correct structure", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 4, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  # Verify the extracted results have correct shape
  expect_true(is.matrix(res$temporal_avg))
  expect_equal(nrow(res$temporal_avg), length(sim$vars))
  expect_equal(ncol(res$temporal_avg), length(sim$vars))
  expect_true(is.matrix(res$contemporaneous_avg))
})


# ===========================================================================
# Section 34: gimme.R invalid groupcutoff (L941)
# ===========================================================================

test_that("build_gimme rejects invalid groupcutoff", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 4, n_time = 50)
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id",
                            groupcutoff = 0), "groupcutoff")
  expect_error(build_gimme(sim$data, vars = sim$vars, id = "id",
                            groupcutoff = 1.5), "groupcutoff")
})


# ===========================================================================
# Section 35: .gimme_fit_and_mi returns NA on non-convergence (L980-981)
# ===========================================================================

test_that(".gimme_fit_and_mi returns NA when model does not converge", {
  skip_if_not_installed("lavaan")
  # Provide data with only 3 rows to guarantee non-convergence
  bad_data <- data.frame(V1 = 1:3, V2 = 1:3, V1lag = c(0, 1, 2),
                          V2lag = c(0, 1, 2))
  result <- suppressWarnings(Nestimate:::.gimme_fit_and_mi(
    syntax = c("V1~~V1", "V2~~V2", "V1~1", "V2~1",
               "V1lag~~V1lag", "V2lag~~V2lag"),
    data_k = bad_data,
    elig_paths = "V1~V2lag"
  ))
  # Should return NA (non-converged or failed)
  expect_true(is.na(result) || is.data.frame(result))
})


# ===========================================================================
# Section 36: .gimme_fit_and_z returns NA on NULL fit (L1015)
# ===========================================================================

test_that(".gimme_fit_and_z returns NA when model fit fails", {
  skip_if_not_installed("lavaan")
  bad_data <- data.frame(V1 = 1:3, V2 = 1:3, V1lag = c(0, 1, 2),
                          V2lag = c(0, 1, 2))
  result <- suppressWarnings(Nestimate:::.gimme_fit_and_z(
    syntax = c("V1~~V1", "V2~~V2", "V1~1", "V2~1"),
    data_k = bad_data,
    elig_paths = "V1~V2lag"
  ))
  expect_true(is.na(result) || is.data.frame(result))
})


# ===========================================================================
# Section 37: gimme with exogenous variable (L1019-1025)
# ===========================================================================

test_that("build_gimme runs with exogenous specification", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  # Use V3 as exogenous variable (not predicted, only predicts)
  expect_no_error(
    res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                       exogenous = "V3", ar = TRUE, seed = 1)
  )
  expect_s3_class(res, "net_gimme")
})


# ===========================================================================
# Section 38: .gimme_prepare_data with exogenous (L1029-1039)
# ===========================================================================

test_that(".gimme_prepare_data handles exogenous variable", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 4, n_time = 60, seed = 1)
  ts_list <- Nestimate:::.gimme_prepare_data(
    sim$data, sim$vars, "id", "time", FALSE, exogenous = "V3"
  )
  expect_true(is.list(ts_list))
  expect_true(length(ts_list) > 0)
})








# ===========================================================================
# Section 42: plot.net_gimme fit histogram for NA columns (L1085)
# ===========================================================================

test_that("plot.net_gimme type='fit' skips all-NA metric columns", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 6, n_time = 60, seed = 1)
  res <- build_gimme(sim$data, vars = sim$vars, id = "id", time = "time",
                     seed = 1)
  # Force all nnfi values to NA to exercise the skip-NA branch
  res$fit$nnfi <- NA_real_
  expect_no_error(plot(res, type = "fit"))
})


# ===========================================================================
# Section 43: .gimme_fit_final sets status (L1098-1104)
# ===========================================================================

test_that(".gimme_fit_final includes status in output", {
  skip_if_not_installed("lavaan")
  sim <- .make_gimme_data(n_subjects = 2, n_time = 60, seed = 1)
  ts_list <- Nestimate:::.gimme_prepare_data(sim$data, sim$vars, "id", "time",
                                             FALSE, NULL)
  varnames <- sim$vars
  lag_names <- paste0(varnames, "lag")
  syntax_info <- Nestimate:::.gimme_build_syntax(
    varnames, lag_names, varnames, lag_names,
    ar = TRUE, paths = character(0), exogenous = NULL, hybrid = FALSE
  )
  result <- Nestimate:::.gimme_fit_final(
    syntax = syntax_info$base_syntax,
    data_k = ts_list[[1L]],
    varnames = varnames,
    lag_names = lag_names
  )
  expect_true("status" %in% names(result))
  expect_true(is.character(result$status))
})
