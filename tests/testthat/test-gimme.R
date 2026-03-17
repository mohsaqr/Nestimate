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
# Section 8: gimme package equivalence
# ===========================================================================
test_that("gimme equivalence: path counts match on simulateVAR data", {
  skip_if_not_installed("gimme")
  library(gimme)

  set.seed(42)
  sim <- gimme::simulateVAR(
    A = matrix(c(.4, .2, 0, .1, .3, .15, 0, .1, .35), 3, 3, byrow = TRUE),
    Phi = matrix(c(.15, .05, 0, 0, .1, 0, 0, 0, .1), 3, 3, byrow = TRUE),
    Psi = diag(.4, 3),
    subAssign = rep(1, 12), N = 12, Obs = 120
  )

  # Run original gimme
  tmpdir <- tempdir()
  datadir <- file.path(tmpdir, "gimme_equiv")
  outdir <- file.path(tmpdir, "gimme_equiv_out")
  dir.create(datadir, showWarnings = FALSE)
  dir.create(outdir, showWarnings = FALSE)
  for (i in seq_along(sim$dataList)) {
    df <- as.data.frame(sim$dataList[[i]])
    colnames(df) <- paste0("V", 1:3)
    write.csv(df, file.path(datadir, sprintf("sub_%02d.csv", i)),
              row.names = FALSE)
  }
  gimme_res <- gimme::gimme(data = datadir, out = outdir, sep = ",",
                             header = TRUE, ar = TRUE, plot = FALSE,
                             subgroup = FALSE)

  # Run our version
  data_list <- lapply(seq_along(sim$dataList), function(i) {
    df <- as.data.frame(sim$dataList[[i]])
    colnames(df) <- paste0("V", 1:3)
    df$id <- i
    df$time <- seq_len(nrow(df))
    df
  })
  long_data <- do.call(rbind, data_list)
  our_res <- build_gimme(long_data, vars = c("V1", "V2", "V3"), id = "id",
                         time = "time", ar = TRUE, seed = 42)

  # Path counts should be identical
  expect_equal(our_res$path_counts, gimme_res$path_counts,
               info = "Path count matrices should match exactly")
})

test_that("gimme equivalence: coefficients are close", {
  skip_if_not_installed("gimme")
  library(gimme)

  set.seed(42)
  sim <- gimme::simulateVAR(
    A = matrix(c(.4, .2, 0, .1, .3, .15, 0, .1, .35), 3, 3, byrow = TRUE),
    Phi = matrix(c(.15, .05, 0, 0, .1, 0, 0, 0, .1), 3, 3, byrow = TRUE),
    Psi = diag(.4, 3),
    subAssign = rep(1, 12), N = 12, Obs = 120
  )

  tmpdir <- tempdir()
  datadir <- file.path(tmpdir, "gimme_equiv2")
  outdir <- file.path(tmpdir, "gimme_equiv_out2")
  dir.create(datadir, showWarnings = FALSE)
  dir.create(outdir, showWarnings = FALSE)
  for (i in seq_along(sim$dataList)) {
    df <- as.data.frame(sim$dataList[[i]])
    colnames(df) <- paste0("V", 1:3)
    write.csv(df, file.path(datadir, sprintf("sub_%02d.csv", i)),
              row.names = FALSE)
  }
  gimme_res <- gimme::gimme(data = datadir, out = outdir, sep = ",",
                             header = TRUE, ar = TRUE, plot = FALSE,
                             subgroup = FALSE)

  data_list <- lapply(seq_along(sim$dataList), function(i) {
    df <- as.data.frame(sim$dataList[[i]])
    colnames(df) <- paste0("V", 1:3)
    df$id <- i
    df$time <- seq_len(nrow(df))
    df
  })
  long_data <- do.call(rbind, data_list)
  our_res <- build_gimme(long_data, vars = c("V1", "V2", "V3"), id = "id",
                         time = "time", ar = TRUE, seed = 42)

  # Coefficients should be very close (same paths, same data, same lavaan)
  max_diffs <- vapply(seq_along(our_res$coefs), function(k) {
    max(abs(our_res$coefs[[k]] - gimme_res$path_est_mats[[k]]))
  }, numeric(1))

  # Coefficients should be identical (same standardized betas, same rounding)
  expect_true(all(max_diffs == 0),
              info = sprintf("Max coefficient diffs: %s",
                             paste(round(max_diffs, 6), collapse = ", ")))
})

test_that("gimme equivalence: data preparation matches", {
  skip_if_not_installed("gimme")
  library(gimme)

  set.seed(42)
  sim <- gimme::simulateVAR(
    A = matrix(c(.3, .1, 0, 0, .3, 0, 0, 0, .3), 3, 3, byrow = TRUE),
    Phi = matrix(c(.1, 0, 0, 0, .1, 0, 0, 0, .1), 3, 3, byrow = TRUE),
    Psi = diag(.5, 3),
    subAssign = rep(1, 8), N = 8, Obs = 100
  )

  tmpdir <- tempdir()
  datadir <- file.path(tmpdir, "gimme_equiv3")
  outdir <- file.path(tmpdir, "gimme_equiv_out3")
  dir.create(datadir, showWarnings = FALSE)
  dir.create(outdir, showWarnings = FALSE)
  for (i in seq_along(sim$dataList)) {
    df <- as.data.frame(sim$dataList[[i]])
    colnames(df) <- paste0("V", 1:3)
    write.csv(df, file.path(datadir, sprintf("sub_%02d.csv", i)),
              row.names = FALSE)
  }
  gimme_res <- gimme::gimme(data = datadir, out = outdir, sep = ",",
                             header = TRUE, ar = TRUE, plot = FALSE,
                             subgroup = FALSE)

  data_list <- lapply(seq_along(sim$dataList), function(i) {
    df <- as.data.frame(sim$dataList[[i]])
    colnames(df) <- paste0("V", 1:3)
    df$id <- i
    df$time <- seq_len(nrow(df))
    df
  })
  long_data <- do.call(rbind, data_list)
  ts <- Nestimate:::.gimme_prepare_data(long_data, c("V1", "V2", "V3"),
                                       "id", "time", FALSE, NULL)

  # Data should be identical
  for (k in seq_along(ts)) {
    g <- gimme_res$data[[k]]
    o <- ts[[k]][, c("V1lag", "V2lag", "V3lag", "V1", "V2", "V3")]
    om <- unname(as.matrix(o))
    gm <- unname(as.matrix(g))
    expect_equal(om, gm, tolerance = 1e-12,
                 info = sprintf("Subject %d data mismatch", k))
  }
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
