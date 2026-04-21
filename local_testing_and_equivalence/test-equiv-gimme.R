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


