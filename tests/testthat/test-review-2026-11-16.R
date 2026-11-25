testthat::skip_on_cran()

test_that("Pearson mosaic tests reject invalid counts before pooling", {
  for (bad in list(-1, 0.5, NA_real_, Inf, "x")) {
    a <- array(2, c(2, 2, 3))
    a[1] <- bad
    for (method in list("chisq", c("chisq", "rao_scott"))) {
      expect_error(mosaic_test(a, method = method),
                   class = "nestimate_bad_input")
    }
  }
})

test_that("unlabelled sequences cannot enter sequence comparison permutations", {
  x <- matrix(c("A", "B", "A", "A", "B", "B", "B", "A"),
              ncol = 2, byrow = TRUE)
  g <- c("x", "x", "y", "y")
  set.seed(14)
  ref <- sequence_compare(x, group = g, sub = 1, min_freq = 1, iter = 39)
  # An unlabelled sequence has arbitrarily many occurrences of a new state.
  # Its presence must change neither the observed counts nor the null draws.
  set.seed(14)
  expect_warning(actual <- sequence_compare(rbind(x, c("Z", "Z")),
    group = c(g, NA), sub = 1, min_freq = 1, iter = 39),
    class = "nestimate_missing_group")
  expect_equal(actual, ref)
})

test_that("sequence comparison validates integer controls before coercion", {
  x <- matrix(rep(c("A", "B"), 8), nrow = 4)
  base <- list(x = x, group = c("x", "x", "y", "y"), sub = 1,
               min_freq = 1, iter = 2)
  for (arg in c("sub", "min_freq", "iter")) {
    for (bad in list(1.5, numeric(), Inf, NA_real_, "2", 2^32)) {
      args <- base
      args[[arg]] <- bad
      expect_error(do.call(sequence_compare, args), class = "nestimate_bad_input")
    }
  }
  for (arg in c("min_freq", "iter")) {
    args <- base
    args[[arg]] <- c(1, 2)
    expect_error(do.call(sequence_compare, args), class = "nestimate_bad_input")
  }
})

test_that("nested GLMM residuals cluster influences at the outermost level", {
  skip_if_not_installed("lme4")
  set.seed(19)
  a <- array(sample(1:12, 2 * 2 * 12, replace = TRUE), c(2, 2, 12))
  lv <- data.frame(team = rep(1:4, each = 3), person = seq_len(12))
  fit <- suppressWarnings(suppressMessages(
    Nestimate:::.nested_association(a, "glmm", levels = lv)))
  # Independently calculate the outer-team influence for cell (1, 1).
  pooled <- apply(a, 1:2, sum)
  N <- sum(a)
  p <- sum(pooled[, 1]) / N
  row_n <- sum(pooled[1, ])
  influence <- vapply(1:4, function(team) {
    z <- apply(a[, , lv$team == team, drop = FALSE], 1:2, sum)
    z[1, 1] - sum(z[1, ]) * p - row_n / N * (sum(z[, 1]) - sum(z) * p)
  }, numeric(1))
  se <- sqrt(4 / 3 * sum((influence - mean(influence))^2))
  expect_equal(unname(fit$residuals[1, 1]),
               (pooled[1, 1] - row_n * p) / se, tolerance = 1e-10)
})
