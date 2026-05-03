# Equivalence test_that() blocks extracted from
# tests/testthat/test-mgm-moderated.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# Smoke tests for moderated MGM implementation.
#
# Tests structure, conditioning at both moderator levels, and symmetry.
# Full equivalence validation is in test-equiv-mgm-mod.R.

skip_if_not_installed("mgm")
skip_if_not_installed("glmnet")

TOL <- 1e-10

test_that(".mgm_estimate_moderated matches mgm on the recon case", {
  set.seed(42)
  n <- 300
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n) + 0.4 * x1
  x3 <- stats::rnorm(n)
  x4 <- sample(0:1, n, TRUE)
  x5 <- stats::rnorm(n) + ifelse(x4 == 1, 0.7 * x1, -0.1 * x1)
  dat <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)
  type  <- c("g", "g", "g", "c", "g")
  level <- c(1, 1, 1, 2, 1)

  fit <- Nestimate:::.mgm_estimate_moderated(
    data = dat, type = type, level = level, moderator = 4L,
    lambdaGam = 0.25, ruleReg = "AND", threshold = "LW", scale = TRUE
  )
  ours_0 <- Nestimate:::condition_moderated(fit, mod_value = 0)
  ours_1 <- Nestimate:::condition_moderated(fit, mod_value = 1)

  ref <- suppressWarnings(suppressMessages(mgm::mgm(
    data       = as.matrix(dat),
    type       = type,
    level      = level,
    moderators = 4L,
    lambdaSel  = "EBIC",
    lambdaGam  = 0.25,
    ruleReg    = "AND",
    threshold  = "LW",
    overparameterize = FALSE,
    scale      = TRUE,
    pbar       = FALSE,
    signInfo   = FALSE,
    warnings   = FALSE
  )))
  ref_0 <- mgm::condition(ref, values = list("4" = 0))$pairwise$wadj
  ref_1 <- mgm::condition(ref, values = list("4" = 1))$pairwise$wadj

  expect_equal(unname(ours_0), unname(ref_0), tolerance = TOL,
               label = "condition(mod = 0)")
  expect_equal(unname(ours_1), unname(ref_1), tolerance = TOL,
               label = "condition(mod = 1)")
})


