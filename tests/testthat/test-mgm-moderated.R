# Smoke tests for moderated MGM implementation.
#
# Tests structure, conditioning at both moderator levels, and symmetry.
# Full equivalence validation is in test-equiv-mgm-mod.R.

skip_if_not_installed("mgm")
skip_if_not_installed("glmnet")

TOL <- 1e-10

test_that(".mgm_estimate_moderated returns expected structure", {
  set.seed(7)
  n <- 200
  dat <- data.frame(
    V1 = stats::rnorm(n),
    V2 = stats::rnorm(n),
    V3 = stats::rnorm(n),
    V4 = sample(0:1, n, TRUE)
  )
  fit <- Nestimate:::.mgm_estimate_moderated(
    data = dat, type = c("g", "g", "g", "c"), level = c(1, 1, 1, 2),
    moderator = 4L
  )
  expect_s3_class(fit, "mgm_moderated")
  expect_equal(fit$moderator, 4L)
  expect_equal(fit$p, 4L)
  expect_length(fit$fits, 4L)
  expect_true(is.matrix(fit$pw_alive))
  expect_true(is.matrix(fit$int_alive))

  cond <- Nestimate:::condition_moderated(fit, mod_value = 0)
  expect_equal(dim(cond), c(4L, 4L))
  expect_true(isSymmetric(cond))
  expect_true(all(cond[fit$moderator, ] == 0))
  expect_true(all(cond[, fit$moderator] == 0))
})


test_that("condition_moderated produces symmetric wadj", {
  set.seed(99)
  n <- 250
  dat <- data.frame(
    V1 = stats::rnorm(n),
    V2 = stats::rnorm(n),
    V3 = stats::rnorm(n),
    V4 = stats::rnorm(n),
    V5 = sample(0:1, n, TRUE)
  )
  dat$V4 <- dat$V4 + ifelse(dat$V5 == 1, 0.5 * dat$V1, 0)
  fit <- Nestimate:::.mgm_estimate_moderated(
    dat, type = c("g", "g", "g", "g", "c"), level = c(1, 1, 1, 1, 2),
    moderator = 5L
  )
  for (val in c(0, 1)) {
    w <- Nestimate:::condition_moderated(fit, mod_value = val)
    expect_true(isSymmetric(w), info = sprintf("mod_value = %d", val))
    expect_true(all(w[5, ] == 0), info = sprintf("mod_value = %d row", val))
    expect_true(all(w[, 5] == 0), info = sprintf("mod_value = %d col", val))
  }
})
