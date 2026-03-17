# ---- Tests for internal utility functions ----

check_val_in_range <- Nestimate:::check_val_in_range
safe_median <- Nestimate:::safe_median
safe_mean <- Nestimate:::safe_mean
safe_sd <- Nestimate:::safe_sd

test_that("check_val_in_range returns TRUE for values inside range", {
  expect_true(check_val_in_range(5, c(1, 10)))
  expect_true(check_val_in_range(0.5, c(0, 1)))
  expect_true(check_val_in_range(-3, c(-5, 0)))
})

test_that("check_val_in_range returns TRUE at boundaries", {
  expect_true(check_val_in_range(1, c(1, 10)))
  expect_true(check_val_in_range(10, c(1, 10)))
  expect_true(check_val_in_range(0, c(0, 0)))
})

test_that("check_val_in_range returns FALSE for values outside range", {
  expect_false(check_val_in_range(0, c(1, 10)))
  expect_false(check_val_in_range(11, c(1, 10)))
  expect_false(check_val_in_range(-1, c(0, 1)))
})

test_that("check_val_in_range returns TRUE when range_val is NULL", {
  expect_true(check_val_in_range(42, NULL))
  expect_true(check_val_in_range(-999, NULL))
})

test_that("check_val_in_range returns FALSE for NA or non-numeric", {
  expect_false(check_val_in_range(NA, c(1, 10)))
  expect_false(check_val_in_range("a", c(1, 10)))
  expect_false(check_val_in_range(NA_real_, c(0, 1)))
})

test_that("safe_median returns correct median", {
  expect_equal(safe_median(c(1, 2, 3)), 2)
  expect_equal(safe_median(c(1, 2, 3, 4)), 2.5)
  expect_equal(safe_median(c(1, NA, 3)), 2)
})

test_that("safe_median returns NA_real_ for empty vector", {
  expect_identical(safe_median(numeric(0)), NA_real_)
})

test_that("safe_mean returns correct mean", {
  expect_equal(safe_mean(c(2, 4, 6)), 4)
  expect_equal(safe_mean(c(1, NA, 3)), 2)
})

test_that("safe_mean returns NA_real_ for empty vector", {
  expect_identical(safe_mean(numeric(0)), NA_real_)
})

test_that("safe_sd returns correct sd", {
  expect_equal(safe_sd(c(1, 2, 3)), sd(c(1, 2, 3)))
  expect_equal(safe_sd(c(10, NA, 30)), sd(c(10, 30), na.rm = TRUE))
})

test_that("safe_sd returns NA_real_ for single or empty vector", {
  expect_identical(safe_sd(numeric(0)), NA_real_)
  expect_identical(safe_sd(5), NA_real_)
  expect_identical(safe_sd(c(NA, NA)), NA_real_)
})
