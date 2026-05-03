# Equivalence test_that() blocks extracted from
# tests/testthat/test-utils.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- Tests for internal utility functions ----

check_val_in_range <- Nestimate:::check_val_in_range
safe_median <- Nestimate:::safe_median
safe_mean <- Nestimate:::safe_mean
safe_sd <- Nestimate:::safe_sd

test_that(".as_netobject returns netobject unchanged (L120 short-circuit)", {
  skip_if_pkg_broken("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  result <- Nestimate:::.as_netobject(net)
  expect_identical(result, net)
})

