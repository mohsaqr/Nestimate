# Equivalence test_that() blocks extracted from
# tests/testthat/test-markov.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- Tests for passage_time() and markov_stability() ----

# Three-state ergodic chain from Kemeny & Snell example
.P3 <- matrix(
  c(0.7, 0.2, 0.1,
    0.3, 0.5, 0.2,
    0.2, 0.3, 0.5),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
)

# ---- .mpt_stationary ----

test_that("passage_time accepts a tna object", {
  skip_if_pkg_broken("tna")
  m  <- tna::tna(tna::group_regulation)
  pt <- passage_time(m)
  expect_s3_class(pt, "net_mpt")
  expect_equal(ncol(pt$matrix), nrow(pt$matrix))
  expect_true(all(pt$matrix > 0))
})

test_that("markov_stability accepts a tna object", {
  skip_if_pkg_broken("tna")
  m  <- tna::tna(tna::group_regulation)
  ms <- markov_stability(m)
  expect_s3_class(ms, "net_markov_stability")
  expect_equal(nrow(ms$stability), ncol(m$weights))
})

test_that("passage_time matches netobject when built from same data as tna", {
  skip_if_pkg_broken("tna")
  m   <- tna::tna(tna::group_regulation)
  pt_tna <- passage_time(m)
  net    <- build_network(tna::group_regulation, method = "relative")
  pt_net <- passage_time(net)
  # Both extract the same weight matrix so results must be identical
  expect_equal(pt_tna$matrix, pt_net$matrix, tolerance = 1e-8)
})

# ---- wide data.frame support ----

