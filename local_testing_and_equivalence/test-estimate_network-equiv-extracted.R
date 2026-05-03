# Equivalence test_that() blocks extracted from
# tests/testthat/test-estimate_network.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- estimate_network() Tests ----
# estimate_network() is deprecated; these tests verify backward compat.

# Helper: generate reproducible frequency-like data
.make_assoc_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}

# Helper: generate wide sequence data
.make_wide_seq <- function(n = 50, t = 10, states = c("A", "B", "C"),
                           seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}


# ---- Deprecation ----

test_that("estimate_network works with tna::group_regulation (wide)", {
  skip_if_not_installed("tna")

  net <- suppressWarnings(
    estimate_network(tna::group_regulation, method = "relative")
  )

  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  expect_equal(net$n_nodes, 9)
  # Rows sum to ~1
  row_sums <- rowSums(net$weights)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("estimate_network frequency matches frequencies() output", {
  skip_if_not_installed("tna")

  net <- suppressWarnings(
    estimate_network(tna::group_regulation, method = "frequency")
  )
  freq_mat <- frequencies(tna::group_regulation, format = "wide")

  expect_equal(net$weights, freq_mat)
})


# ---- Coverage gap tests ----

# estimate_network.R L81: minmax scaling when all non-zero values equal
