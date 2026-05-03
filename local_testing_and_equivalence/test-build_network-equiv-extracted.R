# Equivalence test_that() blocks extracted from
# tests/testthat/test-build_network.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# Helper: generate reproducible frequency-like data
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$rid <- seq_len(n)
  df
}


# ---- Input validation ----

test_that("$data is a data frame for transition methods", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "tna")

  expect_true(is.data.frame(net$data))
  expect_equal(nrow(net$data), 2000)
  expect_equal(ncol(net$data), 26)
})

test_that("print.netobject shows metadata column names", {
  skip_if_not_installed("tna")
  # tna::group_regulation has sequences with metadata-like columns
  # Build a dataset where there are extra non-state numeric columns
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B"), 50, replace = TRUE),
    T2 = sample(c("A", "B"), 50, replace = TRUE),
    T3 = sample(c("A", "B"), 50, replace = TRUE),
    Age = rpois(50, 25),
    stringsAsFactors = FALSE
  )
  net <- build_network(df, method = "relative",
                       params = list(format = "wide"))
  # metadata should be present
  if (!is.null(net$metadata)) {
    out <- capture.output(print(net))
    expect_true(any(grepl("Metadata:", out)))
  } else {
    skip("No metadata column detected in this data")
  }
})

# L569-575: print.netobject_group
test_that("auto-convert: ising does NOT auto-convert sequences", {
  skip_if_not_installed("IsingFit")
  seqs <- .make_seq_data(n = 100, states = c("A", "B"))
  # Ising requires binary data; frequency counts are integers > 1
  # so auto-conversion is skipped and the estimator errors
  expect_error(build_network(seqs, method = "ising"))
})

# ---- 5. Aliases work: ebicglasso, corr ----
test_that("auto-convert: glasso on tna::group_regulation", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  net <- build_network(group_regulation, method = "glasso")
  expect_s3_class(net, "netobject")
  expect_equal(nrow(net$weights), 9)  # 9 states
  expect_true(any(net$weights != 0))
})

# ---- 8. pcor on tna::group_regulation ----
test_that("auto-convert: pcor on tna::group_regulation", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  net <- build_network(group_regulation, method = "pcor")
  expect_s3_class(net, "netobject")
  expect_equal(nrow(net$weights), 9)
})

# ---- 9. cor on tna::group_regulation ----
test_that("auto-convert: cor on tna::group_regulation", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  net <- build_network(group_regulation, method = "cor")
  expect_s3_class(net, "netobject")
  expect_equal(nrow(net$weights), 9)
})

# ---- 10. No conversion for transition methods (still works as before) ----
