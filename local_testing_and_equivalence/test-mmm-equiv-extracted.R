# Equivalence test_that() blocks extracted from
# tests/testthat/test-mmm.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# Tests for mmm.R: Mixed Markov Model

# ============================================
# Synthetic data helper
# ============================================

.make_mmm_data <- function(n_per_group = 30, n_cols = 10, seed = 42) {
  set.seed(seed)
  g1 <- data.frame(matrix(sample(c("A", "B", "C"), n_per_group * n_cols,
                                  replace = TRUE, prob = c(0.7, 0.2, 0.1)),
                           ncol = n_cols))
  g2 <- data.frame(matrix(sample(c("A", "B", "C"), n_per_group * n_cols,
                                  replace = TRUE, prob = c(0.1, 0.2, 0.7)),
                           ncol = n_cols))
  rbind(g1, g2)
}

# ============================================
# build_mmm basic
# ============================================

test_that("build_mmm works with tna input", {
  skip_if_not_installed("tna")
  model <- tna::tna(tna::group_regulation)
  mmm <- build_mmm(model, k = 2, n_starts = 3, seed = 1)

  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 2)
  expect_true(length(mmm$states) > 2)
})

# ============================================
# compare_mmm
# ============================================

