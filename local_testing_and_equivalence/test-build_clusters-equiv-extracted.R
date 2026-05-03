# Equivalence test_that() blocks extracted from
# tests/testthat/test-build_clusters.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ==============================================================================
# Tests for build_clusters()
# ==============================================================================

# ---- Test data ----

make_test_data <- function(n = 50, k = 26, n_states = 4, seed = 42) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  mat <- matrix(sample(states, n * k, replace = TRUE), nrow = n, ncol = k)
  colnames(mat) <- paste0("T", seq_len(k))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

make_data_with_na <- function(n = 30, k = 20, n_states = 3, seed = 42) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  mat <- matrix(sample(states, n * k, replace = TRUE), nrow = n, ncol = k)
  # Add void markers to last few columns
  for (i in seq_len(n)) {
    trail_start <- sample(k - 5, 1) + 5
    if (trail_start <= k) mat[i, trail_start:k] <- "%"
  }
  # Add some mid-sequence NAs
  mat[sample(n * k, 10)] <- "*"
  colnames(mat) <- paste0("T", seq_len(k))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# ==============================================================================
# 1. Input validation
# ==============================================================================

test_that("distance matrices match tna for metrics with matching implementations", {
  skip_if_not_installed("tna")
  data <- tna::group_regulation[1:50, ]

  # Only compare metrics where tna and stringdist agree.
  # tna's own C implementations of osa/lv/dl/lcs/jw differ from
  # stringdist (which is the reference implementation). Confirmed:
  # stringdist::stringdist("acfhicbc", tna3, method="lv") matches
  # our result (21) while tna:::levenshtein_dist gives 18.
  for (metric in c("hamming", "qgram", "cosine", "jaccard")) {
    tna_r <- tna::cluster_data(data, k = 2, dissimilarity = metric, q = 2)
    our_r <- build_clusters(data, k = 2, dissimilarity = metric, q = 2L)
    expect_equal(
      as.matrix(our_r$distance), as.matrix(tna_r$distance),
      tolerance = 1e-10, info = metric
    )
  }
})

test_that("weighted hamming matches tna (lambda = 0.5)", {
  skip_if_not_installed("tna")
  data <- tna::group_regulation[1:50, ]

  tna_r <- tna::cluster_data(data, k = 2, dissimilarity = "hamming",
                             weighted = TRUE, lambda = 0.5)
  our_r <- build_clusters(data, k = 2, dissimilarity = "hamming",
                        weighted = TRUE, lambda = 0.5)
  expect_equal(
    as.matrix(our_r$distance), as.matrix(tna_r$distance),
    tolerance = 1e-10
  )
})

test_that("weighted hamming matches tna (lambda = 2.0)", {
  skip_if_not_installed("tna")
  data <- tna::group_regulation[1:50, ]

  tna_r <- tna::cluster_data(data, k = 2, dissimilarity = "hamming",
                             weighted = TRUE, lambda = 2.0)
  our_r <- build_clusters(data, k = 2, dissimilarity = "hamming",
                        weighted = TRUE, lambda = 2.0)
  expect_equal(
    as.matrix(our_r$distance), as.matrix(tna_r$distance),
    tolerance = 1e-10
  )
})

test_that("cluster assignments match tna for PAM", {
  skip_if_not_installed("tna")
  data <- tna::group_regulation[1:50, ]

  tna_r <- tna::cluster_data(data, k = 3, dissimilarity = "hamming")
  our_r <- build_clusters(data, k = 3, dissimilarity = "hamming")

  # Distance matrices must match exactly
  expect_equal(
    as.matrix(our_r$distance), as.matrix(tna_r$distance),
    tolerance = 1e-10
  )
  # PAM is deterministic on same distance matrix → same assignments
  expect_equal(our_r$assignments, tna_r$assignments)
  expect_equal(our_r$silhouette, tna_r$silhouette, tolerance = 1e-10)
})

test_that("cluster assignments match tna for hclust methods", {
  skip_if_not_installed("tna")
  data <- tna::group_regulation[1:50, ]

  for (m in c("complete", "average")) {
    tna_r <- tna::cluster_data(data, k = 3, dissimilarity = "hamming",
                               method = m)
    our_r <- build_clusters(data, k = 3, dissimilarity = "hamming",
                          method = m)
    expect_equal(
      as.matrix(our_r$distance), as.matrix(tna_r$distance),
      tolerance = 1e-10, info = m
    )
    expect_equal(our_r$assignments, tna_r$assignments, info = m)
  }
})

# ==============================================================================
# 10. R fallback vs stringdist consistency
# ==============================================================================

test_that("build_clusters works on tna model", {
  skip_if_not_installed("tna")
  data <- tna::group_regulation[1:30, ]
  model <- tna::tna(data)
  cl_tna <- build_clusters(model, k = 2)
  cl_df <- build_clusters(data, k = 2)
  expect_s3_class(cl_tna, "net_clustering")
  expect_equal(as.matrix(cl_tna$distance), as.matrix(cl_df$distance),
               tolerance = 1e-10)
})

test_that("build_clusters works on cograph_network", {
  skip_if_not_installed("cograph")
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  # Build a mock cograph_network with $data
  cg <- structure(
    list(data = df, weights = matrix(0, 4, 4), directed = TRUE),
    class = c("cograph_network", "list")
  )
  cl_cg <- build_clusters(cg, k = 3)
  cl_df <- build_clusters(df, k = 3)
  expect_s3_class(cl_cg, "net_clustering")
  expect_equal(as.matrix(cl_cg$distance), as.matrix(cl_df$distance),
               tolerance = 1e-10)
})

# ==============================================================================
# 15. build_network dispatch for net_clustering
# ==============================================================================

test_that("tna input with column-name covariates errors", {
  skip_if_not_installed("tna")
  model <- tna::tna(tna::group_regulation[1:20, ])
  expect_error(
    build_clusters(model, k = 2, covariates = "some_col"),
    "tna/cograph_network"
  )
})

