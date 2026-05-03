# Equivalence test_that() blocks extracted from
# tests/testthat/test-simplicial.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ---- Simplicial Complex Tests ----

# Helper: build a known network
.make_sc_net <- function(seed = 42) {
  set.seed(seed)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:5], 40, TRUE),
    V2 = sample(LETTERS[1:5], 40, TRUE),
    V3 = sample(LETTERS[1:5], 40, TRUE),
    V4 = sample(LETTERS[1:5], 40, TRUE)
  )
  build_network(seqs, method = "relative")
}

# Helper: build a simple known matrix for exact verification
.make_sc_mat <- function() {
  # 4-node complete graph => K4 clique complex
  mat <- matrix(0.5, 4, 4)
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- LETTERS[1:4]
  mat
}

# Helper: triangle matrix
.make_triangle_mat <- function() {
  mat <- matrix(0, 3, 3)
  mat[1, 2] <- mat[2, 1] <- 1
  mat[1, 3] <- mat[3, 1] <- 1
  mat[2, 3] <- mat[3, 2] <- 1
  rownames(mat) <- colnames(mat) <- c("X", "Y", "Z")
  mat
}


# =========================================================================
# build_simplicial — clique type
# =========================================================================

test_that("build_simplicial accepts tna object input", {
  skip_if_not_installed("tna")
  tna_model <- tna::tna(tna::group_regulation)
  sc <- build_simplicial(tna_model, threshold = 0.1)
  expect_s3_class(sc, "simplicial_complex")
})

test_that(".sc_extract_matrix handles different input types", {
  # Matrix
  mat <- .make_sc_mat()
  expect_identical(.sc_extract_matrix(mat), mat)

  # netobject
  net <- .make_sc_net()
  expect_identical(.sc_extract_matrix(net), net$weights)

  # tna object
  skip_if_not_installed("tna")
  tna_model <- tna::tna(tna::group_regulation)
  expect_identical(.sc_extract_matrix(tna_model), tna_model$weights)

  # net_hon object
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  hon <- build_hon(trajs, max_order = 2, min_freq = 1)
  expect_identical(.sc_extract_matrix(hon), hon$matrix)

  # Unsupported
  expect_error(.sc_extract_matrix(data.frame(a = 1)), "Cannot extract")
})


# =========================================================================
# betti_numbers
# =========================================================================

