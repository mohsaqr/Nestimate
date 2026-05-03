# Equivalence test_that() blocks extracted from
# tests/testthat/test-mcml.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# Tests for mcml.R: net_aggregate_weights, cluster_summary, build_mcml

# ============================================
# net_aggregate_weights / wagg
# ============================================

test_that("cluster_summary with tna object extracts weights (L300)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.6, 0.4, 0.7, 0, 0.3, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  cs <- cluster_summary(tna_obj, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml handles tna object with data (tna_data branch, L581-598)", {
  skip_if_not_installed("tna")
  seqs <- data.frame(
    T1 = c("A", "B", "C"),
    T2 = c("B", "C", "A"),
    T3 = c("C", "A", "B")
  )
  # Build a tna object that contains $data (integer-encoded: 1=A, 2=B, 3=C)
  tna_obj <- tna::tna(seqs)
  if (!is.null(tna_obj$data)) {
    # tna encodes states as integers 1,2,3 in $data
    cs <- build_mcml(tna_obj, list(G1 = c("1", "2"), G2 = "3"))
    expect_s3_class(cs, "mcml")
  } else {
    skip("tna object has no $data field in this version")
  }
})

test_that("build_mcml handles tna object without data (tna_matrix branch, L583-585)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.6, 0.4, 0.7, 0, 0.3, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  cs <- build_mcml(tna_obj, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

# ---- build_mcml: netobject_data branch (L586-599) ----

test_that(".detect_mcml_input returns tna_data for tna with data (L620)", {
  skip_if_not_installed("tna")
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  tna_obj <- tna::tna(seqs)
  if (!is.null(tna_obj$data)) {
    result <- Nestimate:::.detect_mcml_input(tna_obj)
    expect_equal(result, "tna_data")
  } else {
    skip("tna object has no $data in this version")
  }
})

test_that(".detect_mcml_input returns tna_matrix for tna without data (L621)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.5, 0.5, 0.3, 0, 0.7, 0.6, 0.4, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  result <- Nestimate:::.detect_mcml_input(tna_obj)
  expect_equal(result, "tna_matrix")
})

test_that("as_tna.default returns tna object unchanged (L1303-1305)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.5, 0.5, 0.3, 0, 0.7, 0.6, 0.4, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  result <- as_tna(tna_obj)
  expect_true(inherits(result, "tna"))
})

