# Equivalence test_that() blocks extracted from
# tests/testthat/test-hon.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ===========================================================================
# Tests for build_hon() — Higher-Order Network construction
# ===========================================================================

# --- Helper: simple trajectories with known structure ---
.make_hon_data <- function() {
  data.frame(
    T1 = c("A", "B", "C", "A", "D"),
    T2 = c("B", "A", "A", "B", "A"),
    T3 = c("C", "B", "B", "C", "B"),
    T4 = c("D", "C", "C", "D", "C"),
    T5 = c("A", "D", "D", "A", NA),
    T6 = c("B", "A", "A", "B", NA),
    T7 = c("C", "B", "B", "C", NA),
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# Section 1: Input validation
# ===========================================================================
test_that("hon+ matches pyHON+ on group_regulation subset", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  gr <- group_regulation[1:100, ]
  trajs <- lapply(seq_len(nrow(gr)), function(i) {
    row <- as.character(gr[i, ])
    row[!is.na(row)]
  })

  r <- build_hon(trajs, max_order = 3L, min_freq = 5L, method = "hon+")
  py <- .run_pyhon_plus(trajs, 3L, 5L)

  py$from <- vapply(py$from, .pipe_to_arrow, character(1L))
  py$to <- vapply(py$to, function(x) {
    arrow <- .pipe_to_arrow(x)
    tail(strsplit(arrow, " -> ", fixed = TRUE)[[1L]], 1L)
  }, character(1L))

  expect_equal(nrow(r$ho_edges), nrow(py),
    info = sprintf("Edge count: R=%d, Python=%d", nrow(r$ho_edges), nrow(py)))
  merged <- merge(
    r$ho_edges[, c("from", "to", "probability")],
    py, by = c("from", "to"), suffixes = c("_r", "_py"))
  expect_equal(nrow(merged), nrow(r$ho_edges),
    info = "Not all edges matched by from/to")
  expect_equal(merged$probability, merged$weight, tolerance = 1e-10)
})

# ===========================================================================
# Section 11: Coverage for previously uncovered paths
# ===========================================================================

# --- .hon_parse_input: all-NA row returns empty trajectory and is dropped ---
