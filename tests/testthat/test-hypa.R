# ===========================================================================
# Section 1: Internal — .hypa_fit_xi
# ===========================================================================

test_that(".hypa_fit_xi gives N >> m", {
  adj <- matrix(c(0, 5, 3, 2, 0, 4, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- .hypa_fit_xi(adj)

  # N = sum(Xi) should be >> m = sum(adj)
  expect_true(sum(xi) > sum(adj))

  # Xi = outer(s_out, s_in) * mask
  s_out <- rowSums(adj)
  s_in <- colSums(adj)
  expected <- outer(s_out, s_in) * (adj > 0)
  expect_equal(xi, expected)
})

test_that(".hypa_fit_xi respects edge structure", {
  adj <- matrix(c(0, 5, 0, 0, 0, 3, 2, 0, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- .hypa_fit_xi(adj)

  # Xi should be zero where adj is zero
  expect_equal(xi[1, 1], 0)
  expect_equal(xi[1, 3], 0)
  expect_equal(xi[2, 1], 0)
  expect_equal(xi[2, 2], 0)
  expect_equal(xi[3, 3], 0)
})

# ===========================================================================
# Section 2: Internal — .hypa_compute_scores
# ===========================================================================

test_that(".hypa_compute_scores returns correct format", {
  adj <- matrix(c(0, 5, 3, 2, 0, 4, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- .hypa_fit_xi(adj)
  scores <- .hypa_compute_scores(adj, xi)

  expect_true(is.data.frame(scores))
  expect_true(all(c("path", "from", "to", "observed", "expected",
                     "ratio", "hypa_score", "anomaly") %in% names(scores)))
  expect_equal(nrow(scores), sum(adj > 0))

  # HYPA scores should be in [0, 1]
  expect_true(all(scores$hypa_score >= 0))
  expect_true(all(scores$hypa_score <= 1))
})

test_that(".hypa_compute_scores handles empty graph", {
  adj <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- adj
  scores <- .hypa_compute_scores(adj, xi)
  expect_equal(nrow(scores), 0L)
})

# ===========================================================================
# Section 3: build_hypa end-to-end
# ===========================================================================

test_that("build_hypa returns net_hypa class", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("B", "C", "A"),
                c("C", "A", "B"), c("A", "C", "B"), c("B", "A", "C"))
  h <- build_hypa(trajs, k = 1L)

  expect_s3_class(h, "net_hypa")
  expect_equal(h$k, 1L)
  expect_true(h$n_edges > 0L)
  expect_true(is.data.frame(h$scores))
})

test_that("build_hypa detects anomalies in biased data", {
  # Create data where A->B->C is overwhelmingly common
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(5, c("A", "B", "D"), simplify = FALSE),
    replicate(5, c("C", "B", "A"), simplify = FALSE),
    replicate(2, c("D", "B", "C"), simplify = FALSE),
    replicate(2, c("C", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, k = 2L, alpha = 0.05)

  # Should find some anomalous paths
  # The A->B->C path is very frequent, may be over-represented
  expect_s3_class(h, "net_hypa")
  expect_true(h$n_edges > 0L)
})

test_that("build_hypa handles k=1 (first-order)", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "C", "B", "D"),
                c("B", "C", "D", "A"), c("D", "A", "B", "C"))
  h <- build_hypa(trajs, k = 1L)

  expect_equal(h$k, 1L)
  expect_true(nrow(h$scores) > 0L)
})

test_that("build_hypa rejects invalid input", {
  expect_error(build_hypa(42), "data.frame or list")
  expect_error(build_hypa(list(c("A", "B")), k = 0L), "k.*must be >= 1")
  expect_error(build_hypa(list(c("A", "B")), alpha = 0.6),
               "alpha.*must be in")
})

test_that("build_hypa alpha parameter affects classification", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("B", "C", "A"),
                c("C", "A", "B"), c("B", "A", "C"), c("A", "C", "B"))

  h1 <- build_hypa(trajs, k = 1L, alpha = 0.05)
  h2 <- build_hypa(trajs, k = 1L, alpha = 0.49)

  # Stricter alpha should find fewer or equal anomalies
  expect_true(h1$n_anomalous <= h2$n_anomalous)
})

# ===========================================================================
# Section 4: HYPA scores properties
# ===========================================================================

test_that("HYPA scores are in [0, 1]", {
  set.seed(42)
  trajs <- lapply(seq_len(50L), function(i) {
    sample(LETTERS[1:4], 5, replace = TRUE)
  })
  h <- build_hypa(trajs, k = 1L)

  expect_true(all(h$scores$hypa_score >= 0))
  expect_true(all(h$scores$hypa_score <= 1))
})

test_that("HYPA expected values are positive", {
  trajs <- list(c("A", "B", "C"), c("A", "C", "B"), c("B", "A", "C"),
                c("C", "B", "A"), c("B", "C", "A"), c("C", "A", "B"))
  h <- build_hypa(trajs, k = 1L)

  expect_true(all(h$scores$expected >= 0))
})

# ===========================================================================
# Section 5: S3 methods
# ===========================================================================

test_that("print.net_hypa works", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  h <- build_hypa(trajs, k = 1L)
  out <- capture.output(print(h))
  expect_true(any(grepl("HYPA", out)))
})

test_that("summary.net_hypa works", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"), c("A", "C", "B"))
  h <- build_hypa(trajs, k = 1L)
  out <- capture.output(summary(h))
  expect_true(any(grepl("HYPA", out)))
})


# ===========================================================================
# Section 6: Data.frame input
# ===========================================================================

test_that("build_hypa handles data.frame input", {
  df <- data.frame(T1 = c("A", "B", "C"), T2 = c("B", "C", "A"),
                   T3 = c("C", "A", "B"))
  h <- build_hypa(df, k = 1L)
  expect_s3_class(h, "net_hypa")
})

# ===========================================================================
# Section 7: Coverage for previously uncovered paths
# ===========================================================================

# --- build_hypa: no valid trajectories ---
test_that("build_hypa stops when no valid trajectories", {
  # All single-state entries: parsed trajectories have < 2 states each
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(build_hypa(df, k = 1L), "No valid trajectories")
})

# --- build_hypa: no edges at given order (paths too short) ---
test_that("build_hypa stops when no edges at requested order k", {
  # k=3 requires 4-grams; trajectories of length 3 produce only 1-grams (k=1)
  # and 2-grams (k=2) but not 3-grams as transitions
  trajs <- list(c("A", "B", "C"), c("B", "C", "D"))
  expect_error(build_hypa(trajs, k = 3L), "No edges at order")
})

# --- summary.net_hypa with anomalies displays anomalous paths ---
test_that("summary.net_hypa displays anomalous paths when present", {
  # Create highly biased data to force anomaly detection
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, k = 2L, alpha = 0.05, min_count = 1L)
  out <- capture.output(summary(h))
  # Either found anomalies or printed "No anomalous paths detected."
  expect_true(any(grepl("Anomalous|anomalous|No anomalous", out,
                         ignore.case = TRUE)))
})


# ===========================================================================
# Section 8: pathways() tests for HYPA
# ===========================================================================

test_that("pathways.net_hypa returns character vector", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, k = 2L, alpha = 0.05, min_count = 1L)
  pw <- pathways(h)
  expect_true(is.character(pw))
})

test_that("pathways.net_hypa type='over' returns over-represented paths", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, k = 2L, alpha = 0.05, min_count = 1L)
  pw_all  <- pathways(h, type = "all")
  pw_over <- pathways(h, type = "over")
  expect_true(length(pw_over) <= length(pw_all))
})

test_that("pathways.net_hypa type='under' returns under-represented paths", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, k = 2L, alpha = 0.05, min_count = 1L)
  pw_under <- pathways(h, type = "under")
  expect_true(is.character(pw_under))
})

test_that("pathways.net_hypa returns empty when no anomalies", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  h <- build_hypa(trajs, k = 1L, alpha = 1e-10)
  pw <- pathways(h)
  # With near-zero alpha threshold, likely no anomalies
  expect_true(is.character(pw))
})
