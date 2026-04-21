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

test_that(".mpt_stationary returns valid stationary distribution", {
  pi <- Nestimate:::.mpt_stationary(.P3)
  expect_equal(sum(pi), 1, tolerance = 1e-10)
  expect_true(all(pi > 0))
  expect_length(pi, 3)
  # Verify π P = π
  expect_equal(as.vector(pi %*% .P3), pi, tolerance = 1e-8)
})

# ---- passage_time() basics ----

test_that("passage_time returns net_mpt with correct structure", {
  pt <- passage_time(.P3)
  expect_s3_class(pt, "net_mpt")
  expect_true(is.matrix(pt$matrix))
  expect_equal(dim(pt$matrix), c(3L, 3L))
  expect_named(pt, c("matrix", "stationary", "return_times", "states"))
  expect_equal(pt$states, c("A", "B", "C"))
})

test_that("passage_time diagonal equals return times (1/pi)", {
  pt <- passage_time(.P3)
  expect_equal(diag(pt$matrix), pt$return_times, tolerance = 1e-8)
  expect_equal(pt$return_times, 1 / pt$stationary, tolerance = 1e-8)
})

test_that("passage_time is asymmetric for directed chain", {
  pt <- passage_time(.P3)
  # M[i,j] != M[j,i] in general for an asymmetric chain
  expect_false(isSymmetric(pt$matrix))
})

test_that("passage_time stationary sums to 1", {
  pt <- passage_time(.P3)
  expect_equal(sum(pt$stationary), 1, tolerance = 1e-10)
})

test_that("passage_time accepts a netobject", {
  net <- build_network(
    data.frame(
      V1 = c("A","B","C","A","B"),
      V2 = c("B","C","A","C","A"),
      stringsAsFactors = FALSE
    ),
    method = "relative"
  )
  pt <- passage_time(net)
  expect_s3_class(pt, "net_mpt")
  expect_equal(ncol(pt$matrix), nrow(pt$matrix))
})

test_that("passage_time states argument subsets correctly", {
  pt <- passage_time(.P3, states = c("A", "B"))
  expect_equal(pt$states, c("A", "B"))
  expect_equal(dim(pt$matrix), c(2L, 2L))
})

test_that("passage_time errors on unknown states", {
  expect_error(passage_time(.P3, states = c("A", "Z")), "Unknown states")
})

test_that("passage_time normalizes non-stochastic matrix with warning", {
  P_bad <- .P3 * 2
  expect_warning(pt <- passage_time(P_bad), "normalizing")
  expect_s3_class(pt, "net_mpt")
})

test_that("passage_time normalize=FALSE errors on non-stochastic input", {
  P_bad <- .P3 * 2
  expect_error(passage_time(P_bad, normalize = FALSE), "must sum to 1")
})

test_that("passage_time errors on unsupported input", {
  expect_error(passage_time("not_a_matrix"), "numeric matrix")
})

# ---- zero-row guard (regression: silent NaN from row_sums == 0) ----

test_that("passage_time errors on zero-sum row and names the dead state", {
  P_dead <- .P3
  P_dead["B", ] <- 0  # B has no outgoing transitions
  expect_error(
    passage_time(P_dead),
    "zero-sum row.*B.*not.*ergodic"
  )
})

test_that("markov_stability errors on zero-sum row and names the dead state", {
  P_dead <- .P3
  P_dead["C", ] <- 0
  expect_error(
    markov_stability(P_dead),
    "zero-sum row.*C.*not.*ergodic"
  )
})

test_that("passage_time zero-row error lists multiple dead states", {
  P_dead <- .P3
  P_dead[c("A", "C"), ] <- 0
  expect_error(
    passage_time(P_dead),
    "zero-sum row.*A, C"
  )
})

test_that("passage_time zero-row guard fires before normalize=FALSE check", {
  P_dead <- .P3
  P_dead["B", ] <- 0
  # Even with normalize=FALSE, the zero-row diagnostic takes priority over
  # the generic "must sum to 1" error — the message must point at B.
  expect_error(
    passage_time(P_dead, normalize = FALSE),
    "zero-sum row.*B"
  )
})

# ---- Numerical equivalence vs linear-system solve ----

test_that("passage_time matches column-by-column linear solve", {
  # Reference: solve column j of MPT independently
  n <- 3L
  P <- .P3
  pi_ref <- Nestimate:::.mpt_stationary(P)

  M_ref <- matrix(0, n, n, dimnames = dimnames(P))
  for (j in seq_len(n)) {
    A        <- diag(n) - P
    A[j, ]   <- 0
    A[j, j]  <- 1
    b        <- rep(1, n); b[j] <- 0
    M_ref[, j] <- solve(A, b)
  }
  diag(M_ref) <- 1 / pi_ref

  pt <- passage_time(P)
  expect_equal(pt$matrix, M_ref, tolerance = 1e-8)
})

# ---- print / summary ----

test_that("print.net_mpt runs without error", {
  pt <- passage_time(.P3)
  expect_output(print(pt), "Mean First Passage Times")
})

test_that("summary.net_mpt returns expected data frame", {
  pt <- passage_time(.P3)
  sm <- summary(pt)
  expect_s3_class(sm, "summary.net_mpt")
  expect_true(is.data.frame(sm$table))
  expect_equal(nrow(sm$table), 3L)
  expect_named(sm$table, c("state", "return_time", "stationary",
                            "mean_out", "mean_in"))
})

# ---- plot.net_mpt ----

test_that("plot.net_mpt returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  pt <- passage_time(.P3)
  g  <- plot(pt)
  expect_s3_class(g, "gg")
})

test_that("plot.net_mpt with log_scale=FALSE runs", {
  skip_if_not_installed("ggplot2")
  pt <- passage_time(.P3)
  expect_s3_class(plot(pt, log_scale = FALSE), "gg")
})

# ---- markov_stability() ----

test_that("markov_stability returns net_markov_stability with correct structure", {
  ms <- markov_stability(.P3)
  expect_s3_class(ms, "net_markov_stability")
  expect_named(ms, c("stability", "mpt"))
  df <- ms$stability
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 3L)
  expect_named(df, c("state", "persistence", "stationary_prob",
                     "return_time", "sojourn_time",
                     "avg_time_to_others", "avg_time_from_others"))
})

test_that("markov_stability persistence equals diagonal of P", {
  ms <- markov_stability(.P3)
  expect_equal(ms$stability$persistence, unname(round(diag(.P3), 4)))
})

test_that("markov_stability stationary_prob sums to 1", {
  ms <- markov_stability(.P3)
  expect_equal(sum(ms$stability$stationary_prob), 1, tolerance = 1e-3)
})

test_that("markov_stability sojourn_time is 1/(1-persistence)", {
  ms <- markov_stability(.P3)
  df <- ms$stability
  expected <- round(1 / (1 - df$persistence), 2)
  expect_equal(df$sojourn_time, expected, tolerance = 1e-4)
})

test_that("markov_stability mpt field is net_mpt", {
  ms <- markov_stability(.P3)
  expect_s3_class(ms$mpt, "net_mpt")
})

test_that("print.net_markov_stability runs without error", {
  ms <- markov_stability(.P3)
  expect_output(print(ms), "Markov Stability")
})

test_that("summary.net_markov_stability prints attractor info", {
  ms <- markov_stability(.P3)
  expect_output(summary(ms), "Most accessible")
})

test_that("plot.net_markov_stability returns a ggplot", {
  skip_if_not_installed("ggplot2")
  ms <- markov_stability(.P3)
  g  <- plot(ms)
  expect_s3_class(g, "gg")
})

test_that("plot.net_markov_stability metrics argument filters panels", {
  skip_if_not_installed("ggplot2")
  ms <- markov_stability(.P3)
  g  <- plot(ms, metrics = c("persistence", "return_time"))
  expect_s3_class(g, "gg")
})

# ---- integration with trajectories dataset ----

test_that("passage_time works on trajectories netobject", {
  net <- build_network(as.data.frame(trajectories), method = "relative")
  pt  <- passage_time(net)
  expect_s3_class(pt, "net_mpt")
  expect_equal(length(pt$states), 3L)
  # Diagonal must equal return times
  expect_equal(diag(pt$matrix), pt$return_times, tolerance = 1e-6)
  # All MFPT values must be positive
  expect_true(all(pt$matrix > 0))
})

test_that("markov_stability works on trajectories netobject", {
  net <- build_network(as.data.frame(trajectories), method = "relative")
  ms  <- markov_stability(net)
  expect_s3_class(ms, "net_markov_stability")
  expect_equal(nrow(ms$stability), 3L)
})

# ---- tna object support ----

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

test_that("passage_time accepts a wide sequence data.frame", {
  seqs <- data.frame(
    V1 = c("A","B","C","A","B"),
    V2 = c("B","C","A","B","C"),
    V3 = c("C","A","B","A","A"),
    stringsAsFactors = FALSE
  )
  pt <- passage_time(seqs)
  expect_s3_class(pt, "net_mpt")
  expect_equal(length(pt$states), 3L)
})

test_that("passage_time on wide data.frame matches netobject route", {
  seqs <- data.frame(
    V1 = c("A","B","C","A","B"),
    V2 = c("B","C","A","B","C"),
    V3 = c("C","A","B","A","A"),
    stringsAsFactors = FALSE
  )
  pt_direct <- passage_time(seqs)
  net       <- build_network(seqs, method = "relative")
  pt_via    <- passage_time(net)
  expect_equal(pt_direct$matrix, pt_via$matrix, tolerance = 1e-8)
})

test_that("markov_stability accepts a wide sequence data.frame", {
  seqs <- data.frame(
    V1 = c("A","B","C","A","B"),
    V2 = c("B","C","A","B","C"),
    stringsAsFactors = FALSE
  )
  ms <- markov_stability(seqs)
  expect_s3_class(ms, "net_markov_stability")
})
