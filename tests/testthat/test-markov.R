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

test_that("markov_stability keeps recurrence times for separate closed classes", {
  P <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 1), 3, 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  result <- markov_stability(P)
  expect_equal(result$stability$return_time, c(2, 2, 1))
  expect_equal(result$stability$return_time, unname(result$mpt$return_times))
  expect_true(all(is.na(result$stability$stationary_prob)))
  expect_equal(result$mpt$matrix["A", "C"], Inf)
})

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
    "zero-sum row.*B.*no outgoing transitions"
  )
})

test_that("markov_stability errors on zero-sum row and names the dead state", {
  P_dead <- .P3
  P_dead["C", ] <- 0
  expect_error(
    markov_stability(P_dead),
    "zero-sum row.*C.*no outgoing transitions"
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

# --- Group dispatch -----------------------------------------------------

test_that("passage_time dispatches over netobject_group", {
  seqs <- data.frame(
    V1 = rep(c("A","B","C"), 6L),
    V2 = rep(c("B","C","A"), 6L),
    V3 = rep(c("C","A","B"), 6L),
    grp = rep(c("g1", "g2"), each = 9L),
    stringsAsFactors = FALSE
  )
  grp_net <- build_network(seqs, method = "relative", group = "grp")
  expect_s3_class(grp_net, "netobject_group")

  pt_grp <- passage_time(grp_net)
  expect_s3_class(pt_grp, "net_mpt_group")
  expect_named(pt_grp, names(grp_net))
  expect_true(all(vapply(pt_grp, inherits, logical(1), "net_mpt")))

  # Group result for each member must match the per-member call.
  for (nm in names(grp_net)) {
    expect_equal(pt_grp[[nm]]$matrix,
                 passage_time(grp_net[[nm]])$matrix,
                 tolerance = 1e-12,
                 info = sprintf("group %s", nm))
  }
  expect_invisible(print(pt_grp))
})

test_that("markov_stability dispatches over netobject_group", {
  seqs <- data.frame(
    V1 = rep(c("A","B","C"), 6L),
    V2 = rep(c("B","C","A"), 6L),
    V3 = rep(c("C","A","B"), 6L),
    grp = rep(c("g1", "g2"), each = 9L),
    stringsAsFactors = FALSE
  )
  grp_net <- build_network(seqs, method = "relative", group = "grp")

  ms_grp <- markov_stability(grp_net)
  expect_s3_class(ms_grp, "net_markov_stability_group")
  expect_named(ms_grp, names(grp_net))
  expect_true(all(vapply(ms_grp, inherits, logical(1),
                         "net_markov_stability")))
  for (nm in names(grp_net)) {
    expect_equal(ms_grp[[nm]]$stability,
                 markov_stability(grp_net[[nm]])$stability,
                 info = sprintf("group %s", nm))
  }
  expect_invisible(print(ms_grp))
})

test_that("passage_time gives Inf for states not reached with certainty", {
  # a <-> b transient, c absorbing. By hand: m_ab = 1 / 0.5 = 2;
  # m_ac = 1 + .5 m_ac + .5 m_bc and m_bc = 1 + .2 m_ac + .3 m_bc give
  # m_ac = 4.8, m_bc = 2.8; c never leaves, so every time out of c is Inf,
  # and b -> a is Inf (b may be absorbed first). The former epsilon floor
  # returned ~1e15 here with a warning.
  P <- rbind(c(.5, .5, 0), c(.2, .3, .5), c(0, 0, 1))
  dimnames(P) <- list(c("a", "b", "c"), c("a", "b", "c"))
  pt <- expect_silent(passage_time(P))
  M <- pt$matrix
  expect_equal(M["a", "b"], 2)
  expect_equal(M["a", "c"], 4.8)
  expect_equal(M["b", "c"], 2.8)
  expect_identical(unname(M["c", c("a", "b")]), c(Inf, Inf))
  expect_identical(M["b", "a"], Inf)
  expect_identical(unname(pt$return_times), c(Inf, Inf, 1))
  expect_equal(unname(pt$stationary), c(0, 0, 1))

  # Two absorbing states: no error, Inf across, NA stationary (not unique).
  P2 <- rbind(c(1, 0, 0), c(.2, .3, .5), c(0, 0, 1))
  dimnames(P2) <- dimnames(P)
  pt2 <- passage_time(P2)
  expect_identical(unname(diag(pt2$matrix)), c(1, Inf, 1))
  expect_true(all(is.infinite(pt2$matrix[row(P2) != col(P2)])))
  expect_true(all(is.na(pt2$stationary)))
})

test_that("hitting-time route equals Kemeny-Snell on an irreducible chain", {
  P <- rbind(c(.2, .5, .3), c(.4, .1, .5), c(.3, .3, .4))
  dimnames(P) <- list(c("a", "b", "c"), c("a", "b", "c"))
  pi <- .mpt_stationary(P)
  expect_equal(.mpt_hitting_times(P), {
    M <- .mpt_full(P, pi); dimnames(M) <- dimnames(P); M
  }, tolerance = 1e-10)
})
