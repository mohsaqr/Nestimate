# Regression tests for the review of 2026-11-20, higher-order area (H1-H17).
# Each expectation failed on the code before the fix, and encodes a value
# from the definition or a hand computation.

# ---- H1: sequence_compare() permutation null ignores the min_freq filter ----

test_that("sequence_compare() p-values of a pattern do not depend on min_freq", {
  # the statistic is computed on every pattern of the length, so the patterns
  # a stricter filter drops cannot change the null of the ones it keeps
  set.seed(3)
  m <- matrix(sample(c("A", "B", "C", "D"), 20 * 8, TRUE,
                     prob = c(.45, .3, .15, .1)), 20, 8)
  g <- rep(c("g1", "g2"), each = 10)
  set.seed(1)
  loose <- summary(sequence_compare(m, group = g, sub = 2, min_freq = 1,
                                    iter = 200, adjust = "none"))
  set.seed(1)
  strict <- summary(sequence_compare(m, group = g, sub = 2, min_freq = 4,
                                     iter = 200, adjust = "none"))
  both <- merge(loose, strict, by = "pattern", suffixes = c("_loose", "_strict"))
  expect_gt(nrow(loose), nrow(strict))
  expect_identical(nrow(both), nrow(strict))
  expect_identical(both$p_value_loose, both$p_value_strict)
  expect_identical(both$effect_size_loose, both$effect_size_strict)
})

test_that("the pattern statistic is the row norm of observed minus expected over all patterns", {
  counts <- rbind(c(4, 0), c(2, 6), c(1, 3))
  tot <- colSums(counts)
  # pattern 1: e = r * c / N = 4 * (7, 9) / 16
  expected_1 <- 4 * tot / sum(tot)
  expect_equal(Nestimate:::.cs_pattern_statistic(counts)[1],
               sqrt(sum((counts[1, ] - expected_1)^2)), tolerance = 1e-12)
})

# ---- H4: numeric states are not padded per group ----

test_that("sequence_compare() reads numeric states the same in every group", {
  d <- data.frame(g = c("x", "x", "y", "y"), T1 = c(1, 10, 1, 1), T2 = "A")
  net <- build_network(d, method = "relative", group = "g",
                       state_cols = c("T1", "T2"))
  res <- summary(sequence_compare(net, sub = 2, min_freq = 1, test = "none"))
  # x holds 1->A once (and 10->A), y holds 1->A twice
  expect_identical(res$pattern, "1->A")
  expect_identical(res$freq_x, 1L)
  expect_identical(res$freq_y, 2L)
  wide <- summary(sequence_compare(d, group = "g", sub = 2, min_freq = 1,
                                   test = "none"))
  expect_identical(wide$freq_y, res$freq_y)
})

# ---- H5: the LRT starts from the first-order null ----

test_that("build_mogen() LRT does not select order 0 for first-order data", {
  x <- rep(list(c("A", "B", "A", "B")), 20)
  m <- build_mogen(x, max_order = 2, criterion = "lrt")
  expect_identical(m$optimal_order, 1L)
  # order 1 is the null of the tests (pathpy estimate_order), never tested
  expect_true(is.na(m$lrt_p[["order_1"]]))
  expect_identical(Nestimate:::.mogen_lrt_order(c(NA, 0.5, 0.5), 0.01), 1L)
  expect_identical(Nestimate:::.mogen_lrt_order(c(NA, 1e-9), 0.01), 1L)
})

test_that("build_mogen() LRT selects order 2 for second-order data", {
  x <- c(rep(list(c("A", "X", "C")), 15), rep(list(c("B", "X", "D")), 15))
  m <- build_mogen(x, max_order = 2, criterion = "lrt")
  expect_identical(m$optimal_order, 2L)
  # layer 2 has 2 path dof (X has two successors, reached from A and B):
  # a chi-square on 2 df has survival exp(-x / 2)
  expect_identical(unname(m$layer_dof[["order_2"]]), 2)
  lr <- -2 * (m$log_likelihood[["order_1"]] - m$log_likelihood[["order_2"]])
  expect_equal(m$lrt_p[["order_2"]], exp(-lr / 2), tolerance = 1e-10)
})

# ---- H6: an order-0 fit has a node table that matches its edges ----

test_that("build_mogen() at order 0 returns the states with no edges", {
  x <- list(c("A", "A"), c("A", "B"), c("B", "A"), c("B", "B"))
  m <- build_mogen(x, max_order = 1)
  expect_identical(m$optimal_order, 0L)
  expect_identical(m$nodes$label, c("A", "B"))
  expect_identical(m$n_nodes, 2L)
  expect_identical(nrow(m$edges), 0L)
  expect_identical(m$n_edges, 0L)
  expect_identical(dim(m$weights), c(2L, 2L))
  expect_true(all(m$weights == 0))
  # the marginal is still there
  expect_equal(m$transition_matrices[[1L]], c(A = 0.5, B = 0.5))
})

# ---- H7: a state containing the node separator is refused ----

test_that("build_hon() refuses a state containing ' -> '", {
  x <- c(rep(list(c("A", "B", "C")), 3), rep(list(c("B", "D")), 12),
         list(c("A -> B", "E")))
  expect_error(build_hon(x, max_order = 2, method = "hon"),
               class = "nestimate_bad_input")
  # "A->B" (no spaces) is not the separator: one row per node, each a
  # distribution (row sums 1 or 0)
  ok <- c(rep(list(c("A", "B", "C")), 3), rep(list(c("B", "D")), 12),
          list(c("A->B", "E")))
  h <- build_hon(ok, max_order = 2, method = "hon")
  expect_false(anyDuplicated(rownames(h$weights)) > 0L)
  rs <- unname(rowSums(h$weights))
  expect_true(all(abs(rs - 1) < 1e-12 | abs(rs) < 1e-12))
})

# ---- H8: entropy_bayes() on one state ----

test_that("entropy_bayes() refuses a one-state chain with a classed error", {
  expect_error(entropy_bayes(matrix(12, 1, 1, dimnames = list("A", "A")),
                             draws = 100, seed = 1),
               class = "nestimate_bad_input")
})

# ---- H10: entropy bases must exceed 1 ----

test_that("entropy verbs refuse a base at or below 1", {
  p <- matrix(.5, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(transition_entropy(p, base = 0.5), class = "nestimate_bad_input")
  expect_error(entropy_network(p, base = 0.5), class = "nestimate_bad_input")
  expect_error(entropy_network(p, base = Inf), class = "nestimate_bad_input")
  expect_error(entropy_bayes(p * 10, base = 0.5, draws = 100, seed = 1),
               class = "nestimate_bad_input")
  d <- data.frame(action = rep(c("A", "B"), 20))
  expect_error(entropy_trajectory(d, action = "action", window = 10,
                                  base = 0.5),
               class = "nestimate_bad_input")
})

test_that("scaling = 'share' gives each edge its percentage of the entropy rate", {
  # uniform 2-state chain: pi = 1/2, P = 1/2, -log2 P = 1, so every edge
  # contributes 1/4 of a rate of 1 bit: 25%
  p <- matrix(.5, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  w <- entropy_network(p, scaling = "share")$weights
  expect_equal(unname(w), matrix(25, 2, 2), tolerance = 1e-12)
  te <- transition_entropy(p)
  expect_equal(te$entropy_rate, 1, tolerance = 1e-12)
  expect_equal(te$entropy_rate_norm, 1, tolerance = 1e-12)
})

# ---- H11: markov_stability() group plot keeps states with undefined pi ----

test_that("markov_stability() group plot keeps every state", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(g = rep(c("x", "y"), each = 2),
                  T1 = rep(c("A", "B"), 2), T2 = rep(c("A", "B"), 2))
  net <- build_network(d, method = "relative", group = "g")
  ms <- markov_stability(net)
  p <- plot(ms, metrics = "return_time")
  expect_false(anyNA(p$data$state))
  expect_setequal(levels(p$data$state), c("A", "B"))
})

# ---- H12: entropy_network() returns a bare cograph_network as a netobject ----

test_that("entropy_network() returns a bare cograph_network as a netobject, as documented", {
  x <- build_network(data.frame(T1 = c("A", "B", "A"), T2 = c("B", "A", "A")),
                     method = "relative")
  class(x) <- "cograph_network"
  y <- entropy_network(x)
  expect_identical(class(y), c("netobject", "cograph_network"))
  expect_identical(y$method, "entropy")
  expect_identical(y$nodes$label, x$nodes$label)
  # the contributions sum to the entropy rate
  expect_equal(sum(y$weights), transition_entropy(x)$entropy_rate,
               tolerance = 1e-12)
})

# ---- H13: entropy_trajectory() with a missing group ----

test_that("entropy_trajectory() refuses a missing group with a classed error", {
  d <- data.frame(action = c("A", "B", "A"), group = c("x", "x", NA))
  expect_error(entropy_trajectory(d, action = "action", group = "group",
                                  window = 10),
               class = "nestimate_bad_input")
  d$a <- 1
  expect_error(entropy_trajectory(d, action = "action", actor = "a",
                                  group = "group", window = 10),
               class = "nestimate_bad_input")
})

# ---- H14: truncation warning at an explicit cap equal to the build cap ----

test_that("persistent_homology() warns on a truncated 0-dimensional VR complex at its own cap", {
  d <- matrix(c(0, 1, 1, 0), 2)
  s <- build_simplicial(d, type = "vr", max_dim = 0, max_scale = 2)
  expect_true(s$truncated)
  expect_warning(persistent_homology(s, max_dim = 0),
                 class = "nestimate_truncated_complex")
  expect_warning(persistent_homology(s, max_dim = 0, max_scale = 2),
                 class = "nestimate_truncated_complex")
})

# ---- H16: a pattern template ignores type, len and gap ----

test_that("sequence_patterns() with a pattern ignores type, len and gap", {
  d <- data.frame(T1 = "A", T2 = "B")
  ref <- sequence_patterns(d, pattern = "A->B", min_freq = 1)
  expect_identical(ref$pattern, "A->B")
  expect_identical(ref$occurrences, 1L)
  expect_identical(sequence_patterns(d, pattern = "A->B", len = integer(),
                                     min_freq = 1), ref)
  expect_identical(sequence_patterns(d, pattern = "A->B", type = "bogus",
                                     min_freq = 1), ref)
})

# ---- H17: max_power is the upper limit of k, the sum runs to D^(L+1) ----

test_that("build_honem() sums exp(-k) D^(k+1) for k = 0..max_power", {
  d <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, byrow = TRUE)
  S <- (d + exp(-1) * d %*% d) / (1 + exp(-1))
  em <- build_honem(d, dim = 1, max_power = 1)
  expect_equal(em$singular_values, svd(S)$d[1], tolerance = 1e-12)
  expect_true(any(grepl("D^2", capture.output(print(em)), fixed = TRUE)))
})
