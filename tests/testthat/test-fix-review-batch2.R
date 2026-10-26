# Regression tests for the second batch of review fixes (2026-10-26):
# permutation column alignment, boot_glasso ties, betweenness ties, Ising
# resampling, failed-resample reporting, factor columns in higher-order input.
testthat::skip_on_cran()

test_that("permutation() on association networks ignores column order", {
  set.seed(5)
  make <- function() {
    d <- as.data.frame(matrix(rnorm(400), 100, 4,
                              dimnames = list(NULL, c("a", "b", "c", "d"))))
    d$b <- d$a + rnorm(100)
    d
  }
  dx <- make()
  dy <- make()
  invisible(lapply(c("cor", "glasso"), function(m) {
    x <- build_network(dx, method = m)
    y <- build_network(dy, method = m)
    y_rev <- build_network(dy[, c("d", "c", "b", "a")], method = m)
    s1 <- summary(permutation(x, y, iter = 60, seed = 1))
    s2 <- summary(permutation(x, y_rev, iter = 60, seed = 1))
    s2 <- s2[match(paste(s1$from, s1$to), paste(s2$from, s2$to)), ]
    expect_equal(s2$p_value, s1$p_value)
  }))
})

test_that("boot_glasso difference test counts ties in both tails", {
  boot <- cbind(e1 = c(0, 0, 0, 0), e2 = c(0, 0, 0, 0),
                e3 = c(1, 1, 1, 1), e4 = c(1, 0, 1, 0))
  p <- .bg_edge_diff_test(boot)
  expect_equal(p["e1", "e2"], 1)      # always tied: no evidence
  expect_equal(p["e1", "e3"], 0)      # e3 above e1 in every draw
  # e4 - e1 is 1, 0, 1, 0: P(>= 0) = 1, P(<= 0) = 0.5 -> p = 1
  expect_equal(p["e1", "e4"], 1)
  expect_true(isSymmetric(p))

  set.seed(2)
  X <- as.data.frame(matrix(rnorm(1500), 300))
  X$V2 <- X$V1 + rnorm(300)
  bg <- boot_glasso(X, iter = 60, seed = 1, cs_iter = 10)
  sig_edges <- sum(bg$edge_ci$ci_lower > 0 | bg$edge_ci$ci_upper < 0)
  n_diff <- sum(bg$edge_diff_p[upper.tri(bg$edge_diff_p)] < 0.05)
  # one present edge: at most it differs from each of the other nine
  expect_identical(sig_edges, 1L)
  expect_lte(n_diff, 9L)
})

test_that("betweenness splits tied shortest paths", {
  # a->b x10, b->c x15, a->c x6: 1/10 + 1/15 ties 1/6 exactly, but not in
  # floating point. Two shortest a->c paths, one through b: b gets 1/2.
  pairs <- c(rep(list(c("a", "b")), 10), rep(list(c("b", "c")), 15),
             rep(list(c("a", "c")), 6))
  d <- do.call(rbind, lapply(pairs, function(v) {
    data.frame(T1 = v[1], T2 = v[2])
  }))
  net <- build_network(d, method = "frequency")
  expect_equal(unname(.betweenness(net$weights)), c(0, 0.5, 0))
  eb <- .edge_betweenness(net$weights)
  expect_equal(c(eb["a", "b"], eb["a", "c"], eb["b", "c"]), c(1.5, 0.5, 1.5))
})

test_that("Ising networks resample: reliability, stability, bootstrap", {
  set.seed(4)
  y <- as.data.frame(matrix(rbinom(1500, 1, 0.5), 300))
  y$V2 <- ifelse(runif(300) < 0.8, y$V1, 1 - y$V1)
  y$V4 <- ifelse(runif(300) < 0.75, y$V3, 1 - y$V3)
  net <- build_network(y, method = "ising")
  rel <- network_reliability(net, iter = 6, seed = 1)
  expect_true(all(is.finite(summary(rel)$mean)))
  boot <- bootstrap_network(net, iter = 10, seed = 1)
  expect_true(all(is.finite(boot$mean)))
  expect_gt(boot$mean["V1", "V2"], 0)
})

test_that("resamplers report draws whose network cannot be estimated", {
  net <- build_network(srl_strategies, method = "glasso")
  broken <- net
  broken$data <- broken$data[, 1L, drop = FALSE]
  expect_warning(network_reliability(broken, iter = 3, seed = 1),
                 class = "nestimate_resample_failed")
  expect_no_warning(network_reliability(net, iter = 3, seed = 1))
})

test_that("factor and character columns mix in higher-order input", {
  df <- data.frame(T1 = factor(c("A", "B")), T2 = c("B", "C"),
                   T3 = factor(c("C", "A")))
  expect_identical(build_hon(df, max_order = 1)$first_order_states,
                   c("A", "B", "C"))
  expect_identical(build_mogen(df, max_order = 1)$states, c("A", "B", "C"))
  expect_setequal(state_frequencies(df)$state, c("A", "B", "C"))
  pc <- path_counts(df)
  expect_identical(pc$count[match(c("A -> B", "B -> C", "C -> A"), pc$path)],
                   c(1L, 2L, 1L))
})

test_that("build_network() refuses a compare_mmm() table", {
  set.seed(1)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C"), 60 * 8, TRUE),
                               60, 8))
  cmp <- suppressWarnings(compare_mmm(seqs, k = 2:3, n_starts = 2, seed = 1))
  expect_error(build_network(cmp, method = "relative"),
               class = "nestimate_bad_input")
  expect_error(build_network(cmp), class = "nestimate_bad_input")
})
