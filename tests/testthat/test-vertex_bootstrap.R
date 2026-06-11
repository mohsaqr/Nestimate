# Tests for vertex_bootstrap() and bootstrap_network(ci_method =)

make_test_net <- function(method = "relative", seed = 7, n = 40) {
  set.seed(seed)
  states <- c("plan", "code", "debug", "test")
  seqs <- data.frame(
    T1 = sample(states, n, TRUE),
    T2 = sample(states, n, TRUE),
    T3 = sample(states, n, TRUE),
    T4 = sample(states, n, TRUE),
    stringsAsFactors = FALSE
  )
  build_network(seqs, method = method)
}

test_that("vertex_bootstrap returns a tidy net_vertex_bootstrap", {
  net <- make_test_net()
  vb <- vertex_bootstrap(net, iter = 100, seed = 1)

  expect_s3_class(vb, "net_vertex_bootstrap")
  expect_s3_class(vb$summary, "data.frame")
  expect_identical(
    names(vb$summary),
    c("statistic", "observed", "boot_mean", "boot_sd", "bias",
      "ci_lower", "ci_upper")
  )
  # directed network: all four built-ins, one row each
  expect_identical(
    vb$summary$statistic,
    c("density", "mean_weight", "centralization", "reciprocity")
  )
  expect_identical(dim(vb$boot_stats), c(100L, 4L))
  expect_true(all(vb$summary$ci_lower <= vb$summary$ci_upper))
})

test_that("vertex_bootstrap is reproducible with a seed", {
  net <- make_test_net()
  vb1 <- vertex_bootstrap(net, iter = 50, seed = 42)
  vb2 <- vertex_bootstrap(net, iter = 50, seed = 42)
  expect_identical(vb1$boot_stats, vb2$boot_stats)
})

test_that("undirected networks drop reciprocity and keep symmetry", {
  net <- make_test_net(method = "co_occurrence")
  expect_false(net$directed)
  vb <- vertex_bootstrap(net, iter = 50, seed = 3)
  expect_false("reciprocity" %in% vb$summary$statistic)

  # every replicate of a symmetric matrix stays symmetric
  set.seed(5)
  S <- matrix(runif(36), 6, 6)
  S <- (S + t(S)) / 2
  syms <- vapply(seq_len(200), function(i) {
    isSymmetric(unname(.vertex_resample(S, directed = FALSE)))
  }, logical(1))
  expect_true(all(syms))
})

test_that("replicate diagonals come only from the original diagonal", {
  set.seed(6)
  D <- matrix(runif(36), 6, 6)
  diag(D) <- 100 + seq_len(6)  # values that cannot occur off-diagonal
  ok <- vapply(seq_len(200), function(i) {
    Wb <- .vertex_resample(D, directed = TRUE)
    all(diag(Wb) %in% diag(D)) &&
      all(Wb[row(Wb) != col(Wb)] %in% D[row(D) != col(D)])
  }, logical(1))
  expect_true(all(ok))
})

test_that("vertex_bootstrap accepts a plain weight matrix", {
  net <- make_test_net()
  vb <- vertex_bootstrap(net$weights, iter = 50, seed = 2)
  expect_s3_class(vb, "net_vertex_bootstrap")
  expect_true(vb$directed)  # transition matrix is asymmetric
})

test_that("vertex_bootstrap works on a data-less netobject", {
  net <- make_test_net()
  net$data <- NULL
  expect_error(bootstrap_network(net, iter = 10), "does not contain \\$data")
  vb <- vertex_bootstrap(net, iter = 50, seed = 1)
  expect_s3_class(vb, "net_vertex_bootstrap")
})

test_that("custom statistic_fn rows are appended", {
  net <- make_test_net()
  vb <- vertex_bootstrap(
    net, iter = 50, seed = 4,
    statistic_fn = list(max_weight = function(W) max(W[row(W) != col(W)]))
  )
  expect_true("max_weight" %in% vb$summary$statistic)
  expect_identical(
    vb$summary$observed[vb$summary$statistic == "max_weight"],
    max(net$weights[row(net$weights) != col(net$weights)])
  )
})

test_that("statistics subset and validation work", {
  net <- make_test_net()
  vb <- vertex_bootstrap(net, iter = 50, seed = 1,
                         statistics = c("density", "reciprocity"))
  expect_identical(vb$summary$statistic, c("density", "reciprocity"))

  expect_error(
    vertex_bootstrap(net, iter = 50, statistics = "bogus"),
    "Unknown statistics"
  )
  expect_error(vertex_bootstrap(net, iter = 1), "iter")
  expect_error(vertex_bootstrap(list(a = 1)), "must be a netobject")
  small <- matrix(1:4, 2, 2)
  expect_error(vertex_bootstrap(small, iter = 10), "at least 3 nodes")
})

test_that("basic CIs are the reflection of percentile CIs", {
  net <- make_test_net()
  vp <- vertex_bootstrap(net, iter = 100, seed = 9, ci_method = "percentile")
  vb <- vertex_bootstrap(net, iter = 100, seed = 9, ci_method = "basic")
  expect_equal(vb$summary$ci_lower,
               2 * vp$summary$observed - vp$summary$ci_upper)
  expect_equal(vb$summary$ci_upper,
               2 * vp$summary$observed - vp$summary$ci_lower)
})

test_that("print, summary, and plot methods work", {
  net <- make_test_net()
  vb <- vertex_bootstrap(net, iter = 50, seed = 1)
  expect_output(print(vb), "Vertex Bootstrap")
  expect_s3_class(summary(vb), "data.frame")
  p <- plot(vb)
  expect_s3_class(p, "ggplot")
})

test_that("bootstrap_network ci_method = 'basic' reflects percentile CIs", {
  net <- make_test_net()
  bp <- bootstrap_network(net, iter = 100, seed = 5)
  bb <- bootstrap_network(net, iter = 100, seed = 5, ci_method = "basic")

  expect_identical(bp$ci_method, "percentile")
  expect_identical(bb$ci_method, "basic")
  expect_equal(bb$ci_lower, 2 * net$weights - bp$ci_upper)
  expect_equal(bb$ci_upper, 2 * net$weights - bp$ci_lower)
  # everything except the CIs must be unchanged
  expect_identical(bp$p_values, bb$p_values)
  expect_identical(bp$mean, bb$mean)

  expect_error(
    bootstrap_network(net, iter = 10, ci_method = "bogus"),
    "should be one of"
  )
})

# ---- vertex_compare() ----

test_that("vertex_compare returns a tidy two-network comparison", {
  net1 <- make_test_net(seed = 1)
  net2 <- make_test_net(seed = 2)
  cmp <- vertex_compare(net1, net2, iter = 100, seed = 5,
                        labels = c("g1", "g2"))

  expect_s3_class(cmp, "net_vertex_comparison")
  expect_identical(
    names(cmp$summary),
    c("statistic", "observed_g1", "observed_g2", "diff", "se_diff",
      "z", "p_value", "ci_lower", "ci_upper")
  )
  expect_identical(nrow(cmp$summary), 4L)
  # z consistent with diff / se where defined
  defined <- !is.na(cmp$summary$z)
  expect_equal(cmp$summary$z[defined],
               (cmp$summary$diff / cmp$summary$se_diff)[defined])
  # normal CI consistent with se
  expect_equal(cmp$summary$ci_upper - cmp$summary$ci_lower,
               2 * qnorm(0.975) * cmp$summary$se_diff)
})

test_that("vertex_compare of a network with itself gives zero differences", {
  net <- make_test_net()
  cmp <- vertex_compare(net, net, iter = 100, seed = 1)
  expect_true(all(cmp$summary$diff == 0))
  defined <- !is.na(cmp$summary$p_value)
  expect_true(all(cmp$summary$p_value[defined] == 1))
})

test_that("vertex_compare accepts precomputed net_vertex_bootstrap objects", {
  net1 <- make_test_net(seed = 1)
  net2 <- make_test_net(seed = 2)
  vb1 <- vertex_bootstrap(net1, iter = 100, seed = 1)
  vb2 <- vertex_bootstrap(net2, iter = 100, seed = 2)
  cmp <- vertex_compare(vb1, vb2)
  expect_s3_class(cmp, "net_vertex_comparison")
  expect_identical(cmp$x, vb1)
  expect_identical(cmp$y, vb2)
})

test_that("vertex_compare guards zero-SE statistics with NA", {
  net <- make_test_net()
  cmp <- vertex_compare(net, net, iter = 50, seed = 1,
                        statistics = "density")
  # density of a complete transition network never varies -> se 0 -> NA
  if (cmp$summary$se_diff == 0) {
    expect_true(is.na(cmp$summary$z))
    expect_true(is.na(cmp$summary$p_value))
  } else {
    succeed()
  }
})

test_that("vertex_compare print and plot methods work", {
  cmp <- vertex_compare(make_test_net(seed = 1), make_test_net(seed = 2),
                        iter = 50, seed = 1)
  expect_output(print(cmp), "Two-Network Vertex Bootstrap Comparison")
  expect_s3_class(summary(cmp), "data.frame")
  expect_s3_class(plot(cmp), "ggplot")
})

test_that("vertex_compare validates labels", {
  net <- make_test_net()
  expect_error(vertex_compare(net, net, iter = 10, labels = "one"),
               "labels")
})

test_that("vertex_compare errors on mismatched statistic sets", {
  net <- make_test_net()
  vb_full <- vertex_bootstrap(net, iter = 50, seed = 1)
  vb_density <- vertex_bootstrap(net, iter = 50, seed = 1,
                                 statistics = "density")
  expect_error(vertex_compare(vb_full, vb_density),
               "different statistics")
  # error message names the remedy with the shared set
  expect_error(vertex_compare(vb_full, vb_density),
               'statistics = c\\("density"\\)')

  # directed vs undirected default sets mismatch on reciprocity
  net_undir <- make_test_net(method = "co_occurrence")
  expect_error(vertex_compare(net, net_undir, iter = 50, seed = 1),
               "reciprocity")
  # aligning the selection fixes it
  cmp <- vertex_compare(net, net_undir, iter = 50, seed = 1,
                        statistics = c("density", "mean_weight"))
  expect_identical(cmp$summary$statistic, c("density", "mean_weight"))
})
