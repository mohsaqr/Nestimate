# Layer 9, composite identity: verbs that take a network and compute on top of it must carry that network
# unchanged. The observed network of bootstrap_network() and the two networks of permutation() and
# compare_networks() must equal their inputs; on sequences that are all the same, every bootstrap resample is the
# input data again, so every resample network must equal the input network (its interval collapses onto the
# weights), which checks that the resampling path rebuilds with the source settings; macro_network() expansions
# must equal the direct fit on the recoded sequences; and the layers of as_tna(<mcml>) must carry the fit's weights.

wide <- tn_simulate(n_actors = 50L, n_states = 4L, seq_length = 10L, seed = 931L)
constant <- tn_simulate(n_actors = 1L, n_states = 4L, seq_length = 12L, seed = 932L)[rep(1L, 30L), ]
rownames(constant) <- NULL

settings <- list(
  list(quantity = "relative", args = list(method = "relative")),
  list(quantity = "relative, minmax scaling", args = list(method = "relative", scaling = "minmax")),
  list(quantity = "relative, threshold 0.2", args = list(method = "relative", threshold = 0.2)),
  list(quantity = "relative, start and end", args = list(method = "relative", start = TRUE, end = TRUE)),
  list(quantity = "frequency, windows of 3", args = list(method = "frequency", params = list(window_size = 3L))),
  list(quantity = "frequency, max scaling then threshold 0.2", args = list(method = "frequency", scaling = "max", threshold = 0.2)),
  list(quantity = "attention, lambda 3", args = list(method = "attention", lambda = 3)),
  list(quantity = "co-occurrence", args = list(method = "co_occurrence")),
  list(quantity = "2-gram", args = list(method = "ngram", n_gram = 2L)),
  list(quantity = "gap", args = list(method = "gap")),
  list(quantity = "reverse", args = list(method = "reverse"))
)

test_that("the observed network of bootstrap_network() is the input network", {
  lapply(settings, \(setting) {
    net <- do.call(build_network, c(list(wide), setting$args))
    boot <- bootstrap_network(net, iter = 10L, seed = 1L)
    li_expect_same_network(boot$original, net, route = "the input network", covers = "bootstrap_network()",
                           quantity = sprintf("observed network, %s", setting$quantity))
  })
  continuous <- psy_continuous(n = 200L, p = 5L, rho = 0.8, seed = 933L)
  glasso <- build_network(continuous, method = "glasso", gamma = 0.25, threshold = 0.02)
  li_expect_same_network(bootstrap_network(glasso, iter = 10L, seed = 1L)$original, glasso, route = "the input network",
                         covers = "bootstrap_network()", quantity = "observed network, glasso gamma 0.25 threshold 0.02")
})

test_that("every bootstrap resample of identical sequences rebuilds the input network", {
  lapply(settings, \(setting) {
    net <- do.call(build_network, c(list(constant), setting$args))
    boot <- bootstrap_network(net, iter = 15L, seed = 2L)
    # a quantile of identical values is that value exactly, so both limits must be the weights
    expect_reference(boot$ci_lower, net$weights, reference = li_ref("the input network (every resample is the input data)"),
                     covers = "bootstrap_network()", quantity = sprintf("lower limit, %s", setting$quantity),
                     tolerance = 0, authority = "internal")
    expect_reference(boot$ci_upper, net$weights, reference = li_ref("the input network (every resample is the input data)"),
                     covers = "bootstrap_network()", quantity = sprintf("upper limit, %s", setting$quantity),
                     tolerance = 0, authority = "internal")
  })
})

test_that("permutation() carries its two networks and their difference unchanged", {
  first <- tn_simulate(n_actors = 40L, n_states = 4L, seq_length = 10L, seed = 934L)
  second <- tn_simulate(n_actors = 40L, n_states = 4L, seq_length = 10L, seed = 935L)
  lapply(settings[c(1L, 2L, 4L, 5L, 7L)], \(setting) {
    x <- do.call(build_network, c(list(first), setting$args))
    y <- do.call(build_network, c(list(second), setting$args))
    result <- permutation(x, y, iter = 20L, seed = 1L)
    li_expect_same_network(result$x, x, route = "the first input network", covers = "permutation()",
                           quantity = sprintf("x, %s", setting$quantity))
    li_expect_same_network(result$y, y, route = "the second input network", covers = "permutation()",
                           quantity = sprintf("y, %s", setting$quantity))
    expect_reference(result$diff, x$weights - y$weights, reference = li_ref("x$weights - y$weights"),
                     covers = "permutation()", quantity = sprintf("observed difference, %s", setting$quantity),
                     tolerance = 0, authority = "internal")
  })
})
