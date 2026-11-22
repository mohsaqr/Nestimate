# Part 2 of test-layer9-identity-resampling.R, split so that one Prova case holds at most 7 references;
# the shared set-up is repeated.

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

test_that("compare_networks() compares the input weights and the centralities net_centrality() returns", {
  first <- tn_simulate(n_actors = 40L, n_states = 4L, seq_length = 10L, seed = 936L)
  second <- tn_simulate(n_actors = 40L, n_states = 4L, seq_length = 10L, seed = 937L)
  continuous_a <- psy_continuous(n = 300L, p = 5L, rho = 0.9, seed = 938L)
  continuous_b <- psy_continuous(n = 300L, p = 5L, rho = 0.9, seed = 939L)
  pairs <- list(
    relative = list(build_network(first, method = "relative"), build_network(second, method = "relative")),
    co_occurrence = list(build_network(first, method = "co_occurrence"), build_network(second, method = "co_occurrence")),
    pcor = list(build_network(continuous_a, method = "pcor"), build_network(continuous_b, method = "pcor"))
  )
  measures <- c("InStrength", "OutStrength", "Betweenness")
  Map(\(nets, label) {
    comparison <- compare_networks(a = nets[[1L]], b = nets[[2L]], measures = measures)
    # the table verbs round to `digits` decimals for reading; the stored tables hold the full values
    edges <- comparison$edges
    # every ordered cell of a directed pair; the upper triangle with the diagonal of an undirected one
    reference_edges <- function(net) {
      nodes <- rownames(net$weights)
      cells <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
      cells <- if (isTRUE(net$directed)) cells else subset(cells, match(from, nodes) <= match(to, nodes))
      transform(cells, weight = net$weights[cbind(from, to)])
    }
    expect_reference(subset(edges, select = c(from, to, weight_a)),
                     setNames(reference_edges(nets[[1L]]), c("from", "to", "weight_a")),
                     reference = li_ref("the first input network's weights"), covers = "compare_networks()",
                     quantity = sprintf("weight_a, %s", label), key = c("from", "to"), tolerance = 0, authority = "internal")
    expect_reference(subset(edges, select = c(from, to, weight_b)),
                     setNames(reference_edges(nets[[2L]]), c("from", "to", "weight_b")),
                     reference = li_ref("the second input network's weights"), covers = "compare_networks()",
                     quantity = sprintf("weight_b, %s", label), key = c("from", "to"), tolerance = 0, authority = "internal")
    nodes <- comparison$nodes
    central <- withCallingHandlers(net_centrality(nets[[1L]], measures = measures),
                                   message = \(condition) invokeRestart("muffleMessage"))
    expected <- do.call(rbind, lapply(measures, \(measure) {
      data.frame(node = as.character(central$state), measure = measure, value_a = central[[measure]])
    }))
    expect_reference(subset(nodes, select = c(node, measure, value_a)), expected,
                     reference = li_ref("net_centrality() of the first network"), covers = c("compare_networks()", "net_centrality()"),
                     quantity = sprintf("value_a, %s", label), key = c("node", "measure"), tolerance = 1e-12,
                     authority = "internal",
                     note = "Two code paths sum the same weights; 1e-12 allows for a different order of summation.")
  }, pairs, names(pairs))
})

test_that("macro_network() collapsed equals the macro layer and expanded equals the direct fit", {
  ragged <- tn_simulate(n_actors = 40L, n_states = 6L, seq_length = 12L, seed = 940L, ragged = TRUE)
  fit <- build_mcml(ragged, clusters = list(g1 = c("A", "B"), g2 = c("C", "D"), g3 = c("E", "F")))
  collapsed <- macro_network(fit)
  expect_reference(collapsed$weights, fit$macro$weights, reference = li_ref("build_mcml() macro layer"),
                   covers = c("macro_network()", "build_mcml()"), quantity = "collapsed weights", tolerance = 0,
                   authority = "internal")
  expect_reference(collapsed$initial, fit$macro$inits, reference = li_ref("build_mcml() macro layer"),
                   covers = c("macro_network()", "build_mcml()"), quantity = "collapsed initial probabilities", tolerance = 0,
                   authority = "internal")
  expansions <- list(
    list(quantity = "relative", args = list(method = "relative")),
    list(quantity = "frequency with start and threshold 2", args = list(method = "frequency", start = TRUE, threshold = 2)),
    list(quantity = "relative with end and minmax scaling", args = list(method = "relative", end = TRUE, scaling = "minmax"))
  )
  lapply(expansions, \(expansion) {
    expanded <- do.call(macro_network, c(list(fit, expand = "all"), expansion$args))
    li_expect_same_network(expanded, do.call(build_network, c(list(ragged), expansion$args)),
                           route = "build_network() on the sequences", covers = "macro_network()",
                           quantity = sprintf("every cluster expanded, %s", expansion$quantity),
                           parts = c("weights", "edges", "initial"))
  })
})

test_that("the layers of as_tna(<mcml>) carry the fit's weights", {
  fit <- build_mcml(wide, clusters = list(first = c("A", "B"), second = c("C", "D")))
  layers <- as_tna(fit)
  expect_reference(layers$macro$weights, fit$macro$weights, reference = li_ref("build_mcml() macro layer"),
                   covers = "build_mcml()", quantity = "as_tna() macro layer", tolerance = 0, authority = "internal")
  Map(\(layer, name) {
    expect_reference(layer$weights, fit$clusters[[name]]$weights, reference = li_ref("build_mcml() cluster layer"),
                     covers = "build_mcml()", quantity = sprintf("as_tna() layer of cluster %s", name), tolerance = 0,
                     authority = "internal")
  }, layers[c("first", "second")], c("first", "second"))
})

test_that("a tna model rebuilt by bootstrap_network() equals the direct fit", {
  skip_if_not_installed("tna")
  model <- tna::tna(wide)
  direct <- build_network(wide, method = "relative")
  li_expect_same_network(bootstrap_network(model, iter = 10L, seed = 1L)$original, direct,
                         route = "build_network() on the sequences", covers = "bootstrap_network()",
                         quantity = "bootstrap_network(<tna>) observed network", parts = c("weights", "edges", "initial"))
})
