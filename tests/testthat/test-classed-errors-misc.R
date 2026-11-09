# User-input errors in the misc verbs carry the class nestimate_bad_input.
testthat::skip_on_cran()

seqs <- data.frame(
  T1 = c("A", "B", "A", "C"), T2 = c("B", "A", "C", "A"),
  T3 = c("A", "C", "B", "B"), T4 = c("C", "B", "A", "A"),
  stringsAsFactors = FALSE
)
net <- build_network(seqs, method = "relative")
p2 <- matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2,
             dimnames = list(c("A", "B"), c("A", "B")))

test_that("entropy verbs class their argument checks", {
  expect_error(transition_entropy(p2, base = 1), class = "nestimate_bad_input")
  expect_error(entropy_network(p2, base = -1), class = "nestimate_bad_input")
  expect_error(entropy_network(p2, weight = "contribution",
                               scaling = "chance"),
               class = "nestimate_bad_input")
  expect_error(entropy_trajectory(1, action = "a"),
               class = "nestimate_bad_input")
  expect_error(entropy_bayes(p2 * 10, prior = -1),
               class = "nestimate_bad_input")
  expect_error(entropy_bayes(p2), class = "nestimate_bad_input")
})

test_that("stability and reliability verbs class their checks", {
  expect_error(centrality_stability(1), class = "nestimate_bad_input")
  expect_error(centrality_stability(net, iter = 1),
               class = "nestimate_bad_input")
  expect_error(casedrop_reliability(1), class = "nestimate_bad_input")
  expect_error(casedrop_reliability(net, iter = 1),
               class = "nestimate_bad_input")
  expect_error(network_reliability(1), class = "nestimate_bad_input")
  expect_error(network_reliability(net, iter = 1),
               class = "nestimate_bad_input")
  expect_error(certainty(net, prior = -1), class = "nestimate_bad_input")
})

test_that("Markov verbs class their checks", {
  expect_error(passage_time(1), class = "nestimate_bad_input")
  expect_error(passage_time(p2, states = "Z"), class = "nestimate_bad_input")
  expect_error(markov_stability(1), class = "nestimate_bad_input")
  expect_error(markov_order_test(seqs, max_order = 0),
               class = "nestimate_bad_input")
  expect_error(path_dependence(seqs, order = 1),
               class = "nestimate_bad_input")
})

test_that("higher-order builders class their checks", {
  expect_error(build_hon(1), class = "nestimate_bad_input")
  expect_error(build_hypa(seqs, p_adjust = "nope"),
               class = "nestimate_bad_input")
  expect_error(build_mogen(1), class = "nestimate_bad_input")
  expect_error(mogen_transitions(1), class = "nestimate_bad_input")
  expect_error(path_counts(seqs, k = 1), class = "nestimate_bad_input")
  expect_error(build_honem(p2, dim = 0), class = "nestimate_bad_input")
  expect_error(build_honem(1), class = "nestimate_bad_input")
})

test_that("hypergraph verbs class their checks", {
  expect_error(build_hypergraph(net, p = 2), class = "nestimate_bad_input")
  expect_error(build_hypergraph(1), class = "nestimate_bad_input")
  expect_error(bipartite_groups(1, "a", "b"), class = "nestimate_bad_input")
  expect_error(clique_expansion(1), class = "nestimate_bad_input")
  expect_error(hypergraph_centrality(1), class = "nestimate_bad_input")
  expect_error(hypergraph_measures(1), class = "nestimate_bad_input")
  expect_error(hypergraph_cluster(1, k = 2), class = "nestimate_bad_input")
  expect_error(hypergraph_transduction(1, labels = "a"),
               class = "nestimate_bad_input")
})

test_that("sequence verbs class their checks", {
  expect_error(sequence_compare(seqs, sub = 0), class = "nestimate_bad_input")
  expect_error(actor_endpoints(1), class = "nestimate_bad_input")
  expect_error(mark_first_state(seqs, state = 1),
               class = "nestimate_bad_input")
  expect_error(mark_terminal_state(seqs, state = ""),
               class = "nestimate_bad_input")
  expect_error(extract_sequences(seqs, meta = NA),
               class = "nestimate_bad_input")
  expect_error(sequence_patterns(seqs, top = -1),
               class = "nestimate_bad_input")
  expect_error(pattern_outcome(seqs, outcome = "T4", n = 0),
               class = "nestimate_bad_input")
  expect_error(sequence_metrics(seqs, omega = NA),
               class = "nestimate_bad_input")
  expect_error(extract_pathways(1, action = "a"),
               class = "nestimate_bad_input")
})

test_that("network accessors and models class their checks", {
  expect_error(extract_transition_matrix(1), class = "nestimate_bad_input")
  expect_error(extract_edges(net, threshold = NA),
               class = "nestimate_bad_input")
  expect_error(extract_edges(net, include_self = NA),
               class = "nestimate_bad_input")
  expect_error(association_rules(seqs, min_support = 2),
               class = "nestimate_bad_input")
  expect_error(association_rules(1), class = "nestimate_bad_input")
  expect_error(predict_links(1), class = "nestimate_bad_input")
  expect_error(evaluate_links(1, data.frame(from = "A", to = "B")),
               class = "nestimate_bad_input")
  expect_error(net_edge_betweenness(1), class = "nestimate_bad_input")
  expect_error(rename_models(1, "a"), class = "nestimate_bad_input")
  expect_error(coefs(1), class = "nestimate_bad_input")
  expect_error(build_mlvar(1, vars = c("a", "b"), id = "id"),
               class = "nestimate_bad_input")
  expect_error(outcome_model(1, outcome = "y", predictors = "x"),
               class = "nestimate_bad_input")
  expect_error(effects_table(1), class = "nestimate_bad_input")
})

test_that("landscape verbs class their argument checks", {
  expect_error(stability_landscape(net, iter = -1),
               class = "nestimate_bad_input")
  expect_error(landscape_surface(net, k = 1), class = "nestimate_bad_input")
  expect_error(landscape_trend(net, window = 1),
               class = "nestimate_bad_input")
})

test_that("summary of a group classes its argument check", {
  grp <- build_network(cbind(seqs, g = c("x", "x", "y", "y")),
                       method = "relative", group = "g")
  expect_error(summary(grp, combined = "yes"), class = "nestimate_bad_input")
})

test_that("build_gimme classes its argument checks", {
  skip_if_not_installed("lavaan")
  expect_error(build_gimme(1, vars = c("a", "b"), id = "id"),
               class = "nestimate_bad_input")
})
