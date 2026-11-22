# frequencies(), convert_sequence_format(), state_frequencies(), state_distribution() and summary() of a network:
# counts of transitions and of states, and the network-level metrics of summary.netobject(), against their
# definitions in base R and against igraph. Re-expressed from
# local_testing_and_equivalence/test-frequencies-equiv-extracted.R and test-estimate_network-equiv-extracted.R
# (which only checked shapes and that estimate_network() equals frequencies()), with the numbers now compared.

tolerance <- 1e-10
wide <- tn_simulate(n_actors = 30L, n_states = 5L, seq_length = 12L, seed = 301L, ragged = TRUE)
sequences <- tn_sequences(wide)
states <- sort(unique(unlist(sequences)))

test_that("frequencies() and its summary() count consecutive pairs", {
  counts <- tn_pair_counts(sequences)
  found <- frequencies(wide, format = "wide")
  expect_reference(found, counts, reference = "independent: table of consecutive pairs, base R",
                   covers = "frequencies()", quantity = "transition counts, wide", tolerance = 0)
  edges <- data.frame(from = rep(states, times = length(states)), to = rep(states, each = length(states)),
                      count = as.vector(counts))
  expect_reference(summary(found), subset(edges, count > 0), reference = "independent: table of consecutive pairs, base R",
                   covers = "summary.nest_transition_counts()", quantity = "non-zero transition counts", tolerance = 0,
                   key = c("from", "to"))
  expect_reference(frequencies(tn_long(wide), action = "action", id = "actor", time = "time", format = "long"), counts,
                   reference = "independent: table of consecutive pairs, base R",
                   covers = "frequencies()", quantity = "transition counts, long", tolerance = 0)
  expect_reference(extract_transition_matrix(withCallingHandlers(estimate_network(wide, method = "frequency"),
                     deprecatedWarning = \(w) invokeRestart("muffleWarning"))), counts,
                   reference = "independent: table of consecutive pairs, base R",
                   covers = "estimate_network()", quantity = "frequency network of the deprecated wrapper", tolerance = 0)
})

test_that("convert_sequence_format() counts, indicates and lists transitions as their definitions do", {
  identified <- cbind(id = seq_len(nrow(wide)), wide)
  per_row <- t(vapply(sequences, \(s) as.integer(table(factor(s, levels = states))), integer(length(states))))
  colnames(per_row) <- states
  expected_counts <- cbind(data.frame(id = seq_len(nrow(wide))), as.data.frame(per_row))
  expect_reference(convert_sequence_format(identified, id_col = "id", format = "frequency"), expected_counts,
                   reference = "independent: table of each sequence's states, base R",
                   covers = "convert_sequence_format()", quantity = "state counts per sequence", tolerance = 0, key = "id")
  expected_onehot <- cbind(data.frame(id = seq_len(nrow(wide))), as.data.frame(1L * (per_row > 0L)))
  expect_reference(convert_sequence_format(identified, id_col = "id", format = "onehot"), expected_onehot,
                   reference = "independent: presence of each state in each sequence, base R",
                   covers = "convert_sequence_format()", quantity = "state presence per sequence", tolerance = 0, key = "id")
  edgelist <- convert_sequence_format(identified, id_col = "id", format = "edgelist")
  expect_reference(unclass(table(factor(edgelist$from, levels = states), factor(edgelist$to, levels = states), dnn = NULL)),
                   tn_pair_counts(sequences), reference = "independent: table of consecutive pairs, base R",
                   covers = "convert_sequence_format()", quantity = "pairs of the edge list format", tolerance = 0)
})

test_that("state_frequencies() and state_distribution() count every event's state", {
  counts <- table(factor(unlist(sequences), levels = states), dnn = NULL)
  expected <- data.frame(state = states, count = as.integer(counts), proportion = as.numeric(counts) / sum(counts))
  # state_frequencies() rounds the proportion to 4 decimals, as its documentation states: the reference does too.
  expect_reference(state_frequencies(wide), transform(expected, proportion = round(proportion, 4)),
                   reference = "independent: table of all events, base R",
                   covers = "state_frequencies()", quantity = "state frequencies", tolerance = 1e-12, key = "state")
  expect_reference(state_distribution(build_tna(wide)), expected, reference = "independent: table of all events, base R",
                   covers = c("state_distribution()", "state_distribution.netobject()"),
                   quantity = "state distribution of a network", tolerance = tolerance, key = "state")
})

test_that("summary() of a transition network equals igraph's metrics of the weighted graph", {
  skip_without_reference("igraph")
  weights <- tn_row_normalise(tn_pair_counts(sequences))
  # Self-loops count as edges, strengths and degrees; centralization is igraph's without loops.
  with_loops <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE, diag = TRUE)
  without_loops <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE, diag = FALSE)
  n <- length(states)
  out_strength <- igraph::strength(with_loops, mode = "out")
  in_strength <- igraph::strength(with_loops, mode = "in")
  out_degree <- igraph::degree(with_loops, mode = "out")
  expected <- data.frame(
    metric = c("Node Count", "Edge Count", "Mean Distance", "Mean Out-Strength", "SD Out-Strength",
               "Mean In-Strength", "SD In-Strength", "Mean Out-Degree", "SD Out-Degree", "Centralization (Out-Degree)",
               "Centralization (In-Degree)", "Reciprocity"),
    value = c(n, igraph::ecount(with_loops), igraph::mean_distance(without_loops), mean(out_strength), stats::sd(out_strength), mean(in_strength),
              stats::sd(in_strength), mean(out_degree), stats::sd(out_degree),
              igraph::centr_degree(without_loops, mode = "out", loops = FALSE)$centralization,
              igraph::centr_degree(without_loops, mode = "in", loops = FALSE)$centralization,
              igraph::reciprocity(with_loops)))
  # Density is compared on its own below (a known difference on this input).
  expect_reference(subset(summary(build_tna(wide)), metric != "Network Density"), expected, reference = "igraph",
                   covers = "summary.netobject()",
                   quantity = "network metrics", tolerance = tolerance, key = "metric",
                   note = "Mean Distance is igraph's mean_distance() with the probabilities as edge lengths, as Nestimate computes it.")
  expect_known_difference(subset(summary(build_tna(wide)), metric == "Network Density"),
                          data.frame(metric = "Network Density", value = igraph::edge_density(with_loops, loops = FALSE)),
                          reference = "igraph", covers = "summary.netobject()", quantity = "network density, complete graph with loops",
                          tolerance = tolerance, key = "metric",
                          reason = "Nestimate counts self-loops as edges over n(n - 1) possible ones, as igraph::edge_density(loops = FALSE) does, but caps the density at 1; this complete graph with loops gives igraph 1.25.")
  grouped <- cbind(group = rep(c("a", "b"), length.out = nrow(wide)), wide)
  summary_group <- summary(build_network(grouped, method = "relative", group = "group"))
  first <- tn_row_normalise(tn_pair_counts(tn_sequences(subset(wide, grouped$group == "a")), states))
  graph_a <- igraph::graph_from_adjacency_matrix(first, mode = "directed", weighted = TRUE, diag = TRUE)
  expect_reference(subset(summary_group, metric %in% c("Node Count", "Edge Count", "Reciprocity")),
                   data.frame(metric = c("Node Count", "Edge Count", "Reciprocity"),
                              a = c(n, igraph::ecount(graph_a), igraph::reciprocity(graph_a))),
                   reference = "igraph", covers = "summary.netobject_group()", quantity = "group a network metrics",
                   tolerance = tolerance, key = "metric")
})

test_that("state distributions of a group, an mcml, an htna and a plotted table count every event's state", {
  counts <- table(factor(unlist(sequences), levels = states), dnn = NULL)
  overall <- data.frame(state = states, count = as.integer(counts), proportion = as.numeric(counts) / sum(counts))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot_state_frequencies(build_tna(wide))
  expect_reference(as.data.frame(plotted), overall, reference = "independent: table of all events, base R",
                   covers = "as.data.frame.state_freq()", quantity = "table of a state frequency plot",
                   tolerance = tolerance, key = "state")
  expect_reference(state_distribution(plotted), overall, reference = "independent: table of all events, base R",
                   covers = "state_distribution.state_freq()", quantity = "state distribution of a state frequency plot",
                   tolerance = tolerance, key = "state")
  # Per group: the events of the rows in that group.
  group <- rep(c("a", "b"), length.out = nrow(wide))
  per_group <- do.call(rbind, lapply(c("a", "b"), \(g) {
    events <- table(factor(unlist(tn_sequences(subset(wide, group == g))), levels = states), dnn = NULL)
    data.frame(group = g, state = states, count = as.integer(events), proportion = as.numeric(events) / sum(events))
  }))
  expect_reference(state_distribution(build_network(cbind(group = group, wide), method = "relative", group = "group")),
                   subset(per_group, count > 0), reference = "independent: table of all events per group, base R",
                   covers = "state_distribution.netobject_group()", quantity = "state distribution per group",
                   tolerance = tolerance, key = c("group", "state"))
  # Per cluster: a state's share among the events of its cluster's states.
  clusters <- list(c1 = states[1:2], c2 = states[3:5])
  per_cluster <- do.call(rbind, lapply(names(clusters), \(cluster) {
    events <- as.integer(counts[clusters[[cluster]]])
    data.frame(group = cluster, state = clusters[[cluster]], count = events, proportion = events / sum(events))
  }))
  fit <- build_mcml(wide, clusters = clusters)
  expect_reference(state_distribution(fit), per_cluster, reference = "independent: table of all events per cluster, base R",
                   covers = "state_distribution.mcml()", quantity = "state distribution per cluster of an mcml",
                   tolerance = tolerance, key = c("group", "state"))
  expect_reference(state_distribution(as_htna(fit)), per_cluster,
                   reference = "independent: table of all events per cluster, base R",
                   covers = c("state_distribution.htna()", "as_htna.mcml()"), quantity = "state distribution per cluster of an htna",
                   tolerance = tolerance, key = c("group", "state"))
})

test_that("as_htna() weights are the row-normalised pairs of the node-level sequences", {
  clusters <- list(c1 = states[1:2], c2 = states[3:5])
  expected <- tn_row_normalise(tn_pair_counts(sequences))
  expect_reference(extract_transition_matrix(as_htna(build_mcml(wide, clusters = clusters))), expected,
                   reference = "independent: table of consecutive pairs, base R", covers = c("as_htna()", "as_htna.mcml()"),
                   quantity = "htna from an mcml", tolerance = tolerance)
  expect_reference(extract_transition_matrix(as_htna(wide, clusters = clusters)), expected,
                   reference = "independent: table of consecutive pairs, base R",
                   covers = c("as_htna()", "as_htna.default()"), quantity = "htna from wide sequences", tolerance = tolerance)
})

test_that("as_htna() of a mixed Markov model fit holds each cluster's transitions and the actor partition", {
  partition <- list(X = states[1:2], Y = states[3:5])
  # The fit's assignment is an input here: the EM itself is compared in test-build_mmm.R. Run to a tight tolerance
  # with smooth = 0 so that the components are the M-step of the final posterior.
  fit <- cluster_mmm(as_htna(wide, clusters = partition), k = 2, seed = 1, smooth = 0, tol = 1e-12, max_iter = 2000L)
  per_sequence <- lapply(sequences, \(s) tn_pair_counts(list(s), states))
  # API gap: the posterior and the hard assignment are fields of the fit (`$posterior`, `$assignments`), with no
  # per-sequence accessor for a fit built from wide data (session_ids() needs long data).
  posterior <- fit$posterior
  assignments <- fit$assignments
  relative <- as_htna(fit)
  frequency <- as_htna(fit, method = "frequency")
  invisible(lapply(seq_len(2L), \(cluster) {
    weighted <- Reduce(`+`, Map(`*`, per_sequence, posterior[, cluster]))
    # The EM stops when the log-likelihood changes by less than 1e-12; the stored posterior is that of the last
    # E-step, so the components and the M-step of that posterior differ by at most the last step's change.
    expect_reference(extract_transition_matrix(relative[[cluster]]), tn_row_normalise(weighted),
                     reference = "independent: posterior-weighted transition counts (M-step), base R",
                     covers = "as_htna.net_mmm()", quantity = sprintf("cluster %d, relative", cluster), tolerance = 1e-6)
    expect_reference(extract_transition_matrix(frequency[[cluster]]),
                     tn_pair_counts(sequences[assignments == cluster], states),
                     reference = "independent: table of consecutive pairs of the assigned sequences, base R",
                     covers = "as_htna.net_mmm()", quantity = sprintf("cluster %d, frequency", cluster), tolerance = 0)
    expect_reference(subset(relative[[cluster]]$node_groups, select = c(node, group)),
                     data.frame(node = states, group = rep(names(partition), lengths(partition))),
                     reference = "independent: the partition given to as_htna()", covers = "as_htna.net_mmm()",
                     quantity = sprintf("cluster %d, actor of each node", cluster), tolerance = 0, key = "node")
  }))
})

test_that("state_distribution() of a tna model counts the events of the sequences it carries", {
  skip_without_reference("tna")
  counts <- table(factor(unlist(sequences), levels = states), dnn = NULL)
  expect_reference(state_distribution(tna::tna(wide)),
                   data.frame(state = states, count = as.integer(counts), proportion = as.numeric(counts) / sum(counts)),
                   reference = "independent: table of all events, base R", covers = "state_distribution.default()",
                   quantity = "state distribution of a tna model", tolerance = tolerance, key = "state")
})

# Outside references for the same counts: TraMineR (Gabadinho et al. 2011) for state frequencies, per-sequence
# state counts and transition counts; tna for the transition networks.
test_that("transition counts equal TraMineR::seqtrate(count = TRUE) and tna::ftna()", {
  skip_without_reference(c("TraMineR", "tna"))
  traminer <- unclass(TraMineR::seqtrate(out_seqdef(wide, states), count = TRUE))
  dimnames(traminer) <- list(states, states)
  frequency <- unclass(tna::ftna(wide)$weights)
  expect_reference(unclass(frequencies(wide, format = "wide")), traminer, reference = "TraMineR",
                   covers = "frequencies()", quantity = "transition counts, wide, TraMineR", tolerance = 0)
  expect_reference(unclass(frequencies(tn_long(wide), action = "action", id = "actor", time = "time", format = "long")),
                   frequency, reference = "tna", covers = "frequencies()", quantity = "transition counts, long, tna",
                   tolerance = 0)
  expect_reference(summary(frequencies(wide, format = "wide")),
                   stats::setNames(out_edges(frequency, nonzero = TRUE), c("from", "to", "count")), reference = "tna",
                   covers = "summary.nest_transition_counts()", quantity = "non-zero transition counts, tna",
                   tolerance = 0, key = c("from", "to"))
  expect_reference(extract_transition_matrix(withCallingHandlers(estimate_network(wide, method = "frequency"),
                     deprecatedWarning = \(w) invokeRestart("muffleWarning"))), frequency, reference = "tna",
                   covers = "estimate_network()", quantity = "frequency network of the deprecated wrapper, tna",
                   tolerance = 0)
  expect_reference(extract_transition_matrix(withCallingHandlers(estimate_network(wide, method = "relative"),
                     deprecatedWarning = \(w) invokeRestart("muffleWarning"))), unclass(tna::tna(wide)$weights),
                   reference = "tna", covers = "estimate_network()",
                   quantity = "relative network of the deprecated wrapper, tna", tolerance = tolerance)
})

test_that("convert_sequence_format() counts and pairs equal TraMineR::seqistatd() and seqtrate()", {
  skip_without_reference("TraMineR")
  identified <- cbind(id = seq_len(nrow(wide)), wide)
  per_sequence <- unclass(TraMineR::seqistatd(out_seqdef(wide, states)))
  counts <- cbind(data.frame(id = seq_len(nrow(wide))), as.data.frame(matrix(as.integer(per_sequence), nrow(wide),
                                                                              dimnames = list(NULL, states))))
  expect_reference(convert_sequence_format(identified, id_col = "id", format = "frequency"), counts,
                   reference = "TraMineR", covers = "convert_sequence_format()",
                   quantity = "state counts per sequence, TraMineR", tolerance = 0, key = "id")
  presence <- cbind(data.frame(id = seq_len(nrow(wide))), as.data.frame(matrix(1L * (per_sequence > 0), nrow(wide),
                                                                                dimnames = list(NULL, states))))
  expect_reference(convert_sequence_format(identified, id_col = "id", format = "onehot"), presence,
                   reference = "TraMineR", covers = "convert_sequence_format()",
                   quantity = "state presence per sequence, TraMineR", tolerance = 0, key = "id")
  edgelist <- convert_sequence_format(identified, id_col = "id", format = "edgelist")
  traminer <- unclass(TraMineR::seqtrate(out_seqdef(wide, states), count = TRUE))
  dimnames(traminer) <- list(states, states)
  expect_reference(unclass(table(factor(edgelist$from, levels = states), factor(edgelist$to, levels = states), dnn = NULL)),
                   traminer, reference = "TraMineR", covers = "convert_sequence_format()",
                   quantity = "pairs of the edge list format, TraMineR", tolerance = 0)
})

test_that("state frequencies and distributions equal TraMineR::seqstatf()", {
  skip_without_reference(c("TraMineR", "tna"))
  overall <- out_state_table(wide, states)
  # state_frequencies() rounds the proportion to 4 decimals, as documented: the reference is rounded too.
  expect_reference(state_frequencies(wide), transform(overall, proportion = round(proportion, 4)), reference = "TraMineR",
                   covers = "state_frequencies()", quantity = "state frequencies, TraMineR", tolerance = 1e-12,
                   key = "state")
  expect_reference(state_distribution(build_tna(wide)), overall, reference = "TraMineR",
                   covers = c("state_distribution()", "state_distribution.netobject()"),
                   quantity = "state distribution of a network, TraMineR", tolerance = tolerance, key = "state")
  expect_reference(state_distribution(tna::tna(wide)), overall, reference = "TraMineR",
                   covers = "state_distribution.default()", quantity = "state distribution of a tna model, TraMineR",
                   tolerance = tolerance, key = "state")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot_state_frequencies(build_tna(wide))
  expect_reference(as.data.frame(plotted), overall, reference = "TraMineR", covers = "as.data.frame.state_freq()",
                   quantity = "table of a state frequency plot, TraMineR", tolerance = tolerance, key = "state")
  expect_reference(state_distribution(plotted), overall, reference = "TraMineR",
                   covers = "state_distribution.state_freq()",
                   quantity = "state distribution of a state frequency plot, TraMineR", tolerance = tolerance,
                   key = "state")
  group <- rep(c("a", "b"), length.out = nrow(wide))
  per_group <- do.call(rbind, lapply(c("a", "b"), \(g) cbind(group = g, out_state_table(subset(wide, group == g), states))))
  expect_reference(state_distribution(build_network(cbind(group = group, wide), method = "relative", group = "group")),
                   subset(per_group, count > 0), reference = "TraMineR", covers = "state_distribution.netobject_group()",
                   quantity = "state distribution per group, TraMineR", tolerance = tolerance, key = c("group", "state"))
  # Per cluster: TraMineR's counts, each state's share among the events of its cluster's states.
  clusters <- list(c1 = states[1:2], c2 = states[3:5])
  per_cluster <- do.call(rbind, lapply(names(clusters), \(cluster) {
    rows <- subset(overall, state %in% clusters[[cluster]])
    data.frame(group = cluster, state = rows$state, count = rows$count, proportion = rows$count / sum(rows$count))
  }))
  fit <- build_mcml(wide, clusters = clusters)
  expect_reference(state_distribution(fit), per_cluster, reference = "TraMineR", covers = "state_distribution.mcml()",
                   quantity = "state distribution per cluster of an mcml, TraMineR counts", tolerance = tolerance,
                   key = c("group", "state"))
  expect_reference(state_distribution(as_htna(fit)), per_cluster, reference = "TraMineR",
                   covers = "state_distribution.htna()",
                   quantity = "state distribution per cluster of an htna, TraMineR counts", tolerance = tolerance,
                   key = c("group", "state"))
})

test_that("as_htna() networks equal tna::tna() and tna::ftna() of the same sequences", {
  skip_without_reference("tna")
  clusters <- list(c1 = states[1:2], c2 = states[3:5])
  expected <- unclass(tna::tna(wide)$weights)
  expect_reference(extract_transition_matrix(as_htna(build_mcml(wide, clusters = clusters))), expected,
                   reference = "tna", covers = c("as_htna()", "as_htna.mcml()"), quantity = "htna from an mcml, tna",
                   tolerance = tolerance)
  expect_reference(extract_transition_matrix(as_htna(wide, clusters = clusters)), expected, reference = "tna",
                   covers = c("as_htna()", "as_htna.default()"), quantity = "htna from wide sequences, tna",
                   tolerance = tolerance)
  # The frequency networks of a mixed Markov fit are tna's count networks of the sequences assigned to each cluster
  # (the assignment is an input here; the EM is compared in test-build_mmm.R).
  fit <- cluster_mmm(as_htna(wide, clusters = clusters), k = 2, seed = 1, smooth = 0, tol = 1e-12, max_iter = 2000L)
  frequency <- as_htna(fit, method = "frequency")
  invisible(lapply(seq_len(2L), \(cluster) {
    assigned <- subset(wide, fit$assignments == cluster)
    expect_reference(out_edges(extract_transition_matrix(frequency[[cluster]]), nonzero = TRUE),
                     out_edges(tna::ftna(assigned)$weights, nonzero = TRUE), reference = "tna",
                     covers = "as_htna.net_mmm()", quantity = sprintf("cluster %d, frequency, tna", cluster),
                     tolerance = 0, key = c("from", "to"))
  }))
})
