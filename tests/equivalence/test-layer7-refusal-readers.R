# Layer 7, refusal parity of the sequence readers, the chain readers and the network readers. A sequence reader takes
# the states stored in a network's $data, so it refuses a network whose data are not positions of states: one-hot
# codes, association variables, an mcml object, a netobject_group where one set of sequences is meant (CLAUDE.md,
# "Column roles are declared"), and a frame whose columns have no declared role. A chain reader refuses signed or
# association weights with nestimate_not_markov. Helpers: helper-layer-refusal.R.

bi <- "nestimate_bad_input"
nm <- "nestimate_not_markov"

event_log <- rf_event_log()
edge_fit <- build_mcml(rf_edges9(), clusters = rf_clusters9, from = "from", to = "to", weight = "w")
inputs <- list(
  latency = build_network(event_log, method = "latency", actor = "chat", action = "sender", time = "at"),
  sum_mcml = build_mcml(rf_wide9(), clusters = rf_clusters9),
  edge_mcml = edge_fit,
  edge_layer = as_tna(edge_fit)$macro,
  onehot_no_unit = build_network(rf_onehot(), method = "cna", codes = c("A", "B", "C"), window_size = 2),
  cor = build_cor(rf_numeric()),
  glasso = build_glasso(rf_numeric()),
  ising = build_ising(as.data.frame(lapply(rf_numeric(), \(v) as.integer(v > 0)))),
  mgm = build_network(rf_numeric(), method = "mgm"),
  bare = rf_bare(),
  group = build_network(rf_wide(group = TRUE), method = "relative", group = "g"),
  mlvar = build_mlvar(rf_panel(), vars = c("V1", "V2", "V3"), id = "id", day = "day", beep = "beep"),
  relative = build_network(rf_wide(), method = "relative")
)

verbs <- list(
  build_hypa = \(x) build_hypa(x, order = 2L),
  build_hon = \(x) build_hon(x, max_order = 2L),
  build_mogen = \(x) build_mogen(x, max_order = 2L),
  markov_order_test = \(x) markov_order_test(x, max_order = 1L, n_perm = 5L),
  sequence_patterns = \(x) sequence_patterns(x),
  sequence_metrics = \(x) sequence_metrics(x),
  extract_sequences = \(x) extract_sequences(x),
  build_mmm = \(x) build_mmm(x, k = 2L, n_starts = 1L, max_iter = 5L, seed = 1),
  build_clusters = \(x) build_clusters(x, k = 2L),
  path_counts = \(x) path_counts(x),
  path_dependence = \(x) path_dependence(x, min_count = 1L),
  state_frequencies = \(x) state_frequencies(x),
  stability_landscape = \(x) stability_landscape(x, iter = 5L, seed = 1),
  transition_entropy = \(x) transition_entropy(x),
  passage_time = \(x) passage_time(x),
  markov_stability = \(x) markov_stability(x),
  chain_structure = \(x) chain_structure(x),
  net_centrality = \(x) net_centrality(x),
  extract_edges = \(x) extract_edges(x),
  net_prune = \(x) net_prune(x),
  predictability = \(x) predictability(x),
  build_network = \(x) build_network(x, method = "relative"),
  sequence_compare = \(x) sequence_compare(x, iter = 5L),
  entropy_bayes = \(x) entropy_bayes(x, draws = 200, seed = 1),
  frequencies = \(x) frequencies(x),
  actor_endpoints = \(x) actor_endpoints(x),
  association_rules = \(x) association_rules(x),
  cooccurrence = \(x) cooccurrence(x),
  build_mcml = \(x) build_mcml(x, clusters = list(g1 = c("A", "B"), g2 = c("C", "D")))
)
covers <- stats::setNames(as.list(paste0(names(verbs), "()")), names(verbs))

sequence_verbs <- c("build_hypa", "build_hon", "build_mogen", "markov_order_test", "sequence_patterns",
                    "sequence_metrics", "extract_sequences", "build_mmm", "build_clusters", "path_counts",
                    "path_dependence", "state_frequencies")
chain_verbs <- c("transition_entropy", "passage_time", "markov_stability", "chain_structure")
latency_results <- c(build_hypa = "net_hypa", build_hon = "net_hon", build_mogen = "net_mogen",
                     markov_order_test = "net_markov_order", sequence_patterns = "data.frame",
                     sequence_metrics = "data.frame", extract_sequences = "data.frame", build_mmm = "net_mmm",
                     build_clusters = "net_clustering", path_counts = "data.frame",
                     path_dependence = "net_path_dependence", state_frequencies = "data.frame")
chain_results <- c(transition_entropy = "net_transition_entropy", passage_time = "net_mpt",
                   markov_stability = "net_markov_stability", chain_structure = "chain_structure")

# One cell table for many verbs on one input, as verb = expected.
by_verb <- function(input, expected, defect = NULL, note = NULL, allow_warning = NULL) {
  pick <- function(values, verb) if (is.null(values)) NULL else stats::setNames(values[verb], input)
  do.call(rbind, Map(\(verb, outcome) {
    rf_cells(verb, stats::setNames(outcome, input), defect = pick(defect, verb), note = pick(note, verb),
             allow_warning = pick(allow_warning, verb))
  }, names(expected), expected))
}
same <- function(verbs, value) stats::setNames(rep(value, length(verbs)), verbs)

test_that("sequence readers refuse one-hot, association and grouped networks, an mcml and a bare network alike", {
  single_verbs <- setdiff(sequence_verbs, c("sequence_patterns", "sequence_metrics", "extract_sequences",
                                            "markov_order_test"))
  cells <- rbind(
    by_verb("onehot_no_unit", same(sequence_verbs, bi)),
    by_verb("cor", same(sequence_verbs, bi)),
    by_verb("glasso", same(sequence_verbs, bi)),
    by_verb("ising", same(sequence_verbs, bi)),
    by_verb("mgm", same(sequence_verbs, bi)),
    by_verb("bare", same(sequence_verbs, bi)),
    by_verb("mlvar", same(sequence_verbs, bi)),
    by_verb("group", same(single_verbs, bi)),
    by_verb("sum_mcml", same(sequence_verbs, bi)),
    by_verb("edge_mcml", same(c("markov_order_test", "build_mmm", "sequence_patterns", "build_hon"), bi))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("an edge-list mcml layer is not read as sequences", {
  # the layer's $data holds from, to and weight, tagged attr(, "source") = "edgelist"; a reader that takes each row
  # for a sequence counts the weight as a third state
  # .coerce_sequence_input() (and .as_sequences() through it), markov_order_test(), build_mmm(), build_clusters() and
  # path_dependence() refuse it
  readers <- sequence_verbs
  cells <- rbind(
    by_verb("edge_layer", same(readers, bi)),
    rf_cells("stability_landscape", c(edge_layer = bi))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("a latency network keeps its senders' sequences, which every sequence reader reads", {
  cells <- rbind(
    by_verb("latency", latency_results,
            defect = c(build_clusters = "refuses (nestimate_bad_input) the sender sequences that build_mmm() and the other sequence readers read: .extract_sequence_data() admits relative/frequency/co_occurrence/attention only (R/cluster_data.R:505-512)"),
            note = same(names(latency_results), "a latency network's $data are the senders of each chat in order, a wide sequence table")),
    by_verb("latency", chain_results, allow_warning = same(names(chain_results), "nestimate_rows_normalized"),
            note = same(names(chain_results), "a chain reader row-normalises non-negative weights and says so with nestimate_rows_normalized"))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("chain readers refuse association weights with nestimate_not_markov", {
  cells <- rbind(
    by_verb("cor", same(chain_verbs, nm)), by_verb("glasso", same(chain_verbs, nm)),
    by_verb("ising", same(chain_verbs, nm)), by_verb("mgm", same(chain_verbs, nm)),
    by_verb("mlvar", same(chain_verbs, nm))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("network readers take any network; predictability() refuses what it cannot read with a class", {
  cells <- rbind(
    by_verb("bare", c(net_centrality = "net_centrality", extract_edges = "data.frame", net_prune = "netobject",
                      predictability = bi),
            allow_warning = c(net_centrality = "nestimate_constant_measure")),
    by_verb("sum_mcml", c(net_centrality = "net_centrality_group", predictability = bi)),
    by_verb("mlvar", c(net_centrality = "net_centrality_group", net_prune = "net_mlvar", predictability = bi),
            allow_warning = c(net_centrality = "nestimate_constant_measure")),
    by_verb("latency", c(net_centrality = "net_centrality", extract_edges = "data.frame", net_prune = "netobject",
                         predictability = bi)),
    by_verb("glasso", c(predictability = "data.frame"),
            note = c(predictability = "a glasso network carries the precision matrix predictability() reads"))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("tna models: predictability() and build_mmm() refuse them with a class, like their siblings", {
  skip_if_not_installed("tna")
  wide <- rf_wide(group = TRUE)
  tna_inputs <- list(tna = tna::tna(rf_wide()),
                     group_tna = tna::group_tna(wide, group = "g", cols = paste0("T", seq_len(10L))))
  cells <- rbind(
    by_verb("tna", c(predictability = bi, build_mmm = "net_mmm", sequence_metrics = "data.frame",
                     markov_order_test = "net_markov_order")),
    by_verb("group_tna", c(build_mmm = bi, build_clusters = bi, sequence_metrics = "data.frame"))
  )
  rf_run(cells, tna_inputs, verbs, covers)
})

test_that("frames without declared roles, with 0/1 columns or repeated names are refused by every reader", {
  dup <- rf_wide(length = 4L)
  names(dup) <- c("T1", "T2", "T2", "T3")
  base <- rf_wide()
  frames <- list(
    undeclared = data.frame(student = paste0("s", 1:30), first = base$T1, second = base$T2, third = base$T3),
    action_column = data.frame(action = base$T1, T1 = base$T2),
    binary = as.data.frame(lapply(rf_numeric(), \(v) as.integer(v > 0))),
    duplicated_names = dup
  )
  frame_verbs <- c("build_network", "build_hypa", "build_hon", "build_mogen", "markov_order_test",
                   "sequence_patterns", "sequence_metrics", "extract_sequences", "build_mmm", "build_clusters",
                   "stability_landscape", "path_counts", "path_dependence", "state_frequencies", "sequence_compare",
                   "transition_entropy", "entropy_bayes", "passage_time", "frequencies", "chain_structure",
                   "actor_endpoints", "association_rules", "cooccurrence", "build_mcml")
  cells <- do.call(rbind, lapply(names(frames), \(frame) by_verb(frame, same(frame_verbs, bi))))
  rf_run(cells, frames, verbs, covers)
})

test_that("a frame of one sequence is refused by the resampling landscape and read by the descriptive verbs", {
  one_row <- list(single_row = data.frame(T1 = "A", T2 = "B", T3 = "C", T4 = "A"))
  cells <- by_verb("single_row",
                   c(stability_landscape = "nestimate_single_sequence", sequence_metrics = "data.frame",
                     build_hon = "net_hon", markov_order_test = "net_markov_order"),
                   note = c(stability_landscape = "a bootstrap of one sequence would give intervals of zero width",
                            sequence_metrics = "describes the one sequence; nothing is resampled",
                            build_hon = "counts the paths of the one sequence; nothing is resampled",
                            markov_order_test = "the within-context permutation is defined on the windows of one sequence"))
  rf_run(cells, one_row, verbs, covers)
})
