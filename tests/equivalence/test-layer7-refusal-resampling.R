# Layer 7, refusal parity of the resampling and posterior verbs: one degenerate or unsupported input through every
# verb that resamples (or draws a posterior from) a network. A latency network is refused by every resampling verb
# but bootstrap_network() with nestimate_latency_resampling, a network resting on one unit (one sequence, or one-hot rows with neither actor
# nor session) with nestimate_single_sequence (CLAUDE.md, "Estimator Registry" and "One-hot data"), and an mcml layer
# whose transitions were combined by a mean or a density is read as counts by no verb but bootstrap_network(), which
# recombines each resample the same way. Helpers: helper-layer-refusal.R.

lr <- "nestimate_latency_resampling"
legit_latency <- paste("bootstrap_network() resamples a latency network by its sequences and tests stability on the",
                       "speed score (validated by simulation against known reply speeds)")
ss <- "nestimate_single_sequence"
bi <- "nestimate_bad_input"
nm <- "nestimate_not_markov"

event_log <- rf_event_log()
mean_fit <- build_mcml(rf_wide9(), clusters = rf_clusters9, method = "mean", type = "frequency")
mean_layers <- as_tna(mean_fit)
uneven <- list(A = c("a1", "a2"), B = c("a3", "b1", "b2", "b3"), C = c("c1", "c2", "c3"))
set.seed(2)
inputs <- list(
  latency = build_network(event_log, method = "latency", actor = "chat", action = "sender", time = "at"),
  latency_group = build_network(event_log, method = "latency", actor = "chat", action = "sender", time = "at",
                                group = "team"),
  single_sequence = rf_quiet_build(build_tna(data.frame(T1 = "A", T2 = "B", T3 = "C", T4 = "A", T5 = "B", T6 = "C",
                                                        T7 = "A", T8 = "A"))),
  onehot_no_unit = build_network(rf_onehot(), method = "cna", codes = c("A", "B", "C"), window_size = 2),
  onehot_frequency = build_network(rf_onehot(), method = "frequency", codes = c("A", "B", "C")),
  onehot_wtna = wtna(rf_onehot()),
  mean_mcml = mean_fit,
  mean_group = mean_layers,
  mean_layer = mean_layers$A,
  density_macro = as_tna(build_mcml(rf_wide9(), clusters = uneven, method = "density"))$macro,
  edge_mcml = build_mcml(rf_edges9(), clusters = rf_clusters9, from = "from", to = "to", weight = "w"),
  bare = rf_bare(),
  mlvar_temporal = build_mlvar(rf_panel(), vars = c("V1", "V2", "V3"), id = "id", day = "day", beep = "beep")$temporal,
  cor = build_cor(rf_numeric()),
  glasso = build_glasso(rf_numeric()),
  mgm = build_network(rf_numeric(), method = "mgm"),
  group = build_network(rf_wide(group = TRUE), method = "relative", group = "g"),
  one_state = build_network(rf_wide(n = 10L, length = 6L, states = "A"), method = "relative"),
  no_transition = build_network(data.frame(T1 = c("A", "B", "C", "A"), T2 = NA_character_), method = "relative")
)
empty <- rf_wide(group = TRUE)
empty[empty$g == "y", paste0("T", seq_len(10L))] <- NA_character_
inputs$empty_group <- build_network(empty, method = "relative", group = "g")

verbs <- list(
  bootstrap_network = \(x) bootstrap_network(x, iter = 5L, seed = 1),
  permutation = \(x) permutation(x, x, iter = 5L, seed = 1),
  centrality_stability = \(x) centrality_stability(x, iter = 5L, seed = 1),
  casedrop_reliability = \(x) casedrop_reliability(x, iter = 5L, seed = 1),
  network_reliability = \(x) network_reliability(x, iter = 5L, seed = 1),
  compare_networks_permutation = \(x) compare_networks(x, x, test = "permutation", iter = 5L),
  compare_networks_bootstrap = \(x) compare_networks(x, x, test = "bootstrap", iter = 5L),
  permutation_diagnostics = \(x) permutation_diagnostics(x, x, actor = "id", iter = 5L, seed = 1),
  vertex_bootstrap = \(x) vertex_bootstrap(x, iter = 5L, seed = 1),
  vertex_compare = \(x) vertex_compare(x, x, iter = 5L, seed = 1),
  stability_landscape = \(x) stability_landscape(x, iter = 5L, seed = 1),
  landscape_trend = \(x) landscape_trend(x, iter = 5L, seed = 1, window = 3L),
  certainty = \(x) certainty(x),
  bayes_compare = \(x) bayes_compare(x, x, draws = 200L, seed = 1),
  entropy_bayes = \(x) entropy_bayes(x, draws = 200, seed = 1),
  boot_glasso = \(x) boot_glasso(x, iter = 5L, cs_iter = 2L, seed = 1),
  mosaic_test = \(x) mosaic_test(x)
)
covers <- list(
  bootstrap_network = "bootstrap_network()", permutation = "permutation()",
  centrality_stability = "centrality_stability()", casedrop_reliability = "casedrop_reliability()",
  network_reliability = "network_reliability()", compare_networks_permutation = "compare_networks()",
  compare_networks_bootstrap = "compare_networks()", permutation_diagnostics = "permutation_diagnostics()",
  vertex_bootstrap = "vertex_bootstrap()", vertex_compare = "vertex_compare()",
  stability_landscape = "stability_landscape()", landscape_trend = "landscape_trend()", certainty = "certainty()",
  bayes_compare = "bayes_compare()", entropy_bayes = "entropy_bayes()", boot_glasso = "boot_glasso()",
  mosaic_test = "mosaic_test()"
)

# Why an input is legitimately supported, by verb family.
legit_aggregated <- "bootstrap_network() recombines each resample's rows by the fit's own aggregation (.bootstrap_aggregated())"
legit_weights <- "resamples the nodes of the weight matrix and re-estimates nothing, so no unit or count is assumed"
legit_posterior <- "an analytic Dirichlet posterior on the observed transitions: nothing is resampled, and one sequence still gives a posterior of non-zero width"
legit_edge <- "the edge list declares its row weights (params$row_weight), so each resample sums the weights of the rows it draws"
legit_tna <- "a tna model is rebuilt as a Nestimate network with its data (.tna_input())"
degenerate <- "every measure is constant on this network; the verb says so with nestimate_constant_measure and reports cs = 0"
one_state_similarity <- "one state: the similarities that need a spread (correlations, CV ratio, RV, rank agreement) are NA, the module's rule for an undefined similarity"
no_unit_resampled <- "silent result: resamples of a network with no observed transition are all zero, intervals of zero width (siblings: sequence readers, transition_entropy() and stability_landscape() refuse it)"

test_that("the five sequence-resampling verbs refuse the same inputs with the same class", {
  # a mean of per-transition ones is 1 for every observed pair, so the mean layer's edge weights are constant
  core <- c(latency = lr, latency_group = lr, single_sequence = ss, onehot_no_unit = ss, onehot_wtna = ss,
            mean_layer = bi, mean_mcml = bi, density_macro = bi, bare = bi, mlvar_temporal = bi)
  cells <- rbind(
    rf_cells("bootstrap_network", c(core[setdiff(names(core), c("latency", "latency_group"))],
                                    latency = "net_bootstrap", latency_group = "net_bootstrap_group",
                                    mean_layer = "net_bootstrap", mean_mcml = "net_bootstrap_group",
                                    density_macro = "net_bootstrap", edge_mcml = "net_bootstrap_group",
                                    one_state = "net_bootstrap", no_transition = bi),
             defect = c(no_transition = no_unit_resampled),
             note = c(latency = legit_latency, latency_group = legit_latency, mean_layer = legit_aggregated, mean_mcml = legit_aggregated, density_macro = legit_aggregated,
                      edge_mcml = legit_edge, one_state = "one state: every resample gives P(A -> A) = 1, which is exact"),
             allow_warning = c(edge_mcml = "nestimate_edgelist_resampling")),
    rf_cells("permutation", c(core, edge_mcml = "net_permutation_group", one_state = "net_permutation",
                              no_transition = bi),
             defect = c(no_transition = no_unit_resampled), note = c(edge_mcml = legit_edge)),
    rf_cells("centrality_stability", c(core, edge_mcml = "net_stability_group", one_state = "net_stability",
                                       no_transition = "net_stability"),
             note = c(edge_mcml = legit_edge, one_state = degenerate, no_transition = degenerate),
             allow_warning = c(one_state = "nestimate_constant_measure", no_transition = "nestimate_constant_measure",
                               mean_mcml = "nestimate_constant_measure", edge_mcml = "nestimate_constant_measure")),
    rf_cells("casedrop_reliability", c(core, edge_mcml = "net_casedrop_reliability_group",
                                       one_state = "net_casedrop_reliability", no_transition = "net_casedrop_reliability"),
             note = c(edge_mcml = legit_edge, one_state = degenerate, no_transition = degenerate,
                      mean_layer = "refused before the constant-edge answer, as by .precompute_per_sequence()"),
             allow_warning = c(one_state = "nestimate_constant_measure", no_transition = "nestimate_constant_measure",
                               edge_mcml = "nestimate_constant_measure")),
    rf_cells("network_reliability", c(core, edge_mcml = "net_reliability", one_state = "net_reliability",
                                      no_transition = bi),
             defect = c(no_transition = no_unit_resampled),
             note = c(edge_mcml = legit_edge, one_state = one_state_similarity))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("compare_networks() refuses what its backends refuse, and the vertex bootstrap agrees with itself", {
  cells <- rbind(
    rf_cells("compare_networks_permutation",
             c(latency = lr, latency_group = lr, single_sequence = ss, onehot_no_unit = ss, mean_layer = bi,
               density_macro = bi, bare = bi, mlvar_temporal = bi, one_state = "net_network_comparison",
               no_transition = bi, cor = "net_network_comparison", glasso = "net_network_comparison",
               mgm = "net_network_comparison"),
             defect = c(no_transition = no_unit_resampled),
             note = c(one_state = paste(one_state_similarity, "(permutation() answers it too)"),
                      cor = "a cor network keeps its data as a numeric matrix, which permutation() reads",
                      glasso = "a glasso network keeps its data as a numeric matrix, which permutation() reads",
                      mgm = "an mgm network keeps its data as a data.frame, so the permutation backend runs")),
    rf_cells("compare_networks_bootstrap",
             c(latency = lr, latency_group = lr, single_sequence = "net_network_comparison",
               onehot_no_unit = "net_network_comparison", mean_layer = "net_network_comparison",
               bare = "net_network_comparison", one_state = bi),
             note = c(one_state = "the vertex bootstrap needs at least 3 nodes, as vertex_bootstrap() says",
                      single_sequence = legit_weights, onehot_no_unit = legit_weights, mean_layer = legit_weights,
                      bare = legit_weights)),
    rf_cells("vertex_bootstrap",
             c(latency = lr, single_sequence = "net_vertex_bootstrap", onehot_no_unit = "net_vertex_bootstrap",
               mean_layer = "net_vertex_bootstrap", bare = "net_vertex_bootstrap",
               mlvar_temporal = "net_vertex_bootstrap", latency_group = bi, group = bi),
             note = c(latency = "latency networks are resampled only by bootstrap_network(), as compare_networks(test = \"bootstrap\") refuses",
                      single_sequence = legit_weights, onehot_no_unit = legit_weights, mean_layer = legit_weights,
                      bare = legit_weights, mlvar_temporal = legit_weights)),
    rf_cells("vertex_compare", c(latency = lr, single_sequence = "net_vertex_comparison"),
             note = c(single_sequence = legit_weights)),
    rf_cells("permutation_diagnostics",
             c(latency = lr, latency_group = lr, single_sequence = ss, onehot_no_unit = ss, mean_layer = bi,
               density_macro = bi, bare = bi))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("the landscape verbs refuse latency networks and single sequences like the other resampling verbs", {
  cells <- rbind(
    rf_cells("stability_landscape",
             c(latency = lr, latency_group = lr, single_sequence = ss, onehot_frequency = bi, onehot_no_unit = nm,
               mean_layer = bi, density_macro = bi, bare = bi, mlvar_temporal = nm, one_state = bi,
               no_transition = "nestimate_landscape_unestimable"),
             note = c(mean_layer = "a mean-combined layer is refused even where its 0/1 weights equal its counts",
                      onehot_frequency = "every one-hot network is refused (its rows are not positions), with or without units",
                      onehot_no_unit = "a co-occurrence network is not a chain, with or without units",
                      density_macro = "its weights cannot be rebuilt from its sequences")),
    rf_cells("landscape_trend", c(latency = lr, single_sequence = ss))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("the posterior verbs read counts only, and refuse what has none", {
  counts_gone <- c(latency = bi, mean_layer = bi, density_macro = bi, bare = bi, mlvar_temporal = bi)
  cells <- rbind(
    rf_cells("certainty", c(counts_gone, single_sequence = "net_certainty", onehot_no_unit = bi, cor = bi,
                            no_transition = "net_certainty"),
             note = c(single_sequence = legit_posterior, no_transition = "no transition observed: the posterior is the prior, a legitimate Bayesian answer")),
    rf_cells("bayes_compare", c(counts_gone, single_sequence = "net_bayes", onehot_no_unit = bi, cor = bi,
                                no_transition = "net_bayes"),
             note = c(single_sequence = legit_posterior, no_transition = "no transition observed: the posterior is the prior")),
    rf_cells("entropy_bayes", c(counts_gone, single_sequence = "net_entropy_bayes", cor = bi, one_state = bi,
                                no_transition = bi, empty_group = bi),
             note = c(single_sequence = legit_posterior)),
    rf_cells("mosaic_test", c(latency = bi, mean_layer = bi, density_macro = bi, cor = bi)),
    rf_cells("boot_glasso", c(glasso = "boot_glasso", latency = bi, latency_group = bi, mean_group = bi, cor = bi,
                              bare = bi, group = bi),
             note = c(glasso = "boot_glasso() resamples the rows of the glasso network's data"))
  )
  rf_run(cells, inputs, verbs, covers)
})

test_that("tna models and group_tna are rebuilt and resampled", {
  skip_if_not_installed("tna")
  wide <- rf_wide(group = TRUE)
  tna_inputs <- list(tna = tna::tna(rf_wide()),
                     group_tna = tna::group_tna(wide, group = "g", cols = paste0("T", seq_len(10L))))
  expected <- list(
    bootstrap_network = c(tna = "net_bootstrap", group_tna = "net_bootstrap_group"),
    permutation = c(tna = "net_permutation", group_tna = "net_permutation_group"),
    centrality_stability = c(tna = "net_stability", group_tna = "net_stability_group"),
    casedrop_reliability = c(tna = "net_casedrop_reliability", group_tna = "net_casedrop_reliability_group"),
    network_reliability = c(tna = "net_reliability", group_tna = "net_reliability"),
    certainty = c(tna = "net_certainty", group_tna = "net_certainty_group"),
    entropy_bayes = c(tna = "net_entropy_bayes", group_tna = "net_entropy_bayes_group"),
    stability_landscape = c(tna = "net_landscape", group_tna = "net_landscape")
  )
  cells <- do.call(rbind, Map(\(verb, outcome) {
    rf_cells(verb, outcome, note = stats::setNames(rep(legit_tna, length(outcome)), names(outcome)))
  }, names(expected), expected))
  rf_run(cells, tna_inputs, verbs, covers)
})
