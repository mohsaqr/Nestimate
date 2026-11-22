# Layer 6, adversarial data: the transition-network verbs on the inputs most likely to break them (a missing first
# state, trailing and internal NA, a state seen once, a state never left and one never entered, a row of length one,
# one actor, identical rows, two states, twenty states with sparse counts, numeric-looking and odd labels, a cycle
# whose centralities tie, a group of two and a group of one, a long log with duplicated time stamps, constant
# weights). Every quantity a verb returns is compared: whole matrices, every edge, every metric row and every
# centrality column, with no filter on the package's side. References: base R from the definitions
# (helper-layer-adversarial.R), tna 1.5 (skip_without_reference()), igraph for the network metrics.
# Tolerances, fixed before comparing: counts are exact integers (0); a probability is one division (1e-12);
# attention, gap and centrality values are sums of a few dozen terms (1e-10).

sets <- adv_sequence_sets()
exact <- 0
division <- 1e-12
sums <- 1e-10

test_that("relative and frequency weights and initial probabilities equal the base-R tables on every adversarial set", {
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    counts <- adv_counts(wide)
    expect_reference(extract_transition_matrix(build_network(wide, method = "frequency")), counts,
                     reference = "independent: adjacent observed pairs counted, base R",
                     covers = c("build_network()", "extract_transition_matrix()"),
                     quantity = sprintf("frequency weights, %s", name), tolerance = exact)
    relative <- build_network(wide, method = "relative")
    expect_reference(extract_transition_matrix(relative), adv_relative(counts),
                     reference = "independent: adjacent observed pairs counted and row-normalised, base R",
                     covers = c("build_network()", "extract_transition_matrix()"),
                     quantity = sprintf("relative weights, %s", name), tolerance = division)
    expect_reference(extract_initial_probs(relative), adv_initial(wide),
                     reference = "independent: share of rows by their first observed state, base R",
                     covers = c("build_network()", "extract_initial_probs()"),
                     quantity = sprintf("initial probabilities, %s", name), tolerance = division)
    expect_reference(summary(extract_initial_probs(relative)),
                     data.frame(state = adv_states(wide), prob = unname(adv_initial(wide))),
                     reference = "independent: share of rows by their first observed state, base R",
                     covers = "summary.nest_initial_probs()", quantity = sprintf("initial probabilities table, %s", name),
                     tolerance = division, key = "state")
  }))
})

test_that("co-occurrence, attention, n-gram, gap and reverse weights equal their positional definitions", {
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    weights <- \(method, params = list()) extract_transition_matrix(build_network(wide, method = method, params = params))
    expect_reference(weights("co_occurrence"), adv_cooccurrence(wide),
                     reference = "independent: every unordered pair of observed positions in a row, base R",
                     covers = "build_network()", quantity = sprintf("co-occurrence weights, %s", name), tolerance = exact)
    expect_reference(weights("attention"), adv_attention(wide, lambda = 1),
                     reference = "independent: exp(-(j - i) / lambda) over ordered pairs of observed positions, base R",
                     covers = "build_network()", quantity = sprintf("attention weights lambda 1, %s", name),
                     tolerance = sums)
    expect_reference(weights("attention", list(lambda = 3)), adv_attention(wide, lambda = 3),
                     reference = "independent: exp(-(j - i) / lambda) over ordered pairs of observed positions, base R",
                     covers = "build_network()", quantity = sprintf("attention weights lambda 3, %s", name),
                     tolerance = sums)
    expect_reference(weights("ngram"), adv_ngram(wide, 2L),
                     reference = "independent: adjacent observed pairs of every window of n positions, base R",
                     covers = "build_network()", quantity = sprintf("n-gram weights n 2, %s", name), tolerance = exact)
    expect_reference(weights("ngram", list(n_gram = 3L)), adv_ngram(wide, 3L),
                     reference = "independent: adjacent observed pairs of every window of n positions, base R",
                     covers = "build_network()", quantity = sprintf("n-gram weights n 3, %s", name), tolerance = exact)
    expect_reference(weights("gap"), adv_gap(wide, 1L),
                     reference = "independent: 1 / distance over observed pairs at most max_gap + 1 apart, base R",
                     covers = "build_network()", quantity = sprintf("gap weights max_gap 1, %s", name), tolerance = sums)
    expect_reference(weights("gap", list(max_gap = 3L)), adv_gap(wide, 3L),
                     reference = "independent: 1 / distance over observed pairs at most max_gap + 1 apart, base R",
                     covers = "build_network()", quantity = sprintf("gap weights max_gap 3, %s", name), tolerance = sums)
    expect_reference(weights("reverse"), t(adv_counts(wide)),
                     reference = "independent: the transposed table of adjacent observed pairs, base R",
                     covers = "build_network()", quantity = sprintf("reverse weights, %s", name), tolerance = exact)
  }))
})

test_that("every transition method equals tna on every adversarial set", {
  skip_without_reference("tna")
  models <- list(relative = \(w) tna::tna(w), frequency = \(w) tna::ftna(w), co_occurrence = \(w) tna::ctna(w),
                 attention = \(w) tna::atna(w), ngram = \(w) tna::build_model(w, type = "n-gram"),
                 gap = \(w) tna::build_model(w, type = "gap"), reverse = \(w) tna::build_model(w, type = "reverse"))
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    invisible(lapply(names(models), \(method) {
      expect_reference(extract_transition_matrix(build_network(wide, method = method)),
                       unclass(models[[method]](wide)$weights), reference = "tna", covers = "build_network()",
                       quantity = sprintf("%s weights, %s", method, name), tolerance = sums)
    }))
  }))
})

test_that("initial probabilities equal tna's, except where the first position is missing", {
  skip_without_reference("tna")
  invisible(lapply(setdiff(names(sets), "leading_na"), \(name) {
    expect_reference(extract_initial_probs(build_network(sets[[name]], method = "relative")), tna::tna(sets[[name]])$inits,
                     reference = "tna", covers = "extract_initial_probs()",
                     quantity = sprintf("initial probabilities, %s", name), tolerance = division)
  }))
  expect_known_difference(extract_initial_probs(build_network(sets$leading_na, method = "relative")),
                          tna::tna(sets$leading_na)$inits, reference = "tna", covers = "extract_initial_probs()",
                          quantity = "initial probabilities, leading_na", tolerance = division,
                          reason = paste("A row whose first position is missing starts, for Nestimate, at its first",
                                         "observed state; tna reads the first column only and leaves such a row out."),
                          agrees_with = "the share of rows by their first observed state (compared above)")
})

test_that("extract_edges() returns every non-zero cell and no zero cell, with and without self-loops", {
  invisible(lapply(names(sets), \(name) {
    net <- build_network(sets[[name]], method = "relative")
    weights <- adv_relative(adv_counts(sets[[name]]))
    expect_reference(extract_edges(net), adv_edge_list(weights),
                     reference = "independent: non-zero off-diagonal cells of the relative table, base R",
                     covers = "extract_edges()", quantity = sprintf("edge list, %s", name), tolerance = division,
                     key = c("from", "to"))
    expect_reference(extract_edges(net, include_self = TRUE), adv_edge_list(weights, include_self = TRUE),
                     reference = "independent: non-zero cells of the relative table, base R",
                     covers = "extract_edges()", quantity = sprintf("edge list with loops, %s", name),
                     tolerance = division, key = c("from", "to"))
    expect_reference(extract_edges(net, threshold = 0.3), adv_edge_list(weights, threshold = 0.3),
                     reference = "independent: off-diagonal cells above the threshold, base R",
                     covers = "extract_edges()", quantity = sprintf("edge list above 0.3, %s", name),
                     tolerance = division, key = c("from", "to"))
    expect_reference(summary(extract_transition_matrix(net)), adv_edge_list(weights, include_self = TRUE),
                     reference = "independent: non-zero cells of the relative table, base R",
                     covers = "summary.nest_transition_matrix()", quantity = sprintf("transition matrix table, %s", name),
                     tolerance = division, key = c("from", "to"))
    frequency <- build_network(sets[[name]], method = "frequency")
    expect_reference(extract_edges(frequency, include_self = TRUE), adv_edge_list(adv_counts(sets[[name]]), include_self = TRUE),
                     reference = "independent: non-zero cells of the count table, base R",
                     covers = "extract_edges()", quantity = sprintf("frequency edge list with loops, %s", name),
                     tolerance = exact, key = c("from", "to"))
  }))
})

test_that("summary() of a network equals igraph's metrics on every adversarial set", {
  skip_without_reference("igraph")
  invisible(lapply(names(sets), \(name) {
    expect_reference(summary(build_network(sets[[name]], method = "relative")),
                     adv_igraph_summary(adv_relative(adv_counts(sets[[name]]))), reference = "igraph",
                     covers = "summary.netobject()", quantity = sprintf("network metrics, %s", name), tolerance = sums,
                     key = "metric")
  }))
})

test_that("net_centrality() equals tna::centralities() on every adversarial set, all nine measures", {
  skip_without_reference("tna")
  # tna >= 1.5 range-normalises Diffusion as net_centrality() does; older tna does not
  skip_if(utils::packageVersion("tna") < "1.5.7", "needs tna >= 1.5.7 (range-normalised Diffusion)")
  measures <- c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut", "Closeness", "Betweenness",
                "BetweennessRSP", "Diffusion", "Clustering")
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    reference <- as.data.frame(tna::centralities(tna::tna(wide), measures = measures))
    expected <- data.frame(state = as.character(reference$state), reference[measures], check.names = FALSE)
    net <- build_network(wide, method = "relative")
    # tna >= 1.5 range-normalises Diffusion, as net_centrality() does by default; a Diffusion equal for every node
    # cannot be normalised: both sides give NaN (not comparable), Nestimate with a classed warning
    found <- if (name %in% c("identical", "ties")) {
      expect_warning(centralities <- net_centrality(net, measures = measures), class = "nestimate_constant_measure")
      centralities
    } else net_centrality(net, measures = measures)
    expect_reference(found, expected, reference = "tna", covers = "net_centrality()",
                     quantity = sprintf("nine centralities, %s", name), tolerance = sums, key = "state")
  }))
})

test_that("groups of very different sizes and a group of one equal the per-group tables and tna::group_tna()", {
  groups <- adv_group_sets(sets$leading_na)
  invisible(lapply(names(groups), \(name) {
    data <- groups[[name]]
    states <- adv_states(sets$leading_na)
    fit <- build_network(data, method = "relative", group = "grp")
    invisible(lapply(unique(data$grp), \(g) {
      rows <- data[data$grp == g, setdiff(names(data), "grp"), drop = FALSE]
      expect_reference(extract_transition_matrix(fit[[g]]), adv_relative(adv_counts(rows, states)),
                       reference = "independent: adjacent observed pairs of the group's rows, base R",
                       covers = "build_network()", quantity = sprintf("%s, group %s weights", name, g), tolerance = division)
      expect_reference(extract_initial_probs(fit[[g]]), adv_initial(rows, states),
                       reference = "independent: share of the group's rows by first observed state, base R",
                       covers = "build_network()", quantity = sprintf("%s, group %s initial", name, g), tolerance = division)
    }))
  }))
})

test_that("summary() of a grouped network equals igraph's metrics for every group and every metric", {
  skip_without_reference("igraph")
  groups <- adv_group_sets(sets$absorbing)
  invisible(lapply(names(groups), \(name) {
    data <- groups[[name]]
    states <- adv_states(sets$absorbing)
    per_group <- lapply(stats::setNames(nm = sort(unique(data$grp))), \(g) {
      rows <- data[data$grp == g, setdiff(names(data), "grp"), drop = FALSE]
      adv_igraph_summary(adv_relative(adv_counts(rows, states)))$value
    })
    expected <- data.frame(metric = adv_igraph_summary(adv_relative(adv_counts(sets$absorbing)))$metric, per_group,
                           check.names = FALSE)
    expect_reference(summary(build_network(data, method = "relative", group = "grp")), expected, reference = "igraph",
                     covers = "summary.netobject_group()", quantity = sprintf("grouped network metrics, %s", name),
                     tolerance = sums, key = "metric")
  }))
})

test_that("a long log with duplicated time stamps keeps tied events in log order", {
  log <- adv_long_duplicated_times()
  # the reference orders each actor's events by time, ties broken by their row in the log
  ordered <- log[order(log$actor, log$time, seq_len(nrow(log))), ]
  wide <- adv_wide(unname(split(ordered$action, ordered$actor)))
  fit <- build_network(log, method = "frequency", action = "action", actor = "actor", time = "time")
  expect_reference(extract_transition_matrix(fit), adv_counts(wide),
                   reference = "independent: pairs of each actor's events ordered by time, ties in log order, base R",
                   covers = "build_network()", quantity = "frequency weights, duplicated time stamps", tolerance = exact)
  expect_reference(extract_initial_probs(build_network(log, method = "relative", action = "action", actor = "actor",
                                                       time = "time")),
                   adv_initial(wide), reference = "independent: first event of each actor by time, base R",
                   covers = "build_network()", quantity = "initial probabilities, duplicated time stamps",
                   tolerance = division)
})

test_that("a matrix of constant weights gives uniform rows, tied centralities and igraph's metrics", {
  constant <- matrix(2, 4L, 4L, dimnames = list(letters[1:4], letters[1:4]))
  net <- build_network(constant, method = "relative")
  expect_reference(extract_transition_matrix(net), constant / rowSums(constant),
                   reference = "independent: constant rows normalised, base R", covers = "build_network()",
                   quantity = "relative weights of a constant matrix", tolerance = division)
  skip_without_reference("igraph")
  expect_reference(summary(net), adv_igraph_summary(constant / rowSums(constant)), reference = "igraph",
                   covers = "summary.netobject()", quantity = "network metrics of a constant matrix", tolerance = sums,
                   key = "metric")
  weights <- constant / rowSums(constant)
  diag(weights) <- 0
  graph <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE)
  distance <- 1 / igraph::E(graph)$weight
  expect_reference(net_centrality(net, measures = c("OutStrength", "InStrength", "Betweenness", "ClosenessIn", "ClosenessOut")),
                   data.frame(state = letters[1:4], OutStrength = igraph::strength(graph, mode = "out"),
                              InStrength = igraph::strength(graph, mode = "in"),
                              Betweenness = igraph::betweenness(graph, weights = distance),
                              ClosenessIn = igraph::closeness(graph, mode = "in", weights = distance),
                              ClosenessOut = igraph::closeness(graph, mode = "out", weights = distance)),
                   reference = "igraph", covers = "net_centrality()", quantity = "centralities of a constant matrix",
                   tolerance = sums, key = "state")
})

test_that("start and end states frame every observed sequence, whatever its missing positions", {
  boundary_defect <- paste("start = TRUE prepends the start state as a column, so a row whose first position is",
                           "missing loses its start -> first-state transition, and a row with no observed state gets",
                           "a start -> end transition (R/estimators.R:99 and :106); tna's begin_state does the same")
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    counts <- adv_boundary_counts(wide)
    # the documented markers: start before each sequence, end in the cell after its last observed state
    expect_reference(extract_transition_matrix(build_network(wide, method = "frequency", start = TRUE, end = TRUE)),
                     counts, reference = "independent: start before the first and end after the last observed state, base R",
                     covers = "build_network()", quantity = sprintf("frequency weights with start and end, %s", name),
                     tolerance = exact, defect = if (anyNA(wide[[1L]])) boundary_defect)
    expect_reference(extract_transition_matrix(build_network(wide, method = "relative", start = TRUE, end = TRUE)),
                     adv_relative(counts),
                     reference = "independent: start before the first and end after the last observed state, base R",
                     covers = "build_network()", quantity = sprintf("relative weights with start and end, %s", name),
                     tolerance = division, defect = if (anyNA(wide[[1L]])) boundary_defect)
  }))
})

test_that("start and end states equal tna's begin_state and end_state", {
  skip_without_reference("tna")
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    expect_reference(extract_transition_matrix(build_network(wide, method = "frequency", start = TRUE, end = TRUE)),
                     unclass(tna::ftna(wide, begin_state = "Start", end_state = "End")$weights), reference = "tna",
                     covers = "build_network()", quantity = sprintf("frequency weights with start and end, %s", name),
                     tolerance = exact)
  }))
})

test_that("a network without a single transition: zero weights, igraph's metrics, an empty edge list", {
  none <- data.frame(T1 = c("A", "B", "C"), T2 = NA_character_)
  net <- build_network(none, method = "relative")
  zero <- adv_relative(adv_counts(none))
  expect_reference(extract_transition_matrix(net), zero, reference = "independent: no adjacent observed pair, base R",
                   covers = "build_network()", quantity = "relative weights, no transition", tolerance = exact)
  expect_identical(nrow(extract_edges(net)), 0L)
  skip_without_reference("igraph")
  # Mean Distance and Reciprocity are undefined without an edge: NaN on both sides, recorded as not comparable
  expect_reference(summary(net), adv_igraph_summary(zero), reference = "igraph", covers = "summary.netobject()",
                   quantity = "network metrics, no transition", tolerance = sums, key = "metric")
})
