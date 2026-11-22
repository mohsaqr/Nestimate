# Layer 8, the argument sweep of net_centrality(), extract_edges(), extract_transition_matrix() and
# extract_initial_probs(). net_centrality() is fitted at every measure, every combination of loops, normalize and
# invert, on relative, frequency and co-occurrence networks, against tna::centralities() at the same options; raw
# Diffusion (normalize_diffusion = FALSE) against its definition, sum over k = 1..n of the row sums of W^k; the
# undirected strengths and shortest-path measures against igraph. The extractors are swept over threshold,
# include_self, sort_by and type against the weight matrix read in base R. test-net_centrality.R checks relative
# networks at the tna defaults; this file adds the other options and network kinds.

ct_wide <- tn_simulate(n_actors = 50L, n_states = 5L, seq_length = 8L, seed = 831L, ragged = TRUE)
ct_all <- c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut", "Closeness", "Betweenness",
            "BetweennessRSP", "Diffusion", "Clustering")
ct_models <- c(relative = "tna", frequency = "ftna", co_occurrence = "ctna")

ct_table <- function(table, columns) {
  table <- as.data.frame(table)
  data.frame(state = as.character(table$state), table[columns], check.names = FALSE)
}

test_that("every measure at every loops / normalize / invert setting equals tna::centralities(), three kinds of network", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  cases <- expand.grid(method = names(ct_models), loops = c(FALSE, TRUE), normalize = c(FALSE, TRUE),
                       invert = c(TRUE, FALSE), stringsAsFactors = FALSE)
  invisible(Map(\(method, loops, normalize, invert) {
    net <- build_network(ct_wide, method = method)
    found <- withCallingHandlers(
      net_centrality(net, measures = ct_all, loops = loops, normalize = normalize, invert = invert),
      message = \(m) invokeRestart("muffleMessage"),
      nestimate_constant_measure = \(w) invokeRestart("muffleWarning"))
    reference <- ct_table(tna::centralities(getExportedValue("tna", ct_models[[method]])(ct_wide), measures = ct_all, loops = loops,
                                            normalize = normalize, invert = invert), ct_all)
    label <- sprintf("%s, loops %s, normalize %s, invert %s", method, loops, normalize, invert)
    undirected <- identical(method, "co_occurrence")
    # Betweenness of an undirected network: tna counts each unordered pair twice (it reads the network as
    # directed); igraph agrees with Nestimate (test below). Range normalisation cancels the factor 2, and with
    # invert = FALSE no node lies on a shortest path here, so only the raw inverted values differ.
    doubled <- undirected && invert && !normalize
    # Closeness ignoring direction with invert = FALSE: the weights are distances, so the undirected edge is the
    # shorter of w_ij and w_ji (R/centrality_measures.R .closeness()), as in tna and igraph. Checked on its own.
    shared <- setdiff(ct_all, c(if (doubled) "Betweenness", if (!undirected && !invert) "Closeness"))
    expect_reference(found[c("state", shared)], reference[c("state", shared)], reference = "tna",
                     covers = "net_centrality()", quantity = label, tolerance = 1e-10, key = "state")
    if (doubled) {
      expect_known_difference(found[c("state", "Betweenness")], reference[c("state", "Betweenness")], reference = "tna",
                              covers = "net_centrality()", quantity = sprintf("Betweenness, %s", label),
                              tolerance = 1e-10, key = "state",
                              reason = "tna counts the shortest paths of an undirected network in both directions; Nestimate counts each unordered pair once, as igraph::betweenness(directed = FALSE) does.",
                              agrees_with = "igraph::betweenness(directed = FALSE)")
    }
    if (!undirected && !invert) {
      expect_reference(found[c("state", "Closeness")], reference[c("state", "Closeness")], reference = "tna",
                       covers = "net_centrality()", quantity = sprintf("Closeness, %s", label), tolerance = 1e-10,
                       key = "state")
    }
  }, cases$method, cases$loops, cases$normalize, cases$invert))
})

test_that("raw Diffusion equals the summed walks of length 1 to n, with and without loops", {
  cases <- expand.grid(method = names(ct_models), loops = c(FALSE, TRUE), stringsAsFactors = FALSE)
  invisible(Map(\(method, loops) {
    net <- build_network(ct_wide, method = method)
    w <- unclass(extract_transition_matrix(net))
    if (!loops) diag(w) <- 0
    powers <- Reduce(\(acc, k) acc %*% w, seq_len(nrow(w)), diag(nrow(w)), accumulate = TRUE)[-1L]
    expected <- data.frame(state = rownames(w), Diffusion = rowSums(Reduce(`+`, powers)))
    found <- withCallingHandlers(
      net_centrality(net, measures = "Diffusion", loops = loops, normalize_diffusion = FALSE),
      message = \(m) invokeRestart("muffleMessage"))
    expect_reference(found, expected, reference = "independent: rowSums(sum_k W^k), k = 1..n, Banerjee et al. 2013",
                     covers = "net_centrality()", quantity = sprintf("raw Diffusion, %s, loops %s", method, loops),
                     tolerance = 1e-9, measure = "relative", key = "state")
    if (identical(method, "relative") && loops) {
      # ?net_centrality: every row of a relative network with loops sums to 1, so Diffusion is constant and its
      # range normalisation is undefined (NaN with a warning).
      expect_warning(net_centrality(net, measures = "Diffusion", loops = TRUE), class = "nestimate_constant_measure")
    } else {
      scaled <- withCallingHandlers(net_centrality(net, measures = "Diffusion", loops = loops),
                                    message = \(m) invokeRestart("muffleMessage"))
      raw <- expected$Diffusion
      expect_reference(scaled, data.frame(state = expected$state, Diffusion = (raw - min(raw)) / (max(raw) - min(raw))),
                       reference = "independent: the raw Diffusion range-normalised to 0-1",
                       covers = "net_centrality()", quantity = sprintf("Diffusion normalised, %s, loops %s", method, loops),
                       tolerance = 1e-10, key = "state")
    }
  }, cases$method, cases$loops))
})

test_that("strength, betweenness and closeness of an undirected network equal igraph, inverted and not", {
  skip_without_reference("igraph")
  net <- build_network(ct_wide, method = "co_occurrence")
  w <- unclass(extract_transition_matrix(net))
  diag(w) <- 0
  graph <- igraph::graph_from_adjacency_matrix(w, mode = "undirected", weighted = TRUE)
  invisible(lapply(c(TRUE, FALSE), \(invert) {
    distance <- if (invert) 1 / igraph::E(graph)$weight else igraph::E(graph)$weight
    expected <- data.frame(state = rownames(w), InStrength = igraph::strength(graph), OutStrength = igraph::strength(graph),
                           Betweenness = igraph::betweenness(graph, weights = distance, directed = FALSE),
                           Closeness = igraph::closeness(graph, weights = distance))
    found <- withCallingHandlers(
      net_centrality(net, measures = c("InStrength", "OutStrength", "Betweenness", "Closeness"), invert = invert),
      message = \(m) invokeRestart("muffleMessage"))
    expect_reference(found, expected, reference = "igraph", covers = "net_centrality()",
                     quantity = sprintf("co-occurrence network, invert %s", invert), tolerance = 1e-10, key = "state")
  }))
})

test_that("measure names, aliases and the defaults of net_centrality() equal the values written out", {
  net <- build_network(ct_wide, method = "relative")
  quiet <- \(...) withCallingHandlers(net_centrality(net, ...), message = \(m) invokeRestart("muffleMessage"))
  expect_reference(quiet(), quiet(measures = c("InStrength", "Betweenness", "Diffusion"), loops = FALSE,
                                  normalize = FALSE, invert = TRUE, normalize_diffusion = TRUE, centrality_fn = NULL),
                   reference = "independent: ?net_centrality, the default measures and options written out",
                   covers = "net_centrality()", quantity = "defaults", tolerance = 0, key = "state")
  expect_reference(quiet(measures = "all"), quiet(measures = ct_all),
                   reference = "independent: ?net_centrality, \"all\" is every built-in measure",
                   covers = "net_centrality()", quantity = "measures = \"all\"", tolerance = 0, key = "state")
  expect_reference(quiet(measures = c("InCloseness", "OutCloseness")),
                   stats::setNames(quiet(measures = c("ClosenessIn", "ClosenessOut")), c("state", "InCloseness", "OutCloseness")),
                   reference = "independent: ?net_centrality, InCloseness and OutCloseness are aliases",
                   covers = "net_centrality()", quantity = "closeness aliases", tolerance = 0, key = "state")
  custom <- \(m) list(RowMax = apply(m, 1L, max))
  expect_reference(quiet(measures = "RowMax", centrality_fn = custom),
                   data.frame(state = rownames(net$weights),
                              RowMax = apply(unclass(net$weights) * (1 - diag(nrow(net$weights))), 1L, max)),
                   reference = "independent: the centrality_fn applied to the weight matrix without loops (loops = FALSE)",
                   covers = "net_centrality()", quantity = "centrality_fn", tolerance = 0, key = "state")
})

test_that("every threshold, include_self and sort_by of extract_edges() equals the weight matrix read in base R", {
  correlations <- psy_continuous(n = 120L, p = 5L, rho = 1, seed = 832L)
  networks <- list(relative = build_network(ct_wide, method = "relative"),
                   frequency = build_network(ct_wide, method = "frequency"),
                   co_occurrence = build_network(ct_wide, method = "co_occurrence"),
                   cor = build_network(correlations, method = "cor"))
  cases <- expand.grid(network = names(networks), threshold = c(0, 0.1, 0.25), include_self = c(FALSE, TRUE),
                       sort_by = c("weight", "from", "to", "none"), stringsAsFactors = FALSE)
  invisible(Map(\(network, threshold, include_self, sort_by) {
    net <- networks[[network]]
    w <- unclass(net$weights)
    cells <- expand.grid(from = rownames(w), to = colnames(w), stringsAsFactors = FALSE)
    cells$weight <- as.vector(w)
    keep <- abs(cells$weight) > threshold & (net$directed | as.vector(row(w) <= col(w))) &
      (include_self | cells$from != cells$to)
    expected <- cells[keep, , drop = FALSE]
    found <- extract_edges(net, threshold = threshold, include_self = include_self,
                           sort_by = if (identical(sort_by, "none")) NULL else sort_by)
    label <- sprintf("%s, threshold %s, include_self %s, sort_by %s", network, threshold, include_self, sort_by)
    expect_reference(found, expected, reference = "independent: cells of the weight matrix above the threshold, base R",
                     covers = "extract_edges()", quantity = label, tolerance = 0, key = c("from", "to"))
    order_key <- switch(sort_by,
      weight = order(-expected$weight),
      from = order(expected$from, expected$to),
      to = order(expected$to, expected$from),
      none = seq_len(nrow(expected)))
    expect_reference(data.frame(rank = seq_len(nrow(found)), from = found$from, to = found$to),
                     data.frame(rank = seq_len(nrow(expected)), from = expected$from[order_key],
                                to = expected$to[order_key]),
                     reference = "independent: ?extract_edges row order", covers = "extract_edges()",
                     quantity = sprintf("row order, %s", label), tolerance = 0, key = "rank")
  }, cases$network, cases$threshold, cases$include_self, cases$sort_by))
  expect_error(extract_edges(networks$relative, threshold = -1), class = "nestimate_bad_input")
  expect_error(extract_edges(networks$relative, include_self = NA), class = "nestimate_bad_input")
  expect_error(extract_edges(networks$relative, sort_by = "size"), class = "nestimate_bad_input")
})

test_that("extract_transition_matrix(type) and the summaries of the extractors follow their definitions", {
  invisible(lapply(c("relative", "frequency", "co_occurrence", "attention", "gap"), \(method) {
    net <- build_network(ct_wide, method = method)
    w <- unclass(net$weights)
    expect_reference(extract_transition_matrix(net, type = "raw"), w,
                     reference = "independent: the stored weights", covers = "extract_transition_matrix()",
                     quantity = sprintf("type raw, %s", method), tolerance = 0)
    expect_reference(extract_transition_matrix(net, type = "scaled"), sw_row_normalise(w),
                     reference = "independent: each row divided by its sum", covers = "extract_transition_matrix()",
                     quantity = sprintf("type scaled, %s", method), tolerance = 1e-12)
    long <- summary(extract_transition_matrix(net))
    cells <- expand.grid(from = rownames(w), to = colnames(w), stringsAsFactors = FALSE)
    cells$weight <- as.vector(w)
    cells <- cells[cells$weight != 0, , drop = FALSE]
    expect_reference(long, cells, reference = "independent: non-zero cells of the weight matrix",
                     covers = "summary.nest_transition_matrix()", quantity = sprintf("non-zero cells, %s", method),
                     tolerance = 0, key = c("from", "to"))
    expect_reference(as.numeric(!is.unsorted(rev(abs(long$weight)))), 1,
                     reference = "independent: ?extract_transition_matrix, sorted by decreasing absolute weight",
                     covers = "summary.nest_transition_matrix()", quantity = sprintf("row order, %s", method),
                     tolerance = 0)
  }))
  cooc <- build_network(ct_wide, method = "co_occurrence")
  expect_warning(uniform <- extract_initial_probs(cooc))
  expect_reference(unclass(uniform), stats::setNames(rep(1 / nrow(cooc$weights), nrow(cooc$weights)), rownames(cooc$weights)),
                   reference = "independent: ?extract_initial_probs, a uniform placeholder for co-occurrence",
                   covers = "extract_initial_probs()", quantity = "co-occurrence placeholder", tolerance = 1e-15)
})
