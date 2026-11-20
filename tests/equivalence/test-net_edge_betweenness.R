# net_edge_betweenness(): the edge-betweenness network of a transition network, against
# tna::betweenness_network() and against igraph::edge_betweenness() on the inverted weights (the definition tna
# follows). Re-expressed from local_testing_and_equivalence/test-equiv-edge-betweenness.R (96 random configs,
# maximum of differences) and the tna comparison of tests/testthat/test-edge-betweenness.R.

inputs <- list(rs_random_sequences(120L, 3L, 6L, 1L), rs_random_sequences(120L, 4L, 6L, 2L),
               rs_random_sequences(120L, 5L, 6L, 3L), rs_random_sequences(120L, 6L, 6L, 4L))

test_that("edge betweenness equals tna::betweenness_network(), inverted and not", {
  skip_without_reference("tna")
  invisible(lapply(seq_along(inputs), \(i) lapply(c(TRUE, FALSE), \(invert) {
    found <- net_edge_betweenness(build_network(inputs[[i]], method = "relative"), invert = invert)
    expected <- tna::betweenness_network(tna::tna(inputs[[i]]), invert = invert)$weights
    expect_reference(extract_transition_matrix(found), expected, reference = "tna",
                     covers = c("net_edge_betweenness()", "net_edge_betweenness.netobject()"),
                     quantity = sprintf("random set %d, invert = %s", i, invert), tolerance = 1e-10)
  })))
})

test_that("edge betweenness equals igraph::edge_betweenness() on inverted weights", {
  skip_without_reference("igraph")
  invisible(lapply(seq_along(inputs), \(i) {
    sequences <- as.matrix(inputs[[i]])
    states <- sort(unique(as.vector(sequences)))
    counts <- table(factor(as.vector(sequences[, -ncol(sequences)]), states), factor(as.vector(sequences[, -1L]), states))
    weights <- unclass(prop.table(counts, margin = 1L))
    graph <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE, diag = TRUE)
    ends <- igraph::ends(graph, igraph::E(graph))
    expected <- matrix(0, length(states), length(states), dimnames = list(states, states))
    expected[ends] <- igraph::edge_betweenness(graph, directed = TRUE, weights = 1 / igraph::E(graph)$weight)
    found <- net_edge_betweenness(build_network(inputs[[i]], method = "relative"))
    expect_reference(extract_transition_matrix(found), expected, reference = "igraph",
                     covers = c("net_edge_betweenness()", "net_edge_betweenness.netobject()"),
                     quantity = sprintf("random set %d, igraph", i), tolerance = 1e-10)
  }))
})

test_that("edge betweenness of a grouped network equals tna::betweenness_network() of each group", {
  skip_without_reference("tna")
  sequences <- inputs[[2L]]
  group <- rep(c("x", "y"), length.out = nrow(sequences))
  found <- net_edge_betweenness(build_network(cbind(sequences, grp = group), method = "relative", group = "grp"))
  # A group of networks has no tidy accessor: each group is read by name (API gap).
  invisible(lapply(c("x", "y"), \(g) {
    expect_reference(extract_transition_matrix(found[[g]]),
                     tna::betweenness_network(tna::tna(sequences[group == g, ]))$weights, reference = "tna",
                     covers = "net_edge_betweenness.netobject_group()", quantity = sprintf("group %s", g),
                     tolerance = 1e-10)
  }))
})
