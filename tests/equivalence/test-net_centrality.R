# net_centrality(): node centralities of a transition network, against tna::centralities() (the nine measures the
# two share, with and without loops and range normalisation) and against igraph from the definitions (strength,
# shortest-path betweenness and closeness on inverted weights). Nestimate computes these itself: igraph is not
# called by net_centrality(). Re-expressed from local_testing_and_equivalence/test-equiv-moved-centrality_stability.R
# and the centrality parts of test-equiv-tna-functions.R.

measures <- c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut", "Closeness", "Betweenness",
              "BetweennessRSP", "Clustering")
inputs <- list(rs_random_sequences(100L, 5L, 6L, 5L), rs_random_sequences(40L, 4L, 12L, 6L),
               rs_random_sequences(60L, 6L, 8L, 7L))

rs_centralities <- function(table, columns) {
  table <- as.data.frame(table)
  data.frame(state = as.character(table$state), table[columns], check.names = FALSE)
}

test_that("eight measures equal tna::centralities() on three random sequence sets", {
  skip_without_reference("tna")
  invisible(lapply(seq_along(inputs), \(i) {
    found <- net_centrality(build_network(inputs[[i]], method = "relative"), measures = measures)
    expected <- rs_centralities(tna::centralities(tna::tna(inputs[[i]]), measures = measures), measures)
    expect_reference(found, expected, reference = "tna", covers = "net_centrality()",
                     quantity = sprintf("random set %d", i), tolerance = 1e-10, key = "state")
  }))
})

test_that("raw Diffusion equals tna 1.3.1's Diffusion; the default is scaled to 0-1", {
  skip_without_reference("tna")
  net <- build_network(inputs[[1L]], method = "relative")
  # tna 1.3.1 numbers the states instead of naming them when it is asked for one measure: ask for two.
  expected <- rs_centralities(tna::centralities(tna::tna(inputs[[1L]]), measures = c("InStrength", "Diffusion")),
                              "Diffusion")
  expect_reference(net_centrality(net, measures = "Diffusion", normalize_diffusion = FALSE), expected,
                   reference = "tna", covers = "net_centrality()", quantity = "Diffusion, normalize_diffusion = FALSE",
                   tolerance = 1e-10, key = "state", note = "tna 1.3.1 returns Diffusion raw")
  expect_known_difference(net_centrality(net, measures = "Diffusion"), expected, reference = "tna",
                          covers = "net_centrality()", quantity = "Diffusion, default", tolerance = 1e-10, key = "state",
                          reason = "net_centrality() range-normalises Diffusion to 0-1 by default (normalize_diffusion = TRUE); tna 1.3.1 returns it raw. tna >= 1.5 scales it as Nestimate does (FINDINGS.md).",
                          agrees_with = "normalize_diffusion = FALSE")
})

test_that("range-normalised centralities equal tna::centralities(normalize = TRUE)", {
  skip_without_reference("tna")
  all_measures <- c(measures, "Diffusion")
  invisible(lapply(seq_along(inputs), \(i) {
    found <- net_centrality(build_network(inputs[[i]], method = "relative"), measures = all_measures, normalize = TRUE)
    expected <- rs_centralities(tna::centralities(tna::tna(inputs[[i]]), measures = all_measures, normalize = TRUE),
                                 all_measures)
    expect_reference(found, expected, reference = "tna", covers = "net_centrality()",
                     quantity = sprintf("random set %d, normalize = TRUE", i), tolerance = 1e-10, key = "state")
  }))
})

test_that("centralities with loops equal tna::centralities(loops = TRUE)", {
  skip_without_reference("tna")
  found <- net_centrality(build_network(inputs[[2L]], method = "relative"), measures = measures, loops = TRUE)
  expected <- rs_centralities(tna::centralities(tna::tna(inputs[[2L]]), measures = measures, loops = TRUE), measures)
  expect_reference(found, expected, reference = "tna", covers = "net_centrality()",
                   quantity = "random set 2, loops = TRUE", tolerance = 1e-10, key = "state")
})

test_that("strength, betweenness and closeness equal igraph on inverted weights", {
  skip_without_reference("igraph")
  invisible(lapply(seq_along(inputs), \(i) {
    # The transition matrix of the sequences, without self-loops: counts of consecutive pairs, row-normalised.
    sequences <- as.matrix(inputs[[i]])
    pairs <- data.frame(from = as.vector(sequences[, -ncol(sequences)]), to = as.vector(sequences[, -1L]))
    counts <- table(factor(pairs$from, sort(unique(as.vector(sequences)))),
                    factor(pairs$to, sort(unique(as.vector(sequences)))))
    weights <- unclass(prop.table(counts, margin = 1L))
    diag(weights) <- 0
    graph <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE)
    distance <- 1 / igraph::E(graph)$weight
    expected <- data.frame(state = rownames(weights),
                           InStrength = igraph::strength(graph, mode = "in"),
                           OutStrength = igraph::strength(graph, mode = "out"),
                           Betweenness = igraph::betweenness(graph, weights = distance),
                           ClosenessIn = igraph::closeness(graph, mode = "in", weights = distance),
                           ClosenessOut = igraph::closeness(graph, mode = "out", weights = distance))
    found <- net_centrality(build_network(inputs[[i]], method = "relative"),
                            measures = c("InStrength", "OutStrength", "Betweenness", "ClosenessIn", "ClosenessOut"))
    expect_reference(found, expected, reference = "igraph", covers = "net_centrality()",
                     quantity = sprintf("random set %d, igraph", i), tolerance = 1e-10, key = "state")
  }))
})
