# vertex_bootstrap(), vertex_compare(): the observed network-level statistics (density, mean non-zero weight,
# Freeman-type strength centralisation and weighted reciprocity, on the off-diagonal cells) against their
# definitions written in base R on the transition matrix of the sequences, and density also against
# igraph::edge_density(). The bootstrap replicates themselves draw vertices inside the function with no argument to
# supply them, so only the observed values are compared.

rs_vertex_statistics <- function(sequences) {
  sequences <- as.matrix(sequences)
  states <- sort(unique(as.vector(sequences)))
  counts <- table(factor(as.vector(sequences[, -ncol(sequences)]), states), factor(as.vector(sequences[, -1L]), states))
  weights <- unclass(prop.table(counts, margin = 1L))
  off <- weights[row(weights) != col(weights)]
  strength <- rowSums(weights) + colSums(weights) - 2 * diag(weights)
  c(density = mean(off != 0), mean_weight = mean(off[off != 0]),
    centralization = sum(max(strength) - strength) / ((length(strength) - 1) * max(strength)),
    reciprocity = sum(pmin(weights, t(weights))[row(weights) != col(weights)]) / sum(off))
}

sparse <- data.frame(T1 = c("A", "B", "C", "A", "D", "B"), T2 = c("B", "C", "A", "C", "A", "D"),
                     T3 = c("C", "C", "B", "A", "B", "A"), T4 = c("A", "B", "B", "D", "B", "C"))
dense <- rs_random_sequences(60L, 5L, 8L, 41L)

test_that("observed statistics equal their definitions", {
  invisible(lapply(list(sparse = sparse, dense = dense), \(sequences) {
    statistics <- rs_vertex_statistics(sequences)
    found <- vertex_bootstrap(build_network(sequences, method = "relative"), iter = 20L, seed = 1)
    expect_reference(summary(found), data.frame(statistic = names(statistics), observed = unname(statistics)),
                     reference = "independent: network statistics from their definitions, base R",
                     covers = c("vertex_bootstrap()", "summary.net_vertex_bootstrap()"),
                     quantity = sprintf("observed, %d sequences", nrow(sequences)), tolerance = 1e-12,
                     key = "statistic")
  }))
})

test_that("observed density equals igraph::edge_density() without loops", {
  skip_without_reference("igraph")
  found <- vertex_bootstrap(build_network(sparse, method = "relative"), iter = 20L, seed = 1)
  weights <- unclass(prop.table(table(factor(as.vector(as.matrix(sparse)[, -4L])), factor(as.vector(as.matrix(sparse)[, -1L]))), 1L))
  diag(weights) <- 0
  graph <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE)
  expect_reference(subset(summary(found), statistic == "density"),
                   data.frame(statistic = "density", observed = igraph::edge_density(graph, loops = FALSE)),
                   reference = "igraph", covers = "vertex_bootstrap()", quantity = "observed density", tolerance = 1e-12,
                   key = "statistic")
})

test_that("vertex_compare() reports both observed statistics and their difference", {
  x <- rs_vertex_statistics(sparse)
  y <- rs_vertex_statistics(dense)
  found <- vertex_compare(build_network(sparse, method = "relative"), build_network(dense, method = "relative"),
                          iter = 20L, seed = 1)
  expect_reference(summary(found), data.frame(statistic = names(x), observed_x = unname(x), observed_y = unname(y),
                                              diff = unname(x - y)),
                   reference = "independent: network statistics from their definitions, base R",
                   covers = c("vertex_compare()", "summary.net_vertex_comparison()"),
                   quantity = "observed, sparse vs dense", tolerance = 1e-12, key = "statistic")
})
