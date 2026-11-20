# hypergraph_laplacian(), hypergraph_transduction(), label_scores(), label_predictions(), hypergraph_cluster(),
# node_clusters(), hypergraph_centrality(): spectral methods on a hypergraph. References: the Zhou, Huang and
# Schoelkopf (2006) normalised Laplacian and transductive closed form written from the definition; HyperNetX
# (PNNL; Hayashi, Aksoy, Park and Park 2020) for the random-walk Laplacian, its stationary distribution and the
# spectral partition; igraph's eigenvector centrality of the clique expansion (Benson 2019).
# Re-expressed from local_testing_and_equivalence/test-equiv-hypergraph-laplacian.R (HyperNetX, 30 random
# hypergraphs by maxima) and test-equiv-hypergraph.R / test-equiv-moved-hypergraph_centrality.R (igraph, cosine).
# The Laplacian is returned as a matrix with the stationary distribution in its attribute "pi" (API gap).

inputs <- list(list(n_people = 8L, n_groups = 7L, seed = 711L), list(n_people = 12L, n_groups = 10L, seed = 712L))

# A connected membership table: a chain of pairs p1-p2, p2-p3, ... guarantees connection.
connected_memberships <- function(input, weighted) {
  random <- hg_memberships(input$n_people, input$n_groups, 4L, input$seed, weighted = weighted)
  chain <- data.frame(person = sprintf("p%02d", c(seq_len(input$n_people - 1L), seq_len(input$n_people - 1L) + 1L)),
                      group = sprintf("c%02d", rep(seq_len(input$n_people - 1L), 2L)), w = 1)
  rbind(random, chain)
}

test_that("the Zhou Laplacian and its stationary distribution equal the definition", {
  invisible(lapply(inputs, \(input) {
    memberships <- connected_memberships(input, weighted = FALSE)
    laplacian <- hypergraph_laplacian(bipartite_groups(memberships, player = "person", group = "group"), type = "zhou")
    zhou <- hg_zhou(hg_incidence(memberships) * 1)
    reference <- "independent: Zhou, Huang and Schoelkopf (2006) normalised Laplacian, base R"
    expect_reference(unclass(laplacian), zhou$laplacian, reference = reference, covers = "hypergraph_laplacian()",
                     quantity = sprintf("Zhou Laplacian, %d people", input$n_people), tolerance = 1e-12)
    expect_reference(attr(laplacian, "pi"), zhou$pi, reference = reference, covers = "hypergraph_laplacian()",
                     quantity = sprintf("Zhou stationary distribution, %d people", input$n_people), tolerance = 1e-12)
  }))
})

test_that("transductive scores and predictions equal (1 - xi) (I - xi Theta)^-1 Y", {
  invisible(lapply(inputs, \(input) {
    memberships <- connected_memberships(input, weighted = FALSE)
    hypergraph <- bipartite_groups(memberships, player = "person", group = "group")
    labels <- c(p01 = "a", p02 = "b", p05 = "a")
    found <- hypergraph_transduction(hypergraph, labels = labels, xi = 0.9)
    zhou <- hg_zhou(hg_incidence(memberships) * 1)
    nodes <- rownames(zhou$theta)
    seeds <- cbind(a = as.numeric(nodes %in% c("p01", "p05")), b = as.numeric(nodes == "p02"))
    scores <- 0.1 * solve(diag(length(nodes)) - 0.9 * zhou$theta, seeds)
    reference <- "independent: Zhou, Huang and Schoelkopf (2006) closed form, base R"
    label <- sprintf("%d people", input$n_people)
    expect_reference(label_scores(found),
                     data.frame(node = rep(nodes, 2L), class = rep(c("a", "b"), each = length(nodes)), score = as.vector(scores)),
                     reference = reference, covers = c("hypergraph_transduction()", "label_scores()"),
                     quantity = sprintf("class scores, %s", label), tolerance = 1e-12, key = c("node", "class"))
    winner <- data.frame(node = nodes, predicted = c("a", "b")[max.col(scores, ties.method = "first")],
                         score = apply(scores, 1L, max))
    expect_reference(label_predictions(found), winner, reference = reference,
                     covers = c("hypergraph_transduction()", "label_predictions()"),
                     quantity = sprintf("predicted class, %s", label), tolerance = 1e-12, key = "node")
    expect_reference(as.data.frame(found), winner, reference = reference,
                     covers = c("hypergraph_transduction()", "as.data.frame.net_hypergraph_transduction()"),
                     quantity = sprintf("predicted class, as.data.frame(), %s", label), tolerance = 1e-12, key = "node")
    expect_reference(summary(found), data.frame(class = c("a", "b"), n_predicted = as.vector(table(factor(winner$predicted, c("a", "b"))))),
                     reference = reference, covers = c("hypergraph_transduction()", "summary.net_hypergraph_transduction()"),
                     quantity = sprintf("nodes predicted per class, %s", label), tolerance = 0, key = "class")
  }))
})

test_that("the random-walk Laplacian and its stationary distribution equal HyperNetX", {
  skip_without_reference("reticulate")
  pandas <- reticulate::import("pandas")
  hypernetx <- reticulate::import("hypernetx")
  clustering <- reticulate::import("hypernetx.algorithms.clustering.laplacians_clustering")
  invisible(lapply(inputs, \(input) {
    memberships <- connected_memberships(input, weighted = TRUE)
    laplacian <- hypergraph_laplacian(bipartite_groups(memberships, player = "person", group = "group", weight = "w"),
                                      type = "random_walk")
    hyper <- hypernetx$Hypergraph(pandas$DataFrame(list(edges = memberships$group, nodes = memberships$person,
                                                        cell_weights = memberships$w)),
                                  edge_col = "edges", node_col = "nodes", cell_weight_col = "cell_weights")
    normalised <- clustering$norm_lap(hyper, weights = TRUE, index = TRUE)
    order <- unlist(normalised[[2L]])
    expected <- as.matrix(normalised[[1L]])
    dimnames(expected) <- list(order, order)
    transitions <- clustering$prob_trans(hyper, weights = TRUE, index = TRUE)
    stationary <- stats::setNames(as.vector(clustering$get_pi(transitions[[1L]])), unlist(transitions[[2L]]))
    expect_reference(unclass(laplacian), expected, reference = "hypernetx", covers = "hypergraph_laplacian()",
                     quantity = sprintf("random-walk Laplacian, %d people", input$n_people), tolerance = 1e-12)
    expect_reference(attr(laplacian, "pi"), stationary, reference = "hypernetx", covers = "hypergraph_laplacian()",
                     quantity = sprintf("random-walk stationary distribution, %d people", input$n_people), tolerance = 1e-12)
  }))
})

test_that("the spectral partition of two planted communities equals HyperNetX spec_clus()", {
  skip_without_reference("reticulate")
  pandas <- reticulate::import("pandas")
  hypernetx <- reticulate::import("hypernetx")
  clustering <- reticulate::import("hypernetx.algorithms.clustering.laplacians_clustering")
  memberships <- data.frame(person = c("a", "b", "c", "a", "b", "c", "d", "e", "f", "d", "e", "f", "c", "d"),
                            group = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3", "m4", "m4", "m4", "m5", "m5"),
                            w = 1)
  found <- hypergraph_cluster(bipartite_groups(memberships, player = "person", group = "group"), k = 2,
                              type = "random_walk", seed = 1)
  hyper <- hypernetx$Hypergraph(pandas$DataFrame(list(edges = memberships$group, nodes = memberships$person,
                                                      cell_weights = memberships$w)),
                                edge_col = "edges", node_col = "nodes", cell_weight_col = "cell_weights")
  parts <- clustering$spec_clus(hyper, k = 2L, weights = FALSE)
  members <- lapply(parts, unlist)
  reference_nodes <- unlist(members)
  reference_cluster <- rep(seq_along(members), lengths(members))
  ours <- node_clusters(found)
  # Labels of clusters are arbitrary on both sides: the partition is compared as "same cluster or not" per pair.
  expect_reference(hg_same_cluster(ours$node, ours$cluster), hg_same_cluster(reference_nodes, reference_cluster),
                   reference = "hypernetx", covers = c("hypergraph_cluster()", "node_clusters()"),
                   quantity = "same cluster, two planted communities", tolerance = 0, key = c("from", "to"))
  transitions <- clustering$prob_trans(hyper, weights = FALSE, index = TRUE)
  stationary <- data.frame(node = unlist(transitions[[2L]]), pi = as.vector(clustering$get_pi(transitions[[1L]])))
  expect_reference(as.data.frame(found), stationary, reference = "hypernetx",
                   covers = c("hypergraph_cluster()", "as.data.frame.net_hypergraph_cluster()"),
                   quantity = "stationary distribution, two planted communities", tolerance = 1e-12, key = "node")
  # Cluster labels are arbitrary: the sizes and shares are compared largest first (rank 1 is the largest cluster).
  sizes <- summary(found)
  ordered <- order(-sizes$size)
  reference_sizes <- sort(lengths(members), decreasing = TRUE)
  expect_reference(data.frame(rank = seq_along(ordered), size = sizes$size[ordered], share = sizes$share[ordered]),
                   data.frame(rank = seq_along(reference_sizes), size = reference_sizes,
                              share = reference_sizes / sum(reference_sizes)),
                   reference = "hypernetx", covers = c("hypergraph_cluster()", "summary.net_hypergraph_cluster()"),
                   quantity = "cluster sizes and shares, largest first", tolerance = 1e-12, key = "rank")
})

test_that("clique-motif centrality equals igraph's eigenvector centrality of the clique expansion, unit length", {
  skip_without_reference("igraph")
  invisible(lapply(inputs, \(input) {
    memberships <- connected_memberships(input, weighted = FALSE)
    incidence <- hg_incidence(memberships) * 1
    graph <- igraph::graph_from_adjacency_matrix(hg_offdiag(tcrossprod(incidence)), mode = "undirected", weighted = TRUE)
    vector <- igraph::eigen_centrality(graph)$vector
    found <- hypergraph_centrality(bipartite_groups(memberships, player = "person", group = "group"), type = "clique")
    # igraph scales the Perron vector to a maximum of 1; hypergraph_centrality() to unit Euclidean length.
    expect_reference(found$clique, vector / sqrt(sum(vector^2)), reference = "igraph", covers = "hypergraph_centrality()",
                     quantity = sprintf("clique-motif centrality, %d people", input$n_people), tolerance = 1e-8,
                     note = "igraph's ARPACK solver converges to about 1e-10")
  }))
})
