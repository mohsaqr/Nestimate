# build_hypergraph(method = "clique"): the hyperedges of a network as its 3-cliques (triangles). Reference: every
# triple of nodes tested from the definition with utils::combn(). build_hypergraph() finds the cliques with
# igraph::cliques(), so igraph is not an independent reference here. Re-expressed from
# local_testing_and_equivalence/test-equiv-hypergraph.R ("matches manual combn enumeration"). The incidence matrix
# is read from the result's field (API gap).

# Triangles of a symmetric 0/1 matrix, each written "a,b,c" (names sorted).
hg_triangles <- function(adjacency) {
  triples <- utils::combn(rownames(adjacency), 3L, simplify = FALSE)
  closed <- Filter(\(t) adjacency[t[1L], t[2L]] > 0 && adjacency[t[1L], t[3L]] > 0 && adjacency[t[2L], t[3L]] > 0, triples)
  vapply(closed, \(t) paste(sort(t), collapse = ","), character(1))
}

# The hyperedges of an incidence matrix, each written "a,b,c".
hg_hyperedges <- function(incidence) {
  apply(incidence, 2L, \(column) paste(sort(rownames(incidence)[column > 0]), collapse = ","))
}

networks <- list(
  two_triangles = rbind(c(1, 2), c(1, 3), c(2, 3), c(4, 5), c(4, 6), c(5, 6), c(3, 4)),
  random = local({
    set.seed(721)
    pairs <- t(utils::combn(8L, 2L))
    pairs[stats::runif(nrow(pairs)) < 0.45, , drop = FALSE]
  })
)

test_that("triangle hyperedges equal the closed triples found with combn()", {
  invisible(lapply(names(networks), \(name) {
    edges <- networks[[name]]
    n <- max(edges)
    adjacency <- matrix(0, n, n, dimnames = list(paste0("n", seq_len(n)), paste0("n", seq_len(n))))
    adjacency[edges] <- 1
    adjacency <- pmax(adjacency, t(adjacency))
    found <- build_hypergraph(adjacency, method = "clique", p = 1, include_pairwise = FALSE)
    triangles <- hg_triangles(adjacency)
    expect_reference(data.frame(triangle = sort(hg_hyperedges(found$incidence)), present = 1),
                     data.frame(triangle = sort(triangles), present = 1),
                     reference = "independent: closed triples by utils::combn()", covers = "build_hypergraph()",
                     quantity = sprintf("triangles, %s network", name), tolerance = 0, key = "triangle")
  }))
})

# Outside reference: igraph::triangles() lists the triangles of a graph with its own enumeration (not the clique
# search build_hypergraph() calls).
test_that("triangle hyperedges equal igraph::triangles()", {
  skip_without_reference("igraph")
  invisible(lapply(names(networks), \(name) {
    edges <- networks[[name]]
    n <- max(edges)
    adjacency <- matrix(0, n, n, dimnames = list(paste0("n", seq_len(n)), paste0("n", seq_len(n))))
    adjacency[edges] <- 1
    adjacency <- pmax(adjacency, t(adjacency))
    graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
    corners <- matrix(igraph::V(graph)$name[as.integer(igraph::triangles(graph))], ncol = 3L, byrow = TRUE)
    triangles <- apply(corners, 1L, \(t) paste(sort(t), collapse = ","))
    found <- build_hypergraph(adjacency, method = "clique", p = 1, include_pairwise = FALSE)
    expect_reference(data.frame(triangle = sort(hg_hyperedges(found$incidence)), present = 1),
                     data.frame(triangle = sort(triangles), present = 1), reference = "igraph",
                     covers = "build_hypergraph()", quantity = sprintf("triangles, %s network, igraph", name),
                     tolerance = 0, key = "triangle")
  }))
})
