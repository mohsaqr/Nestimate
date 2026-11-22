# Layer 6, adversarial data: hypergraphs of membership tables with hyperedges of mixed sizes (1, 2 and 3), a member
# whose weight is 0 (not a member), a membership listed twice, an isolated player, a singleton hyperedge, two
# identical components (a degenerate leading eigenvalue), a star of pairs (a bipartite clique expansion, where power
# iteration oscillates), one hyperedge holding every node, and a connected 3-uniform hypergraph with a repeated
# hyperedge. Verbs: bipartite_groups(), summary.net_hypergraph(), clique_expansion(), hypergraph_measures() (every
# field) and hypergraph_centrality().
# References: the incidence algebra written from the definitions, base R (helper-hypergraph.R; density as documented:
# m / choose(n, k) for a k-uniform hypergraph, the mean share of nodes per hyperedge otherwise); the clique-motif
# centrality as the leading eigenvector of the WHOLE clique expansion by base::eigen() (with a repeated leading
# eigenvalue, the projection of the uniform vector onto the leading eigenspace, the documented rule) and igraph's
# eigen_centrality() where the expansion is connected; the Z- and H-eigenvector centralities through their
# eigen-equations (lambda x = f(x), lambda x^(k-1) = f(x), f(x)_i = sum over e containing i of the product of the
# other members), whose positive solution on a connected uniform hypergraph is unique.
# Tolerances, fixed before comparing: counts exact; ratios 1e-12; eigenvectors 1e-8 (igraph's ARPACK and the
# tensor power iteration run to 1e-8).

hypergraph_sets <- list(
  mixed = data.frame(person = c("a", "b", "a", "b", "c", "c", "d", "e", "f", "g", "g"),
                     group = c("g1", "g1", "g2", "g2", "g2", "g3", "g3", "g4", "g4", "g5", "g5"),
                     w = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 1)),
  two_components = data.frame(person = c("a", "b", "c", "d"), group = c("g1", "g1", "g2", "g2"), w = 1),
  star = data.frame(person = c("a", "b", "a", "c", "a", "d", "a", "e"), group = rep(paste0("g", 1:4), each = 2L), w = 1),
  one_edge = data.frame(person = c("a", "b", "c", "d", "e"), group = "g1", w = 1),
  uniform = data.frame(person = c("a", "b", "c", "b", "c", "d", "c", "d", "e", "d", "e", "a", "a", "b", "c"),
                       group = rep(paste0("g", 1:5), each = 3L), w = 1)
)

# The 0/1 incidence of a weighted membership table: a member has a positive summed weight in the group.
adv_hg_incidence <- function(memberships) {
  weights <- tapply(memberships$w, list(memberships$person, memberships$group), sum)
  weights[is.na(weights)] <- 0
  weights <- weights[sort(rownames(weights)), sort(colnames(weights)), drop = FALSE]
  (weights > 0) * 1
}

# Every field of hypergraph_measures() from the incidence algebra.
adv_hg_measures <- function(incidence) {
  sizes <- colSums(incidence)
  overlap <- hg_offdiag(crossprod(incidence))
  union <- outer(sizes, sizes, `+`) - overlap
  co_degree <- hg_offdiag(tcrossprod(incidence))
  shared <- overlap[upper.tri(overlap)]
  n <- nrow(incidence)
  uniform <- length(unique(sizes)) == 1L
  list(hyperdegree = rowSums(incidence), node_strength = stats::setNames(as.vector(incidence %*% sizes), rownames(incidence)),
       max_edge_size = stats::setNames(apply(incidence * rep(sizes, each = n), 1L, max), rownames(incidence)),
       co_degree = co_degree, edge_sizes = unname(sizes), edge_pairwise_overlap = overlap,
       overlap_coefficient = hg_offdiag(overlap / outer(sizes, sizes, pmin)), jaccard = hg_offdiag(overlap / union),
       # documented: m / choose(n, k) when every hyperedge has size k, else the mean share of nodes per hyperedge
       density = if (uniform) ncol(incidence) / choose(n, sizes[[1L]]) else sum(sizes) / (n * ncol(incidence)),
       avg_edge_size = mean(sizes),
       size_distribution = stats::setNames(as.vector(table(sizes)), paste0("size_", sort(unique(sizes)))),
       intersection_profile = if (length(shared)) stats::setNames(as.vector(table(shared)),
                                                                  paste0("overlap_", sort(unique(shared)))) else integer(),
       pairwise_participation = sum(co_degree[upper.tri(co_degree)] > 0) / choose(n, 2L),
       n_nodes = n, n_hyperedges = ncol(incidence))
}

# The clique-motif centrality: the leading eigenvector of the clique expansion, the uniform vector projected onto
# the leading eigenspace when the leading eigenvalue is repeated; unit length.
adv_clique_centrality <- function(incidence) {
  w <- hg_offdiag(tcrossprod(incidence))
  eig <- eigen(w, symmetric = TRUE)
  top <- abs(eig$values - eig$values[[1L]]) <= 1e-9 * max(1, abs(eig$values[[1L]]))
  space <- eig$vectors[, top, drop = FALSE]
  x <- as.vector(space %*% crossprod(space, rep(1, nrow(w))))
  stats::setNames(x / sqrt(sum(x^2)), rownames(incidence))
}

# The tensor map of a uniform hypergraph at x: f(x)_i = sum over hyperedges e containing i of prod_{j in e, j != i} x_j.
adv_tensor_map <- function(incidence, x) {
  members <- lapply(seq_len(ncol(incidence)), \(e) which(incidence[, e] > 0))
  contributions <- lapply(members, \(m) {
    out <- numeric(length(x))
    out[m] <- vapply(m, \(i) prod(x[setdiff(m, i)]), numeric(1))
    out
  })
  stats::setNames(Reduce(`+`, contributions), rownames(incidence))
}

test_that("incidence, degrees, the clique expansion and every hypergraph measure equal the incidence algebra", {
  invisible(lapply(names(hypergraph_sets), \(name) {
    memberships <- hypergraph_sets[[name]]
    hg <- bipartite_groups(memberships, player = "person", group = "group", weight = "w")
    binary <- adv_hg_incidence(memberships)
    expect_reference(summary(hg), data.frame(node = rownames(binary), degree = rowSums(binary)),
                     reference = "independent: positive summed memberships, base R",
                     covers = c("bipartite_groups()", "summary.net_hypergraph()"),
                     quantity = sprintf("node degree, %s", name), tolerance = 0, key = "node")
    expect_reference(extract_transition_matrix(clique_expansion(hg, weighted = FALSE)), hg_offdiag(tcrossprod(binary)),
                     reference = "independent: tcrossprod() of the 0/1 incidence, base R", covers = "clique_expansion()",
                     quantity = sprintf("shared groups, %s", name), tolerance = 0)
    found <- hypergraph_measures(hg)
    expected <- adv_hg_measures(binary)
    invisible(lapply(names(expected), \(field) {
      # a field with nothing in it (no pair of hyperedges) is asserted empty: there is no number to compare
      if (length(expected[[field]]) == 0L) return(expect_length(found[[field]], 0L))
      expect_reference(found[[field]], expected[[field]], reference = "independent: incidence algebra, base R",
                       covers = "hypergraph_measures()", quantity = sprintf("%s, %s", field, name), tolerance = 1e-12)
    }))
  }))
})

test_that("the clique-motif centrality equals the leading eigenspace of the whole clique expansion", {
  invisible(lapply(names(hypergraph_sets), \(name) {
    memberships <- hypergraph_sets[[name]]
    hg <- bipartite_groups(memberships, player = "person", group = "group", weight = "w")
    expect_reference(hypergraph_centrality(hg, type = "clique")$clique, adv_clique_centrality(adv_hg_incidence(memberships)),
                     reference = "independent: leading eigenspace of the clique expansion by eigen(), base R",
                     covers = "hypergraph_centrality()", quantity = sprintf("clique-motif centrality, %s", name),
                     tolerance = 1e-8)
  }))
})

test_that("the clique-motif centrality equals igraph where the clique expansion is connected", {
  skip_without_reference("igraph")
  invisible(lapply(c("star", "one_edge", "uniform"), \(name) {
    memberships <- hypergraph_sets[[name]]
    graph <- igraph::graph_from_adjacency_matrix(hg_offdiag(tcrossprod(adv_hg_incidence(memberships))),
                                                 mode = "undirected", weighted = TRUE)
    vector <- igraph::eigen_centrality(graph)$vector
    hg <- bipartite_groups(memberships, player = "person", group = "group", weight = "w")
    expect_reference(hypergraph_centrality(hg, type = "clique", normalize = FALSE)$clique, vector, reference = "igraph",
                     covers = "hypergraph_centrality()", quantity = sprintf("clique-motif centrality, max scale, %s", name),
                     tolerance = 1e-8, note = "normalize = FALSE scales the vector to a largest entry of 1, as igraph")
  }))
})

test_that("Z- and H-eigenvector centralities of a uniform hypergraph solve their eigen-equations", {
  memberships <- hypergraph_sets$uniform
  hg <- bipartite_groups(memberships, player = "person", group = "group", weight = "w")
  incidence <- adv_hg_incidence(memberships)
  found <- hypergraph_centrality(hg, type = c("Z", "H"))
  z_image <- adv_tensor_map(incidence, found$Z)
  expect_reference(found$Z, z_image / sqrt(sum(z_image^2)),
                   reference = "independent: the Z eigen-equation lambda x = f(x), base R",
                   covers = "hypergraph_centrality()", quantity = "Z-eigenvector centrality, 3-uniform", tolerance = 1e-8)
  h_image <- adv_tensor_map(incidence, found$H)^(1 / 2)
  expect_reference(found$H, h_image / sqrt(sum(h_image^2)),
                   reference = "independent: the H eigen-equation lambda x^(k - 1) = f(x), base R",
                   covers = "hypergraph_centrality()", quantity = "H-eigenvector centrality, 3-uniform", tolerance = 1e-8)
})

test_that("Z- and H-eigenvector centralities refuse hyperedges of mixed sizes with a classed error", {
  hg <- bipartite_groups(hypergraph_sets$mixed, player = "person", group = "group", weight = "w")
  expect_error(hypergraph_centrality(hg, type = "Z"), class = "nestimate_nonuniform_hypergraph")
  expect_error(hypergraph_centrality(hg, type = "H"), class = "nestimate_nonuniform_hypergraph")
})
