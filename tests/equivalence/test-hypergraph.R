# bipartite_groups(), summary.net_hypergraph(), clique_expansion(), hypergraph_measures(): the hypergraph of a
# membership table and its incidence statistics. References: the incidence algebra written from the definitions
# (Battiston et al. 2020) and HyperG's edge sizes and edge adjacency. Re-expressed from
# local_testing_and_equivalence/test-equiv-hypergraph.R (40 random hypergraphs compared by maxima); here three
# fixed membership tables. hypergraph_measures() returns a list and bipartite_groups() keeps its incidence matrix
# as a field: both are read with `$` (API gap: no tidy accessor).

inputs <- list(list(n_people = 10L, n_groups = 8L, max_size = 3L, seed = 701L),
               list(n_people = 15L, n_groups = 15L, max_size = 4L, seed = 702L),
               list(n_people = 25L, n_groups = 25L, max_size = 4L, seed = 703L))

test_that("the incidence matrix and node degrees equal a table of the memberships", {
  invisible(lapply(inputs, \(input) {
    memberships <- hg_memberships(input$n_people, input$n_groups, input$max_size, input$seed)
    hypergraph <- bipartite_groups(memberships, player = "person", group = "group")
    incidence <- hg_incidence(memberships)
    label <- sprintf("%d people, %d groups", input$n_people, input$n_groups)
    expect_reference(hypergraph$incidence, incidence, reference = "independent: table() of the memberships, base R",
                     covers = "bipartite_groups()", quantity = sprintf("incidence, %s", label), tolerance = 0)
    expect_reference(summary(hypergraph), data.frame(node = rownames(incidence), degree = rowSums(incidence)),
                     reference = "independent: table() of the memberships, base R",
                     covers = c("bipartite_groups()", "summary.net_hypergraph()"),
                     quantity = sprintf("node degree, %s", label), tolerance = 0, key = "node")
  }))
})

test_that("the clique expansion equals H H' with a zero diagonal", {
  invisible(lapply(inputs, \(input) {
    memberships <- hg_memberships(input$n_people, input$n_groups, input$max_size, input$seed)
    incidence <- hg_incidence(memberships)
    expected <- hg_offdiag(tcrossprod(incidence) * 1)
    expect_reference(extract_transition_matrix(clique_expansion(bipartite_groups(memberships, player = "person", group = "group"))),
                     expected, reference = "independent: tcrossprod() of the incidence matrix, base R",
                     covers = "clique_expansion()",
                     quantity = sprintf("shared groups, %d people", input$n_people), tolerance = 0)
  }))
})

test_that("hyperdegree, node strength, edge sizes, co-degree, overlaps and density equal the incidence algebra", {
  invisible(lapply(inputs, \(input) {
    memberships <- hg_memberships(input$n_people, input$n_groups, input$max_size, input$seed)
    measures <- hypergraph_measures(bipartite_groups(memberships, player = "person", group = "group"))
    incidence <- hg_incidence(memberships) * 1
    sizes <- colSums(incidence)
    overlap <- hg_offdiag(crossprod(incidence))
    union <- outer(sizes, sizes, `+`) - overlap
    label <- sprintf("%d people, %d groups", input$n_people, input$n_groups)
    reference <- "independent: incidence algebra, base R"
    expect_reference(measures$hyperdegree, rowSums(incidence), reference = reference, covers = "hypergraph_measures()",
                     quantity = sprintf("hyperdegree, %s", label), tolerance = 0)
    # node strength: the summed sizes of the groups a person belongs to
    expect_reference(measures$node_strength, as.vector(incidence %*% sizes), reference = reference,
                     covers = "hypergraph_measures()", quantity = sprintf("node strength, %s", label), tolerance = 0)
    expect_reference(measures$edge_sizes, unname(sizes), reference = reference, covers = "hypergraph_measures()",
                     quantity = sprintf("edge size, %s", label), tolerance = 0)
    expect_reference(measures$co_degree, hg_offdiag(tcrossprod(incidence)), reference = reference,
                     covers = "hypergraph_measures()", quantity = sprintf("co-degree, %s", label), tolerance = 0)
    expect_reference(measures$jaccard, hg_offdiag(overlap / union), reference = reference,
                     covers = "hypergraph_measures()", quantity = sprintf("Jaccard overlap of groups, %s", label),
                     tolerance = 1e-12)
    expect_reference(measures$overlap_coefficient, hg_offdiag(overlap / outer(sizes, sizes, pmin)), reference = reference,
                     covers = "hypergraph_measures()", quantity = sprintf("overlap coefficient of groups, %s", label),
                     tolerance = 1e-12)
    expect_reference(measures$density, sum(incidence) / length(incidence), reference = reference,
                     covers = "hypergraph_measures()", quantity = sprintf("density, %s", label), tolerance = 1e-12)
  }))
})

test_that("edge sizes and edge adjacency equal HyperG", {
  skip_without_reference("HyperG")
  invisible(lapply(inputs, \(input) {
    memberships <- hg_memberships(input$n_people, input$n_groups, input$max_size, input$seed)
    measures <- hypergraph_measures(bipartite_groups(memberships, player = "person", group = "group"))
    reference <- HyperG::hypergraph_from_incidence_matrix(hg_incidence(memberships))
    label <- sprintf("%d people, %d groups", input$n_people, input$n_groups)
    expect_reference(measures$edge_sizes, as.vector(HyperG::hdegree(reference)), reference = "HyperG",
                     covers = "hypergraph_measures()", quantity = sprintf("edge size, %s, HyperG", label), tolerance = 0)
    expect_reference(unname(measures$edge_pairwise_overlap), unname(as.matrix(HyperG::hadjacency(reference))),
                     reference = "HyperG", covers = "hypergraph_measures()",
                     quantity = sprintf("edge pairwise overlap, %s", label), tolerance = 0)
  }))
})

# Outside reference: igraph's bipartite graph of the memberships (people and groups), its degrees and its
# projection onto the people, with the number of shared groups as the edge multiplicity.
test_that("the incidence matrix, node degrees and clique expansion equal igraph's bipartite graph", {
  skip_without_reference("igraph")
  invisible(lapply(inputs, \(input) {
    memberships <- hg_memberships(input$n_people, input$n_groups, input$max_size, input$seed)
    hypergraph <- bipartite_groups(memberships, player = "person", group = "group")
    pairs <- unique(subset(memberships, select = c(person, group)))
    graph <- igraph::graph_from_data_frame(pairs, directed = FALSE)
    igraph::V(graph)$type <- igraph::V(graph)$name %in% pairs$group
    people <- sort(unique(pairs$person))
    groups <- sort(unique(pairs$group))
    biadjacency <- unclass(as.matrix(igraph::as_biadjacency_matrix(graph)))[people, groups]
    label <- sprintf("%d people, %d groups, igraph", input$n_people, input$n_groups)
    expect_reference(hypergraph$incidence[people, groups], biadjacency, reference = "igraph",
                     covers = "bipartite_groups()", quantity = sprintf("incidence, %s", label), tolerance = 0)
    degree <- igraph::degree(graph)
    expect_reference(summary(hypergraph), data.frame(node = people, degree = unname(degree[people])),
                     reference = "igraph", covers = c("bipartite_groups()", "summary.net_hypergraph()"),
                     quantity = sprintf("node degree, %s", label), tolerance = 0, key = "node")
    projection <- igraph::bipartite_projection(graph, multiplicity = TRUE, which = "false")
    shared <- unclass(as.matrix(igraph::as_adjacency_matrix(projection, attr = "weight")))
    expected <- matrix(0, length(people), length(people), dimnames = list(people, people))
    expected[rownames(shared), colnames(shared)] <- shared
    expect_reference(extract_transition_matrix(clique_expansion(hypergraph))[people, people], expected,
                     reference = "igraph", covers = "clique_expansion()", quantity = sprintf("shared groups, %s", label),
                     tolerance = 0)
  }))
})
