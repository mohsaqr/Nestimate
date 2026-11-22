# build_simplicial(), betti_numbers(), euler_characteristic(), simplicial_degree(), simplicial_features(),
# verify_simplicial(), q_analysis(): the clique complex of a thresholded network and its topology. References: the
# cliques found by testing every subset of nodes with combn(), the Betti numbers as ranks of the boundary
# matrices over the rationals, and Atkin's (1974) q-connectivity, all written from the definitions (build_simplicial()
# finds its cliques with igraph::cliques(), so igraph would only check the hand-over). New, apart from
# verify_simplicial(), which local_testing_and_equivalence/test-simplicial-equiv-extracted.R compared with igraph.
# The f-vector is read from the complex's field `f_vector`, and verify_simplicial() returns a list (API gap).

networks <- lapply(c(731L, 732L, 733L), \(seed) {
  set.seed(seed)
  weights <- matrix(stats::runif(49), 7L, 7L)
  weights <- (weights + t(weights)) / 2
  diag(weights) <- 0
  dimnames(weights) <- list(LETTERS[1:7], LETTERS[1:7])
  weights
})
threshold <- 0.45

test_that("the f-vector, Betti numbers and Euler characteristic equal those of the clique complex", {
  invisible(lapply(seq_along(networks), \(i) {
    adjacency <- (networks[[i]] > threshold) * 1
    cliques <- hg_cliques(adjacency)
    f_vector <- as.vector(table(lengths(cliques)))
    betti <- hg_betti(cliques)
    complex <- build_simplicial(networks[[i]], threshold = threshold)
    reference <- "independent: clique complex by combn() and boundary ranks over the rationals"
    expect_reference(unname(complex$f_vector), f_vector, reference = reference, covers = "build_simplicial()",
                     quantity = sprintf("f-vector, network %d", i), tolerance = 0)
    expect_reference(unname(betti_numbers(complex)), betti, reference = reference, covers = "betti_numbers()",
                     quantity = sprintf("Betti numbers, network %d", i), tolerance = 0)
    expect_reference(euler_characteristic(complex), sum(f_vector * (-1)^(seq_along(f_vector) - 1L)), reference = reference,
                     covers = "euler_characteristic()", quantity = sprintf("Euler characteristic, network %d", i), tolerance = 0)
    checked <- verify_simplicial(networks[[i]], threshold = threshold)
    expect_reference(checked$n_simplices_ours, length(cliques), reference = reference, covers = "verify_simplicial()",
                     quantity = sprintf("number of simplices, network %d", i), tolerance = 0)
    expect_reference(unname(checked$betti), betti, reference = reference, covers = "verify_simplicial()",
                     quantity = sprintf("Betti numbers by verify_simplicial(), network %d", i), tolerance = 0)
  }))
})

test_that("simplicial degrees and the feature table equal counts over the clique complex", {
  invisible(lapply(seq_along(networks), \(i) {
    adjacency <- (networks[[i]] > threshold) * 1
    cliques <- hg_cliques(adjacency)
    top <- max(lengths(cliques)) - 1L
    degrees <- as.data.frame(lapply(stats::setNames(0:top, paste0("d", 0:top)), \(k) {
      vapply(rownames(adjacency), \(node) sum(vapply(cliques, \(s) length(s) == k + 1L && node %in% s, logical(1))), numeric(1))
    }))
    # total: the simplices of dimension 1 and up that contain the node (d1 onwards, as documented)
    expected <- cbind(node = rownames(adjacency), degrees, total = rowSums(degrees) - degrees$d0)
    reference <- "independent: clique complex by combn()"
    expect_reference(simplicial_degree(build_simplicial(networks[[i]], threshold = threshold)), expected,
                     reference = reference, covers = "simplicial_degree()",
                     quantity = sprintf("simplicial degree, network %d", i), tolerance = 0, key = "node")
    betti <- hg_betti(cliques)
    features <- data.frame(threshold = threshold, n_edges = sum(lengths(cliques) == 2L), b0 = betti[[1L]],
                           b1 = if (length(betti) > 1L) betti[[2L]] else 0L,
                           euler = sum(table(lengths(cliques)) * (-1)^(seq_len(top + 1L) - 1L)))
    expect_reference(simplicial_features(networks[[i]], threshold = threshold), features, reference = reference,
                     covers = "simplicial_features()", quantity = sprintf("features, network %d", i), tolerance = 0,
                     key = "threshold")
  }))
})

test_that("the q-vector equals the number of q-connected components of the maximal simplices (Atkin 1974)", {
  invisible(lapply(seq_along(networks), \(i) {
    adjacency <- (networks[[i]] > threshold) * 1
    cliques <- hg_cliques(adjacency)
    maximal <- Filter(\(s) !any(vapply(cliques, \(t) length(t) > length(s) && all(s %in% t), logical(1))), cliques)
    top <- max(lengths(maximal)) - 1L
    # Two maximal simplices are q-near when they share q + 1 vertices; Q_q counts the components among those of
    # dimension >= q.
    components <- \(q) {
      kept <- Filter(\(s) length(s) - 1L >= q, maximal)
      near <- outer(seq_along(kept), seq_along(kept), Vectorize(\(a, b) length(intersect(kept[[a]], kept[[b]])) >= q + 1L))
      reach <- Reduce(\(m, step) (m %*% near > 0) * 1, seq_along(kept), near * 1)
      length(unique(apply(reach, 1L, paste, collapse = "")))
    }
    q_vector <- vapply(top:0, components, integer(1))
    names(q_vector) <- paste0("q_", top:0)
    expect_reference(q_analysis(build_simplicial(networks[[i]], threshold = threshold))$q_vector, q_vector,
                     reference = "independent: Atkin's q-connectivity on the maximal cliques", covers = "q_analysis()",
                     quantity = sprintf("q-vector, network %d", i), tolerance = 0)
  }))
})

# Outside references: TDA (Fasy et al.; GUDHI's Rips filtration and persistence) on the distance matrix that is 1
# for an edge and 2 otherwise, so that the Vietoris-Rips complex at scale 1 is the clique complex (persistence by
# Dionysus: TDA refuses GUDHI above dimension 1 for an arbitrary distance); igraph for the
# connected components (0-connectivity).
out_rips <- function(adjacency) {
  distance <- ifelse(adjacency > 0, 1, 2)
  diag(distance) <- 0
  distance
}

# The simplices of the clique complex from TDA::ripsFiltration() at scale 1, as vectors of node names.
out_tda_simplices <- function(adjacency) {
  filtration <- TDA::ripsFiltration(out_rips(adjacency), maxdimension = nrow(adjacency) - 1L, maxscale = 1,
                                    dist = "arbitrary", library = "GUDHI")
  lapply(filtration$cmplx, \(simplex) rownames(adjacency)[simplex])
}

# Betti numbers at scale 1 of TDA::ripsDiag(): the classes born at or before 1 that die after it.
out_tda_betti <- function(adjacency, top) {
  diagram <- unclass(TDA::ripsDiag(out_rips(adjacency), maxdimension = top, maxscale = 3, dist = "arbitrary",
                                   library = "Dionysus")$diagram)
  vapply(0:top, \(k) sum(diagram[, "dimension"] == k & diagram[, "Birth"] <= 1 & diagram[, "Death"] > 1), numeric(1))
}

test_that("the f-vector, Betti numbers, Euler characteristic and degrees equal TDA's Rips complex at scale 1", {
  skip_without_reference("TDA")
  invisible(lapply(seq_along(networks), \(i) {
    adjacency <- (networks[[i]] > threshold) * 1
    simplices <- out_tda_simplices(adjacency)
    f_vector <- as.vector(table(lengths(simplices)))
    top <- length(f_vector) - 1L
    betti <- out_tda_betti(adjacency, top)
    complex <- build_simplicial(networks[[i]], threshold = threshold)
    expect_reference(unname(complex$f_vector), f_vector, reference = "TDA", covers = "build_simplicial()",
                     quantity = sprintf("f-vector, network %d, TDA ripsFiltration()", i), tolerance = 0)
    expect_reference(unname(betti_numbers(complex)), betti, reference = "TDA", covers = "betti_numbers()",
                     quantity = sprintf("Betti numbers, network %d, TDA ripsDiag()", i), tolerance = 0)
    # Euler-Poincare: the alternating sum of TDA's Betti numbers
    expect_reference(euler_characteristic(complex), sum(betti * (-1)^(0:top)), reference = "TDA",
                     covers = "euler_characteristic()", quantity = sprintf("Euler characteristic, network %d, TDA", i),
                     tolerance = 0)
    checked <- verify_simplicial(networks[[i]], threshold = threshold)
    expect_reference(checked$n_simplices_ours, length(simplices), reference = "TDA", covers = "verify_simplicial()",
                     quantity = sprintf("number of simplices, network %d, TDA", i), tolerance = 0)
    expect_reference(unname(checked$betti), betti, reference = "TDA", covers = "verify_simplicial()",
                     quantity = sprintf("Betti numbers by verify_simplicial(), network %d, TDA", i), tolerance = 0)
    degrees <- as.data.frame(lapply(stats::setNames(0:top, paste0("d", 0:top)), \(k) {
      vapply(rownames(adjacency), \(node) sum(vapply(simplices, \(s) length(s) == k + 1L && node %in% s, logical(1))),
             numeric(1))
    }))
    expect_reference(simplicial_degree(complex), cbind(node = rownames(adjacency), degrees), reference = "TDA",
                     covers = "simplicial_degree()", quantity = sprintf("simplicial degree, network %d, TDA", i),
                     tolerance = 0, key = "node")
    graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
    expect_reference(simplicial_features(networks[[i]], threshold = threshold),
                     data.frame(threshold = threshold, n_edges = igraph::ecount(graph), b0 = igraph::count_components(graph),
                                b1 = betti[[2L]], euler = sum(betti * (-1)^(0:top))),
                     reference = "TDA", covers = "simplicial_features()",
                     quantity = sprintf("features, network %d, TDA and igraph", i), tolerance = 0, key = "threshold")
    # Q_0 counts the classes of 0-near maximal simplices (sharing a vertex): the connected components of the graph.
    expect_reference(unname(q_analysis(complex)$q_vector[["q_0"]]), igraph::count_components(graph), reference = "igraph",
                     covers = "q_analysis()", quantity = sprintf("q_0, network %d, igraph components", i), tolerance = 0)
  }))
})

test_that("simplicial measures of sparser and denser thresholds equal TDA", {
  skip_without_reference("TDA")
  invisible(lapply(c(0.3, 0.6), \(cut) {
    adjacency <- (networks[[1L]] > cut) * 1
    simplices <- out_tda_simplices(adjacency)
    f_vector <- as.vector(table(lengths(simplices)))
    complex <- build_simplicial(networks[[1L]], threshold = cut)
    expect_reference(unname(complex$f_vector), f_vector, reference = "TDA", covers = "build_simplicial()",
                     quantity = sprintf("f-vector, network 1, threshold %.1f, TDA", cut), tolerance = 0)
    expect_reference(unname(betti_numbers(complex)), out_tda_betti(adjacency, length(f_vector) - 1L), reference = "TDA",
                     covers = "betti_numbers()", quantity = sprintf("Betti numbers, network 1, threshold %.1f, TDA", cut),
                     tolerance = 0)
  }))
})
