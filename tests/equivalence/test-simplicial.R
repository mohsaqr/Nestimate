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
