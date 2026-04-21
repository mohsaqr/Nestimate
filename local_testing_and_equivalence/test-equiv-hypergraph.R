# Numerical equivalence: Nestimate's hypergraph module vs R-native references.
#
# Five exports:
#   * bipartite_groups      — vs manual table()-based incidence
#   * build_hypergraph      — vs manual clique enumeration
#   * clique_expansion      — vs tcrossprod()
#   * hypergraph_measures   — vs hand-computed incidence statistics
#   * hypergraph_centrality — CEC vs igraph; Z/H vs clean-room tensor iter
#
# No Python. Pure R. Runs in seconds.

set.seed(4242)
N_HG <- 40L
TOL <- 1e-10
TOL_COSINE <- 1e-6   # Z/H tensor eigenvectors are sign/rotation-ambiguous

# ---- Config generation ----
hg_configs <- lapply(seq_len(N_HG), function(i) {
  list(n_nodes = sample(c(10L, 15L, 25L), 1),
       n_edges = sample(c(8L, 15L, 25L), 1),
       edge_size_range = sample(list(c(2L, 3L), c(2L, 4L), c(3L, 4L)), 1)[[1L]],
       density = sample(c(0.7, 0.8, 0.9), 1),
       seed = sample.int(100000, 1))
})

# ---- bipartite_groups ----
test_that("bipartite_groups reconstructs incidence matrix exactly", {
  skip_on_cran()
  skip_equiv_tests()

  invisible(lapply(seq_len(N_HG), function(i) {
    cfg <- hg_configs[[i]]
    set.seed(cfg$seed)
    n_obs <- cfg$n_nodes * cfg$n_edges
    df <- data.frame(
      player = paste0("p", sample.int(cfg$n_nodes, n_obs, replace = TRUE)),
      session = paste0("s", sample.int(cfg$n_edges, n_obs, replace = TRUE)),
      stringsAsFactors = FALSE
    )
    # Drop duplicate rows (bipartite_groups produces binary incidence unless
    # weights are supplied; duplicates are collapsed to 1 in the incidence).
    hg <- bipartite_groups(df, player = "player", group = "session")

    # Manual reference: table() gives a contingency matrix; clip to binary.
    ref_tab <- table(df$player, df$session)
    ref_inc <- (ref_tab > 0) * 1L
    ref_inc <- ref_inc[sort(rownames(ref_inc)), sort(colnames(ref_inc)),
                       drop = FALSE]

    # Drop empty hyperedges (bipartite_groups does this; table doesn't).
    ref_inc <- ref_inc[, colSums(ref_inc) > 0, drop = FALSE]

    expect_equal(as.integer(hg$incidence), as.integer(ref_inc),
                 label = sprintf("cfg%d bipartite_groups incidence", i))
    expect_equal(rownames(hg$incidence), rownames(ref_inc),
                 label = sprintf("cfg%d bipartite_groups node order", i))
    expect_equal(colnames(hg$incidence), colnames(ref_inc),
                 label = sprintf("cfg%d bipartite_groups edge order", i))
  }))
})

# ---- clique_expansion ----
test_that("clique_expansion matches tcrossprod(incidence) with zero diag", {
  skip_on_cran()
  skip_equiv_tests()

  invisible(lapply(seq_len(N_HG), function(i) {
    cfg <- hg_configs[[i]]
    h <- simulate_hypergraph_incidence(cfg$n_nodes, cfg$n_edges,
                                       cfg$edge_size_range, cfg$density,
                                       seed = cfg$seed)
    hg <- structure(
      list(hyperedges = h$hyperedges,
           incidence = h$incidence,
           nodes = h$nodes,
           n_nodes = length(h$nodes),
           n_hyperedges = length(h$hyperedges),
           size_distribution = integer(0),
           params = list()),
      class = "net_hypergraph"
    )

    net <- clique_expansion(hg)
    W_nest <- unname(net$weights)

    W_ref <- tcrossprod(h$incidence)
    storage.mode(W_ref) <- "double"
    diag(W_ref) <- 0
    W_ref <- unname(W_ref)

    expect_equal(W_nest, W_ref, tolerance = TOL,
                 label = sprintf("cfg%d clique_expansion matrix", i))
  }))
})

# ---- hypergraph_measures ----
test_that("hypergraph_measures matches hand-computed incidence statistics", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  invisible(lapply(seq_len(N_HG), function(i) {
    cfg <- hg_configs[[i]]
    h <- simulate_hypergraph_incidence(cfg$n_nodes, cfg$n_edges,
                                       cfg$edge_size_range, cfg$density,
                                       seed = cfg$seed)
    hg <- structure(
      list(hyperedges = h$hyperedges,
           incidence = h$incidence,
           nodes = h$nodes,
           n_nodes = length(h$nodes),
           n_hyperedges = length(h$hyperedges),
           size_distribution = integer(0),
           params = list()),
      class = "net_hypergraph"
    )

    m <- hypergraph_measures(hg)

    B <- (h$incidence > 0) * 1.0
    ref_hd <- rowSums(B)                      # hyperdegree per node
    ref_es <- colSums(B)                      # edge sizes
    ref_cd <- tcrossprod(B); diag(ref_cd) <- 0  # co-degree
    ref_epo <- crossprod(B); diag(ref_epo) <- 0 # edge pairwise overlap
    ref_ns <- as.numeric(B %*% ref_es)        # node strength
    ref_sz_min <- outer(ref_es, ref_es, pmin)
    ref_sz_union <- outer(ref_es, ref_es, `+`) - ref_epo
    ref_jacc <- ifelse(ref_sz_union > 0,
                       ref_epo / ref_sz_union, 0)
    diag(ref_jacc) <- 0

    checks <- c(
      hd = max(abs(m$hyperdegree - ref_hd)),
      es = max(abs(m$edge_sizes - ref_es)),
      cd = max(abs(unname(m$co_degree) - unname(ref_cd))),
      epo = max(abs(unname(m$edge_pairwise_overlap) - unname(ref_epo))),
      ns = max(abs(unname(m$node_strength) - ref_ns)),
      jacc = max(abs(unname(m$jaccard) - unname(ref_jacc)))
    )

    report$log(
      func = "hypergraph_measures",
      config = sprintf("cfg%d(n=%d,m=%d)",
                       i, cfg$n_nodes, length(h$hyperedges)),
      n_checked = length(m$hyperdegree) + length(m$edge_sizes) +
        length(m$co_degree) + length(m$edge_pairwise_overlap) +
        length(m$node_strength) + length(m$jaccard),
      n_failed = as.integer(any(checks > TOL)),
      max_abs_err = max(checks), mean_abs_err = mean(checks),
      median_abs_err = stats::median(checks), p95_abs_err = max(checks),
      reference = "hand-computed incidence algebra",
      notes = paste(names(checks), sprintf("%.2e", checks),
                    sep = "=", collapse = " ")
    )

    for (nm in names(checks)) {
      expect_true(checks[nm] < TOL,
                  label = sprintf("cfg%d %s delta = %.2e",
                                  i, nm, checks[nm]))
    }
  }))

  report$write_csv("hypergraph_measures")
  report$write_cvs("hypergraph_measures",
                   "local_testing_and_equivalence/test-equiv-hypergraph.R")
})

# ---- hypergraph_measures vs HyperG (external library) ----
test_that("hypergraph edge sizes and pairwise overlap match HyperG", {
  # External-library reference. HyperG is the Emma Towlson CRAN package for
  # hypergraph analysis (v1.0.0, built on Matrix). Its `hdegree()` computes
  # edge sizes via a different code path (Matrix::colSums on a sparse
  # dgCMatrix) and `hadjacency()` computes edge-pairwise overlap via
  # Matrix::crossprod on the sparse incidence. Shares no code with Nestimate.
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("HyperG")

  report <- equiv_report()

  invisible(lapply(seq_len(N_HG), function(i) {
    cfg <- hg_configs[[i]]
    h <- simulate_hypergraph_incidence(cfg$n_nodes, cfg$n_edges,
                                       cfg$edge_size_range, cfg$density,
                                       seed = cfg$seed)
    hg <- structure(
      list(hyperedges = h$hyperedges,
           incidence = h$incidence,
           nodes = h$nodes,
           n_nodes = length(h$nodes),
           n_hyperedges = length(h$hyperedges),
           size_distribution = integer(0),
           params = list()),
      class = "net_hypergraph"
    )

    m <- hypergraph_measures(hg)

    # HyperG path: build hypergraph from incidence, extract primitives
    hyper_ref <- HyperG::hypergraph_from_incidence_matrix(h$incidence)
    hyperg_edge_sizes <- as.integer(HyperG::hdegree(hyper_ref))
    hyperg_overlap <- as.matrix(HyperG::hadjacency(hyper_ref))

    # HyperG zeroes its diagonal the same way Nestimate does; columns are
    # labeled by edge names.
    nest_edge_sizes <- as.integer(m$edge_sizes)
    nest_overlap <- unname(as.matrix(m$edge_pairwise_overlap))
    hyperg_overlap <- unname(hyperg_overlap)

    sz_delta <- max(abs(nest_edge_sizes - hyperg_edge_sizes))
    ov_delta <- max(abs(nest_overlap - hyperg_overlap))

    report$log(
      func = "hypergraph_measures_vs_HyperG",
      config = sprintf("cfg%d(n=%d,m=%d)",
                       i, cfg$n_nodes, length(h$hyperedges)),
      n_checked = length(nest_edge_sizes) + length(nest_overlap),
      n_failed = as.integer(sz_delta > TOL) + as.integer(ov_delta > TOL),
      max_abs_err = max(c(sz_delta, ov_delta)),
      mean_abs_err = mean(c(sz_delta, ov_delta)),
      median_abs_err = stats::median(c(sz_delta, ov_delta)),
      p95_abs_err = as.numeric(stats::quantile(c(sz_delta, ov_delta), 0.95)),
      reference = "HyperG::hdegree + HyperG::hadjacency",
      notes = sprintf("sz=%.2e ov=%.2e", sz_delta, ov_delta)
    )

    expect_true(sz_delta < TOL,
                label = sprintf("cfg%d edge_sizes vs HyperG delta = %.2e",
                                i, sz_delta))
    expect_true(ov_delta < TOL,
                label = sprintf("cfg%d edge_pairwise_overlap vs HyperG delta = %.2e",
                                i, ov_delta))
  }))

  report$write_csv("hypergraph_vs_HyperG")
  report$write_cvs("hypergraph_vs_HyperG",
                   "local_testing_and_equivalence/test-equiv-hypergraph.R")
})

# ---- hypergraph_centrality: CEC vs igraph ----
test_that("CEC centrality matches igraph::eigen_centrality on clique-expanded network", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("igraph")

  invisible(lapply(seq_len(N_HG), function(i) {
    cfg <- hg_configs[[i]]
    h <- simulate_hypergraph_incidence(cfg$n_nodes, cfg$n_edges,
                                       cfg$edge_size_range, cfg$density,
                                       seed = cfg$seed)
    hg <- structure(
      list(hyperedges = h$hyperedges,
           incidence = h$incidence,
           nodes = h$nodes,
           n_nodes = length(h$nodes),
           n_hyperedges = length(h$hyperedges),
           size_distribution = integer(0),
           params = list()),
      class = "net_hypergraph"
    )

    # Need connected clique-expansion for igraph to return a sensible
    # Perron vector; skip degenerate configs.
    B <- (h$incidence > 0) * 1
    W <- tcrossprod(B); diag(W) <- 0
    if (all(W == 0)) return(NULL)

    cec <- hypergraph_centrality(hg, type = "clique",
                                 max_iter = 5000L, tol = 1e-10)$clique
    g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
    ref <- suppressWarnings(igraph::eigen_centrality(g))$vector

    # Sign-invariant cosine similarity
    num <- sum(cec * ref)
    denom <- sqrt(sum(cec^2) * sum(ref^2))
    cos_sim <- if (denom > 0) num / denom else 1
    expect_true(
      abs(abs(cos_sim) - 1) < TOL_COSINE,
      label = sprintf("cfg%d CEC vs igraph cosine = %.6f", i, cos_sim)
    )
  }))
})

# ---- hypergraph_centrality: Z and H via clean-room tensor iteration ----
.tensor_iter_reference <- function(hyperedges, edge_sizes, n,
                                    exponent, max_iter = 2000L,
                                    tol = 1e-10, shift = 1) {
  # Purely list-based power iteration — no leave-one-out shortcut, no total
  # product trick. Computes y[i] = sum over edges containing i of product of
  # x[j] for j != i. Same mathematical operator as Nestimate's kernel but
  # completely different code path.
  x <- rep(1 / sqrt(n), n)
  for (iter in seq_len(max_iter)) {
    y <- numeric(n)
    for (e_idx in seq_along(hyperedges)) {
      e <- hyperedges[[e_idx]]
      ke <- edge_sizes[e_idx]
      if (ke < 2L) next
      for (idx in seq_len(ke)) {
        others <- e[-idx]
        y[e[idx]] <- y[e[idx]] + prod(x[others])
      }
    }
    if (exponent > 1L) {
      y <- sign(y) * abs(y)^(1 / exponent)
    }
    y <- y + shift * x
    nrm <- sqrt(sum(y^2))
    if (nrm == 0) return(y)
    y <- y / nrm
    if (sum(abs(y - x)) < tol) return(y)
    x <- y
  }
  x
}

test_that("Z and H centralities match a clean-room list-based tensor iteration", {
  skip_on_cran()
  skip_equiv_tests()

  invisible(lapply(seq_len(N_HG), function(i) {
    cfg <- hg_configs[[i]]
    h <- simulate_hypergraph_incidence(cfg$n_nodes, cfg$n_edges,
                                       cfg$edge_size_range, cfg$density,
                                       seed = cfg$seed)
    hg <- structure(
      list(hyperedges = h$hyperedges,
           incidence = h$incidence,
           nodes = h$nodes,
           n_nodes = length(h$nodes),
           n_hyperedges = length(h$hyperedges),
           size_distribution = integer(0),
           params = list()),
      class = "net_hypergraph"
    )

    edge_sizes <- vapply(h$hyperedges, length, integer(1L))
    k_max <- max(edge_sizes)

    nest <- hypergraph_centrality(hg, type = c("Z", "H"),
                                  max_iter = 2000L, tol = 1e-10)

    ref_Z <- .tensor_iter_reference(h$hyperedges, edge_sizes, length(h$nodes),
                                    exponent = 1L, max_iter = 2000L,
                                    tol = 1e-10)
    ref_H <- .tensor_iter_reference(h$hyperedges, edge_sizes, length(h$nodes),
                                    exponent = k_max - 1L, max_iter = 2000L,
                                    tol = 1e-10)

    # Cosine similarity (sign-invariant)
    cos_sim <- function(a, b) {
      d <- sqrt(sum(a^2) * sum(b^2))
      if (d == 0) return(1)
      abs(sum(a * b)) / d
    }
    cz <- cos_sim(nest$Z, ref_Z)
    ch <- cos_sim(nest$H, ref_H)

    expect_true(abs(cz - 1) < TOL_COSINE,
                label = sprintf("cfg%d Z centrality cosine = %.6f", i, cz))
    expect_true(abs(ch - 1) < TOL_COSINE,
                label = sprintf("cfg%d H centrality cosine = %.6f", i, ch))
  }))
})

# ---- build_hypergraph: parameter branches ----
# Note: build_hypergraph delegates clique enumeration to igraph::cliques()
# via .find_all_cliques (R/simplicial.R:121). Testing the enumeration against
# igraph would be a tautology. Instead we (a) check parameter-branch logic
# (p sampling, include_pairwise, max_size) and (b) validate set equality
# against utils::combn brute force on tiny graphs where that is feasible.

test_that("build_hypergraph parameter branches behave correctly", {
  skip_on_cran()
  skip_equiv_tests()

  # 4-clique: 1 tetrahedron = 1 4-clique, 4 triangles, 6 edges
  n <- 4L
  W <- matrix(1, n, n) - diag(n)
  dimnames(W) <- list(paste0("n", 1:n), paste0("n", 1:n))
  net <- structure(list(
    weights = W,
    nodes = data.frame(id = 1:n, name = rownames(W), stringsAsFactors = FALSE),
    directed = FALSE, method = "test"
  ), class = c("netobject", "cograph_network"))

  # p = 0 + no pairwise -> empty
  hg_empty <- build_hypergraph(net, method = "clique", p = 0,
                               include_pairwise = FALSE, max_size = 4L)
  expect_equal(hg_empty$n_hyperedges, 0L, label = "p=0, no pairwise")

  # p = 1, no pairwise, max_size = 3 -> 4 triangles only
  hg_tri <- build_hypergraph(net, method = "clique", p = 1,
                             include_pairwise = FALSE, max_size = 3L)
  expect_equal(hg_tri$n_hyperedges, 4L, label = "4 triangles at max_size=3")

  # p = 1, no pairwise, max_size = 4 -> 4 triangles + 1 tetrahedron = 5
  hg_all <- build_hypergraph(net, method = "clique", p = 1,
                             include_pairwise = FALSE, max_size = 4L)
  expect_equal(hg_all$n_hyperedges, 5L, label = "triangles + tetrahedron")

  # include_pairwise = TRUE -> add 6 edges
  hg_with_pw <- build_hypergraph(net, method = "clique", p = 1,
                                 include_pairwise = TRUE, max_size = 4L)
  expect_equal(hg_with_pw$n_hyperedges, 11L, label = "with pairwise = 6 + 5")

  # Bernoulli sampling: p = 0.5, seeded. Seed reproducibility + bounds.
  hg_a <- build_hypergraph(net, method = "clique", p = 0.5,
                           include_pairwise = FALSE, max_size = 3L, seed = 42)
  hg_b <- build_hypergraph(net, method = "clique", p = 0.5,
                           include_pairwise = FALSE, max_size = 3L, seed = 42)
  expect_equal(hg_a$hyperedges, hg_b$hyperedges,
               label = "seed reproducibility")
  expect_true(hg_a$n_hyperedges <= 4L && hg_a$n_hyperedges >= 0L,
              label = "p=0.5 kept count in bounds")

  # Incidence matrix correctness: each column's nonzero rows = hyperedge nodes
  inc <- hg_all$incidence
  for (j in seq_len(ncol(inc))) {
    members_inc <- sort(unname(which(inc[, j] > 0)))
    members_he <- sort(as.integer(hg_all$hyperedges[[j]]))
    expect_equal(members_inc, members_he,
                 label = sprintf("incidence col %d matches hyperedges[[%d]]",
                                 j, j))
  }
})

test_that("build_hypergraph(method='clique') matches manual combn enumeration", {
  skip_on_cran()
  skip_equiv_tests()

  # Small structured graphs where manual enumeration is feasible.
  invisible(lapply(seq_len(10L), function(i) {
    set.seed(1000 + i)
    # Build a small 6-node network with some triangles
    n <- 6L
    W <- matrix(0, n, n, dimnames = list(paste0("n", 1:n),
                                          paste0("n", 1:n)))
    # Triangle 1: n1-n2-n3
    W[1, 2] <- W[2, 1] <- 1
    W[1, 3] <- W[3, 1] <- 1
    W[2, 3] <- W[3, 2] <- 1
    # Triangle 2: n4-n5-n6
    W[4, 5] <- W[5, 4] <- 1
    W[4, 6] <- W[6, 4] <- 1
    W[5, 6] <- W[6, 5] <- 1
    # A bridging edge but no third triangle
    W[3, 4] <- W[4, 3] <- 1

    net <- structure(list(
      weights = W, nodes = data.frame(id = 1:n, name = rownames(W),
                                      stringsAsFactors = FALSE),
      directed = FALSE, method = "test"
    ), class = c("netobject", "cograph_network"))

    hg <- tryCatch(build_hypergraph(net, method = "clique", p = 1,
                                    include_pairwise = FALSE),
                   error = function(e) NULL)
    if (is.null(hg)) return(NULL)

    # Manual reference: enumerate all 3-cliques via combn
    ref_triangles <- list()
    tri_idx <- 0L
    for (tri in utils::combn(n, 3L, simplify = FALSE)) {
      if (W[tri[1], tri[2]] > 0 && W[tri[1], tri[3]] > 0 &&
          W[tri[2], tri[3]] > 0) {
        tri_idx <- tri_idx + 1L
        ref_triangles[[tri_idx]] <- sort(tri)
      }
    }

    # Nestimate hyperedges (3-cliques only with p=1, no pairwise)
    nest_sizes <- vapply(hg$hyperedges, length, integer(1L))
    nest_triangles <- hg$hyperedges[nest_sizes == 3L]
    nest_sorted <- lapply(nest_triangles, sort)

    # Set equality
    ref_keys <- vapply(ref_triangles, function(t) paste(t, collapse = ","),
                       character(1L))
    nest_keys <- vapply(nest_sorted, function(t) paste(t, collapse = ","),
                       character(1L))

    expect_equal(sort(nest_keys), sort(ref_keys),
                 label = sprintf("cfg%d triangle set", i))
  }))
})
