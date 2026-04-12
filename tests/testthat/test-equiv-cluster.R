# ============================================================================
# Equivalence tests — cluster_data() & centrality()
# ============================================================================
# Run with: NESTIMATE_EQUIV_TESTS=true Rscript -e 'testthat::test_file("tests/testthat/test-equiv-cluster.R")'

TOL        <- 1e-10
TOL_EXACT  <- 1e-14
N_CLUST    <- 50L
N_CENT     <- 100L
N_UNDIR    <- 50L
N_BTW      <- 50L

report <- equiv_report()

# ---- Config generators ----

set.seed(5555)
clust_configs <- lapply(seq_len(N_CLUST), function(i) {
  list(n_actors   = sample(c(15L, 20L, 30L), 1L),
       n_states   = sample(3:5, 1L),
       seq_length = sample(c(15L, 20L, 25L), 1L),
       k          = sample(2:4, 1L),
       seed       = sample.int(100000L, 1L))
})

set.seed(6666)
cent_configs <- lapply(seq_len(N_CENT), function(i) {
  list(n_actors   = sample(c(10L, 20L), 1L),
       n_states   = sample(c(3L, 5L, 7L), 1L),
       seq_length = sample(c(10L, 20L), 1L),
       seed       = sample.int(100000L, 1L))
})

set.seed(7777)
undir_configs <- lapply(seq_len(N_UNDIR), function(i) {
  list(n      = sample(c(50L, 100L), 1L),
       p      = sample(c(4L, 5L, 6L), 1L),
       rho    = round(runif(1L, 0.1, 0.5), 2),
       seed   = sample.int(100000L, 1L))
})

set.seed(8888)
btw_configs <- lapply(seq_len(N_BTW), function(i) {
  list(n_actors   = sample(c(10L, 20L), 1L),
       n_states   = sample(c(3L, 5L), 1L),
       seq_length = sample(c(10L, 20L), 1L),
       seed       = sample.int(100000L, 1L))
})


# ============================================================================
# 1. cluster_data cluster assignments are valid (50 configs)
# ============================================================================
test_that("cluster_data produces valid clustering (50 configs)", {
  skip_on_cran()
  skip_equiv_tests()

  results <- vapply(seq_len(N_CLUST), function(i) {
    cfg  <- clust_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states,
                               cfg$seq_length, cfg$seed)
    cl   <- cluster_data(data, k = cfg$k, seed = cfg$seed)

    # k unique cluster labels
    n_unique <- length(unique(cl$assignments))
    ok_k <- identical(n_unique, as.integer(cfg$k))

    # Every actor assigned
    ok_len <- identical(length(cl$assignments), as.integer(cfg$n_actors))

    # Each cluster has >= 1 member
    ok_sizes <- all(cl$sizes >= 1L)

    # Cluster labels are 1..k
    ok_labels <- setequal(unique(cl$assignments), seq_len(cfg$k))

    # Silhouette is numeric scalar in [-1, 1]
    ok_sil <- is.numeric(cl$silhouette) && length(cl$silhouette) == 1L &&
      cl$silhouette >= -1 && cl$silhouette <= 1

    all(ok_k, ok_len, ok_sizes, ok_labels, ok_sil)
  }, logical(1))

  n_pass <- sum(results)
  report$log("cluster_data", "validity_50", N_CLUST,
             N_CLUST - n_pass, 0, 0, 0, 0, "manual", "structural validity")
  expect_true(all(results))
})


# ============================================================================
# 2. centrality InStrength/OutStrength match manual colSums/rowSums (100 configs)
# ============================================================================
test_that("centrality InStrength/OutStrength match colSums/rowSums (100 configs)", {
  skip_on_cran()
  skip_equiv_tests()

  deltas <- vapply(seq_len(N_CENT), function(i) {
    cfg  <- cent_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states,
                               cfg$seq_length, cfg$seed)
    net  <- build_network(data, method = "relative")

    cent <- centrality(net, measures = c("InStrength", "OutStrength"))

    mat <- net$weights
    diag(mat) <- 0  # centrality zeros diagonal by default (loops = FALSE)

    manual_in  <- colSums(mat)
    manual_out <- rowSums(mat)

    c(max(abs(cent$InStrength  - manual_in)),
      max(abs(cent$OutStrength - manual_out)))
  }, numeric(2))

  max_err  <- max(deltas)
  mean_err <- mean(deltas)
  n_fail   <- sum(deltas[1, ] > TOL_EXACT) + sum(deltas[2, ] > TOL_EXACT)

  report$log("centrality", "InOutStrength_100", N_CENT * 2L, n_fail,
             max_err, mean_err, stats::median(deltas),
             stats::quantile(deltas, 0.95), "colSums/rowSums",
             "directed relative networks")
  expect_true(max_err < TOL_EXACT,
              info = sprintf("InStrength/OutStrength max delta = %.2e", max_err))
})


# ============================================================================
# 3. centrality InStrength for undirected networks matches colSums (50 configs)
# ============================================================================
test_that("centrality InStrength for undirected matches colSums (50 configs)", {
  skip_on_cran()
  skip_equiv_tests()

  deltas <- vapply(seq_len(N_UNDIR), function(i) {
    cfg  <- undir_configs[[i]]
    data <- simulate_continuous(cfg$n, cfg$p, cfg$rho, cfg$seed)
    net  <- build_network(data, method = "cor")

    cent <- centrality(net, measures = c("InStrength", "OutStrength"))

    mat <- net$weights
    diag(mat) <- 0

    manual_in  <- colSums(mat)
    manual_out <- rowSums(mat)

    c(max(abs(cent$InStrength  - manual_in)),
      max(abs(cent$OutStrength - manual_out)))
  }, numeric(2))

  max_err  <- max(deltas)
  mean_err <- mean(deltas)
  n_fail   <- sum(deltas[1, ] > TOL_EXACT) + sum(deltas[2, ] > TOL_EXACT)

  report$log("centrality", "Strength_undirected_50", N_UNDIR * 2L, n_fail,
             max_err, mean_err, stats::median(deltas),
             stats::quantile(deltas, 0.95), "colSums/rowSums",
             "undirected cor networks")
  expect_true(max_err < TOL_EXACT,
              info = sprintf("Undirected Strength max delta = %.2e", max_err))
})


# ============================================================================
# 4. centrality betweenness proportional to igraph (50 configs)
# ============================================================================
# Nestimate's Floyd-Warshall betweenness differs from igraph by a constant
# factor of 4 for directed networks (Nestimate counts all ordered (s,t) pairs
# in its summation, yielding 4x igraph's value). We verify:
#   (a) zero-pattern agreement (same nodes have zero betweenness)
#   (b) constant proportionality ratio = 4 across all non-zero nodes
# This confirms identical shortest-path topology and ranking.
test_that("centrality betweenness proportional to igraph (50 configs)", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("igraph")

  EXPECTED_RATIO <- 4

  deltas <- vapply(seq_len(N_BTW), function(i) {
    cfg  <- btw_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states,
                               cfg$seq_length, cfg$seed)
    net  <- build_network(data, method = "relative")

    cent <- centrality(net, measures = "Betweenness")

    mat <- net$weights
    diag(mat) <- 0

    # Build igraph graph with inverted weights (higher weight = shorter path)
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "directed",
                                              weighted = TRUE,
                                              diag = FALSE)
    igraph::E(g)$weight <- 1 / igraph::E(g)$weight
    igraph_btw <- igraph::betweenness(g, directed = TRUE, normalized = TRUE)

    # Check zero-pattern agreement
    nest_zero   <- abs(cent$Betweenness) < 1e-15
    igraph_zero <- abs(igraph_btw) < 1e-15
    if (!identical(nest_zero, setNames(igraph_zero, names(nest_zero)))) {
      return(Inf)
    }

    # For non-zero entries, check ratio = EXPECTED_RATIO
    nz <- !nest_zero
    if (!any(nz)) return(0)  # all zeros, trivially OK

    ratios <- cent$Betweenness[nz] / igraph_btw[nz]
    max(abs(ratios - EXPECTED_RATIO))
  }, numeric(1))

  max_err  <- max(deltas)
  mean_err <- mean(deltas)
  n_fail   <- sum(deltas > TOL)

  report$log("centrality", "betweenness_igraph_50", N_BTW, n_fail,
             max_err, mean_err, stats::median(deltas),
             stats::quantile(deltas, 0.95), "igraph::betweenness",
             "directed, ratio=4 proportionality check")
  expect_true(max_err < TOL,
              info = sprintf("Betweenness ratio deviation max delta = %.2e", max_err))
})


# ============================================================================
# 5. cluster_data per-cluster networks are consistent (50 configs)
# ============================================================================
test_that("cluster_data per-cluster networks match manual build_network (50 configs)", {
  skip_on_cran()
  skip_equiv_tests()

  deltas <- vapply(seq_len(N_CLUST), function(i) {
    cfg  <- clust_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states,
                               cfg$seq_length, cfg$seed)
    cl   <- cluster_data(data, k = cfg$k, seed = cfg$seed)

    # Build per-cluster networks via the package pipeline
    nets <- build_network(cl)

    # For each cluster, manually build network from subset
    max_delta <- vapply(seq_len(cfg$k), function(cl_idx) {
      sub_data <- data[cl$assignments == cl_idx, , drop = FALSE]
      manual_net <- build_network(sub_data, method = "relative")

      pkg_mat    <- nets[[cl_idx]]$weights
      manual_mat <- manual_net$weights

      # Matrices may have different node sets (some states may not appear
      # in a cluster). Align by shared node names.
      shared <- intersect(rownames(pkg_mat), rownames(manual_mat))
      if (length(shared) == 0L) return(0)

      max(abs(pkg_mat[shared, shared] - manual_mat[shared, shared]))
    }, numeric(1))

    max(max_delta)
  }, numeric(1))

  max_err  <- max(deltas)
  mean_err <- mean(deltas)
  n_fail   <- sum(deltas > TOL_EXACT)

  report$log("cluster_data", "per_cluster_networks_50", N_CLUST, n_fail,
             max_err, mean_err, stats::median(deltas),
             stats::quantile(deltas, 0.95), "build_network(subset)",
             "relative method, per-cluster reconstruction")
  expect_true(max_err < TOL_EXACT,
              info = sprintf("Per-cluster network max delta = %.2e", max_err))
})


# ---- Write report ----
test_that("write equivalence report", {
  skip_on_cran()
  skip_equiv_tests()
  report$write_csv("cluster")
  msg <- report$summary()
  message(msg)
  expect_true(TRUE)
})
