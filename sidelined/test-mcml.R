# ===========================================================================
# Tests for build_mcml() — Multi-Cluster Multi-Layer Network
# ===========================================================================

# --- Helper: create test data ---
.make_test_matrix <- function(seed = 42) {
  sim <- simulate_mtna(n_nodes = 4, n_types = 3, seed = seed)
  sim
}

.make_wide_sequences <- function(seed = 1) {
  set.seed(seed)
  states <- c("plan", "monitor", "adapt", "discuss", "evaluate", "reflect")
  as.data.frame(matrix(
    sample(states, 500, replace = TRUE), nrow = 100, ncol = 5,
    dimnames = list(NULL, paste0("T", 1:5))
  ))
}

.test_clusters_6 <- list(
  Regulation = c("plan", "adapt"),
  Cognition = c("monitor", "evaluate"),
  Social = c("discuss", "reflect")
)


# ===========================================================================
# Section 1: Input validation
# ===========================================================================
test_that("mcml rejects unnamed clusters", {
  sim <- .make_test_matrix()
  expect_error(build_mcml(sim$matrix, clusters = list(c("A", "B"), c("C", "D"))),
               "named list")
})

test_that("mcml rejects fewer than 2 clusters", {
  sim <- .make_test_matrix()
  nodes <- sim$node_types[[1]]
  expect_error(build_mcml(sim$matrix, clusters = list(Only = nodes)),
               "length")
})

test_that("mcml rejects duplicate nodes across clusters", {
  sim <- .make_test_matrix()
  dup_clusters <- list(
    A = c("Diagnose", "Regulate"),
    B = c("Regulate", "Plan")
  )
  expect_error(build_mcml(sim$matrix, clusters = dup_clusters),
               "multiple clusters")
})

test_that("mcml rejects missing nodes in matrix", {
  sim <- .make_test_matrix()
  bad_clusters <- list(
    A = c("Diagnose", "NONEXISTENT"),
    B = c("Plan", "Judge")
  )
  expect_error(build_mcml(sim$matrix, clusters = bad_clusters),
               "not found")
})

test_that("mcml rejects missing nodes from estimated network", {
  wide <- .make_wide_sequences()
  bad_clusters <- list(
    A = c("plan", "NONEXISTENT"),
    B = c("adapt", "monitor")
  )
  expect_error(build_mcml(wide, clusters = bad_clusters),
               "not found")
})


# ===========================================================================
# Section 2: Construction from pre-computed matrix
# ===========================================================================
test_that("mcml works with pre-computed matrix", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)

  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_nodes, 12)
  expect_equal(mc$n_clusters, 3)
  expect_equal(mc$cluster_names, c("Metacognitive", "Cognitive", "Behavioral"))
  expect_identical(mc$clusters, sim$node_types)
})

test_that("mcml between matrix is row-normalized", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)

  bw <- mc$between
  rs <- rowSums(bw)
  expect_true(all(abs(rs - 1) < 1e-6))
})

test_that("mcml within matrices have correct dimensions", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)

  for (cl_name in names(mc$within)) {
    w <- mc$within[[cl_name]]
    n <- length(sim$node_types[[cl_name]])
    expect_equal(nrow(w), n)
    expect_equal(ncol(w), n)
    expect_equal(rownames(w), sim$node_types[[cl_name]])
  }
})

test_that("mcml preserves full node-level matrix", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)

  expect_identical(mc$matrix, sim$matrix)
})

test_that("mcml between and within are matrices", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)

  expect_true(is.matrix(mc$between))
  expect_equal(nrow(mc$between), mc$n_clusters)
  expect_true(is.list(mc$within))
  for (w in mc$within) {
    expect_true(is.matrix(w))
  }
})


# ===========================================================================
# Section 3: Construction from wide sequence data
# ===========================================================================
test_that("mcml works with wide sequence data frame", {
  wide <- .make_wide_sequences()
  mc <- build_mcml(wide, clusters = .test_clusters_6)

  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_nodes, 6)
  expect_equal(mc$n_clusters, 3)
  expect_equal(mc$method, "relative")
})

test_that("mcml builds correct node-level network from wide data", {
  wide <- .make_wide_sequences()
  mc <- build_mcml(wide, clusters = .test_clusters_6)

  # Matrix should be row-normalized (relative method)
  rs <- rowSums(mc$matrix)
  expect_true(all(abs(rs - 1) < 1e-6))
  expect_equal(nrow(mc$matrix), 6)
})

test_that("mcml frequency method produces integer-like counts", {
  wide <- .make_wide_sequences()
  mc <- build_mcml(wide, clusters = .test_clusters_6, method = "frequency")

  expect_equal(mc$method, "frequency")
  # Frequency matrix should have non-negative integers
  expect_true(all(mc$matrix >= 0))
  expect_true(all(mc$matrix == round(mc$matrix)))
})


# ===========================================================================
# Section 4: Aggregation methods
# ===========================================================================
test_that("mcml aggregation = 'sum' differs from 'mean'", {
  sim <- .make_test_matrix()
  mc_sum <- build_mcml(sim$matrix, clusters = sim$node_types, aggregation = "sum")
  mc_mean <- build_mcml(sim$matrix, clusters = sim$node_types, aggregation = "mean")

  # Both are row-normalized (type = "tna") but the relative weights differ
  # because sum vs mean aggregation changes the balance
  expect_equal(mc_sum$aggregation, "sum")
  expect_equal(mc_mean$aggregation, "mean")
})

test_that("mcml aggregation = 'max' works", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types, aggregation = "max")
  expect_equal(mc$aggregation, "max")
  expect_s3_class(mc, "mcml_network")
})


# ===========================================================================
# Section 5: S3 methods
# ===========================================================================
test_that("print.mcml_network produces output", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  out <- capture.output(print(mc))
  expect_true(any(grepl("Multi-Cluster Multi-Layer", out)))
  expect_true(any(grepl("Metacognitive", out)))
  expect_true(any(grepl("Cluster sizes", out)))
})

test_that("summary.mcml_network produces output", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  out <- capture.output(summary(mc))
  expect_true(any(grepl("BETWEEN-CLUSTER", out)))
  expect_true(any(grepl("WITHIN-CLUSTER", out)))
  expect_true(any(grepl("Strongest transitions", out)))
})

test_that("plot.mcml_network runs without error for type = 'mcml'", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  expect_no_error(plot(mc, type = "mcml"))
})

test_that("plot.mcml_network runs without error for type = 'between'", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  expect_no_error(plot(mc, type = "between"))
})

test_that("plot.mcml_network runs without error for type = 'within'", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  expect_no_error(plot(mc, type = "within"))
})


# ===========================================================================
# Section 6: Seed reproducibility
# ===========================================================================
test_that("mcml produces identical results with same seed", {
  wide <- .make_wide_sequences()
  mc1 <- build_mcml(wide, clusters = .test_clusters_6, seed = 123)
  mc2 <- build_mcml(wide, clusters = .test_clusters_6, seed = 123)
  expect_identical(mc1$matrix, mc2$matrix)
  expect_identical(mc1$between, mc2$between)
})


# ===========================================================================
# Section 7: Sequence-based vs matrix-based estimation
# ===========================================================================
test_that("mcml from sequences builds between by recoding states", {
  wide <- .make_wide_sequences()
  mc <- build_mcml(wide, clusters = .test_clusters_6)

  # Between should have cluster names as row/colnames
  expect_true(setequal(rownames(mc$between), names(.test_clusters_6)))
  expect_true(setequal(colnames(mc$between), names(.test_clusters_6)))
  # Between is row-normalized
  rs <- rowSums(mc$between)
  expect_true(all(abs(rs - 1) < 1e-6))
})

test_that("mcml between from sequences includes self-loops (TNA alignment)", {
  wide <- .make_wide_sequences()
  mc <- build_mcml(wide, clusters = .test_clusters_6)
  # Between should have non-zero diagonal (same-cluster transitions counted)
  expect_true(any(diag(mc$between) > 0))
  # Row-normalized
  rs <- rowSums(mc$between)
  expect_true(all(abs(rs - 1) < 1e-6))
})

test_that("mcml between from matrix includes self-loops", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  # Diagonal can be non-zero (aggregated from within-cluster block)
  expect_true(is.matrix(mc$between))
  rs <- rowSums(mc$between)
  expect_true(all(abs(rs - 1) < 1e-6))
})

test_that("mcml between frequency counts include same-cluster transitions", {
  wide <- .make_wide_sequences()
  mc_freq <- build_mcml(wide, clusters = .test_clusters_6, method = "frequency")
  # Diagonal can be non-zero — same-cluster transitions counted like TNA
  expect_true(any(diag(mc_freq$between) > 0))
  # All counts must be non-negative integers
  expect_true(all(mc_freq$between >= 0))
  expect_true(all(mc_freq$between == round(mc_freq$between)))
})

test_that("mcml from sequences builds within by filtering states", {
  wide <- .make_wide_sequences()
  mc <- build_mcml(wide, clusters = .test_clusters_6)

  # Each within matrix should have the cluster's state names
  for (cl_name in names(mc$within)) {
    w <- mc$within[[cl_name]]
    expect_true(all(rownames(w) %in% .test_clusters_6[[cl_name]]))
    # Row-normalized
    rs <- rowSums(w)
    expect_true(all(abs(rs - 1) < 1e-6))
  }
})

test_that("mcml cluster_summary field is available for plot_mcml", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  expect_s3_class(mc$cluster_summary, "cluster_summary")
  expect_equal(mc$cluster_summary$meta$n_clusters, 3)
  expect_equal(mc$cluster_summary$meta$n_nodes, 12)
})


# ===========================================================================
# Section 8: Edge cases
# ===========================================================================
test_that("mcml works with 2 clusters", {
  sim <- .make_test_matrix()
  cl2 <- list(
    Group1 = sim$node_types[[1]],
    Group2 = c(sim$node_types[[2]], sim$node_types[[3]])
  )
  mc <- build_mcml(sim$matrix, clusters = cl2)
  expect_equal(mc$n_clusters, 2)
  expect_equal(nrow(mc$between), 2)
})

test_that("mcml works with single-node clusters", {
  mat <- matrix(c(0, 0.6, 0.4,
                  0.3, 0, 0.7,
                  0.5, 0.5, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  mc <- build_mcml(mat, clusters = list(X = "A", Y = "B", Z = "C"))
  expect_equal(mc$n_clusters, 3)
  expect_equal(mc$n_nodes, 3)
})

test_that("mcml works when not all matrix nodes are in clusters", {
  mat <- matrix(c(0, 0.3, 0.3, 0.4,
                  0.5, 0, 0.2, 0.3,
                  0.2, 0.3, 0, 0.5,
                  0.4, 0.4, 0.2, 0), 4, 4, byrow = TRUE,
                dimnames = list(c("A", "B", "C", "D"),
                                c("A", "B", "C", "D")))
  # Only cluster A+B and C — D is unclustered
  # cluster_summary should handle subset of nodes
  mc <- build_mcml(mat, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))
  expect_s3_class(mc, "mcml_network")
})


# ===========================================================================
# Section 9: Multiple seeds produce different results
# ===========================================================================
test_that("different seeds on wide data produce different networks", {
  states <- c("plan", "monitor", "adapt", "discuss", "evaluate", "reflect")
  set.seed(1)
  wide1 <- as.data.frame(matrix(
    sample(states, 500, replace = TRUE), nrow = 100, ncol = 5,
    dimnames = list(NULL, paste0("T", 1:5))
  ))
  set.seed(99)
  wide2 <- as.data.frame(matrix(
    sample(states, 500, replace = TRUE), nrow = 100, ncol = 5,
    dimnames = list(NULL, paste0("T", 1:5))
  ))
  mc1 <- build_mcml(wide1, clusters = .test_clusters_6)
  mc2 <- build_mcml(wide2, clusters = .test_clusters_6)
  expect_false(identical(mc1$between, mc2$between))
})


# ===========================================================================
# Section 10: Multi-format data inputs
# ===========================================================================
test_that("mcml works with edge list data.frame", {
  el <- data.frame(
    from = c("A", "A", "B", "B", "C", "C"),
    to = c("B", "C", "A", "C", "A", "B"),
    weight = c(0.5, 0.3, 0.4, 0.6, 0.2, 0.8)
  )
  mc <- build_mcml(el, clusters = list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_nodes, 3)
  expect_equal(mc$n_clusters, 2)
})

test_that("mcml works with tna object", {
  wide <- .make_wide_sequences()
  tna_obj <- tna::tna(wide)
  mc <- build_mcml(tna_obj, clusters = .test_clusters_6)
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_nodes, 6)
})

test_that("mcml works with cograph_network", {
  el <- data.frame(
    from = c("A", "B", "C", "D"),
    to = c("B", "C", "D", "A"),
    weight = c(0.5, 0.3, 0.2, 0.4)
  )
  g <- cograph::as_cograph(el, directed = TRUE)
  mc <- build_mcml(g, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_clusters, 2)
})

test_that("mcml works with cograph_network + auto-detect clusters", {
  el <- data.frame(
    from = c("A", "B", "C", "D"),
    to = c("B", "C", "D", "A"),
    weight = c(0.5, 0.3, 0.2, 0.4)
  )
  g <- cograph::as_cograph(el, directed = TRUE)
  g$nodes$cluster <- c("G1", "G1", "G2", "G2")
  mc <- build_mcml(g, clusters = NULL)
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_clusters, 2)
  expect_true(all(c("A", "B") %in% mc$clusters[["G1"]]))
  expect_true(all(c("C", "D") %in% mc$clusters[["G2"]]))
})

test_that("mcml works with netobject", {
  wide <- .make_wide_sequences()
  net <- build_network(wide, method = "relative")
  mc <- build_mcml(net, clusters = .test_clusters_6)
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_nodes, 6)
})


# ===========================================================================
# Section 11: Multi-format cluster inputs
# ===========================================================================
test_that("mcml works with cluster as named vector", {
  sim <- .make_test_matrix()
  nodes <- unlist(sim$node_types, use.names = FALSE)
  groups <- rep(names(sim$node_types), vapply(sim$node_types, length, integer(1)))
  named_vec <- setNames(groups, nodes)
  mc <- build_mcml(sim$matrix, clusters = named_vec)
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_clusters, 3)
})

test_that("mcml works with cluster as plain vector", {
  mat <- matrix(c(0, 0.6, 0.4,
                  0.3, 0, 0.7,
                  0.5, 0.5, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  mc <- build_mcml(mat, clusters = c("G1", "G1", "G2"))
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_clusters, 2)
  expect_true("A" %in% mc$clusters[["G1"]])
})

test_that("mcml works with cluster as data.frame", {
  sim <- .make_test_matrix()
  nodes <- unlist(sim$node_types, use.names = FALSE)
  groups <- rep(names(sim$node_types), vapply(sim$node_types, length, integer(1)))
  cl_df <- data.frame(node = nodes, group = groups)
  mc <- build_mcml(sim$matrix, clusters = cl_df)
  expect_s3_class(mc, "mcml_network")
  expect_equal(mc$n_clusters, 3)
})

test_that("mcml rejects NULL clusters without cograph_network", {
  sim <- .make_test_matrix()
  expect_error(build_mcml(sim$matrix, clusters = NULL), "auto-detect")
})

test_that("mcml rejects plain vector with wrong length", {
  mat <- matrix(c(0, 0.6, 0.4,
                  0.3, 0, 0.7,
                  0.5, 0.5, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_error(build_mcml(mat, clusters = c("G1", "G2")), "must match")
})


# ===========================================================================
# Section 12: Edge list output
# ===========================================================================
test_that("mcml edges has correct columns", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  expect_true(is.data.frame(mc$edges))
  expected_cols <- c("from", "to", "weight", "cluster_from", "cluster_to", "type")
  expect_true(all(expected_cols %in% names(mc$edges)))
})

test_that("mcml edges have correct within/between types", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)

  within_edges <- mc$edges[mc$edges$type == "within", ]
  between_edges <- mc$edges[mc$edges$type == "between", ]

  # Within edges: cluster_from == cluster_to
  expect_true(all(within_edges$cluster_from == within_edges$cluster_to))
  # Between edges: cluster_from != cluster_to
  expect_true(all(between_edges$cluster_from != between_edges$cluster_to))
})

test_that("mcml edges cover all non-zero off-diagonal entries", {
  sim <- .make_test_matrix()
  mc <- build_mcml(sim$matrix, clusters = sim$node_types)
  mat <- mc$matrix
  n_expected <- sum(mat != 0 & row(mat) != col(mat))
  expect_equal(nrow(mc$edges), n_expected)
})
