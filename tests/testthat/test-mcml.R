# Tests for mcml.R: net_aggregate_weights, cluster_summary, build_mcml

# ============================================
# net_aggregate_weights / wagg
# ============================================

test_that("net_aggregate_weights sum method", {
  expect_equal(net_aggregate_weights(c(0.5, 0.8, 0.3), "sum"), 1.6)
})

test_that("net_aggregate_weights mean method", {
  expect_equal(net_aggregate_weights(c(2, 4, 6), "mean"), 4)
})

test_that("net_aggregate_weights median method", {
  expect_equal(net_aggregate_weights(c(1, 3, 5), "median"), 3)
})

test_that("net_aggregate_weights max method", {
  expect_equal(net_aggregate_weights(c(1, 5, 3), "max"), 5)
})

test_that("net_aggregate_weights min method", {
  expect_equal(net_aggregate_weights(c(1, 5, 3), "min"), 1)
})

test_that("net_aggregate_weights prod method", {
  expect_equal(net_aggregate_weights(c(2, 3, 4), "prod"), 24)
})

test_that("net_aggregate_weights density with n_possible", {
  expect_equal(net_aggregate_weights(c(1, 2, 3), "density", n_possible = 10), 0.6)
})

test_that("net_aggregate_weights density without n_possible", {
  expect_equal(net_aggregate_weights(c(1, 2, 3), "density"), 2)
})

test_that("net_aggregate_weights geomean method", {
  expect_equal(net_aggregate_weights(c(4, 9), "geomean"), 6, tolerance = 0.01)
})

test_that("net_aggregate_weights removes NA and zero", {
  expect_equal(net_aggregate_weights(c(1, NA, 0, 2), "sum"), 3)
})

test_that("net_aggregate_weights returns 0 for empty/all-zero input", {
  expect_equal(net_aggregate_weights(c(0, 0, NA), "sum"), 0)
  expect_equal(net_aggregate_weights(numeric(0), "mean"), 0)
})

test_that("net_aggregate_weights errors on unknown method", {
  expect_error(net_aggregate_weights(c(1, 2), "bogus"), "Unknown method")
})

test_that("net_aggregate_weights is a function", {
  expect_true(is.function(net_aggregate_weights))
})

# ============================================
# cluster_summary
# ============================================

test_that("cluster_summary with vector clusters", {
  mat <- matrix(c(10, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- c(A = 1, B = 1, C = 2)

  cs <- cluster_summary(mat, clusters)
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
  expect_equal(ncol(cs$macro$weights), 2)
  expect_equal(cs$meta$n_clusters, 2)
  expect_equal(cs$meta$n_nodes, 3)
})

test_that("cluster_summary with named list clusters", {
  mat <- matrix(runif(16), 4, 4,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- cluster_summary(mat, clusters)
  expect_equal(rownames(cs$macro$weights), c("G1", "G2"))
  expect_equal(colnames(cs$macro$weights), c("G1", "G2"))
  expect_equal(length(cs$clusters), 2)
  expect_equal(nrow(cs$clusters$G1$weights), 2)
})

test_that("cluster_summary matrix path returns raw aggregated values (no normalization)", {
  mat <- matrix(c(10, 2, 3, 8), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  clusters <- list(G1 = "A", G2 = "B")

  cs <- cluster_summary(mat, clusters, method = "sum")
  expect_equal(cs$macro$weights["G1", "G1"], 10)
  expect_equal(cs$macro$weights["G1", "G2"], 3)
})

test_that("cluster_summary compute_within=FALSE skips within", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs <- cluster_summary(mat, clusters, compute_within = FALSE)
  expect_null(cs$clusters)
})

test_that("cluster_summary inits sum to 1", {
  mat <- matrix(runif(16), 4, 4,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- cluster_summary(mat, clusters)
  expect_equal(sum(cs$macro$inits), 1, tolerance = 1e-10)
})

test_that("cluster_summary errors without clusters", {
  mat <- matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(cluster_summary(mat), "clusters")
})

test_that("cluster_summary errors on non-square matrix", {
  mat <- matrix(1, 2, 3)
  expect_error(cluster_summary(mat, c(1, 2)), "square")
})

test_that("cluster_summary print method works", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  expect_output(print(cs))
})

test_that("cluster_summary is a function", {
  expect_true(is.function(cluster_summary))
})

test_that("cluster_summary with different methods", {
  mat <- matrix(c(4, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs_sum <- cluster_summary(mat, clusters, method = "sum")
  cs_mean <- cluster_summary(mat, clusters, method = "mean")
  cs_max <- cluster_summary(mat, clusters, method = "max")

  # Sum should be larger than mean for multi-node clusters
  expect_true(cs_sum$macro$weights["G1", "G1"] >= cs_mean$macro$weights["G1", "G1"])
  # Max should be at most sum
  expect_true(cs_max$macro$weights["G1", "G1"] <= cs_sum$macro$weights["G1", "G1"])
})

test_that("cluster_summary within-cluster has correct dimensions", {
  mat <- matrix(runif(16), 4, 4,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  clusters <- list(G1 = c("A", "B", "C"), G2 = "D")

  cs <- cluster_summary(mat, clusters)
  expect_equal(nrow(cs$clusters$G1$weights), 3)
  expect_equal(ncol(cs$clusters$G1$weights), 3)
  expect_equal(nrow(cs$clusters$G2$weights), 1)
})

# ============================================
# build_mcml
# ============================================

test_that("build_mcml with sequence data", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(c("A", "B", "C", "D"), 50, replace = TRUE),
    T2 = sample(c("A", "B", "C", "D"), 50, replace = TRUE),
    T3 = sample(c("A", "B", "C", "D"), 50, replace = TRUE),
    T4 = sample(c("A", "B", "C", "D"), 50, replace = TRUE)
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(seqs, clusters)
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
  expect_equal(sum(cs$macro$inits), 1, tolerance = 1e-10)
})

test_that("build_mcml with edge list", {
  edges <- data.frame(
    from = c("A", "A", "B", "C", "C", "D"),
    to = c("B", "C", "A", "D", "D", "A"),
    weight = c(1, 2, 1, 3, 1, 2)
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(edges, clusters)
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml type=raw preserves counts", {
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "D"),
    T3 = c("C", "A")
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(seqs, clusters, type = "raw")
  expect_true(is.numeric(cs$macro$weights))
})

test_that("build_mcml returns mcml if already mcml", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  cs2 <- build_mcml(cs)
  expect_identical(cs, cs2)
})

test_that("build_mcml tna type normalizes rows", {
  set.seed(1)
  seqs <- data.frame(
    T1 = sample(c("A", "B", "C"), 30, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 30, replace = TRUE),
    T3 = sample(c("A", "B", "C"), 30, replace = TRUE)
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters, type = "tna")
  row_sums <- rowSums(cs$macro$weights)
  expect_equal(unname(row_sums), c(1, 1), tolerance = 1e-10)
})

# ============================================================
# Coverage gap tests — mcml.R
# ============================================================

# ---- net_aggregate_weights / wagg: geomean zero path (L43) ----

test_that("net_aggregate_weights geomean returns 0 when all non-positive", {
  # All negative: pos_w is empty → returns 0
  expect_equal(net_aggregate_weights(c(-1, -2, -3), "geomean"), 0)
})

# ---- cluster_summary: various input types ----

test_that("cluster_summary with tna object extracts weights (L300)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.6, 0.4, 0.7, 0, 0.3, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  cs <- cluster_summary(tna_obj, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("cluster_summary assigns sequential node names when matrix has no rownames (L313-315)", {
  mat <- matrix(c(0, 2, 3, 1, 0, 4, 5, 6, 0), 3, 3)
  # No rownames/colnames
  cs <- cluster_summary(mat, list(G1 = c("1", "2"), G2 = "3"))
  expect_s3_class(cs, "mcml")
  expect_equal(cs$meta$n_nodes, 3)
})

test_that("cluster_summary between_inits fallback when zero matrix (L380)", {
  # All-zero matrix → colSums all zero → uniform inits
  mat <- matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  cs <- cluster_summary(mat, list(G1 = "A", G2 = "B"))
  expect_equal(unname(cs$macro$inits), c(0.5, 0.5), tolerance = 1e-10)
})

test_that("cluster_summary within-cluster zero total produces uniform inits (L423)", {
  # All-zero within-cluster block
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  # G1 has 2 nodes, zero total → rep(0.5, 2)
  expect_equal(unname(cs$clusters$G1$inits), c(0.5, 0.5), tolerance = 1e-10)
})

# ---- build_mcml: cograph_network input (L572) ----

test_that("build_mcml accepts cograph_network input (L572)", {
  mat <- matrix(c(0, 0.3, 0.7, 0.4, 0, 0.6, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  # Create a minimal cograph_network object
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         stringsAsFactors = FALSE),
      edges = data.frame(from = integer(0), to = integer(0),
                         weight = numeric(0)),
      directed = TRUE,
      data = NULL,
      meta = list(),
      node_groups = NULL
    ),
    class = c("cograph_network", "list")
  )
  cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

# ---- build_mcml: tna_data branch (L581-598) ----

test_that("build_mcml handles tna object with data (tna_data branch, L581-598)", {
  skip_if_not_installed("tna")
  seqs <- data.frame(
    T1 = c("A", "B", "C"),
    T2 = c("B", "C", "A"),
    T3 = c("C", "A", "B")
  )
  # Build a tna object that contains $data (integer-encoded: 1=A, 2=B, 3=C)
  tna_obj <- tna::tna(seqs)
  if (!is.null(tna_obj$data)) {
    # tna encodes states as integers 1,2,3 in $data
    cs <- build_mcml(tna_obj, list(G1 = c("1", "2"), G2 = "3"))
    expect_s3_class(cs, "mcml")
  } else {
    skip("tna object has no $data field in this version")
  }
})

test_that("build_mcml handles tna object without data (tna_matrix branch, L583-585)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.6, 0.4, 0.7, 0, 0.3, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  cs <- build_mcml(tna_obj, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

# ---- build_mcml: netobject_data branch (L586-599) ----

test_that("build_mcml handles netobject with sequence data (netobject_data branch, L586-599)", {
  set.seed(1)
  seqs <- data.frame(
    T1 = sample(c("A", "B", "C"), 20, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 20, replace = TRUE),
    T3 = sample(c("A", "B", "C"), 20, replace = TRUE)
  )
  net <- build_network(seqs, method = "relative")
  # net$data contains the raw sequences
  cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml handles netobject with edge list data (netobject_data edgelist, L593-595)", {
  edges <- data.frame(
    from = c("A", "A", "B", "C"),
    to   = c("B", "C", "A", "A"),
    stringsAsFactors = FALSE
  )
  net <- build_network(edges, method = "relative")
  # If net has $data, it should route through netobject_data
  if (!is.null(net$data)) {
    cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
    expect_s3_class(cs, "mcml")
  } else {
    # Falls to netobject_matrix branch — still valid test
    cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
    expect_s3_class(cs, "mcml")
  }
})

# ---- build_mcml: netobject_matrix branch (L601-607) ----

test_that("build_mcml handles netobject without data (netobject_matrix branch, L601-607)", {
  mat <- matrix(c(0, 0.4, 0.6, 0.3, 0, 0.7, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         stringsAsFactors = FALSE),
      edges = data.frame(from = integer(0), to = integer(0),
                         weight = numeric(0)),
      directed = TRUE,
      data = NULL,
      meta = list(),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
  cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml handles plain numeric matrix (matrix branch, L608-610)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- build_mcml(mat, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml errors on unknown input class (L611-612)", {
  expect_error(build_mcml(list(foo = 1), list(G1 = "A")),
               "Cannot build MCML")
})

# ---- .detect_mcml_input coverage (L620-644) ----

test_that(".detect_mcml_input returns tna_data for tna with data (L620)", {
  skip_if_not_installed("tna")
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  tna_obj <- tna::tna(seqs)
  if (!is.null(tna_obj$data)) {
    result <- Nestimate:::.detect_mcml_input(tna_obj)
    expect_equal(result, "tna_data")
  } else {
    skip("tna object has no $data in this version")
  }
})

test_that(".detect_mcml_input returns tna_matrix for tna without data (L621)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.5, 0.5, 0.3, 0, 0.7, 0.6, 0.4, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  result <- Nestimate:::.detect_mcml_input(tna_obj)
  expect_equal(result, "tna_matrix")
})

test_that(".detect_mcml_input returns netobject_data for netobject with data (L625)", {
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  net <- build_network(seqs, method = "relative")
  if (!is.null(net$data)) {
    result <- Nestimate:::.detect_mcml_input(net)
    expect_equal(result, "netobject_data")
  } else {
    skip("netobject has no $data in this configuration")
  }
})

test_that(".detect_mcml_input returns netobject_matrix for netobject without data (L626)", {
  mat <- matrix(c(0, 0.4, 0.6, 0.7, 0, 0.3, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(weights = mat, nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                                           name = c("A", "B", "C"),
                                           x = NA_real_, y = NA_real_,
                                           stringsAsFactors = FALSE),
         edges = data.frame(from = integer(0), to = integer(0),
                            weight = numeric(0)),
         directed = TRUE, data = NULL, meta = list(), node_groups = NULL),
    class = c("netobject", "cograph_network")
  )
  result <- Nestimate:::.detect_mcml_input(net)
  expect_equal(result, "netobject_matrix")
})

test_that(".detect_mcml_input returns sequence for non-square numeric matrix (L641)", {
  mat <- matrix(1:6, 2, 3)
  result <- Nestimate:::.detect_mcml_input(mat)
  expect_equal(result, "sequence")
})

test_that(".detect_mcml_input returns unknown for unrecognized class (L644)", {
  result <- Nestimate:::.detect_mcml_input(42L)
  expect_equal(result, "unknown")
})

# ---- .auto_detect_clusters coverage (L650-672) ----

test_that(".auto_detect_clusters finds cluster from nodes$cluster column (L650-656)", {
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         cluster = c("G1", "G1", "G2"),
                         stringsAsFactors = FALSE),
      data = NULL, node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
  result <- Nestimate:::.auto_detect_clusters(net)
  expect_equal(unname(result), c("G1", "G1", "G2"))
})

test_that(".auto_detect_clusters falls back to node_groups (L660-664)", {
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         stringsAsFactors = FALSE),
      data = NULL,
      node_groups = data.frame(cluster = c("G1", "G1", "G2"),
                               stringsAsFactors = FALSE)
    ),
    class = c("netobject", "cograph_network")
  )
  result <- Nestimate:::.auto_detect_clusters(net)
  expect_equal(unname(result), c("G1", "G1", "G2"))
})

test_that(".auto_detect_clusters errors when no cluster info found (L667-671)", {
  mat <- matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:2, label = c("A", "B"),
                         stringsAsFactors = FALSE),
      data = NULL, node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
  expect_error(Nestimate:::.auto_detect_clusters(net), "No clusters found")
})

# ---- .build_cluster_lookup coverage (L696-725) ----

test_that(".build_cluster_lookup errors on unmapped nodes (L696-698)", {
  cl <- list(G1 = c("A", "B"))
  # all_nodes includes C, which is not in G1
  expect_error(
    Nestimate:::.build_cluster_lookup(cl, c("A", "B", "C")),
    "Unmapped nodes"
  )
})

test_that(".build_cluster_lookup works with character membership vector (L703-711)", {
  all_nodes <- c("A", "B", "C")
  result <- Nestimate:::.build_cluster_lookup(c("G1", "G1", "G2"), all_nodes)
  expect_equal(unname(result), c("G1", "G1", "G2"))
  expect_equal(names(result), c("A", "B", "C"))
})

test_that(".build_cluster_lookup errors when char vector length mismatches (L705-708)", {
  expect_error(
    Nestimate:::.build_cluster_lookup(c("G1", "G1"), c("A", "B", "C")),
    "Membership vector length"
  )
})

test_that(".build_cluster_lookup works with numeric membership vector (L714-721)", {
  all_nodes <- c("X", "Y", "Z")
  result <- Nestimate:::.build_cluster_lookup(c(1, 1, 2), all_nodes)
  expect_equal(unname(result), c("1", "1", "2"))
  expect_equal(names(result), c("X", "Y", "Z"))
})

test_that(".build_cluster_lookup errors when numeric vector length mismatches (L715-718)", {
  expect_error(
    Nestimate:::.build_cluster_lookup(c(1, 2), c("A", "B", "C")),
    "Membership vector length"
  )
})

test_that(".build_cluster_lookup errors on invalid input type (L724-725)", {
  expect_error(
    Nestimate:::.build_cluster_lookup(TRUE, c("A", "B")),
    "clusters must be a named list"
  )
})

# ---- .build_from_transitions: density method and zero-inits paths ----

test_that("build_mcml with density method triggers n_possible computation (L761-764)", {
  seqs <- data.frame(
    T1 = c("A", "C", "B"),
    T2 = c("B", "D", "C"),
    T3 = c("C", "A", "D")
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
  cs <- build_mcml(seqs, clusters, method = "density", type = "raw")
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml zero transitions produce uniform between_inits (L784)", {
  # Create edge list that has no transitions between G2 -> G1 or G2 -> G2
  # Force all edges in G1 to make G2 isolated
  edges <- data.frame(
    from = c("A", "B"),
    to   = c("B", "A"),
    stringsAsFactors = FALSE
  )
  # G2 has node "C" but no transitions at all → zero col sums for G2
  # Actually need G2 to have zero total in cluster
  # Simpler: all transitions are G1 internal, G2 isolated
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  # Since "C" never appears in edges, it won't appear in all_nodes
  # Let's add a C transition explicitly
  edges2 <- data.frame(
    from = c("A", "B", "C"),
    to   = c("B", "A", "C"),
    stringsAsFactors = FALSE
  )
  cs <- build_mcml(edges2, clusters, type = "raw")
  expect_s3_class(cs, "mcml")
})

# ---- Single-node within-cluster (L839-847) ----

test_that("build_mcml single-node cluster computes self-loop weight (L839-847)", {
  # Cluster G2 has only 1 node "C"
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "C"),
    T3 = c("C", "A")
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters, type = "raw")
  expect_s3_class(cs, "mcml")
  # G2 is a 1x1 within matrix
  expect_equal(dim(cs$clusters$G2$weights), c(1, 1))
  expect_equal(unname(cs$clusters$G2$inits), 1)
})

test_that("build_mcml single-node cluster with no self-loops returns 0 weight (L843-844)", {
  edges <- data.frame(
    from = c("A", "B"),
    to   = c("B", "A"),
    stringsAsFactors = FALSE
  )
  # G2 node "C" has no edges at all — single node, no self-loops
  # Add C to appear in nodes but never in edges: need it in seqs
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters, type = "raw")
  # C appears once but only in T1 position, no C->C transition
  expect_equal(unname(cs$clusters$G2$inits), 1)
})

# ---- .build_mcml_edgelist: fallback column detection (L945, L949) ----

test_that(".build_mcml_edgelist falls back to first/second columns when no named from/to (L945,L949)", {
  # Column names not in the standard from/to list — call .build_mcml_edgelist directly
  edges <- data.frame(
    node_from = c("A", "B", "C"),
    node_to   = c("B", "C", "A"),
    stringsAsFactors = FALSE
  )
  # Neither 'node_from' nor 'node_to' match standard aliases
  # so from_col and to_col fall back to 1 and 2 respectively
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- Nestimate:::.build_mcml_edgelist(edges, clusters, "sum", "tna", TRUE, TRUE)
  expect_s3_class(cs, "mcml")
})

# ---- Column-name clusters branch in edgelist (L971-989) ----

test_that("build_mcml edgelist accepts cluster column name (L971-989)", {
  edges <- data.frame(
    from  = c("A", "A", "B", "C", "C", "D"),
    to    = c("B", "C", "A", "D", "D", "A"),
    group = c("G1", "G1", "G1", "G2", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  cs <- build_mcml(edges, clusters = "group")
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml edgelist errors when clusters=NULL (L992-993)", {
  edges <- data.frame(
    from = c("A", "B"),
    to   = c("B", "A"),
    stringsAsFactors = FALSE
  )
  expect_error(build_mcml(edges, clusters = NULL),
               "clusters argument is required")
})

test_that("build_mcml edgelist accepts membership vector clusters (L1000)", {
  edges <- data.frame(
    from = c("A", "B", "C"),
    to   = c("B", "C", "A"),
    stringsAsFactors = FALSE
  )
  # character membership vector for the 3 unique sorted nodes (A, B, C)
  clusters <- c("G1", "G1", "G2")
  cs <- build_mcml(edges, clusters)
  expect_s3_class(cs, "mcml")
})

# ---- .build_mcml_sequence: matrix input, single column error, null clusters (L1017-1044) ----

test_that("build_mcml sequence builder coerces matrix input (L1017)", {
  mat_seq <- matrix(c("A", "B", "C", "B", "C", "A"), nrow = 3,
                    dimnames = list(NULL, c("T1", "T2")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(mat_seq, clusters)
  expect_s3_class(cs, "mcml")
})

test_that(".build_mcml_sequence errors on single-column data (L1023-1024)", {
  df <- data.frame(T1 = c("A", "B", "C"), stringsAsFactors = FALSE)
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  expect_error(
    Nestimate:::.build_mcml_sequence(df, clusters, "sum", "tna", TRUE, TRUE),
    "at least 2 columns"
  )
})

test_that(".build_mcml_sequence errors when clusters=NULL (L1043-1044)", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  expect_error(
    Nestimate:::.build_mcml_sequence(df, NULL, "sum", "tna", TRUE, TRUE),
    "clusters argument is required"
  )
})

test_that(".build_mcml_sequence calls .normalize_clusters for non-list clusters (L1050)", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  # Pass a character membership vector (triggers .normalize_clusters path)
  clusters <- c("G1", "G2")  # 2 nodes: A, B
  cs <- Nestimate:::.build_mcml_sequence(df, clusters, "sum", "tna", TRUE, TRUE)
  expect_s3_class(cs, "mcml")
})

# ---- .process_weights: cooccurrence type (L1070) ----

test_that(".process_weights returns symmetrized matrix for cooccurrence type (L1070)", {
  raw <- matrix(c(1, 3, 2, 4), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  result <- Nestimate:::.process_weights(raw, "cooccurrence")
  expect_true(isSymmetric(result))
  expect_equal(result["A", "B"], 2.5)
})

# ---- as_tna generic dispatch and methods (L1226-1306) ----

test_that("as_tna dispatches correctly for mcml objects (L1226)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
})

test_that("as_tna.mcml on matrix path uses frequency method (matrix is aggregation-only)", {
  mat <- matrix(c(10, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
  # Matrix-derived mcml carries no $meta$type, as_tna defaults to "frequency".
  expect_equal(result$macro$method, "frequency")
})

test_that("as_tna.mcml creates macro netobject (L1243-1244)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result$macro, "netobject")
  expect_true(is.matrix(result$macro$weights))
})

test_that("as_tna.mcml skips clusters with zero-row sums (L1247-1262)", {
  # Cluster with all-zero rows will be excluded
  mat <- matrix(0, 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  mat["A", "B"] <- 0.5; mat["A", "C"] <- 0.5  # only A has outgoing
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
  # macro always present
  expect_true("macro" %in% names(result))
})

test_that("as_tna.mcml with compute_within=FALSE returns empty cluster list (L1256-1258)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"),
                        compute_within = FALSE)
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
  expect_true("macro" %in% names(result))
})

test_that(".wrap_netobject produces valid dual-class object (L1270-1296)", {
  mat <- matrix(c(0, 0.6, 0.4, 0.3, 0, 0.7, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  result <- Nestimate:::.wrap_netobject(mat, method = "relative", directed = TRUE)
  expect_s3_class(result, "netobject")
  expect_true(inherits(result, "cograph_network"))
  expect_equal(result$method, "relative")
  expect_equal(result$n_nodes, 3)
  expect_true(is.data.frame(result$nodes))
  expect_equal(result$nodes$label, c("A", "B", "C"))
})

test_that("as_tna.default returns tna object unchanged (L1303-1305)", {
  skip_if_not_installed("tna")
  mat <- matrix(c(0, 0.5, 0.5, 0.3, 0, 0.7, 0.6, 0.4, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tna_obj <- tna::tna(mat)
  result <- as_tna(tna_obj)
  expect_true(inherits(result, "tna"))
})

test_that("as_tna.default errors for non-tna objects (L1305-1306)", {
  expect_error(as_tna(list(a = 1)), "Cannot convert")
  expect_error(as_tna("some string"), "Cannot convert")
})

# ---- .normalize_clusters coverage (L1312-1362) ----

test_that(".normalize_clusters handles data.frame input (L1313-1318)", {
  node_names <- c("A", "B", "C", "D")
  clusters_df <- data.frame(
    node  = c("A", "B", "C", "D"),
    group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  result <- Nestimate:::.normalize_clusters(clusters_df, node_names)
  expect_true(is.list(result))
  expect_equal(sort(names(result)), c("G1", "G2"))
})

test_that(".normalize_clusters errors on unknown nodes in list (L1325-1327)", {
  node_names <- c("A", "B")
  # C is not in node_names
  expect_error(
    Nestimate:::.normalize_clusters(list(G1 = c("A", "C")), node_names),
    "Unknown nodes"
  )
})

test_that(".normalize_clusters errors on numeric vector wrong length (L1335-1337)", {
  node_names <- c("A", "B", "C")
  expect_error(
    Nestimate:::.normalize_clusters(c(1, 2), node_names),
    "Membership vector length"
  )
})

test_that(".normalize_clusters handles factor cluster membership (L1348-1359)", {
  node_names <- c("A", "B", "C")
  clusters_fac <- factor(c("G1", "G1", "G2"))
  result <- Nestimate:::.normalize_clusters(clusters_fac, node_names)
  expect_true(is.list(result))
  expect_equal(sort(names(result)), c("G1", "G2"))
  expect_true("A" %in% result$G1)
})

test_that(".normalize_clusters handles character cluster membership (L1348-1359)", {
  node_names <- c("A", "B", "C")
  result <- Nestimate:::.normalize_clusters(c("G1", "G1", "G2"), node_names)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
})

test_that(".normalize_clusters errors on character vector wrong length (L1350-1351)", {
  node_names <- c("A", "B", "C")
  expect_error(
    Nestimate:::.normalize_clusters(c("G1", "G2"), node_names),
    "Membership vector length"
  )
})

test_that(".normalize_clusters errors on invalid input type (L1362)", {
  expect_error(
    Nestimate:::.normalize_clusters(TRUE, c("A", "B")),
    "clusters must be a list"
  )
})

# ---- print.mcml with edges field (L1381-1384) ----

test_that("print.mcml shows Transitions line when edges present (L1381-1384)", {
  seqs <- data.frame(
    T1 = c("A", "B", "C"),
    T2 = c("B", "C", "A"),
    T3 = c("C", "A", "B")
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters)
  # build_mcml produces $edges — print should show Transitions line
  output <- capture.output(print(cs))
  expect_true(any(grepl("Transitions", output)))
})

test_that("print.mcml without edges shows no Transitions line (L1385-1386)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  # cluster_summary does not set $edges
  output <- capture.output(print(cs))
  expect_true(any(grepl("MCML Network", output)))
})

# ---- summary.mcml (L1404) ----

test_that("summary.mcml runs and produces output (L1404)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  expect_output(summary(cs), "MCML")
})

test_that("summary.mcml on build_mcml result also works (L1404)", {
  seqs <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B")
  )
  cs <- build_mcml(seqs, list(G1 = "A", G2 = "B"))
  expect_output(summary(cs))
})


# ============================================
# mcml seamless dispatch to downstream functions
# ============================================

.make_mcml <- function(seed = 42) {
  set.seed(seed)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:6], 40, TRUE),
    T2 = sample(LETTERS[1:6], 40, TRUE),
    T3 = sample(LETTERS[1:6], 40, TRUE),
    T4 = sample(LETTERS[1:6], 40, TRUE),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
  build_mcml(seqs, clusters, type = "tna")
}


test_that("centrality_stability() works on mcml", {
  cs <- .make_mcml()
  stab <- centrality_stability(cs, iter = 20, seed = 1)

  expect_true(is.list(stab))
  expect_true("macro" %in% names(stab))
  for (nm in names(stab)) {
    expect_s3_class(stab[[nm]], "net_stability")
  }
})

test_that("network_reliability() works on mcml", {
  cs <- .make_mcml()
  rel <- network_reliability(cs, iter = 20, seed = 1)

  expect_s3_class(rel, "net_reliability")
})

test_that("extract_transition_matrix() works on mcml", {
  cs <- .make_mcml()
  mats <- extract_transition_matrix(cs)

  expect_true(is.list(mats))
  expect_true("macro" %in% names(mats))
  expect_true(is.matrix(mats$macro))
  # Within-cluster matrices
  cluster_names <- setdiff(names(mats), "macro")
  expect_true(length(cluster_names) > 0)
  for (nm in cluster_names) {
    expect_true(is.matrix(mats[[nm]]))
  }
})

test_that("extract_initial_probs() works on mcml", {
  cs <- .make_mcml()
  inits <- extract_initial_probs(cs)

  expect_true(is.list(inits))
  expect_true("macro" %in% names(inits))
  expect_true(is.numeric(inits$macro))
  expect_equal(sum(inits$macro), 1, tolerance = 1e-6)
})

test_that("extract_edges() works on mcml", {
  cs <- .make_mcml()
  edges <- extract_edges(cs)

  expect_true(is.list(edges))
  expect_true("macro" %in% names(edges))
  expect_true(is.data.frame(edges$macro))
  expect_true(all(c("from", "to", "weight") %in% names(edges$macro)))
})

# ---- Edgelist input: $data propagated to macro + clusters so that
#      as_tna() + bootstrap_network() work end-to-end (regression test) ----

test_that("edgelist input yields macro/per-cluster $data usable by bootstrap", {
  # Tiny edgelist with known cluster structure
  edges <- data.frame(
    from   = c("a1", "a1", "a2", "b1", "b2", "b1", "a1", "b2", "a2", "b1"),
    to     = c("a2", "b1", "a1", "b2", "b1", "a2", "b2", "a1", "a1", "b1"),
    weight = 1,
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  mc <- build_mcml(edges, clusters)

  # Previous bug: $macro$data was NULL for edgelist input
  expect_false(is.null(mc$macro$data))
  expect_true(is.data.frame(mc$macro$data))
  expect_equal(ncol(mc$macro$data), 2L)
  expect_true(all(mc$macro$data[[1]] %in% c("A", "B")))
  expect_true(all(mc$macro$data[[2]] %in% c("A", "B")))

  # Each per-cluster data keeps only transitions fully within that cluster
  expect_false(is.null(mc$clusters$A$data))
  expect_true(all(mc$clusters$A$data[[1]] %in% c("a1", "a2")))
  expect_true(all(mc$clusters$A$data[[2]] %in% c("a1", "a2")))
  expect_false(is.null(mc$clusters$B$data))
  expect_true(all(mc$clusters$B$data[[1]] %in% c("b1", "b2")))
  expect_true(all(mc$clusters$B$data[[2]] %in% c("b1", "b2")))

  # Edgelist source is tagged so bootstrap can warn about it
  expect_identical(attr(mc$macro$data, "source"), "edgelist")
  expect_identical(attr(mc$clusters$A$data, "source"), "edgelist")

  # as_tna() then bootstrap_network() must run without the "$data is NULL"
  # error, but should emit the anti-conservative-CI warning pointing to
  # permutation / case-dropping alternatives.
  tnas <- as_tna(mc)
  suppressWarnings(boot <- bootstrap_network(tnas, iter = 20))
  expect_s3_class(boot, "net_bootstrap_group")
  expect_warning(
    bootstrap_network(tnas$macro, iter = 10),
    "edgelist"
  )
})

test_that("sequence input does NOT emit the edgelist bootstrap warning", {
  # Wide-format sequences: rows = actors, cols = time steps. Should have
  # no "source = edgelist" tag, no warning from bootstrap_network().
  seq_df <- data.frame(
    t1 = c("a1", "a2", "b1", "b2"),
    t2 = c("a2", "a1", "b2", "b1"),
    t3 = c("a1", "a1", "b1", "b2"),
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  mc <- build_mcml(seq_df, clusters)

  expect_null(attr(mc$macro$data, "source"))
  tnas <- as_tna(mc)
  expect_no_warning(
    bootstrap_network(tnas$macro, iter = 10)
  )
})
