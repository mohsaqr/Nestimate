# Tests for mcml.R: aggregate_weights, cluster_summary, build_mcml

# ============================================
# aggregate_weights / wagg
# ============================================

test_that("aggregate_weights sum method", {
  expect_equal(aggregate_weights(c(0.5, 0.8, 0.3), "sum"), 1.6)
})

test_that("aggregate_weights mean method", {
  expect_equal(aggregate_weights(c(2, 4, 6), "mean"), 4)
})

test_that("aggregate_weights median method", {
  expect_equal(aggregate_weights(c(1, 3, 5), "median"), 3)
})

test_that("aggregate_weights max method", {
  expect_equal(aggregate_weights(c(1, 5, 3), "max"), 5)
})

test_that("aggregate_weights min method", {
  expect_equal(aggregate_weights(c(1, 5, 3), "min"), 1)
})

test_that("aggregate_weights prod method", {
  expect_equal(aggregate_weights(c(2, 3, 4), "prod"), 24)
})

test_that("aggregate_weights density with n_possible", {
  expect_equal(aggregate_weights(c(1, 2, 3), "density", n_possible = 10), 0.6)
})

test_that("aggregate_weights density without n_possible", {
  expect_equal(aggregate_weights(c(1, 2, 3), "density"), 2)
})

test_that("aggregate_weights geomean method", {
  expect_equal(aggregate_weights(c(4, 9), "geomean"), 6, tolerance = 0.01)
})

test_that("aggregate_weights removes NA and zero", {
  expect_equal(aggregate_weights(c(1, NA, 0, 2), "sum"), 3)
})

test_that("aggregate_weights returns 0 for empty/all-zero input", {
  expect_equal(aggregate_weights(c(0, 0, NA), "sum"), 0)
  expect_equal(aggregate_weights(numeric(0), "mean"), 0)
})

test_that("aggregate_weights errors on unknown method", {
  expect_error(aggregate_weights(c(1, 2), "bogus"), "Unknown method")
})

test_that("wagg is identical to aggregate_weights", {
  expect_identical(wagg, aggregate_weights)
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

test_that("cluster_summary type=tna normalizes rows to 1", {
  mat <- matrix(c(10, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs <- cluster_summary(mat, clusters, type = "tna")
  row_sums <- rowSums(cs$macro$weights)
  expect_equal(unname(row_sums), c(1, 1), tolerance = 1e-10)
})

test_that("cluster_summary type=raw keeps raw values", {
  mat <- matrix(c(10, 2, 3, 8), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  clusters <- list(G1 = "A", G2 = "B")

  cs <- cluster_summary(mat, clusters, type = "raw", method = "sum")
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

test_that("csum is identical to cluster_summary", {
  expect_identical(csum, cluster_summary)
})

test_that("cluster_summary with different methods", {
  mat <- matrix(c(4, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs_sum <- cluster_summary(mat, clusters, method = "sum", type = "raw")
  cs_mean <- cluster_summary(mat, clusters, method = "mean", type = "raw")
  cs_max <- cluster_summary(mat, clusters, method = "max", type = "raw")

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
