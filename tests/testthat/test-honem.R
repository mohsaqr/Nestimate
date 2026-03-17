# ===========================================================================
# Section 1: Internal — .honem_transition_matrix
# ===========================================================================

test_that(".honem_transition_matrix row-normalizes", {
  mat <- matrix(c(0, 3, 1, 2, 0, 2, 0, 0, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  D <- .honem_transition_matrix(mat)

  expect_equal(D["A", "B"], 3 / 4)
  expect_equal(D["A", "C"], 1 / 4)
  expect_equal(D["B", "A"], 2 / 4)
  expect_equal(D["B", "C"], 2 / 4)
  # Row C: all zeros, stays zero

  expect_equal(sum(D["C", ]), 0)
})

# ===========================================================================
# Section 2: Internal — .honem_neighborhood_matrix
# ===========================================================================

test_that(".honem_neighborhood_matrix produces correct shape", {
  D <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3, 3, byrow = TRUE)
  S <- .honem_neighborhood_matrix(D, max_power = 5L)

  expect_equal(dim(S), c(3L, 3L))
  expect_true(all(is.finite(S)))
})

test_that(".honem_neighborhood_matrix respects max_power", {
  D <- diag(3) * 0.5
  S1 <- .honem_neighborhood_matrix(D, max_power = 1L)
  S5 <- .honem_neighborhood_matrix(D, max_power = 5L)

  # With different max_power, results should differ
  expect_false(all(abs(S1 - S5) < 1e-10))
})

# ===========================================================================
# Section 3: Internal — .honem_svd
# ===========================================================================

test_that(".honem_svd returns correct dimensions", {
  S <- matrix(runif(25), 5, 5, dimnames = list(LETTERS[1:5], LETTERS[1:5]))
  result <- .honem_svd(S, dim = 3L)

  expect_equal(nrow(result$embeddings), 5L)
  expect_equal(ncol(result$embeddings), 3L)
  expect_equal(length(result$singular_values), 3L)
  expect_true(result$explained_variance >= 0 && result$explained_variance <= 1)
})

test_that(".honem_svd caps dim at n-1", {
  S <- matrix(runif(9), 3, 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  result <- .honem_svd(S, dim = 100L)

  expect_equal(ncol(result$embeddings), 2L)  # capped at n-1 = 2
})

# ===========================================================================
# Section 4: build_honem end-to-end
# ===========================================================================

test_that("build_honem returns net_honem class from matrix", {
  mat <- matrix(c(0, 2, 0, 0,
                  0, 0, 3, 0,
                  1, 0, 0, 2,
                  0, 1, 0, 0), 4, 4, byrow = TRUE,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  emb <- build_honem(mat, dim = 2L, max_power = 3L)

  expect_s3_class(emb, "net_honem")
  expect_equal(emb$n_nodes, 4L)
  expect_equal(emb$dim, 2L)
  expect_equal(nrow(emb$embeddings), 4L)
  expect_equal(ncol(emb$embeddings), 2L)
})

test_that("build_honem works with net_hon object", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"),
                c("B", "C", "D", "A"), c("C", "D", "A", "B"))
  hon <- build_hon(trajs, max_order = 2L, method = "hon")
  emb <- build_honem(hon, dim = 3L)

  expect_s3_class(emb, "net_honem")
  expect_equal(emb$n_nodes, hon$n_nodes)
})

test_that("build_honem preserves node names", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE,
                dimnames = list(c("X", "Y"), c("X", "Y")))
  emb <- build_honem(mat, dim = 1L)

  expect_equal(emb$nodes, c("X", "Y"))
  expect_equal(rownames(emb$embeddings), c("X", "Y"))
})

test_that("build_honem rejects invalid input", {
  expect_error(build_honem(42), "net_hon object or a square matrix")
  expect_error(build_honem(matrix(1, 1, 1)), "at least 2 nodes")
})

# ===========================================================================
# Section 5: Embedding quality
# ===========================================================================

test_that("build_honem embeddings reflect network structure", {
  # Create a network with two clusters
  mat <- matrix(0, 6, 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
  # Cluster 1: A, B, C (strong connections)
  mat["A", "B"] <- 5; mat["B", "C"] <- 5; mat["C", "A"] <- 5
  # Cluster 2: D, E, F (strong connections)
  mat["D", "E"] <- 5; mat["E", "F"] <- 5; mat["F", "D"] <- 5
  # Weak cross-cluster link
  mat["C", "D"] <- 1

  emb <- build_honem(mat, dim = 2L, max_power = 5L)

  # Nodes within same cluster should be closer than across clusters
  d_within1 <- sqrt(sum((emb$embeddings["A", ] - emb$embeddings["B", ])^2))
  d_across <- sqrt(sum((emb$embeddings["A", ] - emb$embeddings["E", ])^2))

  expect_true(d_within1 < d_across)
})

# ===========================================================================
# Section 6: S3 methods
# ===========================================================================

test_that("print.net_honem works", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- build_honem(mat, dim = 1L)
  out <- capture.output(print(emb))
  expect_true(any(grepl("HONEM", out)))
})

test_that("summary.net_honem works", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- build_honem(mat, dim = 1L)
  out <- capture.output(summary(emb))
  expect_true(any(grepl("Variance", out)))
})

test_that("plot.net_honem works", {
  mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  emb <- build_honem(mat, dim = 2L)
  expect_no_error(plot(emb))
})
