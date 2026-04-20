# ---- hypergraph_centrality() tests ---------------------------------------

# Helpers ------------------------------------------------------------------

.hc_two_overlapping <- function() {
  d <- data.frame(
    player  = c("A", "B", "C",  "A", "B", "D"),
    session = c("S1", "S1", "S1", "S2", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  bipartite_groups(d, "player", "session")
}

# Structure ----------------------------------------------------------------

test_that("default returns named list with all three types", {
  cent <- hypergraph_centrality(.hc_two_overlapping())
  expect_named(cent, c("clique", "Z", "H"))
  for (nm in names(cent)) {
    expect_length(cent[[nm]], 4L)
    expect_named(cent[[nm]], c("A", "B", "C", "D"))
    expect_true(all(is.finite(cent[[nm]])))
  }
})

test_that("single type returns list with just that entry", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = "clique")
  expect_named(cent, "clique")
  expect_length(cent$clique, 4L)
})

# CEC vs igraph validation ------------------------------------------------

test_that("clique centrality matches igraph::eigen_centrality on expansion", {
  skip_if_not_installed("igraph")
  hg  <- .hc_two_overlapping()
  ours <- hypergraph_centrality(hg, type = "clique")$clique
  # Recompute via clique expansion + igraph
  net <- clique_expansion(hg)
  W   <- net$weights
  g   <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
  igr <- igraph::eigen_centrality(g, directed = FALSE)$vector
  # Match up to sign and scale
  ours_n <- ours / sqrt(sum(ours^2))
  igr_n  <- igr  / sqrt(sum(igr^2))
  if (sum(ours_n * igr_n) < 0) ours_n <- -ours_n
  expect_equal(ours_n, igr_n[names(ours_n)],
               tolerance = 1e-5, ignore_attr = TRUE)
})

# Symmetry / structural validation ----------------------------------------

test_that("symmetric nodes get equal centrality (two overlapping triangles)", {
  # A & B are structurally identical (both in S1 and S2).
  # C & D are structurally identical (each in exactly one session).
  cent <- hypergraph_centrality(.hc_two_overlapping())
  for (nm in names(cent)) {
    expect_equal(cent[[nm]][["A"]], cent[[nm]][["B"]], tolerance = 1e-6)
    expect_equal(cent[[nm]][["C"]], cent[[nm]][["D"]], tolerance = 1e-6)
    # A/B should have higher centrality than C/D (in 2 edges vs 1)
    expect_gt(cent[[nm]][["A"]], cent[[nm]][["C"]])
  }
})

# ZEC manual validation on small uniform hypergraph -----------------------

test_that("ZEC power-iteration update matches manual formula on one triangle", {
  # Single triangle (A,B,C), uniform k=3.
  # ZEC eigen-equation: λ x_i = Π_{j≠i} x_j  ⇒ by symmetry x_A=x_B=x_C.
  # Thus the fixed point is (a, a, a), any positive a, normalized.
  d <- data.frame(player  = c("A", "B", "C"),
                  session = c("S1", "S1", "S1"),
                  stringsAsFactors = FALSE)
  hg <- bipartite_groups(d, "player", "session")
  cent <- hypergraph_centrality(hg, type = "Z")$Z
  expected <- rep(1 / sqrt(3), 3)
  names(expected) <- c("A", "B", "C")
  expect_equal(cent, expected, tolerance = 1e-6)
})

test_that("ZEC distinguishes hub-like vs peripheral in a 4-node, 2-edge case", {
  # Edges (A,B,C) and (A,B,D). A and B are in both edges.
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = "Z")$Z
  expect_gt(cent[["A"]], cent[["C"]])
  expect_gt(cent[["B"]], cent[["D"]])
})

# HEC on uniform hypergraph ------------------------------------------------

test_that("HEC reproduces ZEC-like ranking on uniform hypergraph", {
  # For a uniform hypergraph the eigenvector direction is the same;
  # only the scale differs (H takes the k-1 root of the ZEC update).
  cent <- hypergraph_centrality(.hc_two_overlapping())
  # Both should rank A=B above C=D
  expect_equal(rank(cent$Z), rank(cent$H), ignore_attr = TRUE)
})

# Non-negativity, normalization --------------------------------------------

test_that("centralities are non-negative for non-negative hypergraphs", {
  cent <- hypergraph_centrality(.hc_two_overlapping())
  for (nm in names(cent)) {
    expect_true(all(cent[[nm]] >= -1e-10))
  }
})

test_that("normalize=TRUE gives unit-L2-norm vectors (clique at least)", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), normalize = TRUE)
  expect_equal(sqrt(sum(cent$clique^2)), 1, tolerance = 1e-6)
})

test_that("normalize=FALSE gives max-abs = 1", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), normalize = FALSE)
  expect_equal(max(abs(cent$clique)), 1, tolerance = 1e-6)
})

# Empty & degenerate -----------------------------------------------------

test_that("empty hypergraph returns all-zero centralities", {
  inc <- matrix(0L, 3L, 0L, dimnames = list(c("A", "B", "C"), NULL))
  hg <- structure(
    list(hyperedges = list(), incidence = inc, nodes = c("A", "B", "C"),
         n_nodes = 3L, n_hyperedges = 0L,
         size_distribution = integer(0), params = list()),
    class = "net_hypergraph")
  cent <- hypergraph_centrality(hg)
  for (nm in names(cent)) {
    expect_true(all(cent[[nm]] == 0))
  }
})

test_that("single size-1 hyperedge contributes nothing (k < 2)", {
  # Hyperedge of size 1 shouldn't participate (no "other member" to pair with)
  d <- data.frame(player  = c("A"),
                  session = c("S1"),
                  stringsAsFactors = FALSE)
  hg <- bipartite_groups(d, "player", "session")
  expect_no_error(hypergraph_centrality(hg, type = c("Z", "H")))
})

# Input validation --------------------------------------------------------

test_that("rejects non-net_hypergraph input", {
  expect_error(hypergraph_centrality(matrix(0, 3, 3)), "net_hypergraph")
})

test_that("unknown type is rejected by match.arg", {
  expect_error(hypergraph_centrality(.hc_two_overlapping(), type = "bogus"))
})

# Reproducibility ---------------------------------------------------------

test_that("power iteration is deterministic given same input", {
  hg <- .hc_two_overlapping()
  c1 <- hypergraph_centrality(hg)
  c2 <- hypergraph_centrality(hg)
  expect_identical(c1, c2)
})

# Integration with bipartite_groups + bundled data ------------------------

test_that("runs on bundled human_long dataset without error", {
  data("human_long", package = "Nestimate")
  # Use a small subset to keep the test fast
  sub <- head(human_long, 500)
  hg  <- bipartite_groups(sub, "code", "session_id")
  cent <- hypergraph_centrality(hg, type = "clique")
  expect_length(cent$clique, hg$n_nodes)
  expect_true(all(is.finite(cent$clique)))
})
