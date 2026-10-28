# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

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

.hc_all <- c("clique", "Z", "H")

# Hypergraph from a list of hyperedges (character vectors of members)
.hc_from_edges <- function(edges) {
  d <- data.frame(
    player  = unlist(edges),
    session = rep(sprintf("S%02d", seq_along(edges)), lengths(edges)),
    stringsAsFactors = FALSE
  )
  bipartite_groups(d, "player", "session")
}

# Structure ----------------------------------------------------------------

test_that("default returns the clique centrality only", {
  expect_named(hypergraph_centrality(.hc_two_overlapping()), "clique")
})

test_that("all three types return named vectors on a uniform hypergraph", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = .hc_all)
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

# Symmetry / structural validation ----------------------------------------

test_that("symmetric nodes get equal centrality (two overlapping triangles)", {
  # A & B are structurally identical (both in S1 and S2).
  # C & D are structurally identical (each in exactly one session).
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = .hc_all)
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
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = .hc_all)
  # Both should rank A=B above C=D
  expect_equal(rank(cent$Z), rank(cent$H), ignore_attr = TRUE)
})

# Non-negativity, normalization --------------------------------------------

test_that("centralities are non-negative for non-negative hypergraphs", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = .hc_all)
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
  cent <- hypergraph_centrality(hg, type = .hc_all)
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
  c1 <- hypergraph_centrality(hg, type = .hc_all)
  c2 <- hypergraph_centrality(hg, type = .hc_all)
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

# Oracles (2026-10-27 review) ----------------------------------------------

test_that("CEC on a star of pairs is the Perron vector, whatever max_iter", {
  # Clique expansion of {h,a},{h,b},{h,c} is the bipartite star K_{1,3}:
  # power iteration oscillates there; the Perron vector is (sqrt 3, 1, 1, 1).
  hg <- .hc_from_edges(list(c("h", "a"), c("h", "b"), c("h", "c")))
  truth <- c(a = 1, b = 1, c = 1, h = sqrt(3)) / sqrt(6)
  c_even <- hypergraph_centrality(hg, max_iter = 1000L)$clique
  c_odd  <- hypergraph_centrality(hg, max_iter = 1001L)$clique
  expect_equal(c_even[names(truth)], truth, tolerance = 1e-12)
  expect_identical(c_even, c_odd)
})

test_that("CEC equals the leading eigenvector of W on random connected hypergraphs", {
  set.seed(20261027)
  invisible(lapply(seq_len(25), \(r) {
    nodes <- sprintf("n%02d", seq_len(12))
    # a spanning path keeps the clique graph connected
    edges <- c(Map(c, nodes[-12], nodes[-1]),
               replicate(10, sample(nodes, sample(2:5, 1)), simplify = FALSE))
    hg <- .hc_from_edges(edges)
    W <- tcrossprod((hg$incidence != 0) * 1); diag(W) <- 0
    v <- abs(eigen(W, symmetric = TRUE)$vectors[, 1])
    got <- hypergraph_centrality(hg)$clique
    expect_equal(unname(got), v, tolerance = 1e-10)
  }))
})

test_that("Z and H reproduce Benson (2019) Table 2 on the 3-petal sunflower", {
  hg <- .hc_from_edges(list(c("c", "p1", "q1"), c("c", "p2", "q2"),
                            c("c", "p3", "q3")))
  cent <- hypergraph_centrality(hg, type = .hc_all, tol = 1e-13,
                                max_iter = 10000L)
  expect_equal(cent$clique[["c"]] / cent$clique[["p1"]], 2, tolerance = 1e-10)
  expect_equal(cent$Z[["c"]] / cent$Z[["p1"]], sqrt(3), tolerance = 1e-8)
  expect_equal(cent$H[["c"]] / cent$H[["p1"]], 3^(1 / 3), tolerance = 1e-8)
})

test_that("Z and H on mixed hyperedge sizes raise a classed error", {
  hg <- .hc_from_edges(list(c("a", "b"), c("b", "c", "d")))
  expect_error(hypergraph_centrality(hg, type = "Z"),
               class = "nestimate_nonuniform_hypergraph")
  expect_error(hypergraph_centrality(hg, type = "H"),
               class = "nestimate_nonuniform_hypergraph")
  expect_no_error(hypergraph_centrality(hg, type = "clique"))
})

test_that("singletons are ignored by the uniformity check and score zero alone", {
  hg <- .hc_from_edges(list("z", c("a", "b", "c"), c("b", "c", "d")))
  expect_no_error(hypergraph_centrality(hg, type = c("Z", "H")))
  lone <- .hc_from_edges(list("a", "b"))
  cent <- hypergraph_centrality(lone, type = .hc_all)
  for (nm in .hc_all) expect_true(all(cent[[nm]] == 0))
})

test_that("an unconverged tensor iteration warns with a class", {
  expect_warning(
    hypergraph_centrality(.hc_two_overlapping(), type = "Z", max_iter = 1L),
    class = "nestimate_no_converge"
  )
})

test_that("Z solves the eigen-equation on large uniform hyperedges", {
  # 3 petals of 20 members around a hub: Benson's ratios are sqrt(3) for Z
  # and 3^(1/20) for H. An unnormalised shift used to swamp the tiny tensor
  # update and return a ratio of 1 without warning.
  edges <- lapply(0:2, function(j) c("hub", paste0("p", j, "_", 1:19)))
  hg <- .hc_from_edges(edges)
  expect_no_warning(cent <- hypergraph_centrality(hg, type = c("Z", "H")))
  expect_equal(cent$Z[["hub"]] / cent$Z[["p0_1"]], sqrt(3), tolerance = 1e-6)
  expect_equal(cent$H[["hub"]] / cent$H[["p0_1"]], 3^(1 / 20), tolerance = 1e-6)
})
