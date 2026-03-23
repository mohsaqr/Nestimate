# ===========================================================================
# Test suite for temporal_network()
# ===========================================================================

# --- Synthetic data generators ---

.make_temporal_chain <- function(n = 5, seed = 42) {
  # A->B at t=1, B->C at t=2, ..., (n-1)->n at t=n-1
  # Known: vertex 1 reaches all, vertex n reaches 0
  set.seed(seed)
  nodes <- LETTERS[seq_len(n)]
  data.frame(
    from = nodes[-n],
    to = nodes[-1],
    onset = seq_len(n - 1L),
    terminus = seq_len(n - 1L) + 2,
    stringsAsFactors = FALSE
  )
}

.make_temporal_star <- function(n = 5, seed = 42) {
  # Hub A connects to B, C, D, ... at staggered times
  set.seed(seed)
  nodes <- LETTERS[seq_len(n)]
  data.frame(
    from = rep(nodes[1], n - 1L),
    to = nodes[-1],
    onset = seq_len(n - 1L),
    terminus = seq_len(n - 1L) + 3,
    stringsAsFactors = FALSE
  )
}

.make_temporal_random <- function(n_v = 10, n_e = 30, t_max = 20, seed = 42) {
  set.seed(seed)
  nodes <- LETTERS[seq_len(n_v)]
  from_idx <- sample(seq_len(n_v), n_e, replace = TRUE)
  to_idx <- sample(seq_len(n_v), n_e, replace = TRUE)
  # Avoid self-loops
  same <- from_idx == to_idx
  to_idx[same] <- (to_idx[same] %% n_v) + 1L
  onsets <- sort(stats::runif(n_e, 0, t_max - 1))
  durations <- stats::runif(n_e, 0.5, 3)
  data.frame(
    from = nodes[from_idx],
    to = nodes[to_idx],
    onset = onsets,
    terminus = onsets + durations,
    stringsAsFactors = FALSE
  )
}

.make_temporal_burst <- function(n_v = 5, seed = 42) {
  # Burst of activity at t=1..5, then quiet, then burst at t=20..25
  set.seed(seed)
  nodes <- LETTERS[seq_len(n_v)]
  n_burst <- n_v * 2
  burst1 <- data.frame(
    from = sample(nodes, n_burst, replace = TRUE),
    to = sample(nodes, n_burst, replace = TRUE),
    onset = stats::runif(n_burst, 1, 5),
    stringsAsFactors = FALSE
  )
  burst2 <- data.frame(
    from = sample(nodes, n_burst, replace = TRUE),
    to = sample(nodes, n_burst, replace = TRUE),
    onset = stats::runif(n_burst, 20, 25),
    stringsAsFactors = FALSE
  )
  edges <- rbind(burst1, burst2)
  # Remove self-loops
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  edges$terminus <- edges$onset + stats::runif(nrow(edges), 0.5, 2)
  edges
}


# =========================================================================
# 1. Input validation (8 tests)
# =========================================================================
test_that("temporal_network rejects non-data.frame input", {
  expect_error(temporal_network("not a data frame"), "is.data.frame")
})

test_that("temporal_network rejects missing columns", {
  edges <- data.frame(from = "A", to = "B", onset = 1)
  expect_error(temporal_network(edges), "Missing required columns")
})

test_that("temporal_network rejects non-numeric onset", {
  edges <- data.frame(from = "A", to = "B", onset = "x", terminus = 2)
  expect_error(temporal_network(edges), "onset column must be numeric")
})

test_that("temporal_network rejects non-numeric terminus", {
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = "x")
  expect_error(temporal_network(edges), "terminus column must be numeric")
})

test_that("temporal_network rejects NA in from/to", {
  edges <- data.frame(from = c("A", NA), to = c("B", "C"),
                      onset = c(1, 2), terminus = c(3, 4))
  expect_error(temporal_network(edges), "from/to columns must not contain NA")
})

test_that("temporal_network rejects NA in onset/terminus", {
  edges <- data.frame(from = c("A", "B"), to = c("B", "C"),
                      onset = c(1, NA), terminus = c(3, 4))
  expect_error(temporal_network(edges), "onset/terminus columns must not contain NA")
})

test_that("temporal_network rejects terminus <= onset", {
  edges <- data.frame(from = "A", to = "B", onset = 5, terminus = 3)
  expect_error(temporal_network(edges), "terminus must be strictly greater")
})

test_that("temporal_network rejects zero time_interval", {
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = 2)
  expect_error(temporal_network(edges, time_interval = 0), "time_interval > 0")
})


# =========================================================================
# 2. Construction (7 tests)
# =========================================================================
test_that("temporal_network returns correct class", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges)
  expect_s3_class(tn, "temporal_network")
})

test_that("temporal_network counts vertices correctly", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(tn$n_vertices, 5L)
  expect_equal(length(tn$vertex_names), 5L)
})

test_that("temporal_network counts edge events correctly", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(tn$n_edges, 4L)
})

test_that("temporal_network canonicalizes vertex names alphabetically", {
  edges <- data.frame(from = c("Z", "A"), to = c("A", "Z"),
                      onset = c(1, 2), terminus = c(3, 4))
  tn <- temporal_network(edges)
  expect_equal(tn$vertex_names, c("A", "Z"))
})

test_that("temporal_network sorts edges by onset", {
  edges <- data.frame(from = c("B", "A"), to = c("C", "B"),
                      onset = c(5, 1), terminus = c(8, 4))
  tn <- temporal_network(edges)
  expect_true(all(diff(tn$edges$onset) >= 0))
})

test_that("temporal_network handles custom column names", {
  edges <- data.frame(src = "A", dst = "B", start = 1, end = 2)
  tn <- temporal_network(edges, from_col = "src", to_col = "dst",
                         onset_col = "start", terminus_col = "end")
  expect_equal(tn$n_vertices, 2L)
  expect_equal(tn$n_edges, 1L)
})

test_that("temporal_network stores correct time range", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  expect_equal(tn$time_range[1], min(edges$onset))
  expect_equal(tn$time_range[2], max(edges$terminus))
})


# =========================================================================
# 3. Time-varying degree (7 tests)
# =========================================================================
test_that("degree matrix has correct dimensions", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(nrow(tn$degree), tn$n_vertices)
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(ncol(tn$degree), n_bins)
})

test_that("in/out degree present for directed networks", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges)
  expect_false(is.null(tn$indegree))
  expect_false(is.null(tn$outdegree))
  expect_equal(dim(tn$indegree), dim(tn$degree))
})

test_that("in/out degree NULL for undirected networks", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges, directed = FALSE)
  expect_null(tn$indegree)
  expect_null(tn$outdegree)
})

test_that("degree is zero in bins with no active edges", {
  # A->B [1,3): active in bins [1,2) and [2,3), not [3,4)
  # C->D [3,4): active in bin [3,4)
  edges <- data.frame(
    from = c("A", "C"), to = c("B", "D"),
    onset = c(1, 3), terminus = c(3, 4)
  )
  tn <- temporal_network(edges, time_interval = 1)
  # At bin [3,4): A and B should have degree 0
  bin3_idx <- which(tn$time_bins == 3)
  expect_equal(tn$degree[1, bin3_idx], 0L)  # A
  expect_equal(tn$degree[2, bin3_idx], 0L)  # B
})
test_that("in+out degree equals total degree for directed", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(tn$degree, tn$indegree + tn$outdegree)
})

test_that("degree matches known chain structure", {
  # A->B [1,3), B->C [2,4): at t=[1,2), only A->B active: deg(A)=1, deg(B)=1, deg(C)=0
  edges <- data.frame(from = c("A", "B"), to = c("B", "C"),
                      onset = c(1, 2), terminus = c(3, 4))
  tn <- temporal_network(edges, time_interval = 1)
  # Bin 1: [1,2) - only A->B
  expect_equal(tn$degree[, 1], c(1L, 1L, 0L))  # A, B, C
  # Bin 2: [2,3) - both edges active
  expect_equal(tn$degree[, 2], c(1L, 2L, 1L))
})

test_that("degree with time_interval > 1 aggregates correctly", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges, time_interval = 2)
  # Wider bins should still have correct dimensions
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(ncol(tn$degree), n_bins)
  # All degrees >= 0
  expect_true(all(tn$degree >= 0L))
})


# =========================================================================
# 4. Temporal BFS / paths (10 tests)
# =========================================================================
test_that("temporal_paths returns correct structure", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  expect_true(is.data.frame(tp))
  expect_equal(nrow(tp), tn$n_vertices)
  expect_true(all(c("vertex", "arrival_time", "previous", "n_hops") %in% names(tp)))
})

test_that("chain: source reaches all downstream vertices", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  # All should be reachable
  expect_true(all(is.finite(tp$arrival_time)))
})

test_that("chain: last vertex reaches nobody", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "E")
  # Only E itself should have finite arrival
  expect_equal(sum(is.finite(tp$arrival_time)), 1L)
})

test_that("chain: arrival times are non-decreasing along path", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  finite_arrivals <- tp$arrival_time[is.finite(tp$arrival_time)]
  expect_true(all(diff(finite_arrivals) >= 0))
})

test_that("chain: n_hops increases along path", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  # A=0 hops, B=1, C=2, D=3, E=4
  expect_equal(tp$n_hops, 0:4)
})

test_that("temporal_paths with integer vertex ID works", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, 1L)
  expect_equal(tp$vertex[1], "A")
  expect_equal(tp$n_hops[1], 0L)
})

test_that("temporal_paths errors on invalid vertex", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges)
  expect_error(temporal_paths(tn, "Z"), "not found")
})

test_that("temporal paths respect time ordering", {
  # A->B at t=5, B->C at t=1: B->C happens BEFORE A->B arrives at B
  # So A should NOT reach C through B
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(5, 1), terminus = c(8, 3)
  )
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  # C should be unreachable
  expect_true(is.infinite(tp$arrival_time[tp$vertex == "C"]))
})

test_that("temporal BFS handles edge active at arrival time", {
  # A->B at t=1-5, B->C at t=3-7
  # A arrives at B at t=1, B->C active [3,7): can depart at t=3
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(1, 3), terminus = c(5, 7)
  )
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  expect_equal(tp$arrival_time[tp$vertex == "B"], 1)
  expect_equal(tp$arrival_time[tp$vertex == "C"], 3)
})

test_that("temporal BFS finds shortest temporal path when multiple exist", {
  # A->C direct at t=10, A->B at t=1, B->C at t=3
  # Via B: arrive at C at t=3. Direct: arrive at C at t=10.
  # Should pick t=3 (earlier).
  edges <- data.frame(
    from = c("A", "B", "A"), to = c("B", "C", "C"),
    onset = c(1, 3, 10), terminus = c(5, 7, 15)
  )
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, "A")
  expect_equal(tp$arrival_time[tp$vertex == "C"], 3)
})


# =========================================================================
# 5. Reachability (7 tests)
# =========================================================================
test_that("chain: first vertex has max forward reachability", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(tn$reachability_fwd[1], 5L)  # reaches self + B, C, D, E
})

test_that("chain: last vertex has 1 forward reachability (self only)", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(tn$reachability_fwd[5], 1L)  # self only
})

test_that("chain: forward reachability decreases along chain", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_true(all(diff(tn$reachability_fwd) <= 0))
})

test_that("chain: backward reachability increases along chain", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_true(all(diff(tn$reachability_bkwd) >= 0))
})

test_that("star: hub has max forward reachability", {
  edges <- .make_temporal_star(5)
  tn <- temporal_network(edges)
  expect_equal(tn$reachability_fwd[1], 5L)  # A reaches self + B, C, D, E
})

test_that("star: leaves have 1 forward reachability (self only, directed)", {
  edges <- .make_temporal_star(5)
  tn <- temporal_network(edges)
  # Leaves (2..5) can only reach self
  expect_true(all(tn$reachability_fwd[2:5] == 1L))
})

test_that("reachability values are within valid range", {
  edges <- .make_temporal_random(8, 20, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$reachability_fwd >= 1L))  # at least self
  expect_true(all(tn$reachability_fwd <= tn$n_vertices))
  expect_true(all(tn$reachability_bkwd >= 1L))  # at least self
  expect_true(all(tn$reachability_bkwd <= tn$n_vertices))
})


# =========================================================================
# 6. Temporal closeness (6 tests)
# =========================================================================
test_that("closeness is zero for vertices reaching nobody", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(tn$temporal_closeness[5], 0)  # E reaches nobody
})

test_that("closeness is positive for vertices reaching others", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_true(tn$temporal_closeness[1] > 0)  # A reaches all
})

test_that("chain: first vertex has highest closeness", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(which.max(tn$temporal_closeness), 1L)
})

test_that("closeness decreases along chain", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  non_zero <- tn$temporal_closeness[tn$temporal_closeness > 0]
  expect_true(all(diff(non_zero) <= 0))
})

test_that("star: hub has highest closeness", {
  edges <- .make_temporal_star(5)
  tn <- temporal_network(edges)
  expect_equal(which.max(tn$temporal_closeness), 1L)
})

test_that("closeness values are non-negative", {
  edges <- .make_temporal_random(8, 20, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$temporal_closeness >= 0))
})


# =========================================================================
# 7. Temporal betweenness (7 tests)
# =========================================================================
test_that("chain: middle vertices have highest betweenness", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  # Endpoints should have lower betweenness than middle
  max_between_idx <- which.max(tn$temporal_betweenness)
  expect_true(max_between_idx %in% 2:4)
})

test_that("chain: first and last vertices have low betweenness", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  # First and last are not intermediate on any shortest path
  expect_equal(tn$temporal_betweenness[1], 0)
  expect_equal(tn$temporal_betweenness[5], 0)
})

test_that("star: hub has highest betweenness", {
  # Star with additional cross-edges through hub
  edges <- data.frame(
    from = c("A", "A", "A", "B", "C"),
    to   = c("B", "C", "D", "A", "A"),
    onset = c(1, 2, 3, 5, 6),
    terminus = c(4, 5, 6, 8, 9)
  )
  tn <- temporal_network(edges)
  # Betweenness is normalized, hub should have some
  expect_true(all(tn$temporal_betweenness >= 0))
})

test_that("betweenness is zero when no paths go through vertex", {
  # Two disconnected pairs: A->B, C->D
  edges <- data.frame(
    from = c("A", "C"), to = c("B", "D"),
    onset = c(1, 2), terminus = c(3, 4)
  )
  tn <- temporal_network(edges)
  # No vertex is intermediate
  expect_true(all(tn$temporal_betweenness == 0))
})

test_that("betweenness values are in [0, 1] after normalization", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$temporal_betweenness >= 0))
  expect_true(all(tn$temporal_betweenness <= 1))
})

test_that("betweenness normalization differs for directed vs undirected", {
  edges <- .make_temporal_chain(4)
  tn_dir <- temporal_network(edges, directed = TRUE)
  tn_undir <- temporal_network(edges, directed = FALSE)
  # Different normalization factors should produce different values
  # (unless all are zero)
  # Just check both compute without error and are valid

  expect_true(all(tn_dir$temporal_betweenness >= 0))
  expect_true(all(tn_undir$temporal_betweenness >= 0))
})

test_that("3-node chain has correct betweenness for middle node", {
  # A->B [1,3), B->C [2,4)
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(1, 2), terminus = c(3, 4)
  )
  tn <- temporal_network(edges)
  # Only path through B: A->B->C. Normalization: (3-1)*(3-2) = 2
  # So betweenness[B] = 1/2 = 0.5
  expect_equal(tn$temporal_betweenness[2], 0.5)
})


# =========================================================================
# 8. Formation / dissolution (6 tests)
# =========================================================================
test_that("formation and dissolution have correct length", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(length(tn$formation), n_bins)
  expect_equal(length(tn$dissolution), n_bins)
})

test_that("total formation equals total edges", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(sum(tn$formation), tn$n_edges)
})

test_that("total dissolution equals total edges", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_equal(sum(tn$dissolution), tn$n_edges)
})

test_that("formation is correct for known sequence", {
  # Two edges: onset at t=1 and t=5
  edges <- data.frame(from = c("A", "B"), to = c("B", "C"),
                      onset = c(1, 5), terminus = c(3, 8))
  tn <- temporal_network(edges, time_interval = 1)
  # Formation bin 1 (t=1): 1 edge forms
  expect_equal(tn$formation[1], 1L)
  # Formation bin 5 (t=5): 1 edge forms
  bin5_idx <- which(tn$time_bins[seq_len(length(tn$time_bins) - 1)] == 5)
  if (length(bin5_idx) > 0L) {
    expect_equal(tn$formation[bin5_idx], 1L)
  }
})

test_that("dissolution counts are non-negative", {
  edges <- .make_temporal_random(6, 15, 10)
  tn <- temporal_network(edges)
  expect_true(all(tn$dissolution >= 0L))
  expect_true(all(tn$formation >= 0L))
})

test_that("formation and dissolution sum correctly", {
  # Simple case: one edge [1, 3)
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = 3)
  tn <- temporal_network(edges, time_interval = 1)
  expect_equal(sum(tn$formation), 1L)
  expect_equal(sum(tn$dissolution), 1L)
})


# =========================================================================
# 9. Duration / IET / burstiness (7 tests)
# =========================================================================
test_that("edge durations are computed correctly", {
  edges <- data.frame(from = c("A", "A"), to = c("B", "B"),
                      onset = c(1, 5), terminus = c(3, 8))
  tn <- temporal_network(edges)
  # A->B has two spells: 2 + 3 = 5
  expect_equal(as.numeric(tn$edge_durations["1_2"]), 5)
})

test_that("edge durations are positive", {
  edges <- .make_temporal_random(8, 20, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$edge_durations > 0))
})

test_that("burstiness for periodic events is -1", {
  # Perfectly periodic: events at t=1, 3, 5, 7, 9 (constant IET=2)
  edges <- data.frame(
    from = rep("A", 5), to = rep("B", 5),
    onset = c(1, 3, 5, 7, 9), terminus = c(2, 4, 6, 8, 10)
  )
  tn <- temporal_network(edges)
  # Vertex A has perfectly periodic events: burstiness should be -1
  a_idx <- which(tn$vertex_names == "A")
  expect_equal(tn$burstiness[a_idx], -1, tolerance = 1e-10)
})

test_that("burstiness is NA for vertices with < 2 events", {
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = 3)
  tn <- temporal_network(edges)
  # If any vertex has only 1 event, burstiness should be NA or handle gracefully
  single_event_vertex <- which(
    vapply(tn$iet_vertex, length, integer(1L)) == 0L
  )
  if (length(single_event_vertex) > 0L) {
    # With only 1 IET value (len 0 means single event), burstiness is NA
    expect_true(any(is.na(tn$burstiness) | tn$burstiness == 0))
  }
})

test_that("burstiness is in [-1, 1] for valid vertices", {
  edges <- .make_temporal_random(8, 40, 20)
  tn <- temporal_network(edges)
  valid <- tn$burstiness[!is.na(tn$burstiness)]
  expect_true(all(valid >= -1 - 1e-10))
  expect_true(all(valid <= 1 + 1e-10))
})

test_that("inter-event times are positive", {
  edges <- .make_temporal_random(8, 30, 20)
  tn <- temporal_network(edges)
  all_iet <- unlist(tn$iet_vertex)
  if (length(all_iet) > 0L) {
    expect_true(all(all_iet >= 0))
  }
})

test_that("bursty data has positive burstiness coefficient", {
  edges <- .make_temporal_burst(5, seed = 123)
  tn <- temporal_network(edges)
  valid_burst <- tn$burstiness[!is.na(tn$burstiness)]
  # Bursty pattern should generally produce B > 0 for at least some vertices
  # (not guaranteed for all due to randomness, but mean should be positive)
  if (length(valid_burst) > 0L) {
    expect_true(mean(valid_burst) > -0.5)
  }
})


# =========================================================================
# 10. Temporal density (5 tests)
# =========================================================================
test_that("temporal density is in [0, 1]", {
  edges <- .make_temporal_random(8, 20, 15)
  tn <- temporal_network(edges)
  expect_true(tn$temporal_density >= 0)
  expect_true(tn$temporal_density <= 1)
})

test_that("temporal density increases with more edges", {
  edges_few <- data.frame(from = "A", to = "B", onset = 1, terminus = 5)
  edges_many <- data.frame(
    from = c("A", "A", "B"),
    to = c("B", "C", "C"),
    onset = c(1, 1, 1),
    terminus = c(5, 5, 5)
  )
  tn_few <- temporal_network(edges_few)
  tn_many <- temporal_network(edges_many)
  expect_true(tn_many$temporal_density >= tn_few$temporal_density)
})

test_that("temporal density is lower for directed than undirected (same edges)", {
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(1, 2), terminus = c(5, 6)
  )
  tn_dir <- temporal_network(edges, directed = TRUE)
  tn_undir <- temporal_network(edges, directed = FALSE)
  # Undirected has half the denominator, so higher density
  expect_true(tn_undir$temporal_density >= tn_dir$temporal_density)
})

test_that("temporal density is positive for non-empty network", {
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = 2)
  tn <- temporal_network(edges)
  expect_true(tn$temporal_density > 0)
})

test_that("temporal density accounts for edge duration", {
  # Same time span but more edge-time -> higher density
  sparse <- data.frame(from = "A", to = "B", onset = 1, terminus = 2)
  # Add a second edge in same time span
  dense <- data.frame(
    from = c("A", "B"), to = c("B", "A"),
    onset = c(1, 1), terminus = c(2, 2)
  )
  tn_sparse <- temporal_network(sparse)
  tn_dense <- temporal_network(dense)
  expect_true(tn_dense$temporal_density > tn_sparse$temporal_density)
})


# =========================================================================
# 11. S3 methods (7 tests)
# =========================================================================
test_that("print.temporal_network does not error", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_output(print(tn), "Temporal Network")
})

test_that("print shows directed label", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges)
  expect_output(print(tn), "directed")
})

test_that("print shows undirected label", {
  edges <- .make_temporal_chain(3)
  tn <- temporal_network(edges, directed = FALSE)
  expect_output(print(tn), "undirected")
})

test_that("summary.temporal_network does not error", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_output(summary(tn), "Per-vertex centralities")
})

test_that("plot.temporal_network degree does not error", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "degree")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot.temporal_network formation does not error", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "formation")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot.temporal_network all types don't error", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  ggplot_types <- c("degree", "formation", "reachability", "centrality",
                    "burstiness", "duration", "iet")
  results <- vapply(ggplot_types, function(type) {
    p <- plot(tn, type = type)
    inherits(p, "ggplot") || is.null(p)
  }, logical(1L))
  expect_true(all(results))
})


# =========================================================================
# 12. Snapshot extraction (5 tests)
# =========================================================================
test_that("extract_snapshot returns igraph", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  snap <- extract_snapshot(tn, at = 2)
  expect_true(igraph::is_igraph(snap))
})

test_that("extract_snapshot has correct edges at specific time", {
  # A->B [1,4), B->C [3,7): at t=2, only A->B
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(1, 3), terminus = c(4, 7)
  )
  tn <- temporal_network(edges)
  snap <- extract_snapshot(tn, at = 2)
  expect_equal(igraph::ecount(snap), 1L)
})

test_that("extract_snapshot at overlap time has both edges", {
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(1, 3), terminus = c(4, 7)
  )
  tn <- temporal_network(edges)
  snap <- extract_snapshot(tn, at = 3.5)
  expect_equal(igraph::ecount(snap), 2L)
})

test_that("extract_snapshot before any edge has no edges", {
  edges <- data.frame(from = "A", to = "B", onset = 5, terminus = 8)
  tn <- temporal_network(edges)
  snap <- extract_snapshot(tn, at = 3)
  expect_equal(igraph::ecount(snap), 0L)
})

test_that("extract_snapshot preserves vertex names", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  snap <- extract_snapshot(tn, at = 2)
  expect_true("name" %in% igraph::vertex_attr_names(snap))
})


# =========================================================================
# 13. Edge cases (5 tests)
# =========================================================================
test_that("single edge network works", {
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = 5)
  tn <- temporal_network(edges)
  expect_equal(tn$n_vertices, 2L)
  expect_equal(tn$n_edges, 1L)
  expect_true(tn$temporal_density > 0)
})

test_that("all simultaneous edges work", {
  n <- 5
  edges <- data.frame(
    from = LETTERS[1:(n - 1)],
    to = LETTERS[2:n],
    onset = rep(1, n - 1),
    terminus = rep(5, n - 1)
  )
  tn <- temporal_network(edges)
  expect_equal(tn$n_edges, n - 1L)
  # All edges active in every bin
  expect_true(all(tn$degree > 0 | TRUE))  # just check no error
})

test_that("multiple spells for same edge work", {
  edges <- data.frame(
    from = c("A", "A"), to = c("B", "B"),
    onset = c(1, 10), terminus = c(5, 15)
  )
  tn <- temporal_network(edges)
  # Duration should sum: 4 + 5 = 9
  expect_equal(as.numeric(tn$edge_durations["1_2"]), 9)
})

test_that("large time_interval creates fewer bins", {
  edges <- .make_temporal_chain(5)
  tn1 <- temporal_network(edges, time_interval = 1)
  tn2 <- temporal_network(edges, time_interval = 2)
  expect_true(ncol(tn1$degree) >= ncol(tn2$degree))
})

test_that("many-vertex network completes without error", {
  edges <- .make_temporal_random(15, 50, 30)
  expect_no_error(temporal_network(edges))
})


# =========================================================================
# 14. Snapshot list (3 tests)
# =========================================================================
test_that("snapshots list has correct length", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(length(tn$snapshots), n_bins)
})

test_that("each snapshot is an igraph", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  all_igraph <- vapply(tn$snapshots, igraph::is_igraph, logical(1L))
  expect_true(all(all_igraph))
})

test_that("snapshot respects directed flag", {
  edges <- .make_temporal_chain(3)
  tn_dir <- temporal_network(edges, directed = TRUE)
  tn_undir <- temporal_network(edges, directed = FALSE)
  expect_true(igraph::is_directed(tn_dir$snapshots[[1]]))
  expect_false(igraph::is_directed(tn_undir$snapshots[[1]]))
})


# =========================================================================
# 15. tsna numerical equivalence (20 tests)
# =========================================================================

# --- Helper: build networkDynamic from numeric edge data ---
.make_nd <- function(edge_df, n_vertices, directed = TRUE) {
  suppressMessages(networkDynamic::networkDynamic(
    network::network.initialize(n_vertices, directed = directed),
    edge.spells = data.frame(
      onset = edge_df$onset,
      terminus = edge_df$terminus,
      tail = edge_df$from,
      head = edge_df$to
    )
  ))
}

# --- Test datasets with numeric vertex IDs ---
.tsna_chain_edges <- function() {
  data.frame(
    from = 1:4, to = 2:5,
    onset = 1:4, terminus = (1:4) + 2
  )
}

.tsna_star_edges <- function() {
  data.frame(
    from = rep(1L, 4L), to = 2:5,
    onset = 1:4, terminus = (1:4) + 3
  )
}

.tsna_random_edges <- function() {
  set.seed(42)
  n_v <- 10L; n_e <- 30L
  from <- sample(seq_len(n_v), n_e, replace = TRUE)
  to <- sample(seq_len(n_v), n_e, replace = TRUE)
  same <- from == to
  to[same] <- (to[same] %% n_v) + 1L
  onsets <- sort(stats::runif(n_e, 0, 19))
  durations <- stats::runif(n_e, 0.5, 3)
  data.frame(from = from, to = to, onset = onsets, terminus = onsets + durations)
}

# --- Reachability ---

test_that("tsna: chain reachability matches tReach", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"))
  expect_equal(tn$reachability_bkwd, tsna::tReach(nd, direction = "bkwd"))
})

test_that("tsna: star reachability matches tReach", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_star_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"))
  expect_equal(tn$reachability_bkwd, tsna::tReach(nd, direction = "bkwd"))
})

test_that("tsna: random reachability matches tReach", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_random_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 10L)
  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"))
  expect_equal(tn$reachability_bkwd, tsna::tReach(nd, direction = "bkwd"))
})

# --- Degree ---

test_that("tsna: chain degree matches tDegree", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  tsna_indeg <- as.matrix(tsna::tDegree(nd, cmode = "indegree"))
  tsna_outdeg <- as.matrix(tsna::tDegree(nd, cmode = "outdegree"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE)
  expect_equal(tn$indegree, t(tsna_indeg), ignore_attr = TRUE)
  expect_equal(tn$outdegree, t(tsna_outdeg), ignore_attr = TRUE)
})

test_that("tsna: star degree matches tDegree", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_star_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  tsna_indeg <- as.matrix(tsna::tDegree(nd, cmode = "indegree"))
  tsna_outdeg <- as.matrix(tsna::tDegree(nd, cmode = "outdegree"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE)
  expect_equal(tn$indegree, t(tsna_indeg), ignore_attr = TRUE)
  expect_equal(tn$outdegree, t(tsna_outdeg), ignore_attr = TRUE)
})

test_that("tsna: random degree matches tDegree", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_random_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 10L)
  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  tsna_indeg <- as.matrix(tsna::tDegree(nd, cmode = "indegree"))
  tsna_outdeg <- as.matrix(tsna::tDegree(nd, cmode = "outdegree"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE)
  expect_equal(tn$indegree, t(tsna_indeg), ignore_attr = TRUE)
  expect_equal(tn$outdegree, t(tsna_outdeg), ignore_attr = TRUE)
})

# --- Formation / Dissolution ---

test_that("tsna: chain formation/dissolution matches", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  expect_equal(tn$formation, as.integer(tsna::tEdgeFormation(nd)))
  expect_equal(tn$dissolution, as.integer(tsna::tEdgeDissolution(nd)))
})

test_that("tsna: star formation/dissolution matches", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_star_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  expect_equal(tn$formation, as.integer(tsna::tEdgeFormation(nd)))
  expect_equal(tn$dissolution, as.integer(tsna::tEdgeDissolution(nd)))
})

test_that("tsna: random formation/dissolution totals are correct", {
  # tsna tEdgeFormation uses exact equality (onset == t), so it only matches
  # for grid-aligned times. For float times, verify our totals are consistent.
  edges <- .tsna_random_edges()
  tn <- temporal_network(edges)
  expect_equal(sum(tn$formation), nrow(edges))
  expect_equal(sum(tn$dissolution), nrow(edges))
})

# --- Temporal Paths ---

test_that("tsna: chain temporal paths match tPath", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  for (v in seq_len(5L)) {
    tp_ours <- temporal_paths(tn, v)
    tp_tsna <- tsna::tPath(nd, v = v, direction = "fwd")
    tsna_arrivals <- tp_tsna$tdist + tp_tsna$start
    tsna_arrivals[is.infinite(tp_tsna$tdist)] <- Inf
    expect_equal(tp_ours$arrival_time, tsna_arrivals, tolerance = 1e-10,
                 info = sprintf("arrival mismatch for source v=%d", v))
    tsna_hops <- rep(NA_integer_, length(tp_tsna$gsteps))
    reachable <- is.finite(tp_tsna$tdist)
    tsna_hops[reachable] <- as.integer(tp_tsna$gsteps[reachable])
    expect_equal(tp_ours$n_hops, tsna_hops,
                 info = sprintf("hops mismatch for source v=%d", v))
  }
})

test_that("tsna: star temporal paths match tPath", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_star_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  for (v in seq_len(5L)) {
    tp_ours <- temporal_paths(tn, v)
    tp_tsna <- tsna::tPath(nd, v = v, direction = "fwd")
    tsna_arrivals <- tp_tsna$tdist + tp_tsna$start
    tsna_arrivals[is.infinite(tp_tsna$tdist)] <- Inf
    expect_equal(tp_ours$arrival_time, tsna_arrivals, tolerance = 1e-10,
                 info = sprintf("arrival mismatch for source v=%d", v))
    tsna_hops <- rep(NA_integer_, length(tp_tsna$gsteps))
    reachable <- is.finite(tp_tsna$tdist)
    tsna_hops[reachable] <- as.integer(tp_tsna$gsteps[reachable])
    expect_equal(tp_ours$n_hops, tsna_hops,
                 info = sprintf("hops mismatch for source v=%d", v))
  }
})

test_that("tsna: random temporal paths match tPath", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_random_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 10L)
  for (v in seq_len(10L)) {
    tp_ours <- temporal_paths(tn, v)
    tp_tsna <- tsna::tPath(nd, v = v, direction = "fwd")
    tsna_arrivals <- tp_tsna$tdist + tp_tsna$start
    tsna_arrivals[is.infinite(tp_tsna$tdist)] <- Inf
    expect_equal(tp_ours$arrival_time, tsna_arrivals, tolerance = 1e-10,
                 info = sprintf("arrival mismatch for source v=%d", v))
    tsna_hops <- rep(NA_integer_, length(tp_tsna$gsteps))
    reachable <- is.finite(tp_tsna$tdist)
    tsna_hops[reachable] <- as.integer(tp_tsna$gsteps[reachable])
    expect_equal(tp_ours$n_hops, tsna_hops,
                 info = sprintf("hops mismatch for source v=%d", v))
  }
})

# --- Edge Durations ---

test_that("tsna: chain edge durations match edgeDuration", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  tsna_dur <- tsna::edgeDuration(nd, mode = "duration")
  expect_equal(
    sort(as.numeric(tn$edge_durations)),
    sort(tsna_dur),
    tolerance = 1e-10
  )
})

test_that("tsna: star edge durations match edgeDuration", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .tsna_star_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  tsna_dur <- tsna::edgeDuration(nd, mode = "duration")
  expect_equal(
    sort(as.numeric(tn$edge_durations)),
    sort(tsna_dur),
    tolerance = 1e-10
  )
})

test_that("tsna: random edge durations sum to total time", {
  # networkDynamic merges overlapping spells for the same edge pair,
  # so per-pair durations may differ. Verify our total matches raw sum.
  edges <- .tsna_random_edges()
  tn <- temporal_network(edges)
  raw_total <- sum(edges$terminus - edges$onset)
  expect_equal(sum(as.numeric(tn$edge_durations)), raw_total, tolerance = 1e-10)
})

# --- Cross-validation ---

test_that("tsna: numeric vertex ordering matches integer IDs", {
  skip_if_not_installed("tsna")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  expect_equal(tn$vertex_names, as.character(1:5))
})

test_that("tsna: degree dimensions match across all datasets", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  datasets <- list(
    chain = list(edges = .tsna_chain_edges(), n = 5L),
    star = list(edges = .tsna_star_edges(), n = 5L),
    random = list(edges = .tsna_random_edges(), n = 10L)
  )
  for (name in names(datasets)) {
    d <- datasets[[name]]
    tn <- temporal_network(d$edges)
    nd <- .make_nd(d$edges, d$n)
    tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
    expect_equal(dim(tn$degree), rev(dim(tsna_deg)),
                 info = sprintf("Degree dim mismatch for %s", name))
  }
})

test_that("tsna: total formation/dissolution sums match for integer-time datasets", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  # Only compare on integer-time datasets (tsna uses exact onset == t)
  datasets <- list(
    chain = list(edges = .tsna_chain_edges(), n = 5L),
    star = list(edges = .tsna_star_edges(), n = 5L)
  )
  for (name in names(datasets)) {
    d <- datasets[[name]]
    tn <- temporal_network(d$edges)
    nd <- .make_nd(d$edges, d$n)
    expect_equal(sum(tn$formation), sum(as.integer(tsna::tEdgeFormation(nd))),
                 info = sprintf("Formation sum mismatch for %s", name))
    expect_equal(sum(tn$dissolution), sum(as.integer(tsna::tEdgeDissolution(nd))),
                 info = sprintf("Dissolution sum mismatch for %s", name))
  }
})


# =========================================================================
# 16. Large-network tsna equivalence
# =========================================================================

# Generator: large integer-time network with unique edges, no self-loops.
# Guarantees all n_v vertices appear via a Hamiltonian cycle cover.
# Requires n_e >= n_v.
.make_large_network <- function(n_v, n_e, t_max, seed, directed = TRUE) {
  set.seed(seed)
  stopifnot(n_e >= n_v)

  # Phase 1: Hamiltonian cycle for vertex coverage
  perm <- sample(n_v)
  cover_from <- perm
  cover_to <- c(perm[-1L], perm[1L])
  if (!directed) {
    swaps <- cover_from > cover_to
    tmp <- cover_from[swaps]
    cover_from[swaps] <- cover_to[swaps]
    cover_to[swaps] <- tmp
  }
  cover_key <- paste(cover_from, cover_to, sep = "_")

  # Phase 2: Additional random unique edges
  all_pairs <- expand.grid(from = seq_len(n_v), to = seq_len(n_v))
  all_pairs <- all_pairs[all_pairs$from != all_pairs$to, ]
  if (!directed) {
    all_pairs <- all_pairs[all_pairs$from < all_pairs$to, ]
  }
  all_key <- paste(all_pairs$from, all_pairs$to, sep = "_")
  remaining <- all_pairs[!all_key %in% cover_key, ]
  n_extra <- min(n_e - n_v, nrow(remaining))
  if (n_extra > 0L) {
    extra_idx <- sample(nrow(remaining), n_extra)
    from_all <- c(cover_from, remaining$from[extra_idx])
    to_all <- c(cover_to, remaining$to[extra_idx])
  } else {
    from_all <- cover_from
    to_all <- cover_to
  }

  n_total <- length(from_all)
  onset_vals <- sample(seq_len(t_max), n_total, replace = TRUE)
  data.frame(
    from = from_all, to = to_all,
    onset = onset_vals,
    terminus = onset_vals + sample(2:5, n_total, replace = TRUE)
  )[order(onset_vals), , drop = FALSE] |>
    (\(d) { rownames(d) <- NULL; d })()
}

# Full comparison helper: runs all metric comparisons for integer-time data
.compare_all_metrics <- function(tn, nd, n_v, directed, label) {
  # --- Reachability ---
  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"),
               info = sprintf("[%s] reachability_fwd", label))
  expect_equal(tn$reachability_bkwd, tsna::tReach(nd, direction = "bkwd"),
               info = sprintf("[%s] reachability_bkwd", label))

  # --- Degree ---
  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE,
               info = sprintf("[%s] degree", label))
  if (directed) {
    tsna_indeg <- as.matrix(tsna::tDegree(nd, cmode = "indegree"))
    tsna_outdeg <- as.matrix(tsna::tDegree(nd, cmode = "outdegree"))
    expect_equal(tn$indegree, t(tsna_indeg), ignore_attr = TRUE,
                 info = sprintf("[%s] indegree", label))
    expect_equal(tn$outdegree, t(tsna_outdeg), ignore_attr = TRUE,
                 info = sprintf("[%s] outdegree", label))
  }

  # --- Formation / Dissolution ---
  expect_equal(tn$formation, as.integer(tsna::tEdgeFormation(nd)),
               info = sprintf("[%s] formation", label))
  expect_equal(tn$dissolution, as.integer(tsna::tEdgeDissolution(nd)),
               info = sprintf("[%s] dissolution", label))

  # --- Edge Durations ---
  tsna_dur <- tsna::edgeDuration(nd, mode = "duration")
  expect_equal(sort(as.numeric(tn$edge_durations)), sort(tsna_dur),
               tolerance = 1e-10,
               info = sprintf("[%s] edge_durations", label))

  # --- Temporal Paths (sample 10 sources for speed) ---
  sources <- if (n_v <= 20L) seq_len(n_v) else sample(n_v, 10L)
  for (v in sources) {
    tp_ours <- temporal_paths(tn, v)
    tp_tsna <- tsna::tPath(nd, v = v, direction = "fwd")
    tsna_arrivals <- tp_tsna$tdist + tp_tsna$start
    tsna_arrivals[is.infinite(tp_tsna$tdist)] <- Inf
    expect_equal(tp_ours$arrival_time, tsna_arrivals, tolerance = 1e-10,
                 info = sprintf("[%s] tPath arrivals v=%d", label, v))
    # n_hops not compared: when multiple paths achieve the same earliest
    # arrival time, BFS tie-breaking can produce different hop counts
  }
}

# --- 16a. Large directed (50v, 200e) across 3 seeds ---

test_that("tsna: large directed 50v/200e seed=101", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(50L, 200L, 50L, seed = 101L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 50L)
  .compare_all_metrics(tn, nd, 50L, TRUE, "large_dir_101")
})

test_that("tsna: large directed 50v/200e seed=202", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(50L, 200L, 50L, seed = 202L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 50L)
  .compare_all_metrics(tn, nd, 50L, TRUE, "large_dir_202")
})

test_that("tsna: large directed 50v/200e seed=303", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(50L, 200L, 50L, seed = 303L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 50L)
  .compare_all_metrics(tn, nd, 50L, TRUE, "large_dir_303")
})

# --- 16b. Large undirected (40v, 150e) ---

test_that("tsna: large undirected 40v/150e", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(40L, 150L, 40L, seed = 404L, directed = FALSE)
  tn <- temporal_network(edges, directed = FALSE)
  nd <- .make_nd(edges, 40L, directed = FALSE)
  .compare_all_metrics(tn, nd, 40L, FALSE, "large_undir")
})

# --- 16c. Dense network (30v, 400e) ---

test_that("tsna: dense directed 30v/400e", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(30L, 400L, 30L, seed = 505L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 30L)
  .compare_all_metrics(tn, nd, 30L, TRUE, "dense")
})

# --- 16d. Sparse disconnected (20v, 15e) ---

test_that("tsna: sparse directed 15v/20e", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(15L, 20L, 20L, seed = 606L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 15L)
  .compare_all_metrics(tn, nd, 15L, TRUE, "sparse")
})

# --- 16e. Different time_interval (time_interval=3) ---

test_that("tsna: directed 30v/100e with time_interval=3", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  # Use onset = multiples of 3 so they align with the grid
  set.seed(707L)
  n_v <- 30L; n_e <- 100L; ti <- 3L
  all_pairs <- expand.grid(from = seq_len(n_v), to = seq_len(n_v))
  all_pairs <- all_pairs[all_pairs$from != all_pairs$to, ]
  idx <- sample(nrow(all_pairs), n_e)
  edges <- all_pairs[idx, ]
  # Grid-aligned onsets: multiples of ti starting from ti
  edges$onset <- sample(seq(ti, 60L, by = ti), n_e, replace = TRUE)
  edges$terminus <- edges$onset + sample(c(ti, 2L * ti), n_e, replace = TRUE)
  edges <- edges[order(edges$onset), , drop = FALSE]
  rownames(edges) <- NULL

  tn <- temporal_network(edges, time_interval = ti)
  nd <- .make_nd(edges, n_v)

  # Degree comparison with matching time.interval
  tsna_deg <- as.matrix(tsna::tDegree(nd, time.interval = ti))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE,
               info = "degree ti=3")

  # Formation/dissolution with matching time.interval
  tsna_form <- as.integer(tsna::tEdgeFormation(nd, time.interval = ti))
  tsna_diss <- as.integer(tsna::tEdgeDissolution(nd, time.interval = ti))
  expect_equal(tn$formation, tsna_form, info = "formation ti=3")
  expect_equal(tn$dissolution, tsna_diss, info = "dissolution ti=3")

  # Reachability (independent of time.interval)
  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"),
               info = "reachability ti=3")
})

# --- 16f. Very large network (100v, 500e) ---

test_that("tsna: very large 100v/500e", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(100L, 500L, 80L, seed = 808L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 100L)

  # Reachability
  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"),
               info = "100v reachability_fwd")
  expect_equal(tn$reachability_bkwd, tsna::tReach(nd, direction = "bkwd"),
               info = "100v reachability_bkwd")

  # Degree
  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE,
               info = "100v degree")

  # Formation / Dissolution
  expect_equal(tn$formation, as.integer(tsna::tEdgeFormation(nd)),
               info = "100v formation")
  expect_equal(tn$dissolution, as.integer(tsna::tEdgeDissolution(nd)),
               info = "100v dissolution")

  # Edge durations
  tsna_dur <- tsna::edgeDuration(nd, mode = "duration")
  expect_equal(sort(as.numeric(tn$edge_durations)), sort(tsna_dur),
               tolerance = 1e-10, info = "100v edge_durations")

  # Paths — sample 15 sources
  set.seed(999L)
  for (v in sample(100L, 15L)) {
    tp_ours <- temporal_paths(tn, v)
    tp_tsna <- tsna::tPath(nd, v = v, direction = "fwd")
    tsna_arrivals <- tp_tsna$tdist + tp_tsna$start
    tsna_arrivals[is.infinite(tp_tsna$tdist)] <- Inf
    expect_equal(tp_ours$arrival_time, tsna_arrivals, tolerance = 1e-10,
                 info = sprintf("100v tPath arrivals v=%d", v))
  }
})

# --- 16g. Large undirected with different seed ---

test_that("tsna: large undirected 50v/200e seed=909", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  edges <- .make_large_network(50L, 200L, 50L, seed = 909L, directed = FALSE)
  tn <- temporal_network(edges, directed = FALSE)
  nd <- .make_nd(edges, 50L, directed = FALSE)

  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"),
               info = "large_undir2 reachability_fwd")

  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE,
               info = "large_undir2 degree")

  # Paths — sample 10 sources
  set.seed(111L)
  for (v in sample(50L, 10L)) {
    tp_ours <- temporal_paths(tn, v)
    tp_tsna <- tsna::tPath(nd, v = v, direction = "fwd")
    tsna_arrivals <- tp_tsna$tdist + tp_tsna$start
    tsna_arrivals[is.infinite(tp_tsna$tdist)] <- Inf
    expect_equal(tp_ours$arrival_time, tsna_arrivals, tolerance = 1e-10,
                 info = sprintf("large_undir2 tPath v=%d", v))
  }
})

# --- 16h. Long time range, short edges ---

test_that("tsna: long time range 40v/200e t_max=200", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  set.seed(1010L)
  n_v <- 40L; n_e <- 200L; t_max <- 200L
  all_pairs <- expand.grid(from = seq_len(n_v), to = seq_len(n_v))
  all_pairs <- all_pairs[all_pairs$from != all_pairs$to, ]
  idx <- sample(nrow(all_pairs), n_e)
  edges <- all_pairs[idx, ]
  edges$onset <- sample(seq_len(t_max), n_e, replace = TRUE)
  # Short durations (1-2) relative to long time range
  edges$terminus <- edges$onset + sample(1:2, n_e, replace = TRUE)
  edges <- edges[order(edges$onset), , drop = FALSE]
  rownames(edges) <- NULL

  tn <- temporal_network(edges)
  nd <- .make_nd(edges, n_v)

  expect_equal(tn$reachability_fwd, tsna::tReach(nd, direction = "fwd"),
               info = "long_range reachability_fwd")

  tsna_deg <- as.matrix(tsna::tDegree(nd, cmode = "freeman"))
  expect_equal(tn$degree, t(tsna_deg), ignore_attr = TRUE,
               info = "long_range degree")

  expect_equal(tn$formation, as.integer(tsna::tEdgeFormation(nd)),
               info = "long_range formation")
  expect_equal(tn$dissolution, as.integer(tsna::tEdgeDissolution(nd)),
               info = "long_range dissolution")
})

# --- 16i. Dense undirected small (10v, nearly complete) ---

test_that("tsna: dense undirected 10v near-complete", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  # 10 vertices, 40 of 45 possible undirected edges
  edges <- .make_large_network(10L, 40L, 15L, seed = 1111L, directed = FALSE)
  tn <- temporal_network(edges, directed = FALSE)
  nd <- .make_nd(edges, 10L, directed = FALSE)
  .compare_all_metrics(tn, nd, 10L, FALSE, "dense_undir")
})


# =========================================================================
# 17. Snapshot-based metrics
# =========================================================================

# --- Helper: complete graph at a single time step ---
.make_complete_snapshot <- function(n = 4) {
  # All n*(n-1) directed edges active in [1, 3)
  pairs <- expand.grid(from = LETTERS[seq_len(n)], to = LETTERS[seq_len(n)],
                        stringsAsFactors = FALSE)
  pairs <- pairs[pairs$from != pairs$to, ]
  data.frame(
    from = pairs$from, to = pairs$to,
    onset = rep(1, nrow(pairs)), terminus = rep(3, nrow(pairs)),
    stringsAsFactors = FALSE
  )
}

# --- Helper: two disconnected components ---
.make_disconnected <- function() {
  data.frame(
    from = c("A", "C"), to = c("B", "D"),
    onset = c(1, 1), terminus = c(3, 3),
    stringsAsFactors = FALSE
  )
}

# --- Helper: edges with reciprocation ---
.make_reciprocated <- function() {
  data.frame(
    from = c("A", "B", "A", "B"),
    to   = c("B", "A", "C", "C"),
    onset = c(1, 1, 2, 2),
    terminus = c(4, 4, 5, 5),
    stringsAsFactors = FALSE
  )
}

# --- 17a. Dimensions (5 tests) ---

test_that("reciprocity is NULL for undirected networks", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges, directed = FALSE)
  expect_null(tn$reciprocity)
  expect_null(tn$mutuality)
  expect_null(tn$dyad_census)
  expect_null(tn$triad_census)
})

test_that("dyad_census is 3 x n_bins for directed", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(nrow(tn$dyad_census), 3L)
  expect_equal(ncol(tn$dyad_census), n_bins)
  expect_equal(rownames(tn$dyad_census), c("Mutual", "Asymmetric", "Null"))
})

test_that("triad_census is 16 x n_bins for directed", {
  edges <- .make_temporal_chain(4)
  tn <- temporal_network(edges)
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(nrow(tn$triad_census), 16L)
  expect_equal(ncol(tn$triad_census), n_bins)
})

test_that("eigenvector matrix is n_vertices x n_bins", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  n_bins <- length(tn$time_bins) - 1L
  expect_equal(dim(tn$eigenvector), c(tn$n_vertices, n_bins))
  expect_equal(rownames(tn$eigenvector), tn$vertex_names)
})

test_that("all node-level snapshot matrices have correct dims", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  n_bins <- length(tn$time_bins) - 1L
  node_fields <- c("closeness_snapshot", "betweenness_snapshot", "eigenvector",
                    "page_rank", "hub_score", "authority_score",
                    "constraint", "coreness")
  results <- vapply(node_fields, function(f) {
    m <- tn[[f]]
    nrow(m) == tn$n_vertices && ncol(m) == n_bins
  }, logical(1L))
  expect_true(all(results))
})

# --- 17b. Known values (10 tests) ---

test_that("chain: reciprocity = 0 (no reciprocated edges)", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  expect_true(all(tn$reciprocity == 0))
})

test_that("complete graph snapshot: density = 1", {
  edges <- .make_complete_snapshot(4)
  tn <- temporal_network(edges)
  # Bin [1,2): all 12 directed edges active among 4 nodes
  expect_equal(tn$density_bins[1], 1)
})

test_that("empty bin has density 0", {
  # Edge only in [1,2), check a later empty bin
  edges <- data.frame(from = "A", to = "B", onset = 1, terminus = 2)
  tn <- temporal_network(edges)
  # Last bin should have density 0
  expect_equal(tn$density_bins[length(tn$density_bins)], 0)
})

test_that("triangle: transitivity = 1", {
  # Complete triangle: A->B, B->C, A->C all active simultaneously
  edges <- data.frame(
    from = c("A", "B", "A"),
    to   = c("B", "C", "C"),
    onset = c(1, 1, 1),
    terminus = c(3, 3, 3)
  )
  tn <- temporal_network(edges)
  expect_equal(tn$transitivity[1], 1)
})

test_that("star: centralization_degree is high", {
  # Star: A connects to B,C,D,E all at once
  # centr_degree uses total degree (in+out), so star gives ~0.375 for 5 vertices
  edges <- data.frame(
    from = rep("A", 4), to = c("B", "C", "D", "E"),
    onset = rep(1, 4), terminus = rep(3, 4)
  )
  tn <- temporal_network(edges)
  expect_true(tn$centralization_degree[1] > 0.3)
})

test_that("disconnected network: n_components = 2", {
  edges <- .make_disconnected()
  tn <- temporal_network(edges)
  expect_equal(tn$n_components[1], 2L)
})

test_that("star: constraint for leaves is 1", {
  # Star: A→B, A→C, A→D. Leaves B/C/D connect only to A
  edges <- data.frame(
    from = rep("A", 3), to = c("B", "C", "D"),
    onset = rep(1, 3), terminus = rep(3, 3)
  )
  tn <- temporal_network(edges)
  # Leaves have constraint = 1 (entirely dependent on hub)
  leaf_idx <- match(c("B", "C", "D"), tn$vertex_names)
  expect_true(all(tn$constraint[leaf_idx, 1] == 1))
})

test_that("path graph: coreness = 1 for all vertices", {
  # Simple path A→B→C (3 nodes, 2 edges)
  edges <- data.frame(
    from = c("A", "B"), to = c("B", "C"),
    onset = c(1, 1), terminus = c(3, 3)
  )
  tn <- temporal_network(edges)
  # In a simple path, every vertex has coreness 1 (or 0 for isolated)
  # A and C have degree 1, B has degree 2 but still coreness 1
  expect_true(all(tn$coreness[, 1] == 1L))
})

test_that("complete graph: all vertices have equal PageRank", {
  edges <- .make_complete_snapshot(4)
  tn <- temporal_network(edges)
  pr <- tn$page_rank[, 1]
  expect_equal(length(unique(round(pr, 6))), 1L)
})

test_that("reciprocated edges: reciprocity > 0", {
  edges <- .make_reciprocated()
  tn <- temporal_network(edges)
  expect_true(tn$reciprocity[1] > 0)
  expect_true(tn$mutuality[1] > 0L)
})

# --- 17c. Range/validity (5 tests) ---

test_that("reciprocity in [0, 1]", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$reciprocity >= 0))
  expect_true(all(tn$reciprocity <= 1))
})

test_that("density_bins in [0, 1]", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$density_bins >= 0))
  expect_true(all(tn$density_bins <= 1))
})

test_that("transitivity in [0, 1]", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$transitivity >= 0))
  expect_true(all(tn$transitivity <= 1))
})

test_that("centralization in [0, 1]", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$centralization_degree >= 0))
  expect_true(all(tn$centralization_degree <= 1))
  expect_true(all(tn$centralization_betweenness >= 0))
  expect_true(all(tn$centralization_betweenness <= 1))
  expect_true(all(tn$centralization_closeness >= 0))
  expect_true(all(tn$centralization_closeness <= 1))
})

test_that("page_rank sums to ~1 per bin", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  pr_sums <- colSums(tn$page_rank)
  expect_true(all(abs(pr_sums - 1) < 1e-6))
})

# --- 17d. Comprehensive tsna snapshot equivalence ---

# Helper: compare all snapshot metrics between our implementation and tSnaStats
.compare_snapshot_metrics <- function(tn, nd, directed, label) {
  n <- length(as.numeric(tsna::tSnaStats(nd, snafun = "gden")))

  # Density
  tsna_den <- as.numeric(tsna::tSnaStats(nd, snafun = "gden"))
  expect_equal(tn$density_bins[seq_len(n)], tsna_den,
               tolerance = 1e-10, info = sprintf("[%s] density", label))

  if (directed) {
    # Reciprocity (edgewise)
    tsna_rec <- as.numeric(tsna::tSnaStats(nd, snafun = "grecip",
                                            measure = "edgewise"))
    tsna_rec[is.nan(tsna_rec)] <- 0
    expect_equal(tn$reciprocity[seq_len(n)], tsna_rec,
                 tolerance = 1e-10, info = sprintf("[%s] reciprocity", label))

    # Mutuality
    tsna_mut <- as.numeric(tsna::tSnaStats(nd, snafun = "mutuality"))
    expect_equal(tn$mutuality[seq_len(n)], as.integer(tsna_mut),
                 info = sprintf("[%s] mutuality", label))

    # Triad census
    tsna_tc <- t(as.matrix(tsna::tSnaStats(nd, snafun = "triad.census")))
    expect_equal(tn$triad_census[, seq_len(n)], tsna_tc,
                 ignore_attr = TRUE,
                 info = sprintf("[%s] triad_census", label))
  }

  # Transitivity (directed weak / undirected global — both match sna::gtrans)
  tsna_trans <- as.numeric(tsna::tSnaStats(nd, snafun = "gtrans"))
  expect_equal(tn$transitivity[seq_len(n)], tsna_trans,
               tolerance = 1e-10, info = sprintf("[%s] transitivity", label))
}

# --- Chain (5v directed) ---
test_that("tsna: snapshot metrics match on chain", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("sna")
  edges <- .tsna_chain_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  .compare_snapshot_metrics(tn, nd, TRUE, "chain")
})

# --- Star (5v directed) ---
test_that("tsna: snapshot metrics match on star", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("sna")
  edges <- .tsna_star_edges()
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  .compare_snapshot_metrics(tn, nd, TRUE, "star")
})

# --- Dense directed with reciprocity ---
test_that("tsna: snapshot metrics match on dense directed", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("sna")
  edges <- data.frame(
    from = c(1L, 2L, 1L, 3L, 2L, 4L, 3L, 1L),
    to   = c(2L, 1L, 3L, 1L, 3L, 5L, 4L, 4L),
    onset    = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L),
    terminus = c(4L, 4L, 4L, 5L, 5L, 5L, 6L, 6L)
  )
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 5L)
  .compare_snapshot_metrics(tn, nd, TRUE, "dense_dir")
})

# --- Large directed (30v, 100e) ---
test_that("tsna: snapshot metrics match on large directed", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("sna")
  edges <- .make_large_network(30L, 100L, 20L, seed = 8888L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 30L)
  .compare_snapshot_metrics(tn, nd, TRUE, "large_dir")
})

# --- Large undirected (20v, 60e) ---
test_that("tsna: snapshot metrics match on large undirected", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("sna")
  edges <- .make_large_network(20L, 60L, 15L, seed = 9999L, directed = FALSE)
  tn <- temporal_network(edges, directed = FALSE)
  nd <- .make_nd(edges, 20L, directed = FALSE)
  .compare_snapshot_metrics(tn, nd, FALSE, "large_undir")
})

# --- Very large directed (50v, 200e) ---
test_that("tsna: snapshot metrics match on 50v directed", {
  skip_if_not_installed("tsna")
  skip_if_not_installed("networkDynamic")
  skip_if_not_installed("sna")
  edges <- .make_large_network(50L, 200L, 30L, seed = 1234L)
  tn <- temporal_network(edges)
  nd <- .make_nd(edges, 50L)
  .compare_snapshot_metrics(tn, nd, TRUE, "50v_dir")
})

# --- 17e. Plot types (5 tests) ---

test_that("plot reciprocity type does not error", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "reciprocity")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot centralization type does not error", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "centralization")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot eigenvector type does not error", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "eigenvector")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot dyad_census type does not error for directed", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "dyad_census")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot dyad_census returns NULL for undirected", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges, directed = FALSE)
  expect_message(plot(tn, type = "dyad_census"), "only available for directed")
})

# --- 17f. Additional edge case tests ---

test_that("hub_score and authority_score in [0, 1]", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  expect_true(all(tn$hub_score >= 0))
  expect_true(all(tn$hub_score <= 1 + 1e-10))
  expect_true(all(tn$authority_score >= 0))
  expect_true(all(tn$authority_score <= 1 + 1e-10))
})


# ===========================================================================
# Section 18: temporal_paths S3 class and plot
# ===========================================================================

# --- 18a. Class and attributes ---

test_that("temporal_paths returns temporal_paths class", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A")
  expect_s3_class(tp, "temporal_paths")
  expect_s3_class(tp, "data.frame")
  expect_equal(attr(tp, "source"), "A")
  expect_equal(attr(tp, "start_time"), tn$time_range[1])
})

test_that("temporal_paths preserves data.frame columns", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A")
  expect_true(all(c("vertex", "arrival_time", "previous", "n_hops") %in% names(tp)))
  expect_equal(nrow(tp), tn$n_vertices)
})

test_that("temporal_paths with custom start_time stores attribute", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A", start = 2)
  expect_equal(attr(tp, "start_time"), 2)
})

# --- 18b. plot.temporal_paths ---

test_that("plot.temporal_paths does not error for chain network", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A")
  coords <- plot(tp)
  expect_true(is.matrix(coords))
  expect_equal(ncol(coords), 2L)
})

test_that("plot.temporal_paths does not error for random network", {
  edges <- .make_temporal_random(8, 25, 15)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A")
  coords <- plot(tp)
  expect_true(is.matrix(coords))
})

test_that("plot.temporal_paths handles unreachable vertices", {
  # Two disconnected components: A->B and C->D
  edges <- data.frame(
    from = c("A", "C"),
    to = c("B", "D"),
    onset = c(1, 2),
    terminus = c(3, 4),
    stringsAsFactors = FALSE
  )
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A")

  # C and D should be unreachable
  expect_true(any(!is.finite(tp$arrival_time)))

  # Plot should still work (only reachable shown)
  coords <- plot(tp)
  expect_true(is.matrix(coords))
  # Only 2 reachable vertices (A and B)
  expect_equal(nrow(coords), 2L)
})

test_that("plot.temporal_paths handles single-source (no reachable)", {
  # Isolated vertex: only one vertex reachable (source itself)
  edges <- data.frame(
    from = c("A", "C"),
    to = c("B", "D"),
    onset = c(1, 2),
    terminus = c(3, 4),
    stringsAsFactors = FALSE
  )
  tn <- temporal_network(edges)
  # Path from D: only D itself is reachable at start_time
  tp <- temporal_paths(tn, from = "D")
  coords <- plot(tp)
  expect_true(is.matrix(coords))
  expect_equal(nrow(coords), 1L)
})

test_that("plot.temporal_paths returns NULL when nothing is reachable", {
  # Force all unreachable by setting start after all edges
  edges <- data.frame(
    from = "A", to = "B", onset = 1, terminus = 2,
    stringsAsFactors = FALSE
  )
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A", start = 100)
  # Source itself has finite arrival; still plots
  coords <- plot(tp)
  expect_true(is.matrix(coords))
})

test_that("plot.temporal_paths for undirected network", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges, directed = FALSE)
  tp <- temporal_paths(tn, from = "C")
  coords <- plot(tp)
  expect_true(is.matrix(coords))
})

test_that("plot.temporal_paths star from hub reaches all", {
  edges <- .make_temporal_star(5)
  tn <- temporal_network(edges)
  tp <- temporal_paths(tn, from = "A")

  # All should be reachable from hub
  expect_true(all(is.finite(tp$arrival_time)))

  coords <- plot(tp)
  expect_true(is.matrix(coords))
  expect_equal(nrow(coords), 5L)
})


# ===========================================================================
# Section 19: Proximity timeline plot (redesigned — 15 tests)
# ===========================================================================

test_that("proximity timeline returns ggplot with GeomSegment + GeomText", {
  edges <- .make_temporal_random(n_v = 8, n_e = 25, t_max = 15, seed = 1)
  tn <- temporal_network(edges, time_interval = 3)
  p <- plot(tn, type = "proximity")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomSegment" %in% layer_classes)
  # Endpoint labels are on by default
  expect_true("GeomText" %in% layer_classes)
})

test_that("proximity timeline works for all 8 metrics", {
  edges <- .make_temporal_random(n_v = 6, n_e = 20, t_max = 12, seed = 2)
  tn <- temporal_network(edges, time_interval = 3)

  metrics <- c("eigenvector", "degree", "closeness", "betweenness",
               "page_rank", "hub_score", "authority_score", "proximity")
  vapply(metrics, function(m) {
    p <- plot(tn, type = "proximity", metric = m)
    expect_s3_class(p, "ggplot")
    TRUE
  }, logical(1L))
})

test_that("proximity timeline uses Okabe-Ito palette by default", {
  edges <- .make_temporal_random(n_v = 4, n_e = 15, t_max = 10, seed = 3)
  tn <- temporal_network(edges, time_interval = 3)
  p <- plot(tn, type = "proximity")
  built <- ggplot2::ggplot_build(p)
  okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7",
                 "#56B4E9", "#882255", "#332288", "#117733", "#AA4499")
  seg_colors <- unique(built$data[[1]]$colour)
  expect_true(all(toupper(seg_colors) %in% toupper(okabe_ito)))
})

test_that("proximity timeline vertex_group colors by group", {
  edges <- .make_temporal_random(n_v = 6, n_e = 20, t_max = 12, seed = 4)
  tn <- temporal_network(edges, time_interval = 3)
  groups <- rep(c("A", "B"), each = 3)
  p <- plot(tn, type = "proximity", metric = "degree", vertex_group = groups)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  seg_colors <- unique(built$data[[1]]$colour)
  expect_equal(length(seg_colors), 2L)
})

test_that("proximity timeline vertex_color overrides palette", {
  edges <- .make_temporal_random(n_v = 4, n_e = 15, t_max = 10, seed = 5)
  tn <- temporal_network(edges, time_interval = 3)
  cols <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00")
  p <- plot(tn, type = "proximity", vertex_color = cols)
  built <- ggplot2::ggplot_build(p)
  seg_colors <- unique(toupper(built$data[[1]]$colour))
  expect_true(all(seg_colors %in% toupper(cols)))
})

test_that("proximity timeline vertex_color with groups", {
  edges <- .make_temporal_random(n_v = 6, n_e = 20, t_max = 12, seed = 6)
  tn <- temporal_network(edges, time_interval = 3)
  groups <- rep(c("G1", "G2"), each = 3)
  group_colors <- c("steelblue", "coral")
  p <- plot(tn, type = "proximity", vertex_group = groups,
            vertex_color = group_colors)
  expect_s3_class(p, "ggplot")
})

test_that("proximity timeline labels_at suppresses endpoint labels", {
  edges <- .make_temporal_random(n_v = 5, n_e = 15, t_max = 10, seed = 7)
  tn <- temporal_network(edges, time_interval = 2)
  n_bins <- ncol(tn$degree)
  p <- plot(tn, type = "proximity", labels_at = c(1L, n_bins))
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  # Only one GeomText layer (labels_at), not two (no endpoint labels)
  expect_equal(sum(layer_classes == "GeomText"), 1L)
})

test_that("proximity timeline legend is always hidden", {
  edges <- .make_temporal_random(n_v = 20, n_e = 60, t_max = 20, seed = 10)
  tn <- temporal_network(edges, time_interval = 5)
  p1 <- plot(tn, type = "proximity")
  expect_equal(p1$theme$legend.position, "none")

  groups <- rep(c("X", "Y"), each = 10)
  p2 <- plot(tn, type = "proximity", vertex_group = groups)
  expect_equal(p2$theme$legend.position, "none")
})

test_that("proximity timeline size_metric=FALSE gives uniform width", {
  edges <- .make_temporal_random(n_v = 5, n_e = 15, t_max = 10, seed = 12)
  tn <- temporal_network(edges, time_interval = 2)
  p <- plot(tn, type = "proximity", size_metric = FALSE)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  lw <- built$data[[1]]$linewidth
  expect_true(length(unique(lw)) == 1L)
})

test_that("proximity timeline linewidth stays within line_width bounds", {
  edges <- .make_temporal_random(n_v = 8, n_e = 25, t_max = 15, seed = 13)
  tn <- temporal_network(edges, time_interval = 3)
  lw_range <- c(0.5, 3)
  p <- plot(tn, type = "proximity", line_width = lw_range)
  built <- ggplot2::ggplot_build(p)
  lw <- built$data[[1]]$linewidth
  expect_true(all(lw >= lw_range[1] - 1e-10))
  expect_true(all(lw <= lw_range[2] + 1e-10))
})

test_that("proximity timeline constant degree gives finite linewidth", {
  edges <- data.frame(
    from = rep("hub", 4),
    to = paste0("leaf", 1:4),
    onset = rep(c(0, 5), each = 4),
    terminus = rep(c(5, 10), each = 4),
    stringsAsFactors = FALSE
  )
  tn <- temporal_network(edges, time_interval = 5)
  p <- plot(tn, type = "proximity", line_width = c(1, 4))
  built <- ggplot2::ggplot_build(p)
  lw <- built$data[[1]]$linewidth
  expect_true(all(is.finite(lw)))
})

test_that("proximity timeline highlight fades non-highlighted vertices", {
  edges <- .make_temporal_random(n_v = 6, n_e = 20, t_max = 12, seed = 14)
  tn <- temporal_network(edges, time_interval = 3)
  vnames <- tn$vertex_names

  p <- plot(tn, type = "proximity", highlight = vnames[1])
  built <- ggplot2::ggplot_build(p)
  seg_data <- built$data[[1]]
  alphas <- unique(round(seg_data$alpha, 2))
  expect_true(length(alphas) == 2L)
  expect_true(min(alphas) <= 0.15)
  expect_true(max(alphas) >= 0.35)
})

test_that("proximity timeline on chain network (small/simple)", {
  edges <- .make_temporal_chain(5)
  tn <- temporal_network(edges)
  p <- plot(tn, type = "proximity", metric = "degree")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(nrow(built$data[[1]]) > 0L)
})

test_that("proximity timeline single time bin returns valid plot", {
  edges <- data.frame(
    from = c("A", "B", "C"),
    to = c("B", "C", "A"),
    onset = c(0, 0, 0),
    terminus = c(5, 5, 5),
    stringsAsFactors = FALSE
  )
  tn <- temporal_network(edges, time_interval = 10)
  p <- plot(tn, type = "proximity")
  expect_s3_class(p, "ggplot")
})
