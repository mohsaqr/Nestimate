# House-rule fixes in the simplicial, hypergraph and link-prediction files:
# classed conditions, restored random state, corrected labels, and loops
# replaced by vectorised code that returns the same values.

hr_g_weights <- function() {
  matrix(c(0, .6, .5, 0,
           .6, 0, .4, .3,
           .5, .4, 0, .2,
           0, .3, .2, 0), 4, 4,
         dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D")))
}

hr_g_events <- function() {
  data.frame(
    person  = c("a", "b", "c", "a", "c", "d", "d", "e", "f", "e"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3", "m3"),
    hours   = c(1, 2, 0.5, 1.5, 1, 2, 1, 1, 3, 0.25),
    stringsAsFactors = FALSE
  )
}

hr_g_net <- function() {
  build_network(group_regulation_long, method = "relative",
                action = "Action", actor = "Actor", time = "Time")
}

# ---- simplicial.R -----------------------------------------------------------

test_that("build_simplicial() warns when `...` is not read", {
  w <- hr_g_weights()
  expect_warning(build_simplicial(w, foo = 1),
                 class = "nestimate_ignored_argument")
  expect_warning(build_simplicial(1 - w, type = "vr", foo = 1),
                 class = "nestimate_ignored_argument")
  expect_no_warning(build_simplicial(w))

  hon <- build_hon(build_network(trajectories, method = "relative"),
                   max_order = 2)
  expect_warning(sc <- build_simplicial(hon, type = "pathway", method = "hypa"),
                 class = "nestimate_ignored_argument")
  expect_identical(sc, build_simplicial(hon, type = "pathway"))
})

test_that("build_simplicial() refuses a non-numeric matrix with a class", {
  expect_error(build_simplicial(matrix("a", 2, 2)),
               class = "nestimate_bad_input")
  expect_error(persistent_homology(matrix("a", 2, 2)),
               class = "nestimate_bad_input")
})

test_that("simplex bookkeeping matches a count by hand", {
  # triangle A-B-C with a pendant edge B-D
  sc <- build_simplicial(hr_g_weights())
  expect_identical(unname(sc$f_vector), c(4L, 5L, 2L))
  expect_identical(betti_numbers(sc), c(b0 = 1L, b1 = 0L, b2 = 0L))

  deg <- simplicial_degree(sc)
  by_node <- deg[match(c("A", "B", "C", "D"), deg$node), ]
  expect_identical(by_node$d0, c(1L, 1L, 1L, 1L))
  expect_identical(by_node$d1, c(2L, 3L, 3L, 2L))
  expect_identical(by_node$d2, c(1L, 2L, 2L, 1L))
  expect_equal(by_node$total, c(3, 5, 5, 3))

  norm <- simplicial_degree(sc, normalized = TRUE)
  norm <- norm[match(c("A", "B", "C", "D"), norm$node), ]
  expect_equal(norm$d1, c(2, 3, 3, 2) / choose(3, 1))
  expect_equal(norm$d2, c(1, 2, 2, 1) / choose(3, 2))

  # two maximal triangles sharing the edge B-C: one component at q <= 1
  expect_identical(unname(q_analysis(sc)$q_vector), c(2L, 1L, 1L))
})

test_that("a hollow square has one cycle and no isolated node is lost", {
  ring <- matrix(0, 5, 5, dimnames = list(letters[1:5], letters[1:5]))
  ring[cbind(1:4, c(2:4, 1))] <- 1
  sc <- build_simplicial(ring + t(ring))
  expect_identical(betti_numbers(sc), c(b0 = 2L, b1 = 1L))
  expect_identical(sc$n_simplices, 9L)
  expect_true(any(vapply(sc$simplices, identical, logical(1), 5L)))
})

test_that("simplicial_degree() checks `normalized`", {
  sc <- build_simplicial(hr_g_weights())
  expect_error(simplicial_degree(sc, normalized = "yes"),
               class = "nestimate_bad_input")
})

test_that("persistent_homology() says why a clique complex cannot be read", {
  sc <- build_simplicial(hr_g_weights())
  expect_error(persistent_homology(sc), "carries no filtration",
               class = "nestimate_bad_input")
})

test_that("print.persistent_homology() lists the longest-lived features", {
  out <- capture.output(print(persistent_homology(hr_g_weights())))
  lived <- grep("life:", out, value = TRUE)
  expect_length(lived, 3L)
  expect_match(lived[1L], "b0: 0.6000 .* 0.0000 \\(life: 0.6000\\)")
})

test_that("plotting a complex without edges raises a classed error", {
  skip_if_not_installed("ggplot2")
  empty <- build_simplicial(matrix(0, 3, 3,
                                   dimnames = list(letters[1:3], letters[1:3])))
  expect_error(plot(empty), class = "nestimate_bad_input")
  expect_error(plot(empty, combined = FALSE), class = "nestimate_bad_input")
})

test_that("the Betti curve is labelled by the scale of its filtration", {
  skip_if_not_installed("ggplot2")
  w <- hr_g_weights()
  rips <- plot(persistent_homology(1 - w, type = "vr"), combined = FALSE)
  expect_identical(rips$betti_curve$labels$x, "Distance Threshold")
  expect_match(rips$betti_curve$labels$subtitle, "distance thresholds")
  clique <- plot(persistent_homology(w), combined = FALSE)
  expect_identical(clique$betti_curve$labels$x, "Weight Threshold")
})

# ---- simplicial_features.R --------------------------------------------------

test_that("simplicial_features() validates max_dim and threshold up front", {
  w <- hr_g_weights()
  expect_error(simplicial_features(w, max_dim = 0), "`max_dim` must be",
               class = "nestimate_bad_input")
  expect_error(simplicial_features(w, threshold = -1), "`threshold` must be",
               class = "nestimate_bad_input")
  expect_equal(simplicial_features(w, max_dim = 2)$d2, 2)
})

# ---- ph_distances.R ---------------------------------------------------------

test_that("diagram verbs name the argument that broke the contract", {
  ph <- persistent_homology(hr_g_weights())
  expect_error(persistence_landscape(ph, k_max = 0), "`k_max` must be",
               class = "nestimate_bad_input")
  expect_error(bottleneck_distance(ph, ph, dimension = 0.5),
               "`dimension` must be", class = "nestimate_bad_input")
  expect_equal(unname(bottleneck_distance(ph, ph)), 0)
})

# ---- hypergraph.R -----------------------------------------------------------

test_that("build_hypergraph() leaves the caller's random stream alone", {
  w <- hr_g_weights()
  set.seed(1)
  expected <- runif(1)
  set.seed(1)
  sampled <- build_hypergraph(w, p = 0.5, seed = 5)
  expect_identical(runif(1), expected)
  expect_identical(sampled, build_hypergraph(w, p = 0.5, seed = 5))
})

test_that("build_hypergraph() raises classed errors on bad arguments", {
  w <- hr_g_weights()
  expect_error(build_hypergraph(w, p = 2), "`p` must be",
               class = "nestimate_bad_input")
  expect_error(suppressWarnings(build_hypergraph(w, p = 0.5, seed = "a")),
               "`seed` must be", class = "nestimate_bad_input")
})

test_that("the incidence matrix has one column per hyperedge", {
  hg <- build_hypergraph(hr_g_weights())
  expect_identical(dim(hg$incidence), c(4L, 7L))
  members <- lapply(seq_len(hg$n_hyperedges),
                    function(j) unname(which(hg$incidence[, j] == 1L)))
  expect_identical(members, lapply(hg$hyperedges, as.integer))
  expect_identical(build_hypergraph(build_simplicial(hr_g_weights()))$hyperedges,
                   hg$hyperedges)

  out <- capture.output(print(hg))
  expect_identical(out[3:4], c("  size_2   : 5", "  size_3   : 2"))
})

# ---- bipartite_groups.R -----------------------------------------------------

test_that("bipartite_groups() adds the weights of each player-group cell", {
  events <- rbind(hr_g_events(),
                  data.frame(person = "a", meeting = "m1", hours = 0.25))
  hg <- bipartite_groups(events, player = "person", group = "meeting",
                         weight = "hours")
  expect_equal(hg$incidence["a", "m1"], 1.25)
  expect_equal(hg$incidence["f", "m3"], 3)
  expect_equal(hg$incidence["f", "m1"], 0)
  expect_equal(sum(hg$incidence), sum(events$hours))
})

test_that("bipartite_groups() names the argument in its classed errors", {
  events <- hr_g_events()
  expect_error(bipartite_groups(events, player = "who", group = "meeting"),
               "`player` must name", class = "nestimate_bad_input")
  expect_error(bipartite_groups(events, player = "person", group = "meeting",
                                weight = "nope"),
               "`weight` must be", class = "nestimate_bad_input")
})

# ---- clique_expansion.R, hypergraph_measures.R, hypergraph_centrality.R -----

test_that("hypergraph verbs say what `hg` must be", {
  expect_error(clique_expansion(1), "`hg` must be a net_hypergraph",
               class = "nestimate_bad_input")
  expect_error(hypergraph_measures(1), "`hg` must be a net_hypergraph",
               class = "nestimate_bad_input")
  expect_error(hypergraph_centrality(1), "`hg` must be a net_hypergraph",
               class = "nestimate_bad_input")
  hg <- bipartite_groups(hr_g_events(), player = "person", group = "meeting")
  expect_error(hypergraph_centrality(hg, tol = -1), "`tol` must be",
               class = "nestimate_bad_input")
})

test_that("print.hypergraph_measures() lists the intersection profile", {
  hg <- bipartite_groups(hr_g_events(), player = "person", group = "meeting")
  out <- capture.output(print(hypergraph_measures(hg)))
  profile <- out[seq(grep("Intersection profile", out) + 1L, length(out))]
  expect_identical(profile, c("  overlap_0   : 1", "  overlap_1   : 1",
                              "  overlap_2   : 1"))
})

# ---- hypergraph_laplacian.R -------------------------------------------------

test_that("hypergraph_cluster() leaves the caller's random stream alone", {
  hg <- bipartite_groups(hr_g_events(), player = "person", group = "meeting")
  set.seed(1)
  expected <- runif(1)
  set.seed(1)
  cl <- hypergraph_cluster(hg, k = 2, seed = 5)
  expect_identical(runif(1), expected)
  expect_error(suppressWarnings(hypergraph_cluster(hg, k = 2, seed = "a")),
               "`seed` must be", class = "nestimate_bad_input")
})

test_that("the embedding plot draws more than nine clusters", {
  skip_if_not_installed("ggplot2")
  i <- seq_len(23L)
  chain <- data.frame(person  = sprintf("n%02d", c(i, i + 1L)),
                      meeting = sprintf("m%02d", c(i, i)))
  hg <- bipartite_groups(chain, player = "person", group = "meeting")
  cl <- hypergraph_cluster(hg, k = 11, seed = 1)
  built <- ggplot2::ggplot_build(plot(cl, what = "embedding"))
  glyphs <- unique(built$data[[1L]][, c("colour", "shape")])
  expect_identical(nrow(glyphs), 11L)
  expect_error(plot(cl, what = "spectrum", n_values = "many"),
               class = "nestimate_bad_input")
})

# ---- link_prediction.R ------------------------------------------------------

test_that("the Katz damping adjustment is a classed warning", {
  w <- hr_g_weights()
  expect_warning(
    adjusted <- predict_links(w, methods = "katz", katz_damping = 5),
    "auto-adjusted", class = "nestimate_ignored_argument"
  )
  expect_identical(adjusted$scores,
                   predict_links(w, methods = "katz")$scores)
})

test_that("the truncated Katz series is the sum of its six terms", {
  a <- (hr_g_weights() > 0) * 1
  damping <- 0.1
  by_hand <- damping * a +
    damping^2 * (a %*% a) +
    damping^3 * (a %*% a %*% a) +
    damping^4 * (a %*% a %*% a %*% a) +
    damping^5 * (a %*% a %*% a %*% a %*% a) +
    damping^6 * (a %*% a %*% a %*% a %*% a %*% a)
  expect_equal(.lp_katz_series(a, damping, terms = 6L), by_hand)
})

test_that("the Katz index refuses missing weights with a class", {
  w <- hr_g_weights()
  w["A", "B"] <- NA
  expect_error(predict_links(w, methods = "katz"),
               class = "nestimate_bad_input")
})

test_that("evaluate_links() warns when a matrix's column names are not read", {
  pred <- predict_links(hr_g_weights())
  truth <- matrix(0, 4, 4, dimnames = list(c("A", "B", "C", "D"),
                                           c("D", "C", "B", "A")))
  truth["A", "D"] <- 1
  expect_warning(evaluate_links(pred, truth, k = 1),
                 class = "nestimate_ignored_argument")
  colnames(truth) <- rownames(truth)
  expect_no_warning(evaluate_links(pred, truth, k = 1))
  expect_error(evaluate_links(pred, list(1)), "`true_edges` must be",
               class = "nestimate_bad_input")
})

test_that("print.net_link_prediction() numbers the top links", {
  pred <- predict_links(hr_g_net())
  out <- capture.output(print(pred))
  ranked <- grep("^    [0-9]+\\. ", out, value = TRUE)
  expect_length(ranked, min(10L, nrow(pred$consensus)))
  top <- pred$consensus[1L, ]
  expect_identical(
    ranked[1L],
    sprintf("    1. %s -> %s  (avg rank: %.1f, agreed: %d/%d)",
            top$from, top$to, top$avg_rank, top$n_methods,
            length(pred$methods))
  )

  single <- predict_links(hr_g_weights(), methods = "jaccard")
  out <- capture.output(print(single))
  expect_identical(
    out[5L],
    sprintf("    1. %s -- %s  (score: %.4f)", single$predictions$from[1L],
            single$predictions$to[1L], single$predictions$score[1L])
  )
})
