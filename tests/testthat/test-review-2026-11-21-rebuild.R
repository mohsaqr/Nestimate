testthat::skip_on_cran()

# Review 2026-11-21: the networks rebuilt for the clusters of a network are
# built afresh with its method; `inherit = TRUE` adds its threshold, scaling
# and $params; start/end on one-hot data, build_mmm() on an mcml or group_tna
# and predictability() on a tna model, bare network or mcml are refused with
# a class; build_mmm()'s empty-cluster warning is classed.

rb_rows <- function(assignments, cluster) {
  trajectories[assignments == cluster, , drop = FALSE]
}

rb_expect_direct <- function(networks, ..., assignments = attr(networks, "clustering")$assignments) {
  invisible(Map(function(net, cluster) {
    direct <- build_network(rb_rows(assignments, cluster), ...)
    expect_equal(net$weights, direct$weights, tolerance = 0, info = sprintf("cluster %d", cluster))
    expect_identical(net$threshold, direct$threshold)
    expect_identical(net$scaling, direct$scaling)
  }, networks, seq_along(networks)))
}

test_that("cluster networks of a network are built afresh unless inherit = TRUE", {
  net <- build_network(trajectories, method = "attention", lambda = 3, threshold = 0.05, scaling = "minmax")
  clustering <- build_clusters(net, k = 2, seed = 1)
  expect_identical(clustering$network_settings$threshold, 0.05)
  expect_identical(clustering$network_settings$params$lambda, 3)
  # default: the network's method only
  rb_expect_direct(build_network(clustering), method = "attention")
  rb_expect_direct(cluster_network(net, k = 2, seed = 1), method = "attention")
  # inherit = TRUE: its threshold, scaling and params
  rb_expect_direct(build_network(clustering, inherit = TRUE), method = "attention",
                   lambda = 3, threshold = 0.05, scaling = "minmax")
  rb_expect_direct(cluster_network(net, k = 2, seed = 1, inherit = TRUE), method = "attention",
                   lambda = 3, threshold = 0.05, scaling = "minmax")
  # a setting given in the call wins; scaling = NULL keeps the networks unscaled
  rb_expect_direct(build_network(clustering, inherit = TRUE, scaling = NULL, threshold = 0.1),
                   method = "attention", lambda = 3, threshold = 0.1)
  rb_expect_direct(build_network(clustering, inherit = TRUE, params = list(lambda = 1)),
                   method = "attention", lambda = 1, threshold = 0.05, scaling = "minmax")
  # another method inherits nothing
  rb_expect_direct(build_network(clustering, method = "frequency", inherit = TRUE), method = "frequency")
})

test_that("start and end states are inherited only with inherit = TRUE", {
  started <- build_network(trajectories, method = "relative", start = TRUE, end = TRUE)
  clustering <- build_clusters(started, k = 2, seed = 1)
  fresh <- build_network(clustering)
  expect_false("Start" %in% rownames(fresh[[1]]$weights))
  rb_expect_direct(fresh, method = "relative")
  kept <- build_network(clustering, inherit = TRUE)
  expect_true(all(c("Start", "End") %in% rownames(kept[[1]]$weights)))
  rb_expect_direct(kept, method = "relative", start = TRUE, end = TRUE)
})

test_that("a network built with an actor rebuilds its clusters from the positions", {
  wide <- data.frame(learner = paste0("s", seq_len(nrow(trajectories))),
                     as.data.frame(trajectories, stringsAsFactors = FALSE),
                     stringsAsFactors = FALSE)
  net <- build_network(wide, method = "relative", actor = "learner")
  networks <- build_network(build_clusters(net, k = 2), inherit = TRUE)
  expect_s3_class(networks, "netobject_group")
  rb_expect_direct(networks, method = "relative")
})

test_that("the networks of a mixture are built afresh unless inherit = TRUE", {
  net <- build_network(trajectories, method = "relative", threshold = 0.2, scaling = "minmax")
  fit <- build_mmm(net, k = 2, seed = 1, n_starts = 3)
  expect_identical(fit$network_settings$scaling, "minmax")
  given <- build_network(fit, threshold = 0.2, scaling = "minmax")
  invisible(Map(function(fresh, one_step, kept, one_step_kept, model, reshaped) {
    expect_identical(fresh$weights, model$weights)
    expect_identical(one_step$weights, model$weights)
    expect_identical(fresh$threshold, 0)
    expect_identical(kept$weights, reshaped$weights)
    expect_identical(one_step_kept$weights, reshaped$weights)
  }, build_network(fit), cluster_network(net, k = 2, cluster_by = "mmm", seed = 1, n_starts = 3),
  build_network(fit, inherit = TRUE),
  cluster_network(net, k = 2, cluster_by = "mmm", seed = 1, n_starts = 3, inherit = TRUE),
  fit$models, given))
  counts <- build_network(trajectories, method = "frequency", threshold = 5)
  fit_counts <- build_mmm(counts, k = 2, seed = 1, n_starts = 3)
  rb_expect_direct(build_network(fit_counts), method = "frequency", assignments = fit_counts$assignments)
  rb_expect_direct(build_network(fit_counts, inherit = TRUE), method = "frequency", threshold = 5,
                   assignments = fit_counts$assignments)
  expect_identical(build_network(fit_counts, method = "relative", inherit = TRUE)[[1]]$threshold, 0)
})

test_that("L8-3: start/end on one-hot data are refused with nestimate_bad_input", {
  onehot <- data.frame(A = c(1, 0, 1, 0, 1, 1), B = c(0, 1, 1, 0, 0, 1),
                       actor = rep(c("x", "y"), each = 3))
  invisible(lapply(c("relative", "frequency", "co_occurrence"), function(method) {
    expect_error(build_network(onehot, method = method, codes = c("A", "B"), actor = "actor",
                               window_size = 2, start = TRUE), class = "nestimate_bad_input")
    expect_error(build_network(onehot, method = method, codes = c("A", "B"), actor = "actor",
                               window_size = 2, end = "Stop"), class = "nestimate_bad_input")
    net <- build_network(onehot, method = method, codes = c("A", "B"), actor = "actor", window_size = 2)
    expect_identical(rownames(net$weights), c("A", "B"))
  }))
})

test_that("L7-12: build_mmm() refuses an mcml and a group_tna with nestimate_bad_input", {
  wide <- data.frame(T1 = c("a", "b", "c", "a", "b", "c"), T2 = c("b", "c", "a", "c", "a", "b"),
                     T3 = c("c", "a", "b", "b", "c", "a"))
  mcml <- build_mcml(wide, clusters = list(X = c("a", "b"), Y = "c"))
  expect_error(build_mmm(mcml, k = 2), class = "nestimate_bad_input")
  expect_error(build_mmm(structure(list(g1 = list(), g2 = list()), class = "group_tna"), k = 2),
               class = "nestimate_bad_input")
})

test_that("L7-13: predictability() refuses a tna model, a bare network and an mcml with a class", {
  bare <- structure(list(weights = matrix(c(0, .6, .4, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
                         nodes = data.frame(id = 1:2, label = c("A", "B"), name = c("A", "B")),
                         directed = TRUE), class = "cograph_network")
  expect_error(predictability(bare), class = "nestimate_bad_input")
  wide <- data.frame(T1 = c("a", "b", "c", "a", "b", "c"), T2 = c("b", "c", "a", "c", "a", "b"),
                     T3 = c("c", "a", "b", "b", "c", "a"))
  expect_error(predictability(build_mcml(wide, clusters = list(X = c("a", "b"), Y = "c"))),
               class = "nestimate_bad_input")
  expect_error(predictability(list(a = 1)), class = "nestimate_bad_input")
  skip_if_not_installed("tna")
  model <- tna::tna(data.frame(T1 = c("A", "B", "A", "C"), T2 = c("B", "C", "C", "A"),
                               T3 = c("C", "A", "B", "B")))
  expect_error(predictability(model), class = "nestimate_bad_input")
})

test_that("L6: build_mmm()'s empty-cluster warning has class nestimate_empty_cluster", {
  constant <- data.frame(V1 = rep("A", 20), V2 = rep("B", 20), V3 = rep("A", 20), V4 = rep("B", 20))
  expect_warning(build_mmm(constant, k = 2, n_starts = 5, max_iter = 50, seed = 1),
                 class = "nestimate_empty_cluster")
})

test_that("build_mmm() and build_clusters() refuse an edge-list mcml layer with nestimate_bad_input", {
  edges <- data.frame(from = c("a1", "a2", "b1", "b2", "a1", "b1"), to = c("a2", "b1", "b2", "a1", "b1", "a2"),
                      w = c(5, 3, 4, 2, 1, 6))
  fit <- build_mcml(edges, clusters = list(A = c("a1", "a2"), B = c("b1", "b2")),
                    from = "from", to = "to", weight = "w")
  layer <- as_tna(fit)$macro
  expect_identical(attr(layer$data, "source"), "edgelist")
  expect_error(build_mmm(layer, k = 2), class = "nestimate_bad_input")
  expect_error(build_clusters(layer, k = 2), class = "nestimate_bad_input")
  expect_error(cluster_network(layer, k = 2), class = "nestimate_bad_input")
})
