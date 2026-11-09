# Named verbs for the tables that used to be reached with as.data.frame().

.tv_hg <- function() {
  events <- data.frame(
    person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
               "d", "e", "f", "c", "d"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
                "m4", "m4", "m4", "m5", "m5")
  )
  bipartite_groups(events, player = "person", group = "meeting")
}

test_that("node_clusters() returns one row per node", {
  cl <- hypergraph_cluster(.tv_hg(), k = 2, seed = 1)
  out <- node_clusters(cl)
  expect_identical(names(out), c("node", "cluster", "pi", "dim1", "dim2"))
  expect_identical(out$node, c("a", "b", "c", "d", "e", "f"))
  expect_equal(sum(out$pi), 1)
  expect_identical(as.data.frame(cl), out)
  expect_error(node_clusters(1), class = "nestimate_bad_input")
})

test_that("label_predictions() and label_scores() agree with the scores", {
  tr <- hypergraph_transduction(.tv_hg(), labels = c(a = "x", d = "y"),
                                xi = 0.9)
  pred <- label_predictions(tr)
  sc <- label_scores(tr)
  expect_identical(pred$predicted, c("x", "x", "x", "y", "y", "y"))
  expect_identical(nrow(sc), nrow(pred) * 2L)
  # hand check: the winning score is the row maximum of the score matrix
  expect_equal(pred$score, unname(apply(tr$scores, 1, max)))
  expect_equal(sc$score, as.vector(tr$scores))
  expect_identical(as.data.frame(tr), pred)
  expect_identical(as.data.frame(tr, what = "scores"), sc)
  expect_error(label_predictions(1), class = "nestimate_bad_input")
  expect_error(label_scores(1), class = "nestimate_bad_input")
})

test_that("the default xi over-smooths the small example, as documented", {
  tr <- hypergraph_transduction(.tv_hg(), labels = c(a = "x", d = "y"))
  expect_identical(label_predictions(tr)$predicted,
                   c("x", "y", "y", "y", "y", "y"))
})

test_that("state_distribution() of a state_freq result is its table", {
  skip_if_not_installed("ggplot2")
  net <- build_network(group_regulation_long, method = "relative",
                       format = "long", actor = "Actor", action = "Action",
                       order = "Time", group = "Course")
  res <- plot_state_frequencies(net)
  out <- state_distribution(res)
  expect_s3_class(out, "data.frame")
  expect_identical(out, res$table)
  expect_identical(as.data.frame(res), out)
})

test_that("summary() of cluster diagnostics is the per-cluster table", {
  cl <- build_clusters(trajectories, k = 2, seed = 1)
  dg <- cluster_diagnostics(cl)
  out <- summary(dg)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 2L)
  expect_identical(out, dg$per_cluster)
  expect_identical(as.data.frame(dg), out)
  expect_error(summary(dg, typo_arg = 1), class = "nestimate_bad_input")
})
