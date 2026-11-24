testthat::skip_on_cran()

# A palette set on a network is carried to every object derived from it.

.carry_net <- function() {
  set_state_colors(build_network(trajectories, method = "relative"),
                   c(Active = "#D55E00", Average = "#CC79A7"))
}

.all_set <- function(obj, expected) {
  tb <- state_colors(obj)
  identical(stats::setNames(tb$color, tb$state)[names(expected)], expected)
}

test_that("clusterings, mixtures and their networks carry the palette", {
  net <- .carry_net()
  expected <- stats::setNames(state_colors(net)$color, state_colors(net)$state)
  clusters <- build_clusters(net, k = 3, seed = 1)
  expect_identical(.stored_state_colors(clusters)[names(expected)], expected)
  expect_true(.all_set(build_network(clusters), expected))
  expect_true(.all_set(cluster_network(net, k = 3), expected))
  mixture <- build_mmm(net, k = 2, seed = 1, n_starts = 2)
  expect_identical(.stored_state_colors(mixture)[names(expected)], expected)
  expect_true(.all_set(build_network(mixture), expected))
  # each member network also carries the cograph contract
  member <- build_network(clusters)[[1]]
  expect_equal(member$meta$splot$defaults$node_fill,
               unname(expected[member$nodes$label]))
})

test_that("edge-betweenness and entropy networks carry the palette", {
  net <- .carry_net()
  expected <- stats::setNames(state_colors(net)$color, state_colors(net)$state)
  eb <- net_edge_betweenness(net)
  expect_true(.all_set(eb, expected))
  expect_equal(eb$meta$splot$defaults$node_fill,
               unname(expected[eb$nodes$label]))
  ent <- entropy_network(net)
  expect_equal(ent$meta$splot$defaults$node_fill,
               unname(expected[ent$nodes$label]))
  # the entropy styling is kept beside the colours
  expect_true(isTRUE(ent$meta$splot$defaults$tna_styling))
})

test_that("a network without a palette passes on its Okabe-Ito defaults", {
  net <- build_network(trajectories, method = "relative")
  defaults <- stats::setNames(state_colors(net)$color, state_colors(net)$state)
  clusters <- build_clusters(net, k = 3, seed = 1)
  # a clustering resolves the same defaults itself and carries nothing
  expect_null(.stored_state_colors(clusters))
  member <- build_network(clusters)[[1]]
  expect_equal(member$meta$splot$defaults$node_fill,
               unname(defaults[member$nodes$label]))
  eb <- net_edge_betweenness(net)
  expect_equal(eb$meta$splot$defaults$node_fill, unname(defaults[eb$nodes$label]))
  # nothing is reported as set
  expect_true(all(state_colors(build_network(clusters))$source == "default"))
})

test_that("bootstrap, certainty, permutation and Bayes results declare the colours", {
  net <- .carry_net()
  fill <- net$meta$splot$defaults$node_fill
  boot <- bootstrap_network(net, iter = 10, seed = 1)
  expect_equal(boot$meta$splot$defaults$node_fill, fill)
  expect_equal(boot$model$meta$splot$defaults$node_fill, fill)
  expect_equal(certainty(net)$meta$splot$defaults$node_fill, fill)
  grp <- set_state_colors(
    build_network(cbind(trajectories,
                        g = rep(c("a", "b"), length.out = nrow(trajectories))),
                  method = "relative", group = "g"),
    c(Active = "#D55E00", Average = "#CC79A7"))
  pair_fill <- grp[[1]]$meta$splot$defaults$node_fill
  expect_equal(permutation(grp, iter = 10, seed = 1)[[1]]$meta$splot$defaults$node_fill,
               pair_fill)
  bayes <- bayes_compare(grp, seed = 1, draws = 100)[[1]]
  expect_equal(bayes$meta$splot$defaults$node_fill, pair_fill)
  # the Bayesian result keeps its renderer beside the colours
  expect_identical(bayes$meta$splot$renderer, "permutation")
})

test_that("cograph draws a bootstrap in the colours of its network", {
  skip_if_not_installed("cograph")
  net <- .carry_net()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  drawn <- cograph::splot(bootstrap_network(net, iter = 10, seed = 1))
  expect_equal(unname(drawn$nodes$node_fill),
               net$meta$splot$defaults$node_fill)
})
