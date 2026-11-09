# A `seed =` argument makes the call reproducible and leaves the caller's own
# random stream where it was (.local_seed(), as withr::local_seed()).

expect_rng_untouched <- function(expr) {
  set.seed(99)
  untouched <- stats::runif(3)
  set.seed(99)
  force(expr)
  expect_identical(stats::runif(3), untouched)
}

test_that(".local_seed() seeds the calling function and restores on exit", {
  draw <- function(seed = NULL) {
    if (!is.null(seed)) .local_seed(seed)
    stats::runif(2)
  }
  set.seed(5)
  seeded <- stats::runif(2)
  expect_identical(draw(seed = 5), seeded)
  expect_rng_untouched(draw(seed = 1))
  # no stray variable left in the caller's frame
  expect_false(exists("restore", envir = environment(), inherits = FALSE))
})

test_that(".local_seed() leaves a session with no random state without one", {
  set.seed(1)
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  rm(".Random.seed", envir = globalenv())
  draw <- function() {
    .local_seed(3)
    stats::runif(1)
  }
  draw()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("seed = in the resampling verbs does not move the caller's stream", {
  testthat::skip_on_cran()
  net_a <- build_network(utils::head(trajectories, 70), method = "relative")
  net_b <- build_network(utils::tail(trajectories, 68), method = "relative")
  freq_a <- build_network(utils::head(trajectories, 70), method = "frequency")
  freq_b <- build_network(utils::tail(trajectories, 68), method = "frequency")

  expect_rng_untouched(bootstrap_network(net_a, iter = 5, seed = 1))
  expect_rng_untouched(permutation(net_a, net_b, iter = 5, seed = 1))
  expect_rng_untouched(bayes_compare(freq_a, freq_b, draws = 100, seed = 1))
  expect_rng_untouched(entropy_bayes(freq_a, draws = 100, seed = 1))
  expect_rng_untouched(centrality_stability(net_a, iter = 5, seed = 1))
  expect_rng_untouched(network_reliability(net_a, iter = 5, seed = 1))
  expect_rng_untouched(vertex_bootstrap(net_a, iter = 5, seed = 1))
  expect_rng_untouched(build_clusters(net_a, k = 2, seed = 1))
  expect_rng_untouched(build_mmm(net_a, k = 2, n_starts = 1, seed = 1))
  expect_rng_untouched(markov_order_test(net_a, max_order = 1, n_perm = 5,
                                         seed = 1))
  expect_rng_untouched(stability_landscape(trajectories, iter = 5, seed = 1))
})

test_that("seed = still reproduces the same result", {
  testthat::skip_on_cran()
  net <- build_network(utils::head(trajectories, 70), method = "frequency")
  first <- entropy_bayes(net, draws = 100, seed = 7)
  second <- entropy_bayes(net, draws = 100, seed = 7)
  expect_identical(first, second)
})
