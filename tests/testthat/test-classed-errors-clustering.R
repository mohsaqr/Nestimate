# Input errors from clustering / estimator / mcml_pc verbs carry nestimate_bad_input.
testthat::skip_on_cran()

seqs <- trajectories[1:30, ]

test_that("build_clusters() input errors are classed", {
  expect_error(build_clusters(seqs, k = 1.5), class = "nestimate_bad_input")
  expect_error(build_clusters(seqs, k = 2, lambda = -1),
               class = "nestimate_bad_input")
  expect_error(build_clusters(seqs, k = 2, weighted = "a"),
               class = "nestimate_bad_input")
  expect_error(build_clusters(list(1, 2), k = 2),
               class = "nestimate_bad_input")
})

test_that("net_clustering methods and cluster_network() errors are classed", {
  cl <- build_clusters(seqs, k = 2)
  expect_error(plot(cl, type = "predictors"), class = "nestimate_bad_input")
  expect_error(plot(cl, combined = NA_character_),
               class = "nestimate_bad_input")
  expect_error(print(cl, digits = -1), class = "nestimate_bad_input")
  expect_error(cluster_network(seqs, k = 2, bogus = 3),
               class = "nestimate_bad_input")
})

test_that("cluster_choice() and its plot method errors are classed", {
  expect_error(cluster_choice(seqs, k = 2.5), class = "nestimate_bad_input")
  cc <- cluster_choice(seqs, k = 2:3)
  expect_error(plot(cc, abbrev = NA), class = "nestimate_bad_input")
  expect_error(plot(cc, type = "heatmap"), class = "nestimate_bad_input")
})

test_that("cluster_diagnostics() errors are classed", {
  expect_error(cluster_diagnostics(1), class = "nestimate_bad_input")
  cd <- cluster_diagnostics(build_clusters(seqs, k = 2))
  expect_error(print(cd, foo = 1), class = "nestimate_bad_input")
})

test_that("build_network() estimator input errors are classed", {
  expect_error(build_network(trajectories, method = "relative",
                             params = list(concat = 0)),
               class = "nestimate_bad_input")
  expect_error(build_network(trajectories, method = "attention",
                             params = list(lambda = -1)),
               class = "nestimate_bad_input")
  expect_error(build_network(data.frame(a = 1:5), method = "cor"),
               class = "nestimate_bad_input")
  expect_error(build_network(srl_strategies, method = "glasso",
                             params = list(gamma = -1)),
               class = "nestimate_bad_input")
  expect_error(build_network(data.frame(a = c(0, 1, 2, 1), b = c(1, 0, 1, 0)),
                             method = "ising"),
               class = "nestimate_bad_input")
  expect_error(build_network(srl_strategies[, 1:3], method = "mgm",
                             params = list(type = c("g", "x", "g"))),
               class = "nestimate_bad_input")
})

test_that("ngram / gap / reverse estimator errors are classed", {
  expect_error(build_network(trajectories, method = "ngram",
                             params = list(n_gram = 1)),
               class = "nestimate_bad_input")
  expect_error(build_network(trajectories, method = "gap",
                             params = list(max_gap = -1)),
               class = "nestimate_bad_input")
  expect_error(build_network(trajectories, method = "reverse",
                             params = list(weighted = NA)),
               class = "nestimate_bad_input")
})

test_that("build_mcml_pc() and its helpers' errors are classed", {
  expect_error(build_mcml_pc("x", clusters = list(A = "a")),
               class = "nestimate_bad_input")
  expect_error(build_mcml_pc(srl_strategies, clusters = 1:3),
               class = "nestimate_bad_input")
  expect_error(loading_stability(1), class = "nestimate_bad_input")
  expect_error(as_networks(1), class = "nestimate_bad_input")
})
