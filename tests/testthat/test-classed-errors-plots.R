# Bad input to the plotting, MMM and simplicial verbs raises class nestimate_bad_input.
testthat::skip_on_cran()

seqs <- data.frame(V1 = c("A", "B", "A"), V2 = c("B", "A", "B"),
                   V3 = c("A", "B", "B"), stringsAsFactors = FALSE)
sq <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3,
             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))

test_that("mosaic_plot() input checks are classed", {
  expect_error(mosaic_plot(1:3), class = "nestimate_bad_input")
  expect_error(mosaic_plot(data.frame(a = "x", w = 1), x_var = "a",
                           y_var = "b", weight = "w"),
               class = "nestimate_bad_input")
  expect_error(mosaic_plot(matrix(c(1, 2, 3, 4), 2), n_perm = 1L),
               class = "nestimate_bad_input")
})

test_that("plot_state_frequencies() and state_distribution() reject other classes", {
  expect_error(plot_state_frequencies(1), class = "nestimate_bad_input")
  expect_error(state_distribution(1), class = "nestimate_bad_input")
})

test_that("mosaic_analysis() input checks are classed", {
  expect_error(mosaic_analysis(1, "a", "b"), class = "nestimate_bad_input")
  expect_error(mosaic_analysis(data.frame(a = "x"), "a", "zz"),
               class = "nestimate_bad_input")
})

test_that("build_mmm(), compare_mmm() and cluster_mmm() input checks are classed", {
  expect_error(build_mmm(seqs, k = 1), class = "nestimate_bad_input")
  expect_error(build_mmm(seqs, k = "2"), class = "nestimate_bad_input")
  expect_error(compare_mmm(seqs, return_fits = NA),
               class = "nestimate_bad_input")
  expect_error(cluster_mmm(seqs, k = 2, cluster_by = "pam"),
               class = "nestimate_bad_input")
})

test_that("simplicial verbs input checks are classed", {
  expect_error(build_simplicial(sq, threshold = -1),
               class = "nestimate_bad_input")
  expect_error(betti_numbers(1), class = "nestimate_bad_input")
  expect_error(euler_characteristic(1), class = "nestimate_bad_input")
  expect_error(simplicial_degree(1), class = "nestimate_bad_input")
  expect_error(q_analysis(1), class = "nestimate_bad_input")
  expect_error(persistent_homology(sq, n_steps = 0),
               class = "nestimate_bad_input")
  expect_error(simplicial_features(sq, threshold = NA),
               class = "nestimate_bad_input")
  expect_error(simplicial_features("x"), class = "nestimate_bad_input")
})

test_that("persistence diagram verbs input checks are classed", {
  expect_error(bottleneck_distance(1, 1), class = "nestimate_bad_input")
  expect_error(bottleneck_distance(data.frame(birth = 0), data.frame(birth = 0)),
               class = "nestimate_bad_input")
  expect_error(persistence_landscape(1, k_max = 0),
               class = "nestimate_bad_input")
})

test_that("sequence_plot() and distribution_plot() input checks are classed", {
  expect_error(sequence_plot(list(1)), class = "nestimate_bad_input")
  expect_error(sequence_plot(seqs, combined = "yes"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(seqs, trim = c(1, 2, 3)),
               class = "nestimate_bad_input")
  expect_error(distribution_plot(seqs, na = "yes"),
               class = "nestimate_bad_input")
})

test_that("state colour verbs input checks are classed", {
  expect_error(set_state_colors(1, c(A = "#000000")),
               class = "nestimate_bad_input")
  expect_error(state_colors(1), class = "nestimate_bad_input")
})
