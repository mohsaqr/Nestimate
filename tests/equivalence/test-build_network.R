# build_network() and its wrappers on wide sequence data: relative (build_tna()), frequency (build_ftna()),
# co-occurrence (build_cna()) and attention (build_atna()) weights and initial probabilities, against tna's
# tna(), ftna(), ctna() and atna() (a separate implementation by the same authors: open-source, to be confirmed by
# the owner) and against the definitions in base R (helper-transitions.R). Re-expressed from
# local_testing_and_equivalence/test-equiv-tna.R and test-build_network-equiv-extracted.R, whose sweeps of 30
# random configurations stay there as a standing study; here five fixed, seeded inputs.
# The weights are read with extract_transition_matrix() and the initial probabilities with
# extract_initial_probs(). A grouped result (netobject_group) has no accessor that returns all groups' weights,
# so one group is taken with `[[` (API gap).

inputs <- list(
  small = tn_simulate(n_actors = 15L, n_states = 3L, seq_length = 12L, seed = 101L),
  medium = tn_simulate(n_actors = 20L, n_states = 5L, seq_length = 20L, seed = 102L),
  large = tn_simulate(n_actors = 30L, n_states = 7L, seq_length = 30L, seed = 103L),
  ragged = tn_simulate(n_actors = 25L, n_states = 4L, seq_length = 15L, seed = 104L, ragged = TRUE),
  ragged_wide = tn_simulate(n_actors = 30L, n_states = 6L, seq_length = 25L, seed = 105L, ragged = TRUE)
)

# Counts are integers and a probability is one division of two integers: rounding of the last bits only.
tolerance <- 1e-10

test_that("relative and frequency weights equal the table of consecutive pairs", {
  invisible(lapply(names(inputs), \(name) {
    sequences <- tn_sequences(inputs[[name]])
    counts <- tn_pair_counts(sequences)
    expect_reference(extract_transition_matrix(build_network(inputs[[name]], method = "relative")),
                     tn_row_normalise(counts), reference = "independent: table of consecutive pairs, base R",
                     covers = c("build_network()", "extract_transition_matrix()"),
                     quantity = sprintf("relative weights, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_tna(inputs[[name]])),
                     tn_row_normalise(counts), reference = "independent: table of consecutive pairs, base R",
                     covers = "build_tna()", quantity = sprintf("build_tna() weights, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_ftna(inputs[[name]])), counts,
                     reference = "independent: table of consecutive pairs, base R",
                     covers = c("build_ftna()", "build_network()"),
                     quantity = sprintf("frequency weights, %s", name), tolerance = tolerance)
    expect_reference(extract_initial_probs(build_tna(inputs[[name]])), tn_initial(sequences),
                     reference = "independent: proportion of sequences starting in each state, base R",
                     covers = c("build_tna()", "extract_initial_probs()"),
                     quantity = sprintf("initial probabilities, %s", name), tolerance = tolerance)
  }))
})

test_that("co-occurrence and attention weights equal their definitions", {
  invisible(lapply(c("small", "ragged", "ragged_wide"), \(name) {
    sequences <- tn_sequences(inputs[[name]])
    expect_reference(extract_transition_matrix(build_cna(inputs[[name]])), tn_cooccurrence(sequences),
                     reference = "independent: co-occurrence of every pair of positions, base R",
                     covers = c("build_cna()", "build_network()"),
                     quantity = sprintf("co-occurrence weights, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_atna(inputs[[name]])), tn_attention(sequences, lambda = 1),
                     reference = "independent: exp(-(j - i) / lambda) summed over ordered pairs, base R",
                     covers = c("build_atna()", "build_network()"),
                     quantity = sprintf("attention weights lambda 1, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_network(inputs[[name]], method = "attention", params = list(lambda = 2))),
                     tn_attention(sequences, lambda = 2),
                     reference = "independent: exp(-(j - i) / lambda) summed over ordered pairs, base R",
                     covers = "build_network()", quantity = sprintf("attention weights lambda 2, %s", name),
                     tolerance = tolerance)
  }))
})

test_that("relative, frequency, co-occurrence and attention weights equal tna", {
  skip_without_reference("tna")
  invisible(lapply(names(inputs), \(name) {
    wide <- inputs[[name]]
    expect_reference(extract_transition_matrix(build_tna(wide)), unclass(tna::tna(wide)$weights), reference = "tna",
                     covers = c("build_tna()", "build_network()", "extract_transition_matrix()"),
                     quantity = sprintf("relative weights, %s", name), tolerance = tolerance)
    expect_reference(extract_initial_probs(build_tna(wide)), tna::tna(wide)$inits, reference = "tna",
                     covers = c("build_tna()", "extract_initial_probs()"),
                     quantity = sprintf("initial probabilities, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_ftna(wide)), unclass(tna::ftna(wide)$weights), reference = "tna",
                     covers = c("build_ftna()", "build_network()"),
                     quantity = sprintf("frequency weights, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_cna(wide)), unclass(tna::ctna(wide)$weights), reference = "tna",
                     covers = c("build_cna()", "build_network()"),
                     quantity = sprintf("co-occurrence weights, %s", name), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_atna(wide)), unclass(tna::atna(wide)$weights), reference = "tna",
                     covers = c("build_atna()", "build_network()"),
                     quantity = sprintf("attention weights, %s", name), tolerance = tolerance)
  }))
})

test_that("attention weights at several decay rates equal tna::atna()", {
  skip_without_reference("tna")
  invisible(lapply(c(0.5, 2, 5), \(lambda) {
    expect_reference(extract_transition_matrix(build_network(inputs$medium, method = "attention", params = list(lambda = lambda))),
                     unclass(tna::atna(inputs$medium, params = list(lambda = lambda))$weights), reference = "tna",
                     covers = "build_network()", quantity = sprintf("attention weights lambda %s, medium", lambda),
                     tolerance = tolerance)
  }))
})

test_that("grouped relative and frequency networks equal tna's group_tna() and group_ftna()", {
  skip_without_reference("tna")
  grouped <- cbind(group = rep(c("G1", "G2", "G3"), length.out = nrow(inputs$large)), inputs$large)
  fit_relative <- build_network(grouped, method = "relative", group = "group")
  fit_frequency <- build_network(grouped, method = "frequency", group = "group")
  ref_relative <- tna::group_tna(grouped, group = "group")
  ref_frequency <- tna::group_ftna(grouped, group = "group")
  invisible(lapply(c("G1", "G2", "G3"), \(group) {
    expect_reference(extract_transition_matrix(fit_relative[[group]]), unclass(ref_relative[[group]]$weights),
                     reference = "tna", covers = "build_network()",
                     quantity = sprintf("grouped relative weights, %s", group), tolerance = tolerance)
    expect_reference(extract_initial_probs(fit_relative[[group]]), ref_relative[[group]]$inits,
                     reference = "tna", covers = "build_network()",
                     quantity = sprintf("grouped initial probabilities, %s", group), tolerance = tolerance)
    expect_reference(extract_transition_matrix(fit_frequency[[group]]), unclass(ref_frequency[[group]]$weights),
                     reference = "tna", covers = "build_network()",
                     quantity = sprintf("grouped frequency weights, %s", group), tolerance = tolerance)
  }))
})

test_that("a group that lacks a state keeps the global alphabet, as in tna::group_tna()", {
  skip_without_reference("tna")
  data <- data.frame(V1 = c("A", "A", "B", "C"), V2 = c("B", "B", "C", NA), group = c("g1", "g2", "g2", "g1"))
  fit <- build_network(data, method = "relative", group = "group")
  ref <- tna::group_tna(data, group = "group")
  invisible(lapply(c("g1", "g2"), \(group) {
    expect_reference(extract_transition_matrix(fit[[group]]), unclass(ref[[group]]$weights), reference = "tna",
                     covers = "build_network()", quantity = sprintf("global alphabet weights, %s", group),
                     tolerance = tolerance)
    expect_reference(extract_initial_probs(fit[[group]]), ref[[group]]$inits, reference = "tna",
                     covers = "build_network()", quantity = sprintf("global alphabet initial probabilities, %s", group),
                     tolerance = tolerance)
  }))
})

test_that("edge cases equal tna: no transitions, weighted counts, a matrix input, min-max scaling", {
  skip_without_reference("tna")
  empty <- data.frame(V1 = c("A", "B"), V2 = c(NA, NA))
  expect_reference(extract_transition_matrix(build_tna(empty)), unclass(tna::tna(empty)$weights), reference = "tna",
                   covers = "build_tna()", quantity = "no transitions, weights", tolerance = tolerance)
  expect_reference(extract_initial_probs(build_tna(empty)), tna::tna(empty)$inits, reference = "tna",
                   covers = "build_tna()", quantity = "no transitions, initial probabilities", tolerance = tolerance)
  data <- data.frame(V1 = c("A", "A", "B"), V2 = c("B", "C", "C"), V3 = c("C", "C", NA))
  expect_reference(extract_transition_matrix(build_network(data, method = "frequency", params = list(weighted = TRUE))),
                   unclass(tna::ftna(data, params = list(weighted = TRUE))$weights), reference = "tna",
                   covers = "build_network()", quantity = "weighted frequency", tolerance = tolerance)
  expect_reference(extract_transition_matrix(build_network(data, method = "co_occurrence", params = list(weighted = TRUE))),
                   unclass(tna::ctna(data, params = list(weighted = TRUE))$weights), reference = "tna",
                   covers = "build_network()", quantity = "weighted co-occurrence", tolerance = tolerance)
  counts <- matrix(c(1, 2, 3, 4), 2L, 2L, dimnames = list(c("A", "B"), c("A", "B")))
  expect_reference(extract_transition_matrix(build_network(counts, method = "relative")),
                   unclass(tna::tna(counts)$weights), reference = "tna",
                   covers = "build_network()", quantity = "matrix input", tolerance = tolerance)
  expect_reference(extract_transition_matrix(build_network(data, method = "frequency", scaling = "minmax")),
                   unclass(tna::ftna(data, scaling = "minmax")$weights), reference = "tna",
                   covers = "build_network()", quantity = "min-max scaling", tolerance = tolerance)
})

test_that("rank scaling ranks only the edges that exist, where tna ranks every cell", {
  skip_without_reference("tna")
  data <- data.frame(V1 = c("A", "A", "B"), V2 = c("B", "C", "C"), V3 = c("C", "C", NA))
  expect_known_difference(extract_transition_matrix(build_network(data, method = "frequency", scaling = "rank")),
                          unclass(tna::ftna(data, scaling = "rank")$weights), reference = "tna",
                          covers = "build_network()", quantity = "rank scaling", tolerance = tolerance,
                          reason = "Nestimate keeps a zero cell at zero and ranks the non-zero weights only; tna ranks every cell, zeros included.")
  counts <- tn_pair_counts(tn_sequences(data))
  expected <- counts
  expected[counts > 0] <- rank(counts[counts > 0])
  expect_reference(extract_transition_matrix(build_network(data, method = "frequency", scaling = "rank")),
                   expected, reference = "independent: base::rank() of the non-zero counts",
                   covers = "build_network()", quantity = "rank scaling of non-zero counts", tolerance = tolerance)
})

test_that("summary() of a transition matrix and of initial probabilities equal tna's weights and inits", {
  skip_without_reference("tna")
  wide <- inputs$medium
  reference <- tna::tna(wide)
  weights <- unclass(reference$weights)
  edges <- data.frame(from = rep(rownames(weights), times = ncol(weights)),
                      to = rep(colnames(weights), each = nrow(weights)), weight = as.vector(weights))
  expect_reference(summary(extract_transition_matrix(build_tna(wide))), subset(edges, weight > 0),
                   reference = "tna", covers = "summary.nest_transition_matrix()",
                   quantity = "edges of the transition matrix", tolerance = tolerance, key = c("from", "to"))
  expect_reference(summary(extract_initial_probs(build_tna(wide))),
                   data.frame(state = names(reference$inits), prob = unname(reference$inits)),
                   reference = "tna", covers = "summary.nest_initial_probs()",
                   quantity = "initial probabilities", tolerance = tolerance, key = "state")
})

test_that("extract_edges() lists the non-zero, non-loop weights of tna's model", {
  skip_without_reference("tna")
  weights <- unclass(tna::tna(inputs$medium)$weights)
  edges <- data.frame(from = rep(rownames(weights), times = ncol(weights)),
                      to = rep(colnames(weights), each = nrow(weights)), weight = as.vector(weights))
  expect_reference(extract_edges(build_tna(inputs$medium)), subset(edges, weight > 0 & from != to),
                   reference = "tna", covers = "extract_edges()", quantity = "edge list", tolerance = tolerance,
                   key = c("from", "to"))
})

test_that("partial correlations of a well-conditioned table equal the scaled inverse of cor()", {
  # Moved from tests/testthat/test-build_network.R ("pcor output unchanged for well-conditioned data"), which
  # compared to 1e-14 relative; the values are of order one.
  data <- data.frame(V1 = c(1, 3, 2, 5, 4, 6, 8, 7), V2 = c(2, 1, 4, 3, 6, 5, 7, 9),
                     V3 = c(5, 3, 4, 1, 2, 0, 1, 2), V4 = c(1, 2, 2, 3, 5, 4, 3, 6))
  precision <- solve(stats::cor(as.matrix(data)))
  scale <- diag(1 / sqrt(diag(precision)))
  partial <- -scale %*% precision %*% scale
  diag(partial) <- 0
  dimnames(partial) <- dimnames(precision)
  expect_reference(extract_transition_matrix(build_network(data, method = "pcor")), partial,
                   reference = "independent: -D solve(cor) D with a zero diagonal, base R",
                   covers = "build_network()", quantity = "partial correlations, method pcor", tolerance = 1e-14,
                   measure = "relative")
})
