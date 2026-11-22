# The sequence estimators of build_network() beyond the four classical ones: n-gram, gap and reverse
# (tna::build_model() types "n-gram", "gap", "reverse"), one-hot data with windows, actors and sessions
# (tna::import_onehot() then ftna()/tna()/ctna()), and action_to_onehot() against its definition. tna is the
# owner's co-authored package, a separate implementation (open-source; to be confirmed by the owner).
# Re-expressed from local_testing_and_equivalence/test-equiv-moved-estimators-tna-parity.R and
# test-equiv-moved-import_onehot.R (100 random one-hot data sets through an internal helper; here a few fixed
# inputs through the exported build_network(codes = )). Sliding windows are compared in test-wtna.R.

# Sequences of three states with about 10% missing cells, the input of the old parity test.
tn_parity_sequences <- function() {
  set.seed(20260829)
  as.data.frame(matrix(sample(c("a", "b", "c", NA), 60, TRUE, prob = c(.3, .3, .3, .1)), nrow = 12, ncol = 5))
}

# One-hot data: rows are time points, columns A..D are 0/1 indicators, two actors.
tn_onehot_data <- function(seed, n = 24L) {
  set.seed(seed)
  codes <- c("A", "B", "C", "D")
  data <- as.data.frame(matrix(sample(0:1, n * 4L, TRUE, prob = c(.55, .45)), n, dimnames = list(NULL, codes)))
  data$actor <- rep(c("x", "y"), each = n / 2L)
  data$session <- rep(c("s1", "s2"), times = n / 2L)
  data
}

test_that("n-gram, gap and reverse weights equal tna::build_model()", {
  skip_without_reference("tna")
  sequences <- tn_parity_sequences()
  cases <- list(
    list(method = "ngram", type = "n-gram", params = list(), label = "n-gram, n 2"),
    list(method = "ngram", type = "n-gram", params = list(n_gram = 3), label = "n-gram, n 3"),
    list(method = "ngram", type = "n-gram", params = list(n_gram = 4), label = "n-gram, n 4"),
    list(method = "gap", type = "gap", params = list(), label = "gap, default"),
    list(method = "gap", type = "gap", params = list(max_gap = 2), label = "gap, max_gap 2"),
    list(method = "gap", type = "gap", params = list(max_gap = 0), label = "gap, max_gap 0"),
    list(method = "reverse", type = "reverse", params = list(), label = "reverse"),
    list(method = "reverse", type = "reverse", params = list(weighted = TRUE), label = "reverse, weighted")
  )
  aliases <- list(ngram = "build_ngram", gap = "build_gap", reverse = "build_reverse")
  invisible(lapply(cases, \(case) {
    expected <- unclass(tna::build_model(sequences, type = case$type, params = case$params)$weights)
    # The old test used testthat's default tolerance, sqrt(machine epsilon), relative; the weights are of order
    # one or counts, so 1e-10 relative is at least as strict and is what a one-division computation allows.
    expect_reference(extract_transition_matrix(build_network(sequences, method = case$method, params = case$params)),
                     expected, reference = "tna", covers = "build_network()", quantity = case$label, tolerance = 1e-10,
                     measure = "relative")
    # the method's alias, with its setting given by name
    alias <- aliases[[case$method]]
    expect_reference(extract_transition_matrix(do.call(alias, c(list(sequences), case$params))),
                     expected, reference = "tna", covers = paste0(alias, "()"), quantity = case$label,
                     tolerance = 1e-10, measure = "relative")
  }))
  regulation <- utils::head(tna::group_regulation, 200)
  expect_reference(extract_transition_matrix(build_network(regulation, method = "ngram", params = list(n_gram = 3))),
                   unclass(tna::build_model(regulation, type = "n-gram", params = list(n_gram = 3))$weights),
                   reference = "tna", covers = "build_network()", quantity = "n-gram, n 3, group_regulation rows 1-200",
                   tolerance = 1e-10, measure = "relative")
})

test_that("one-hot data with windows and actors gives the networks of tna::import_onehot()", {
  skip_without_reference("tna")
  codes <- c("A", "B", "C", "D")
  cases <- list(
    list(seed = 1L, window_size = 2L, mode = "non-overlapping", type = "tumbling", actor = "actor", session = NULL),
    list(seed = 2L, window_size = 3L, mode = "non-overlapping", type = "tumbling", actor = "actor", session = NULL),
    list(seed = 3L, window_size = 1L, mode = "non-overlapping", type = "tumbling", actor = "actor", session = NULL)
  )
  invisible(lapply(cases, \(case) {
    data <- tn_onehot_data(case$seed)
    # tna 1.3.1, the version on the ci machine, has no `session` argument in import_onehot(): sessions are not
    # compared.
    imported <- tna::import_onehot(data, cols = codes, window_size = case$window_size, window_type = case$type,
                                   actor = case$actor)
    label <- sprintf("seed %d, window %d %s, actor", case$seed, case$window_size, case$mode)
    # Window 1 without sessions draws Nestimate's advisory warning that the network will be sparse.
    fit <- \(method) withCallingHandlers(
      build_network(data, method = method, codes = codes, window_size = case$window_size, mode = case$mode,
                    actor = case$actor),
      nestimate_sparse_onehot = \(w) invokeRestart("muffleWarning"))
    expect_reference(extract_transition_matrix(fit("frequency")), unclass(tna::ftna(imported)$weights),
                     reference = "tna", covers = "build_network()", quantity = sprintf("one-hot frequency, %s", label),
                     tolerance = 1e-10)
    expect_reference(extract_transition_matrix(fit("relative")), unclass(tna::tna(imported)$weights),
                     reference = "tna", covers = "build_network()", quantity = sprintf("one-hot relative, %s", label),
                     tolerance = 1e-10)
    expect_reference(extract_transition_matrix(fit("co_occurrence")), unclass(tna::ctna(imported)$weights),
                     reference = "tna", covers = "build_network()", quantity = sprintf("one-hot co-occurrence, %s", label),
                     tolerance = 1e-10)
  }))
})

test_that("action_to_onehot() writes one 0/1 indicator column per state", {
  events <- data.frame(Actor = rep(1:4, each = 3), Action = c("b", "a", "c", "a", "a", "b", "c", "b", "a", "c", "c", "a"))
  states <- sort(unique(events$Action))
  indicators <- vapply(states, \(state) as.integer(events$Action == state), integer(nrow(events)))
  expected <- cbind(data.frame(row = seq_len(nrow(events))), as.data.frame(indicators))
  found <- cbind(data.frame(row = seq_len(nrow(events))), action_to_onehot(events, sort_states = TRUE))
  expect_reference(found, expected, reference = "independent: one indicator per state, base R",
                   covers = "action_to_onehot()", quantity = "indicators of 12 events", tolerance = 0, key = "row")
})
