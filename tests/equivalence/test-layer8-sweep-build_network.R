# Layer 8, the argument sweep of build_network() on sequence data. Every transition method is fitted at every
# boundary setting (start, end: FALSE / TRUE / a label), every scaling (alone and chained) with thresholds, the
# estimator settings each method documents (lambda, direction, n_gram, max_gap, type, weighted, concat, windows),
# every input format (wide, event log with actor/time, event log with session, grouped), and every alias, and each
# fit is compared with a reference that does not use the code under test: the definitions written in base R
# (helper-layer-sweep.R) and tna::build_model() (a separate implementation; tna is the reference package of the
# transition methods). A default left out must equal the same default written out. A combination the package
# refuses is asserted by its condition class. test-build_network.R and test-estimators.R check the same methods
# at their defaults; nothing here repeats one of their settings without a new option beside it.

sweep_inputs <- list(
  ragged = tn_simulate(n_actors = 30L, n_states = 4L, seq_length = 9L, seed = 801L, ragged = TRUE),
  full = tn_simulate(n_actors = 18L, n_states = 3L, seq_length = 7L, seed = 802L)
)

# Every documented setting of every sequence method (?build_network, Estimators).
method_grid <- list(
  list(method = "relative", params = list(), label = "relative"),
  list(method = "frequency", params = list(), label = "frequency"),
  list(method = "co_occurrence", params = list(), label = "co_occurrence"),
  list(method = "co_occurrence", params = list(type = "relative"), label = "co_occurrence type relative"),
  list(method = "attention", params = list(), label = "attention"),
  list(method = "attention", params = list(lambda = 2.5), label = "attention lambda 2.5"),
  list(method = "attention", params = list(direction = "backward"), label = "attention backward"),
  list(method = "attention", params = list(direction = "both", lambda = 0.5), label = "attention both lambda 0.5"),
  list(method = "ngram", params = list(), label = "ngram"),
  list(method = "ngram", params = list(n_gram = 3L), label = "ngram 3"),
  list(method = "ngram", params = list(n_gram = 5L), label = "ngram 5"),
  list(method = "gap", params = list(), label = "gap"),
  list(method = "gap", params = list(max_gap = 0L), label = "gap 0"),
  list(method = "gap", params = list(max_gap = 3L), label = "gap 3"),
  list(method = "reverse", params = list(), label = "reverse"),
  list(method = "reverse", params = list(weighted = TRUE), label = "reverse weighted")
)

# start / end: FALSE, TRUE (default label) and a custom label.
boundary_grid <- expand.grid(start = c("FALSE", "TRUE", "Begin"), end = c("FALSE", "TRUE", "Stop"),
                             stringsAsFactors = FALSE)
sw_flag <- function(value) switch(value, `FALSE` = FALSE, `TRUE` = TRUE, value)
sw_label <- function(value, default) switch(value, `FALSE` = NULL, `TRUE` = default, value)
sw_directed <- function(config) !identical(config$method, "co_occurrence") || identical(config$params$type, "relative")
# Counts are exact; one division or one exp() per cell leaves the last bits only.
sw_tolerance <- 1e-10

# The reference weights of one configuration on one input, with its boundary labels.
sw_expected <- function(wide, config, start = NULL, end = NULL) {
  mat <- sw_positions(wide, start, end)
  states <- sw_states(mat, start, end)
  list(weights = sw_reference(mat, states, config$method, config$params), initial = sw_initial(mat, states),
       mat = mat, states = states)
}

test_that("every sequence method at every start/end setting equals its definition (wide input)", {
  cases <- expand.grid(input = names(sweep_inputs), config = seq_along(method_grid), boundary = seq_len(nrow(boundary_grid)),
                       stringsAsFactors = FALSE)
  invisible(Map(\(input, config_index, boundary) {
    config <- method_grid[[config_index]]
    start <- boundary_grid$start[[boundary]]
    end <- boundary_grid$end[[boundary]]
    wide <- sweep_inputs[[input]]
    expected <- sw_expected(wide, config, sw_label(start, "Start"), sw_label(end, "End"))
    fit <- build_network(wide, method = config$method, start = sw_flag(start), end = sw_flag(end),
                         params = config$params)
    label <- sprintf("%s, start %s, end %s, %s", config$label, start, end, input)
    expect_reference(extract_transition_matrix(fit), expected$weights,
                     reference = "independent: pair counts of the positions with the boundary states, base R",
                     covers = c("build_network()", "extract_transition_matrix()"), quantity = sprintf("weights, %s", label),
                     tolerance = sw_tolerance)
    if (!identical(config$method, "co_occurrence")) {
      expect_reference(extract_initial_probs(fit), expected$initial,
                       reference = "independent: share of sequences by first observed cell, base R",
                       covers = c("build_network()", "extract_initial_probs()"),
                       quantity = sprintf("initial probabilities, %s", label), tolerance = sw_tolerance)
    }
    expect_reference(as.numeric(fit$directed), as.numeric(sw_directed(config)),
                     reference = "independent: ?build_network, Estimators", covers = "build_network()",
                     quantity = sprintf("directed, %s", label), tolerance = 0)
  }, cases$input, cases$config, cases$boundary))
})

test_that("every tna-parity method at every start/end setting equals tna::build_model(begin_state, end_state)", {
  skip_without_reference("tna")
  tna_types <- c(relative = "relative", frequency = "frequency", co_occurrence = "co-occurrence",
                 attention = "attention", ngram = "n-gram", gap = "gap", reverse = "reverse")
  # tna has no `direction` for attention and no `type` for co-occurrence: those settings have the base-R reference only.
  comparable <- Filter(\(config) is.null(config$params$direction) && is.null(config$params$type), method_grid)
  cases <- expand.grid(input = names(sweep_inputs), config = seq_along(comparable), boundary = seq_len(nrow(boundary_grid)),
                       stringsAsFactors = FALSE)
  invisible(Map(\(input, config_index, boundary) {
    config <- comparable[[config_index]]
    start <- sw_label(boundary_grid$start[[boundary]], "Start")
    end <- sw_label(boundary_grid$end[[boundary]], "End")
    wide <- sweep_inputs[[input]]
    tna_args <- c(list(wide, type = tna_types[[config$method]], params = config$params),
                  if (!is.null(start)) list(begin_state = start), if (!is.null(end)) list(end_state = end))
    reference <- do.call(tna::build_model, tna_args)
    fit <- build_network(wide, method = config$method, start = sw_flag(boundary_grid$start[[boundary]]),
                         end = sw_flag(boundary_grid$end[[boundary]]), params = config$params)
    label <- sprintf("%s, start %s, end %s, %s", config$label, boundary_grid$start[[boundary]],
                     boundary_grid$end[[boundary]], input)
    expect_reference(extract_transition_matrix(fit), unclass(reference$weights), reference = "tna",
                     covers = "build_network()", quantity = sprintf("weights, %s", label), tolerance = sw_tolerance)
    if (!identical(config$method, "co_occurrence")) {
      expect_reference(extract_initial_probs(fit), reference$inits, reference = "tna", covers = "build_network()",
                       quantity = sprintf("initial probabilities, %s", label), tolerance = sw_tolerance)
    }
  }, cases$input, cases$config, cases$boundary))
})

test_that("an end state on sequences that all stop before the last column: n-gram equals tna", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  # Every sequence leaves the last column empty: tna puts the end state in the empty cell and adds no column.
  wide <- tn_simulate(n_actors = 20L, n_states = 3L, seq_length = 8L, seed = 803L, ragged = TRUE)
  wide$T9 <- NA_character_
  invisible(lapply(c(2L, 3L, 4L), \(n) {
    reference <- tna::build_model(wide, type = "n-gram", params = list(n_gram = n), end_state = "End")
    fit <- build_network(wide, method = "ngram", end = TRUE, params = list(n_gram = n))
    # Like tna, the n-gram estimator adds a column for the end state only when a sequence fills the last column
    # (.prepare_sequence_wide(end_column = "needed")), so the windows are tna's.
    expect_reference(extract_transition_matrix(fit), unclass(reference$weights), reference = "tna",
                     covers = "build_network()", quantity = sprintf("n-gram %d with end, last column empty", n),
                     tolerance = sw_tolerance)
  }))
})

test_that("an empty trailing position column leaves every method's weights unchanged", {
  # A wide table padded with one more all-missing position describes the same sequences.
  wide <- sweep_inputs$ragged
  padded <- cbind(wide, T10 = NA_character_)
  invisible(lapply(method_grid, \(config) {
    expect_reference(extract_transition_matrix(build_network(padded, method = config$method, params = config$params)),
                     extract_transition_matrix(build_network(wide, method = config$method, params = config$params)),
                     reference = "independent: the same sequences without the padding column",
                     covers = "build_network()", quantity = sprintf("padding invariance, %s", config$label),
                     tolerance = sw_tolerance,
                     # n-gram windows run over the table's columns (tna's definition), so an empty trailing column
                     # adds windows that contain the last observed pairs.
                     defect = if (identical(config$method, "ngram") && (config$params$n_gram %||% 2L) >= 3L)
                       "ngram (n_gram >= 3) counts change with trailing all-missing padding columns: windows run over columns, not over each sequence (inherited from tna)")
  }))
})

test_that("an event log with actor and time, and one with session, gives the networks of its wide table", {
  wide <- sweep_inputs$ragged
  long <- sw_long(wide)
  long$time <- as.POSIXct(long$time, origin = "1970-01-01", tz = "UTC")
  cases <- expand.grid(config = seq_along(method_grid), boundary = c(1L, 5L, 9L, 3L), route = c("actor", "session"),
                       stringsAsFactors = FALSE)
  invisible(Map(\(config_index, boundary, route) {
    config <- method_grid[[config_index]]
    start <- boundary_grid$start[[boundary]]
    end <- boundary_grid$end[[boundary]]
    expected <- sw_expected(wide, config, sw_label(start, "Start"), sw_label(end, "End"))
    fit <- if (identical(route, "actor")) {
      build_network(long, method = config$method, action = "action", actor = "actor", time = "time",
                    start = sw_flag(start), end = sw_flag(end), params = config$params)
    } else {
      build_network(long, method = config$method, action = "action", session = "actor", time = "time",
                    start = sw_flag(start), end = sw_flag(end), params = config$params)
    }
    label <- sprintf("%s, start %s, end %s, event log by %s", config$label, start, end, route)
    expect_reference(extract_transition_matrix(fit), expected$weights,
                     reference = "independent: pair counts of the positions with the boundary states, base R",
                     covers = "build_network()", quantity = sprintf("weights, %s", label), tolerance = sw_tolerance)
    if (!identical(config$method, "co_occurrence")) {
      expect_reference(extract_initial_probs(fit), expected$initial,
                       reference = "independent: share of sequences by first observed cell, base R",
                       covers = c("build_network()", "extract_initial_probs()"),
                       quantity = sprintf("initial probabilities, %s", label), tolerance = sw_tolerance)
    }
  }, cases$config, cases$boundary, cases$route))
})

test_that("every scaling, alone and chained, with every threshold, equals the documented rescaling", {
  wide <- sweep_inputs$ragged
  scalings <- list("minmax", "max", "rank", "normalize", c("rank", "minmax"), c("normalize", "max"),
                   c("max", "rank"), c("minmax", "normalize"))
  thresholds <- c(0, 0.05, 0.3)
  configs <- method_grid[c(1L, 2L, 3L, 4L, 5L, 9L, 12L, 15L)]
  cases <- expand.grid(config = seq_along(configs), scaling = seq_along(scalings), threshold = thresholds,
                       stringsAsFactors = FALSE)
  invisible(Map(\(config_index, scaling_index, threshold) {
    config <- configs[[config_index]]
    scaling <- scalings[[scaling_index]]
    raw <- sw_expected(wide, config)$weights
    expected <- sw_threshold(sw_scale(raw, scaling, include_zeros = TRUE), threshold)
    fit <- build_network(wide, method = config$method, params = config$params, scaling = scaling,
                         threshold = threshold)
    label <- sprintf("%s, scaling %s, threshold %s", config$label, paste(scaling, collapse = "+"), threshold)
    expect_reference(extract_transition_matrix(fit), expected,
                     reference = "independent: ?build_network scaling then threshold, base R",
                     covers = "build_network()", quantity = label, tolerance = sw_tolerance)
    expect_reference(nrow(fit$edges), sum(expected != 0 & (sw_directed(config) | row(expected) <= col(expected))),
                     reference = "independent: non-zero cells (upper triangle with diagonal when undirected)",
                     covers = "build_network()", quantity = sprintf("edge count, %s", label), tolerance = 0)
  }, cases$config, cases$scaling, cases$threshold))
})

test_that("min-max scaling of a transition network spans every cell, as documented", {
  # ?build_network: for the transition methods minmax rescales every cell, zeros included (tna's rule), so on a
  # network where every transition occurs the smallest weight becomes 0 and the others stay positive.
  dense <- tn_simulate(n_actors = 40L, n_states = 3L, seq_length = 30L, seed = 804L)
  raw <- extract_transition_matrix(build_network(dense, method = "frequency"))
  scaled <- extract_transition_matrix(build_network(dense, method = "frequency", scaling = "minmax"))
  expect_reference(as.numeric(sum(raw == 0)), 0, reference = "independent: every pair occurs in this input",
                   covers = "build_network()", quantity = "dense input has no zero transition", tolerance = 0)
  expect_reference(as.numeric(sum(scaled > 0)), as.numeric(sum(raw > min(raw))),
                   reference = "independent: ?build_network, scaling = \"minmax\" on a transition network spans every cell",
                   covers = "build_network()", quantity = "non-zero weights after min-max, dense frequency network",
                   tolerance = 0)
})

test_that("weighted and concat settings of every tna-parity method equal tna::build_model()", {
  skip_without_reference("tna")
  wide <- sweep_inputs$ragged
  tna_types <- c(relative = "relative", frequency = "frequency", co_occurrence = "co-occurrence",
                 attention = "attention", ngram = "n-gram", gap = "gap", reverse = "reverse")
  settings <- list(list(weighted = TRUE), list(concat = 2L), list(concat = 3L), list(weighted = TRUE, concat = 2L))
  cases <- expand.grid(method = names(tna_types), setting = seq_along(settings), boundary = c(1L, 9L),
                       stringsAsFactors = FALSE)
  invisible(Map(\(method, setting_index, boundary) {
    params <- settings[[setting_index]]
    start <- sw_label(boundary_grid$start[[boundary]], "Start")
    end <- sw_label(boundary_grid$end[[boundary]], "End")
    # tna takes concat as an argument of build_model(), not in params.
    tna_args <- c(list(wide, type = tna_types[[method]], params = params[setdiff(names(params), "concat")]),
                  params["concat"][!is.null(params$concat)],
                  if (!is.null(start)) list(begin_state = start), if (!is.null(end)) list(end_state = end))
    reference <- do.call(tna::build_model, tna_args)
    fit <- build_network(wide, method = method, params = params, start = sw_flag(boundary_grid$start[[boundary]]),
                         end = sw_flag(boundary_grid$end[[boundary]]))
    label <- sprintf("%s, %s, start %s, end %s", method,
                     paste(names(params), unlist(params), sep = " ", collapse = ", "),
                     boundary_grid$start[[boundary]], boundary_grid$end[[boundary]])
    if (identical(method, "relative") && isTRUE(params$weighted)) {
      expect_known_difference(extract_transition_matrix(fit), unclass(reference$weights), reference = "tna",
                              covers = "build_network()", quantity = sprintf("weights, %s", label),
                              tolerance = sw_tolerance,
                              reason = "tna's compute_transitions() applies `weighted` to frequency only and ignores it for relative; Nestimate weights each pair by 1 / sequence length before row-normalising (checked against that definition below).",
                              agrees_with = "weighted = FALSE")
    } else {
      expect_reference(extract_transition_matrix(fit), unclass(reference$weights), reference = "tna",
                       covers = "build_network()", quantity = sprintf("weights, %s", label), tolerance = sw_tolerance)
    }
  }, cases$method, cases$setting, cases$boundary))
})

test_that("weighted relative, frequency and co-occurrence weights equal the inverse-length definition", {
  wide <- sweep_inputs$ragged
  configs <- list(list(method = "relative", params = list(weighted = TRUE)),
                  list(method = "frequency", params = list(weighted = TRUE)),
                  list(method = "co_occurrence", params = list(weighted = TRUE)),
                  list(method = "reverse", params = list(weighted = TRUE)))
  cases <- expand.grid(config = seq_along(configs), boundary = c(1L, 5L, 9L), stringsAsFactors = FALSE)
  invisible(Map(\(config_index, boundary) {
    config <- configs[[config_index]]
    start <- boundary_grid$start[[boundary]]
    end <- boundary_grid$end[[boundary]]
    expected <- sw_expected(wide, config, sw_label(start, "Start"), sw_label(end, "End"))
    fit <- build_network(wide, method = config$method, params = config$params, start = sw_flag(start),
                         end = sw_flag(end))
    expect_reference(extract_transition_matrix(fit), expected$weights,
                     reference = "independent: each pair weighted 1 / observed cells of its row, base R",
                     covers = "build_network()",
                     quantity = sprintf("weighted %s, start %s, end %s", config$method, start, end),
                     tolerance = sw_tolerance)
  }, cases$config, cases$boundary))
})

test_that("windows on sequence data (window_size, window_span, windowed) equal tna::build_model()", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  wide <- tn_simulate(n_actors = 24L, n_states = 3L, seq_length = 12L, seed = 805L, ragged = TRUE)
  tna_types <- c(relative = "relative", frequency = "frequency", co_occurrence = "co-occurrence")
  windows <- list(list(window_size = 2L), list(window_size = 3L), list(window_size = 2L, window_span = 2L),
                  list(window_span = 3L), list(windowed = TRUE), list(window_size = 4L, weighted = TRUE))
  cases <- expand.grid(method = names(tna_types), window = seq_along(windows), stringsAsFactors = FALSE)
  invisible(Map(\(method, window_index) {
    params <- windows[[window_index]]
    tna_params <- params
    # tna reads a window only with windowed = TRUE.
    tna_params$windowed <- TRUE
    reference <- tna::build_model(wide, type = tna_types[[method]], params = tna_params)
    fit <- build_network(wide, method = method, params = params)
    label <- sprintf("%s, %s", method, paste(names(params), unlist(params), sep = " ", collapse = ", "))
    if (isTRUE(params$weighted)) {
      expect_known_difference(extract_transition_matrix(fit), unclass(reference$weights), reference = "tna",
                              covers = "build_network()", quantity = label, tolerance = sw_tolerance,
                              reason = "tna's windowed counter (compute_transitions_windowed()) ignores `weighted`; Nestimate weights each windowed pair by 1 / sequence length, as in its unwindowed count.",
                              agrees_with = "weighted = FALSE")
    } else {
      expect_reference(extract_transition_matrix(fit), unclass(reference$weights), reference = "tna",
                       covers = "build_network()", quantity = label, tolerance = sw_tolerance)
    }
  }, cases$method, cases$window))
})

test_that("windows refused or ignored where the documentation says so", {
  wide <- sweep_inputs$ragged
  invisible(lapply(c("attention", "ngram", "gap", "reverse"), \(method) {
    expect_error(build_network(wide, method = method, params = list(window_size = 2L)), class = "nestimate_bad_input")
  }))
  expect_warning(fit <- build_network(wide, method = "frequency", params = list(window_size = 2L),
                                      mode = "overlapping"), class = "nestimate_ignored_mode")
  expect_reference(extract_transition_matrix(fit),
                   extract_transition_matrix(build_network(wide, method = "frequency", params = list(window_size = 2L))),
                   reference = "independent: ?build_network, overlapping mode ignored on sequence windows",
                   covers = "build_network()", quantity = "overlapping mode on a sequence window", tolerance = 0)
  expect_reference(as.numeric("mode" %in% names(fit$params)), 0,
                   reference = "independent: ?build_network, the ignored mode is dropped from params",
                   covers = "build_network()", quantity = "mode left out of params", tolerance = 0)
})

test_that("a grouped network gives each group the network of its rows on the shared alphabet", {
  wide <- tn_simulate(n_actors = 36L, n_states = 4L, seq_length = 8L, seed = 806L, ragged = TRUE)
  # group g3 never visits D: the shared alphabet keeps D as a row of zeros.
  grp <- rep(c("g1", "g2", "g3"), length.out = nrow(wide))
  wide[grp == "g3", ] <- lapply(wide[grp == "g3", ], \(column) ifelse(column == "D", "A", column))
  grouped <- cbind(wide, grp = grp, stringsAsFactors = FALSE)
  configs <- method_grid[c(1L, 2L, 4L, 6L, 10L, 13L, 16L)]
  cases <- expand.grid(config = seq_along(configs), boundary = c(1L, 9L), scaling = c("none", "max"),
                       stringsAsFactors = FALSE)
  invisible(Map(\(config_index, boundary, scaling) {
    config <- configs[[config_index]]
    start <- boundary_grid$start[[boundary]]
    end <- boundary_grid$end[[boundary]]
    scale <- if (identical(scaling, "none")) NULL else scaling
    fit <- build_network(grouped, method = config$method, group = "grp", params = config$params,
                         start = sw_flag(start), end = sw_flag(end), scaling = scale)
    all_states <- sw_states(sw_positions(wide, sw_label(start, "Start"), sw_label(end, "End")),
                            sw_label(start, "Start"), sw_label(end, "End"))
    lapply(c("g1", "g2", "g3"), \(g) {
      mat <- sw_positions(wide[grp == g, , drop = FALSE], sw_label(start, "Start"), sw_label(end, "End"))
      expected <- sw_reference(mat, all_states, config$method, config$params)
      if (!is.null(scale)) expected <- sw_scale(expected, scale)
      # A grouped result has no accessor for one group's weights: `[[` (API gap, as in test-build_network.R).
      expect_reference(extract_transition_matrix(fit[[g]]), expected,
                       reference = "independent: pair counts of the group's rows on the shared alphabet, base R",
                       covers = "build_network()",
                       quantity = sprintf("%s, group %s, start %s, end %s, scaling %s", config$label, g, start, end,
                                          scaling),
                       tolerance = sw_tolerance)
    })
  }, cases$config, cases$boundary, cases$scaling))
})

test_that("every method alias gives the network of its canonical method", {
  wide <- sweep_inputs$ragged
  aliases <- c(tna = "relative", transition = "relative", ftna = "frequency", counts = "frequency",
               cna = "co_occurrence", wcna = "co_occurrence", `co-occurrence` = "co_occurrence", atna = "attention",
               `n-gram` = "ngram", n_gram = "ngram")
  invisible(Map(\(alias, canonical) {
    expected <- sw_expected(wide, list(method = canonical, params = list()))$weights
    fit <- build_network(wide, method = alias)
    expect_reference(extract_transition_matrix(fit), expected,
                     reference = "independent: pair counts of the canonical method, base R",
                     covers = "build_network()", quantity = sprintf("alias %s", alias), tolerance = sw_tolerance)
    expect_reference(as.numeric(identical(fit$method, canonical)), 1,
                     reference = "independent: ?build_network, aliases resolve to the canonical name",
                     covers = "build_network()", quantity = sprintf("method recorded for alias %s", alias), tolerance = 0)
  }, names(aliases), aliases))
})

test_that("the transition wrappers pass every argument, positional start and end included, to build_network()", {
  wide <- sweep_inputs$ragged
  wrappers <- list(build_tna = list(fn = build_tna, method = "relative"),
                   build_ftna = list(fn = build_ftna, method = "frequency"),
                   build_atna = list(fn = build_atna, method = "attention"),
                   build_cna = list(fn = build_cna, method = "co_occurrence"))
  cases <- expand.grid(wrapper = names(wrappers), boundary = seq_len(nrow(boundary_grid)), stringsAsFactors = FALSE)
  invisible(Map(\(wrapper, boundary) {
    w <- wrappers[[wrapper]]
    start <- boundary_grid$start[[boundary]]
    end <- boundary_grid$end[[boundary]]
    expected <- sw_expected(wide, list(method = w$method, params = list()), sw_label(start, "Start"),
                            sw_label(end, "End"))$weights
    # start and end are the second and third formals of each wrapper: passed by position.
    positional <- w$fn(wide, sw_flag(start), sw_flag(end))
    expect_reference(extract_transition_matrix(positional), expected,
                     reference = "independent: pair counts of the positions with the boundary states, base R",
                     covers = sprintf("%s()", wrapper),
                     quantity = sprintf("positional start %s, end %s", start, end), tolerance = sw_tolerance)
  }, cases$wrapper, cases$boundary))
  settings <- list(list(scaling = "max", threshold = 0.1), list(scaling = c("rank", "minmax"), threshold = 0),
                   list(scaling = NULL, threshold = 0.2))
  cases <- expand.grid(wrapper = names(wrappers), setting = seq_along(settings), stringsAsFactors = FALSE)
  invisible(Map(\(wrapper, setting_index) {
    w <- wrappers[[wrapper]]
    s <- settings[[setting_index]]
    params <- if (identical(w$method, "attention")) list(lambda = 3) else list(weighted = TRUE)
    expected <- sw_threshold(sw_scale(sw_expected(wide, list(method = w$method, params = params))$weights,
                                      s$scaling %||% character()), s$threshold)
    fit <- w$fn(wide, scaling = s$scaling, threshold = s$threshold, params = params)
    expect_reference(extract_transition_matrix(fit), expected,
                     reference = "independent: ?build_network weights, scaling then threshold, base R",
                     covers = sprintf("%s()", wrapper),
                     quantity = sprintf("scaling %s, threshold %s, params %s", paste(s$scaling, collapse = "+"),
                                        s$threshold, names(params)), tolerance = sw_tolerance)
  }, cases$wrapper, cases$setting))
})

test_that("start and end are refused outside the transition methods and for malformed labels", {
  wide <- sweep_inputs$ragged
  numeric_data <- data.frame(x = c(1, 3, 2, 5, 4, 6), y = c(2, 1, 4, 3, 6, 5), z = c(5, 3, 4, 1, 2, 0))
  invisible(lapply(c("cor", "pcor", "glasso"), \(method) {
    expect_error(build_network(numeric_data, method = method, start = TRUE), class = "nestimate_bad_input")
    expect_error(build_network(numeric_data, method = method, end = "Stop"), class = "nestimate_bad_input")
  }))
  invisible(lapply(list(NA, c("a", "b"), "", 1), \(bad) {
    expect_error(build_network(wide, method = "relative", start = bad), class = "nestimate_bad_input")
    expect_error(build_network(wide, method = "relative", end = bad), class = "nestimate_bad_input")
  }))
  expect_error(build_network(wide, method = "relative", latency_floor = 0.2), class = "nestimate_bad_input")
  expect_error(build_network(wide, method = "frequency", level = "between", actor = "T1"),
               class = "nestimate_bad_input")
  expect_error(build_network(wide, method = "relative", scaling = "zscore"), class = "nestimate_bad_input")
  expect_error(build_network(wide, method = "co_occurrence", params = list(type = "rank")),
               class = "nestimate_bad_input")
})

test_that("every default of build_network() and of the sequence estimators equals the value written out", {
  wide <- sweep_inputs$ragged
  # ?build_network gives these defaults; the estimators' own (lambda 1, direction forward, n_gram 2, max_gap 1,
  # type frequency, weighted FALSE, concat 1) are stated in the Estimators section and in tna.
  explicit <- list(
    relative = list(params = list(weighted = FALSE, concat = 1L)),
    frequency = list(params = list(weighted = FALSE, concat = 1L)),
    co_occurrence = list(params = list(type = "frequency", weighted = FALSE)),
    attention = list(params = list(lambda = 1, direction = "forward")),
    ngram = list(params = list(n_gram = 2L)),
    gap = list(params = list(max_gap = 1L)),
    reverse = list(params = list(weighted = FALSE))
  )
  invisible(Map(\(method, spec) {
    omitted <- build_network(wide, method = method)
    written <- build_network(wide, method = method, start = FALSE, end = FALSE, scaling = NULL, threshold = 0,
                             format = "auto", time_threshold = 900, timezone = "UTC", params = spec$params)
    expect_reference(extract_transition_matrix(written), extract_transition_matrix(omitted),
                     reference = "independent: the documented defaults written out", covers = "build_network()",
                     quantity = sprintf("defaults, %s", method), tolerance = 0)
    if (!identical(method, "co_occurrence")) {
      expect_reference(extract_initial_probs(written), extract_initial_probs(omitted),
                       reference = "independent: the documented defaults written out", covers = "build_network()",
                       quantity = sprintf("defaults, initial probabilities, %s", method), tolerance = 0)
    }
  }, names(explicit), explicit))
  # The wrappers' defaults are build_network()'s.
  invisible(Map(\(fn, method) {
    expect_reference(extract_transition_matrix(fn(wide)), extract_transition_matrix(build_network(wide, method = method)),
                     reference = "independent: the wrapper is build_network() with its method", covers = "build_network()",
                     quantity = sprintf("wrapper defaults, %s", method), tolerance = 0)
  }, list(build_tna, build_ftna, build_atna, build_cna), c("relative", "frequency", "attention", "co_occurrence")))
})

test_that("one-hot data at every window size, mode and grouping gives the networks of tna::import_onehot()", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  codes <- c("A", "B", "C", "D")
  set.seed(807)
  n <- 48L
  onehot <- as.data.frame(matrix(stats::rbinom(4L * n, 1L, 0.4), n, dimnames = list(NULL, codes)))
  onehot$actor <- rep(c("x", "y", "z"), each = n / 3L)
  onehot$session <- rep(rep(c("s1", "s2"), each = n / 6L), times = 3L)
  cases <- expand.grid(window = c(2L, 3L, 4L), mode = c("non-overlapping", "overlapping"),
                       grouping = c("actor", "actor+session"), stringsAsFactors = FALSE)
  # tna's sliding rule differs from the tumbling one only for windows wider than 2 in tna < 1.5.5; tna 1.5.7 is
  # the reference here.
  invisible(Map(\(window, mode, grouping) {
    session <- if (identical(grouping, "actor+session")) "session" else NULL
    imported <- do.call(tna::import_onehot, c(list(onehot, cols = codes, actor = "actor", window_size = window,
                                                   window_type = if (identical(mode, "overlapping")) "sliding" else "tumbling"),
                                              if (!is.null(session)) list(session = session)))
    fit <- \(method, params = list()) build_network(onehot, method = method, codes = codes, actor = "actor",
                                                    session = session, window_size = window, mode = mode,
                                                    params = params)
    label <- sprintf("window %d, %s, %s", window, mode, grouping)
    cooccurrence <- unclass(tna::build_model(imported, type = "co-occurrence", scaling = character(0))$weights)
    expect_reference(extract_transition_matrix(fit("frequency")), unclass(tna::ftna(imported)$weights),
                     reference = "tna", covers = "build_network()", quantity = sprintf("one-hot frequency, %s", label),
                     tolerance = sw_tolerance)
    expect_reference(extract_transition_matrix(fit("relative")), unclass(tna::tna(imported)$weights),
                     reference = "tna", covers = "build_network()", quantity = sprintf("one-hot relative, %s", label),
                     tolerance = sw_tolerance)
    expect_reference(extract_transition_matrix(fit("co_occurrence")), cooccurrence,
                     reference = "tna", covers = "build_network()", quantity = sprintf("one-hot co-occurrence, %s", label),
                     tolerance = sw_tolerance)
    expect_reference(extract_transition_matrix(fit("co_occurrence", list(type = "relative"))),
                     sw_row_normalise(cooccurrence), reference = "tna", covers = "build_network()",
                     quantity = sprintf("one-hot co-occurrence, type relative (tna's counts row-normalised), %s", label),
                     tolerance = sw_tolerance)
  }, cases$window, cases$mode, cases$grouping))
  expect_error(build_network(onehot, method = "attention", codes = codes, actor = "actor"), class = "nestimate_bad_input")
  # ?build_network: start and end are for sequences; the one-hot route (.estimator_wtna_core()) never reads
  # begin_state / end_state. Either a refusal or a start node would honour the request.
  honoured <- vapply(c("relative", "frequency", "co_occurrence"), \(method) {
    withCallingHandlers(
      tryCatch("Start" %in% rownames(build_network(onehot, method = method, codes = codes, actor = "actor",
                                                   window_size = 2L, start = TRUE)$weights),
               nestimate_bad_input = \(e) TRUE),
      nestimate_sparse_onehot = \(w) invokeRestart("muffleWarning"))
  }, logical(1))
  expect_reference(as.numeric(honoured), rep(1, 3L),
                   reference = "independent: ?build_network, start adds a Start state or is refused",
                   covers = "build_network()", quantity = "start = TRUE on one-hot data (relative, frequency, co_occurrence)",
                   tolerance = 0, key = c("relative", "frequency", "co_occurrence"))
})
