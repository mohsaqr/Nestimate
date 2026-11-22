# prepare(), extract_sequences(), build_network() on a long event log, long_to_wide() and wide_to_long():
# sequences cut from an event log by actor, session and a time gap, against tna::prepare_data() (open-source; the
# owner's co-authored package, to be confirmed), the networks built from them against tna::tna() and tna::ftna(), and
# the reshaping against base R. Re-expressed from local_testing_and_equivalence/test-equiv-tna-long.R (151
# combinations of the bundled long datasets; the full grid stays there as a standing study) and
# test-equiv-prepare-grouping-tna.R.
# A prepared object (nestimate_data) has no accessor for its session table beyond extract_sequences(meta = TRUE),
# which is used here; tna's prepared object is read with `$` (its API).

tolerance <- 1e-10
regulation <- subset(group_regulation_long, Actor <= 150)

# tna's prepared sequences as one row per session: the session id and T1..Tn.
tn_tna_sequences <- function(prepared) {
  sequences <- as.data.frame(prepared$sequence_data)
  names(sequences) <- sub("^Action_", "", names(sequences))
  cbind(.session_id = as.data.frame(prepared$meta_data)$.session_id, sequences)
}

test_that("prepare() cuts the same sessions as tna::prepare_data()", {
  skip_without_reference("tna")
  prepared <- prepare(regulation, actor = "Actor", action = "Action", time = "Time")
  reference <- suppressMessages(tna::prepare_data(regulation, actor = "Actor", action = "Action", time = "Time"))
  found <- extract_sequences(prepared, meta = TRUE)
  expected <- tn_tna_sequences(reference)
  expect_reference(found, expected, reference = "tna", covers = c("prepare()", "extract_sequences()"),
                   quantity = "sessions of group_regulation_long, actors 1-150", tolerance = 0, key = ".session_id")
})

test_that("a long event log gives the networks of tna on the sequences of tna::prepare_data()", {
  skip_without_reference("tna")
  reference <- suppressMessages(tna::prepare_data(regulation, actor = "Actor", action = "Action", time = "Time"))
  expect_reference(extract_transition_matrix(build_network(regulation, method = "relative", action = "Action",
                                                           actor = "Actor", time = "Time")),
                   unclass(tna::tna(reference)$weights), reference = "tna", covers = "build_network()",
                   quantity = "relative weights, group_regulation_long actor and time", tolerance = tolerance)
  expect_reference(extract_transition_matrix(build_network(regulation, method = "frequency", action = "Action",
                                                           actor = "Actor", time = "Time")),
                   unclass(tna::ftna(reference)$weights), reference = "tna", covers = "build_network()",
                   quantity = "frequency weights, group_regulation_long actor and time", tolerance = tolerance)
  actor_only <- suppressMessages(tna::prepare_data(regulation, actor = "Actor", action = "Action"))
  expect_reference(extract_transition_matrix(build_network(regulation, method = "relative", action = "Action",
                                                           actor = "Actor")),
                   unclass(tna::tna(actor_only)$weights), reference = "tna", covers = "build_network()",
                   quantity = "relative weights, group_regulation_long actor only", tolerance = tolerance)
})

test_that("an event log with projects, sessions and timestamps gives the networks of tna", {
  skip_without_reference("tna")
  invisible(lapply(c("code", "cluster"), \(action) {
    reference <- suppressMessages(tna::prepare_data(human_long, actor = c("project", "session_id"), action = action,
                                                    time = "timestamp"))
    expect_reference(extract_transition_matrix(build_network(human_long, method = "relative", action = action,
                                                             actor = "project", session = "session_id",
                                                             time = "timestamp")),
                     unclass(tna::tna(reference)$weights), reference = "tna", covers = "build_network()",
                     quantity = sprintf("relative weights, human_long %s", action), tolerance = tolerance)
    expect_reference(extract_transition_matrix(build_network(human_long, method = "frequency", action = action,
                                                             actor = "project", session = "session_id",
                                                             time = "timestamp")),
                     unclass(tna::ftna(reference)$weights), reference = "tna", covers = "build_network()",
                     quantity = sprintf("frequency weights, human_long %s", action), tolerance = tolerance)
  }))
  reference <- suppressMessages(tna::prepare_data(ai_long, actor = c("project", "session_id"), action = "code",
                                                  time = "timestamp"))
  expect_reference(extract_transition_matrix(build_network(ai_long, method = "relative", action = "code",
                                                           actor = "project", session = "session_id",
                                                           time = "timestamp")),
                   unclass(tna::tna(reference)$weights), reference = "tna", covers = "build_network()",
                   quantity = "relative weights, ai_long code", tolerance = tolerance)
})

test_that("long_to_wide() orders each actor's actions by time, as split() and order() do", {
  wide <- long_to_wide(regulation, id_col = "Actor", time_col = "Time", action_col = "Action")
  ordered <- regulation[order(regulation$Actor, regulation$Time), ]
  sequences <- split(ordered$Action, ordered$Actor)
  width <- max(lengths(sequences))
  expected <- as.data.frame(do.call(rbind, lapply(sequences, \(s) c(s, rep(NA_character_, width - length(s))))))
  names(expected) <- paste0("V", seq_len(width))
  expected <- cbind(Actor = as.integer(names(sequences)), expected)
  expect_reference(wide, expected, reference = "independent: split() of the actions by actor in time order, base R",
                   covers = "long_to_wide()", quantity = "wide sequences of group_regulation_long, actors 1-150",
                   tolerance = 0, key = "Actor")
})

test_that("wide_to_long() stacks the non-missing cells with their position, as base R does", {
  wide <- tn_simulate(n_actors = 12L, n_states = 4L, seq_length = 10L, seed = 201L, ragged = TRUE)
  names(wide) <- paste0("V", seq_len(ncol(wide)))
  long <- wide_to_long(wide)
  sequences <- tn_sequences(wide)
  expected <- do.call(rbind, lapply(seq_along(sequences), \(i) {
    data.frame(id = i, Time = seq_along(sequences[[i]]), Action = sequences[[i]])
  }))
  expect_reference(long, expected, reference = "independent: the non-missing cells of each row with their position, base R",
                   covers = "wide_to_long()", quantity = "long table of 12 ragged sequences", tolerance = 0,
                   key = c("id", "Time"))
})

test_that("prepare_for_tna() reshapes a long log into each actor's actions in time order", {
  log <- subset(regulation, Actor <= 40, select = c(Actor, Time, Action))
  ordered <- log[order(log$Actor, log$Time), ]
  sequences <- split(ordered$Action, ordered$Actor)
  width <- max(lengths(sequences))
  expected <- as.data.frame(do.call(rbind, lapply(sequences, \(s) c(s, rep(NA_character_, width - length(s))))))
  names(expected) <- paste0("V", seq_len(width))
  found <- cbind(row = seq_len(length(sequences)), prepare_for_tna(log, type = "long"))
  expect_reference(found, cbind(row = seq_len(length(sequences)), expected),
                   reference = "independent: split() of the actions by actor in time order, base R",
                   covers = "prepare_for_tna()", quantity = "wide sequences of actors 1-40", tolerance = 0, key = "row")
})

test_that("long_to_wide() and prepare_for_tna() give the sequences of tna::prepare_data() ordered by time", {
  skip_without_reference("tna")
  reference <- suppressMessages(tna::prepare_data(regulation, actor = "Actor", action = "Action", order = "Time"))
  sequences <- as.data.frame(reference$sequence_data)
  expected <- cbind(Actor = as.integer(as.data.frame(reference$meta_data)$Actor),
                    stats::setNames(sequences, paste0("V", seq_len(ncol(sequences)))))
  expect_reference(long_to_wide(regulation, id_col = "Actor", time_col = "Time", action_col = "Action"), expected,
                   reference = "tna", covers = "long_to_wide()",
                   quantity = "wide sequences of group_regulation_long, actors 1-150, tna", tolerance = 0, key = "Actor")
  log <- subset(regulation, Actor <= 40, select = c(Actor, Time, Action))
  reference <- suppressMessages(tna::prepare_data(log, actor = "Actor", action = "Action", order = "Time"))
  sequences <- as.data.frame(reference$sequence_data)
  expect_reference(cbind(row = seq_len(nrow(sequences)), prepare_for_tna(log, type = "long")),
                   cbind(row = seq_len(nrow(sequences)), stats::setNames(sequences, paste0("V", seq_len(ncol(sequences))))),
                   reference = "tna", covers = "prepare_for_tna()", quantity = "wide sequences of actors 1-40, tna",
                   tolerance = 0, key = "row")
})

test_that("wide_to_long() equals stats::reshape() of the wide table without its missing cells", {
  invisible(lapply(c(201L, 202L), \(seed) {
    wide <- tn_simulate(n_actors = 12L, n_states = 4L, seq_length = 10L, seed = seed, ragged = TRUE)
    names(wide) <- paste0("V", seq_len(ncol(wide)))
    stacked <- stats::reshape(cbind(id = seq_len(nrow(wide)), wide), direction = "long", varying = names(wide),
                              v.names = "Action", timevar = "Time", idvar = "id")
    expected <- subset(stacked, !is.na(Action), select = c(id, Time, Action))
    expect_reference(wide_to_long(wide), expected, reference = "stats", covers = "wide_to_long()",
                     quantity = sprintf("long table of 12 ragged sequences, seed %d, stats reshape()", seed),
                     tolerance = 0, key = c("id", "Time"))
  }))
})
