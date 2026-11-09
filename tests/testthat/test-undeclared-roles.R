# Names give no role: a column called `action` that nobody declared is refused by
# every verb that would otherwise read the log as positions, and a column declared
# under any role is never called unused.

.undeclared_log <- function() {
  data.frame(
    Actor = rep(c("p1", "p2"), each = 4L),
    Time = rep(1:4, 2L),
    Action = c("A", "B", "A", "C", "B", "B", "C", "A")
  )
}

test_that("every verb refuses an undeclared action column", {
  log <- .undeclared_log()
  expect_error(build_network(log, method = "frequency"), class = "nestimate_bad_input")
  expect_error(frequencies(log), class = "nestimate_bad_input")
  expect_error(frequencies(log, id = "Actor"), class = "nestimate_bad_input")
  expect_error(magnitude_difference(log), class = "nestimate_bad_input")
  expect_error(magnitude_difference(log, actor = "Actor"), class = "nestimate_bad_input")
  expect_error(convert_sequence_format(log), class = "nestimate_bad_input")
  expect_error(extract_sequences(log), class = "nestimate_bad_input")
  expect_error(sequence_patterns(log), class = "nestimate_bad_input")
  expect_error(sequence_metrics(log), class = "nestimate_bad_input")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(sequence_plot(log), class = "nestimate_bad_input")
})

test_that("declaring the action gives the sequences, counted by hand", {
  log <- .undeclared_log()
  # p1: A B A C; p2: B B C A
  counts <- unclass(frequencies(log, action = "Action", id = "Actor", time = "Time"))
  expect_equal(counts[["A", "B"]], 1)
  expect_equal(counts[["B", "A"]], 1)
  expect_equal(counts[["A", "C"]], 1)
  expect_equal(counts[["B", "B"]], 1)
  expect_equal(counts[["B", "C"]], 1)
  expect_equal(counts[["C", "A"]], 1)
  expect_equal(sum(counts), 6)
  expect_identical(nrow(extract_sequences(log, action = "Action", actor = "Actor")), 2L)
  expect_silent(magnitude_difference(log, actor = "Actor", action = "Action",
                                     time = "Time"))
})

test_that("a wide table may hold positions, and one-hot a column called action", {
  wide <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"), T3 = c("C", "A"))
  expect_equal(sum(unclass(frequencies(wide))), 4)
  explicit <- data.frame(action = c("A", "B"), T2 = c("B", "C"))
  expect_equal(sum(unclass(frequencies(explicit, cols = c("action", "T2")))), 2)
  # 0/1 columns with nothing declared could be either reading
  expect_error(.refuse_undeclared_roles(data.frame(action = c(0, 1), b = c(1, 0)), "verb"),
               class = "nestimate_bad_input")
})

test_that("only an undeclared time or session column is called unused", {
  log <- .undeclared_log()
  log$session_id <- rep(c("s1", "s2"), 4L)
  says_unused <- function(...) {
    said <- character()
    withCallingHandlers(
      build_network(log, method = "frequency", action = "Action", ...),
      nestimate_column_role = function(m) {
        said <<- c(said, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      message = function(m) invokeRestart("muffleMessage")
    )
    said
  }
  expect_match(says_unused(actor = "Actor"), "`Time`, `session_id` not used")
  # declared under another role, the column is in use
  expect_match(says_unused(actor = "session_id"), "`Time` not used")
  expect_false(any(grepl("session_id", says_unused(actor = "session_id"))))
  expect_false(any(grepl("Time", says_unused(actor = "Actor", order = "Time"))))
  expect_length(says_unused(actor = "session_id", time = "Time"), 0L)
})

test_that("a public mode = 'overlapping' warns on sequence windows too", {
  d <- data.frame(T1 = c("A", "A"), T2 = "B", T3 = "C", T4 = "D")
  expect_warning(
    build_network(d, method = "frequency", mode = "overlapping",
                  params = list(window_size = 2L)),
    class = "nestimate_ignored_mode"
  )
})

test_that("a boot_glasso draw without a correlation counts against the certainty", {
  # 19 of 20 reach it with the NA in the denominator, .95 > .95 is false
  cors <- list(c(rep(1, 19L), NA))
  expect_equal(.bg_cs_coefficient(cors, 0.5, 0.7, 0.95), 0)
  expect_equal(.bg_cs_coefficient(list(rep(1, 20L)), 0.5, 0.7, 0.95), 0.5)
})

test_that("a frame with no role declared must be readable, on every sequence verb", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  log <- head(human_long, 20L)
  verbs <- list(
    function(d) build_network(d, method = "frequency"),
    function(d) frequencies(d),
    function(d) magnitude_difference(d),
    function(d) convert_sequence_format(d),
    function(d) extract_sequences(d),
    function(d) sequence_patterns(d),
    function(d) sequence_metrics(d),
    function(d) stability_landscape(d, iter = 2L),
    function(d) sequence_plot(d),
    function(d) distribution_plot(d),
    function(d) build_mcml(d, list(G = "A")),
    function(d) build_hon(d),
    function(d) build_hypa(d),
    function(d) build_mogen(d),
    function(d) markov_order_test(d, n_perm = 2L),
    function(d) path_counts(d),
    function(d) path_dependence(d),
    function(d) state_frequencies(d),
    function(d) sequence_compare(d, group = rep(c("a", "b"), 10L)),
    function(d) association_rules(d),
    function(d) cooccurrence(d),
    function(d) actor_endpoints(d),
    function(d) mark_first_state(d),
    function(d) mark_terminal_state(d),
    function(d) build_clusters(d, k = 2L),
    function(d) build_mmm(d, k = 2L, n_starts = 1L)
  )
  invisible(lapply(verbs, function(verb) {
    expect_error(verb(log), class = "nestimate_bad_input")
  }))
})

test_that("a declared wide table reads as the same table with numbered names", {
  wide <- data.frame(week_a = c("A", "B", "A"), week_b = c("B", "C", "C"),
                     week_c = c("C", "A", "B"))
  numbered <- wide
  names(numbered) <- paste0("T", 1:3)
  # A->B, B->C, B->C, C->A, A->C, C->B: six transitions
  declared <- build_network(wide, method = "frequency", state_cols = names(wide))
  expect_equal(sum(declared$weights), 6)
  expect_equal(unname(declared$weights),
               unname(build_network(numbered, method = "frequency")$weights))
  expect_equal(sum(unclass(frequencies(wide, cols = names(wide)))), 6)
  expect_identical(nrow(extract_sequences(wide, cols = names(wide))), 3L)
  # the network is the way to declare for verbs with no `action`
  expect_equal(build_hon(declared, max_order = 2L)$n_edges,
               build_hon(numbered, max_order = 2L)$n_edges)
  # resampling and clustering keep working on what was declared
  expect_s3_class(bootstrap_network(declared, iter = 3L, seed = 1L), "net_bootstrap")
  expect_s3_class(build_clusters(declared, k = 2L), "net_clustering")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(sequence_plot(wide, cols = names(wide)))
  expect_no_error(sequence_plot(wide, cols = names(wide), type = "distribution"))
  expect_no_error(distribution_plot(wide, cols = names(wide)))
})

test_that("verbs that read every column refuse a run with something beside it", {
  beside <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), T3 = c("A", "A"),
                       outcome = c("pass", "fail"))
  expect_error(frequencies(beside), class = "nestimate_bad_input")
  expect_error(build_hon(beside), class = "nestimate_bad_input")
  expect_error(path_counts(beside), class = "nestimate_bad_input")
  expect_error(build_mcml(beside, list(G = c("A", "B"))),
               class = "nestimate_bad_input")
  # build_network() sets it aside, and says so
  expect_message(net <- build_network(beside, method = "frequency"),
                 class = "nestimate_column_role")
  expect_identical(sort(rownames(net$weights)), c("A", "B"))
})

test_that("a role naming a column that is not there is refused", {
  log <- .undeclared_log()
  expect_error(frequencies(log, action = "absent"), class = "nestimate_bad_input")
  expect_error(frequencies(log, action = "Action", id = "absent"),
               class = "nestimate_bad_input")
  expect_error(magnitude_difference(log, action = "absent"),
               class = "nestimate_bad_input")
  expect_error(convert_sequence_format(log, action = "absent"),
               class = "nestimate_bad_input")
  expect_error(extract_sequences(log, action = "Action", cols = "absent"),
               class = "nestimate_bad_input")
  expect_error(build_mcml(log, list(G = "A"), from = "absent", to = "Action"),
               class = "nestimate_bad_input")
})

test_that("association methods are not refused for a column called action", {
  scores <- data.frame(action = c(1, 3, 2, 5, 4, 6, 8, 7), b = c(2, 1, 4, 3, 6, 5, 7, 8),
                       c = c(1, 2, 2, 4, 3, 5, 5, 9))
  expect_identical(dim(build_network(scores, method = "cor")$weights), c(3L, 3L))
  # sequence methods still ask
  expect_error(build_network(scores, method = "frequency"),
               class = "nestimate_bad_input")
})

test_that("a column declared as the time is never called unused", {
  log <- .undeclared_log()
  log$session_id <- rep(1:4, 2L)
  said <- character()
  withCallingHandlers(
    build_network(log, method = "frequency", action = "Action", actor = "Actor",
                  time = "session_id"),
    nestimate_column_role = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_false(any(grepl("session_id", said, fixed = TRUE)))
})

test_that("the one-sequence note follows what was built", {
  log <- data.frame(session_id = c("s1", "s1", "s2", "s2"), at = 1:4,
                    code = c("A", "B", "B", "A"))
  notes <- character()
  net <- withCallingHandlers(
    build_network(log, method = "frequency", action = "code",
                  session = "session_id", order = "at"),
    message = function(m) {
      notes <<- c(notes, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_false(any(grepl("one long sequence", notes, fixed = TRUE)))
  expect_identical(nrow(net$data), 2L)
  # no session and no actor: one sequence, and the note is true
  expect_message(build_network(log, method = "frequency", action = "code"),
                 "one long sequence")
})

test_that("one column has nothing beside it to refuse", {
  counted <- state_frequencies(data.frame(code = c("A", "B", "A")))
  expect_equal(counted$count, c(2L, 1L))
  expect_error(state_frequencies(data.frame(action = c("A", "B"))),
               class = "nestimate_bad_input")
})

test_that("positions a verb returned read back in, whatever they are called", {
  months <- data.frame(Sep.93 = c("A", "B"), Oct.93 = c("B", "A"), Nov.93 = c("A", "A"))
  expect_error(frequencies(months), class = "nestimate_bad_input")
  # A->B, B->A, B->A, A->A: four transitions
  routes <- list(extract_sequences(months, cols = names(months)),
                 mark_first_state(months, cols = names(months)),
                 mark_terminal_state(months, cols = names(months)))
  invisible(lapply(routes, function(read) {
    expect_equal(sum(unclass(frequencies(read))), 4)
    expect_equal(sum(build_network(read, method = "frequency")$weights), 4)
    expect_identical(nrow(sequence_metrics(read)), 2L)
    expect_s3_class(build_hon(read, max_order = 2L), "net_hon")
  }))
  # keeping only some columns is a new frame: declared again
  expect_error(frequencies(subset(routes[[1L]], select = c(Sep.93, Oct.93))),
               class = "nestimate_bad_input")
})

test_that("declared positions keep their roles when they are read back in", {
  named <- data.frame(from = c("A", "B"), to = c("B", "C"), T3 = c("C", "A"),
                      T4 = c("A", "B"))
  # A->B B->C C->A, B->C C->A A->B: six transitions
  marked <- extract_sequences(named, cols = names(named))
  expect_equal(sum(build_network(marked, method = "frequency")$weights), 6)
  expect_identical(dim(extract_sequences(marked)), c(2L, 4L))
  clusters <- list(G = c("A", "B"), H = "C")
  expect_equal(sum(build_mcml(marked, clusters, type = "raw")$macro$weights), 6)
  expect_equal(sum(build_mcml(mark_first_state(named, cols = names(named)), clusters,
                              type = "raw")$macro$weights), 6)
})

test_that("a column added to resolved positions is not trusted", {
  weeks <- data.frame(week_a = c("A", "B"), week_b = c("B", "A"), week_c = c("A", "A"))
  marked <- extract_sequences(weeks, cols = names(weeks))
  marked$grade <- c("pass", "fail")
  expect_error(frequencies(marked), class = "nestimate_bad_input")
  expect_error(build_network(marked, method = "frequency"),
               class = "nestimate_bad_input")
})

test_that("a numeric covariate beside numeric positions is set aside and named", {
  coded <- data.frame(T1 = c(1, 2, 1), T2 = c(2, 3, 3), T3 = c(3, 1, 2), grade = c(77, 88, 99))
  said <- character()
  net <- withCallingHandlers(
    build_network(coded, method = "frequency"),
    nestimate_column_role = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  # 1->2 2->3, 2->3 3->1, 1->3 3->2: six transitions among three states
  expect_identical(sort(rownames(net$weights)), c("1", "2", "3"))
  expect_equal(sum(net$weights), 6)
  expect_match(said, "`grade`", fixed = TRUE)
  # one-hot columns stay codes
  hot <- data.frame(code1 = c(1, 0), code2 = c(0, 1), other = c(1, 1))
  expect_no_error(build_network(hot, method = "co_occurrence", format = "onehot"))
})

test_that("each refusal's advice works when it is followed", {
  weeks <- data.frame(week_a = c("A", "B"), week_b = c("B", "A"), week_c = c("A", "A"))
  # metadata_cols, as the message says
  with_cohort <- transform(weeks, cohort = "a")
  expect_error(build_network(with_cohort, method = "relative"),
               class = "nestimate_bad_input")
  expect_identical(build_network(with_cohort, method = "relative",
                                 metadata_cols = "cohort")$n_nodes, 2L)
  # the network, as the message says, on every verb that advises it
  net <- build_network(weeks, method = "relative", state_cols = names(weeks))
  freq_net <- build_network(weeks, method = "frequency", state_cols = names(weeks))
  expect_s3_class(entropy_bayes(net, draws = 100L), "net_entropy_bayes")
  expect_s3_class(transition_entropy(net), "net_transition_entropy")
  expect_s3_class(build_hon(net, max_order = 2L), "net_hon")
  expect_s3_class(build_mogen(net, max_order = 1L), "net_mogen")
  expect_s3_class(build_hypa(net), "net_hypa")
  expect_s3_class(path_counts(net), "data.frame")
  expect_s3_class(state_frequencies(net), "data.frame")
  expect_s3_class(entropy_bayes(freq_net, draws = 100L), "net_entropy_bayes")
  # `cols`, as the message says, on verbs that have no `action`
  expect_s3_class(actor_endpoints(weeks, cols = names(weeks)), "data.frame")
  expect_no_error(stability_landscape(weeks, iter = 2L, cols = names(weeks)))
})

test_that("cooccurrence() refuses a role it would ignore, and its advice works", {
  d <- data.frame(T1 = c("a", "b"), T2 = c("b", "c"))
  expect_error(cooccurrence(d, field = "TYPO"), class = "nestimate_bad_input")
  expect_error(cooccurrence(d, by = "TYPO"), class = "nestimate_bad_input")
  expect_error(cooccurrence(d, field = "T1"), class = "nestimate_bad_input")
  expect_error(cooccurrence(d, by = "T1"), class = "nestimate_bad_input")
  long <- data.frame(item = c("a", "b", "a", "c"), doc = c(1, 1, 2, 2))
  expect_error(cooccurrence(long), class = "nestimate_bad_input")
  # a and b together once, a and c together once
  co <- cooccurrence(long, field = "item", by = "doc")
  expect_equal(co$weights[["a", "b"]], 1)
  expect_equal(co$weights[["a", "c"]], 1)
  expect_equal(co$weights[["b", "c"]], 0)
  delim <- data.frame(tags = c("a;b", "a;c"))
  expect_s3_class(cooccurrence(delim, field = "tags", sep = ";"), "netobject")
})

test_that("the unused-role message says what holds for the call", {
  log <- data.frame(actor = c("p", "p", "q", "q"), code = c("B", "A", "B", "A"),
                    time = c(2, 1, 2, 1), session_id = c("s1", "s1", "s2", "s2"))
  says <- function(...) {
    said <- character()
    withCallingHandlers(
      build_network(log, method = "frequency", action = "code", actor = "actor", ...),
      nestimate_column_role = function(m) {
        said <<- c(said, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    said
  }
  # time declared: events follow it, so the message must not say row order
  with_time <- says(time = "time")
  expect_match(with_time, "`session_id` not used", fixed = TRUE)
  expect_false(grepl("row order", with_time, fixed = TRUE))
  # session declared: the message must not say that no session column splits them
  with_session <- says(session = "session_id")
  expect_match(with_session, "`time` not used", fixed = TRUE)
  expect_match(with_session, "row order", fixed = TRUE)
  expect_false(grepl("no session column", with_session, fixed = TRUE))
  expect_length(says(time = "time", session = "session_id"), 0L)
})

test_that("an mcml from declared positions keeps them for later verbs", {
  weeks <- data.frame(week_a = c("A", "B", "C"), week_b = c("B", "C", "A"),
                      week_c = c("C", "A", "B"))
  clusters <- list(G = c("A", "B"), H = "C")
  fit <- build_mcml(weeks, clusters, type = "raw", cols = names(weeks))
  # A->B B->C, B->C C->A, C->A A->B: six transitions
  expect_equal(sum(fit$macro$weights), 6)
  expect_s3_class(as_htna(fit), "netobject")
  expect_s3_class(macro_network(fit, expand = "G"), "netobject")
  expect_s3_class(build_mcml(build_mcml(weeks, list(G = "A", H = "B", I = "C"),
                                        cols = names(weeks)), combine = c("G", "H")),
                  "mcml")
})

test_that("a role declared in the call comes before the mark", {
  weeks <- data.frame(week_a = c("A", "B"), week_b = c("B", "A"), week_c = c("A", "A"),
                      who = c("p", "q"))
  marked <- extract_sequences(weeks, cols = c("week_a", "week_b", "week_c"))
  # all three: A->B B->A, B->A A->A = 4; without week_a: B->A, A->A = 2
  expect_equal(sum(build_network(marked, method = "frequency")$weights), 4)
  dropped <- build_network(marked, method = "frequency", metadata_cols = "week_a")
  expect_equal(sum(dropped$weights), 2)
  expect_equal(dropped$weights[["B", "A"]], 1)
  expect_equal(dropped$weights[["A", "A"]], 1)
  as_actor <- build_network(marked, method = "frequency", actor = "week_a")
  expect_false("week_a" %in% names(as_actor$data))
})

test_that("positions are read in the order the frame holds them", {
  skip_if_not_installed("tibble")
  weeks <- data.frame(week_a = c("A", "A"), week_b = c("B", "B"), week_c = c("C", "C"))
  marked <- extract_sequences(weeks, cols = names(weeks))
  # reversed: C->B, B->A in each row
  backwards <- tibble::as_tibble(marked)[rev(names(marked))]
  net <- build_network(backwards, method = "frequency")
  expect_equal(net$weights[["C", "B"]], 2)
  expect_equal(net$weights[["A", "B"]], 0)
  expect_equal(net$weights, build_network(backwards, method = "frequency",
                                          state_cols = names(backwards))$weights)
  expect_equal(net$weights[["C", "B"]], unclass(frequencies(backwards))[["C", "B"]])
})

test_that("declared positions are not read as one-hot codes", {
  coded <- data.frame(week_a = c(1, 0), week_b = c(0, 0), week_c = c(1, 1))
  # 1->0 0->1, 0->0 0->1: four transitions among the states 0 and 1
  declared <- build_network(coded, method = "frequency", state_cols = names(coded))
  expect_identical(sort(rownames(declared$weights)), c("0", "1"))
  expect_equal(sum(declared$weights), 4)
  marked <- extract_sequences(data.frame(week_a = c("x", "y"), week_b = c("y", "y"),
                                         week_c = c("x", "x")),
                              cols = c("week_a", "week_b", "week_c"))
  marked[] <- lapply(marked, function(v) as.numeric(v == "x"))
  read <- build_network(marked, method = "frequency")
  expect_identical(sort(rownames(read$weights)), c("0", "1"))
  expect_equal(sum(read$weights), 4)
})

test_that("one-hot data is declared, never guessed from the values", {
  coded <- data.frame(week_a = c(1, 0), week_b = c(0, 0), week_c = c(1, 1))
  numbered <- data.frame(T1 = c(1, 0), T2 = c(0, 0), T3 = c(1, 1))
  # nothing declared: refused, with numbered names too
  expect_error(build_network(coded, method = "frequency"), class = "nestimate_bad_input")
  expect_error(build_network(numbered, method = "relative"), class = "nestimate_bad_input")
  expect_error(frequencies(numbered), class = "nestimate_bad_input")
  expect_error(build_hon(numbered), class = "nestimate_bad_input")
  # a binary column added to positions does not reach a one-hot reading
  months <- data.frame(week_a = c("A", "B"), week_b = c("B", "A"), week_c = c("A", "A"))
  marked <- extract_sequences(months, cols = names(months))
  marked[] <- lapply(marked, function(v) as.numeric(v == "A"))
  marked$grade <- 0
  expect_error(build_network(marked, method = "frequency", metadata_cols = "grade"),
               class = "nestimate_bad_input")
  # the advice works: declared as one-hot, three codes at two time points
  hot <- build_network(coded, method = "co_occurrence", format = "onehot")
  expect_identical(sort(rownames(hot$weights)), c("week_a", "week_b", "week_c"))
  by_codes <- build_network(coded, method = "co_occurrence", codes = names(coded))
  expect_equal(hot$weights, by_codes$weights)
  # declared as positions: 1->0 0->1, 0->0 0->1, four transitions
  expect_equal(sum(build_network(coded, method = "frequency",
                                 state_cols = names(coded))$weights), 4)
  # association methods read 0/1 columns as variables, and wtna methods name one-hot
  vars <- data.frame(a = c(1, 0, 1, 0, 1, 1), b = c(0, 0, 1, 1, 1, 0), c = c(1, 1, 0, 0, 1, 0))
  expect_identical(dim(build_network(vars, method = "cor")$weights), c(3L, 3L))
  expect_no_error(build_network(learning_activities, method = "wtna", actor = "student"))
})

test_that("cooccurrence() does not count a missing cell", {
  b <- data.frame(x = c(1, 1, NA), y = c(1, NA, 1), z = c(0, 1, 1))
  co <- cooccurrence(b, format = "binary")
  # row 1: x y; row 2: x z; row 3: y z
  expect_equal(co$weights[["x", "y"]], 1)
  expect_equal(co$weights[["x", "z"]], 1)
  expect_equal(co$weights[["y", "z"]], 1)
})

test_that("a declared order is not called row order", {
  events <- data.frame(who = "p", what = c("A", "C", "B"), pos = c(1, 3, 2),
                       time = c(100, 200, 300))
  said <- character()
  net <- withCallingHandlers(
    build_network(events, method = "frequency", action = "what", actor = "who",
                  order = "pos"),
    nestimate_column_role = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_identical(unlist(extract_sequences(net)[1L, ], use.names = FALSE),
                   c("A", "B", "C"))
  expect_match(said, "`time` not used", fixed = TRUE)
  expect_false(grepl("row order", said, fixed = TRUE))
})

test_that("cooccurrence() on a list refuses roles it has no columns for", {
  transactions <- list(c("a", "b"), c("b", "c"))
  expect_error(cooccurrence(transactions, field = "TYPO"), class = "nestimate_bad_input")
  expect_error(cooccurrence(transactions, by = "x"), class = "nestimate_bad_input")
  expect_error(cooccurrence(transactions, sep = ";"), class = "nestimate_bad_input")
  co <- cooccurrence(transactions)
  expect_equal(co$weights[["a", "b"]], 1)
  expect_equal(co$weights[["a", "c"]], 0)
})

test_that("binary items are declared, never guessed from the values", {
  items <- data.frame(x = c(1, 1, 0), y = c(1, 0, 1), z = c(0, 1, 1))
  expect_error(cooccurrence(items), class = "nestimate_bad_input")
  expect_error(cooccurrence(as.matrix(items)), class = "nestimate_bad_input")
  expect_error(association_rules(items), class = "nestimate_bad_input")
  expect_error(association_rules(as.matrix(items)), class = "nestimate_bad_input")
  # x and y once, x and z once, y and z once
  co <- cooccurrence(items, format = "binary")
  expect_equal(co$weights[["x", "y"]], 1)
  expect_equal(co$weights[["y", "z"]], 1)
  rules <- association_rules(items, format = "binary", min_support = 0.1,
                             min_confidence = 0.1)
  expect_s3_class(rules, "net_association_rules")
  # the argument is not a way to read anything else as items
  expect_error(cooccurrence(data.frame(a = c("p", "q"), b = c("q", "p")), format = "binary"),
               class = "nestimate_bad_input")
  expect_error(association_rules(data.frame(a = c("p", "q"), b = c("q", "p")),
                                 format = "binary"),
               class = "nestimate_bad_input")
})

test_that("one-hot declarations keep metadata out and refuse positional methods", {
  q <- data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
                  C = rep(c(0, 0, 1, 1), 6L), actor = rep(paste0("p", 1:6), each = 4L))
  p <- q
  p$flag <- rep(c(0, 1), 12L)
  by_format <- build_network(p, method = "frequency", format = "onehot", actor = "actor",
                             metadata_cols = "flag", window_size = 2)
  by_codes <- build_network(p, method = "frequency", codes = c("A", "B", "C"),
                            actor = "actor", window_size = 2)
  expect_identical(sort(rownames(by_format$weights)), c("A", "B", "C"))
  expect_equal(by_format$weights, by_codes$weights)
  invisible(lapply(c("attention", "ngram", "gap", "reverse"), function(m) {
    expect_error(build_network(q, method = m, format = "onehot", actor = "actor"),
                 class = "nestimate_bad_input")
    expect_error(build_network(q, method = m, codes = c("A", "B", "C"), actor = "actor"),
                 class = "nestimate_bad_input")
  }))
})

test_that("a one-hot network is not read as positions of 0 and 1", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  q <- data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
                  C = rep(c(0, 0, 1, 1), 6L), actor = rep(paste0("p", 1:6), each = 4L))
  net <- build_network(q, method = "frequency", format = "onehot", actor = "actor",
                       window_size = 2)
  readers <- list(build_hon, build_hypa, build_mogen, path_counts, path_dependence,
                  state_frequencies, extract_sequences, sequence_metrics, sequence_patterns,
                  sequence_plot, distribution_plot,
                  function(x) markov_order_test(x, n_perm = 2L),
                  function(x) sequence_compare(x, group = rep(c("a", "b"), 12L)))
  invisible(lapply(readers, function(read) {
    expect_error(read(net), class = "nestimate_bad_input")
  }))
  # build_mcml aggregates the weights: nodes A and B into G, C into H
  fit <- build_mcml(net, list(G = c("A", "B"), H = "C"))
  expect_identical(sort(rownames(fit$macro$weights)), c("G", "H"))
})

test_that("params$cols declares positions, as state_cols does", {
  coded <- data.frame(T1 = c(1, 0), T2 = c(0, 0), T3 = c(1, 1))
  net <- build_network(coded, method = "frequency", params = list(cols = names(coded)))
  expect_identical(sort(rownames(net$weights)), c("0", "1"))
  expect_equal(net$weights, build_network(coded, method = "frequency",
                                          state_cols = names(coded))$weights)
})

test_that("entropy_bayes() keeps a network's positions, and says when it has no counts", {
  coded <- data.frame(T1 = c(1, 0, 1, 1), T2 = c(0, 0, 1, 0), T3 = c(1, 1, 0, 1))
  rel <- build_network(coded, method = "relative", state_cols = names(coded))
  freq <- build_network(coded, method = "frequency", state_cols = names(coded))
  expect_s3_class(entropy_bayes(rel, draws = 100L, seed = 1), "net_entropy_bayes")
  # the same counts under the same seed give the same draws
  expect_equal(entropy_bayes(rel, draws = 100L, seed = 1)$draws_entropy_rate,
               entropy_bayes(freq, draws = 100L, seed = 1)$draws_entropy_rate)
  q <- data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
                  C = rep(c(0, 0, 1, 1), 6L), actor = rep(paste0("p", 1:6), each = 4L))
  hot <- build_network(q, method = "relative", format = "onehot", actor = "actor",
                       window_size = 2)
  # a relative one-hot network is recounted per unit, as its frequency twin
  counts <- build_network(q, method = "frequency", format = "onehot",
                          actor = "actor", window_size = 2)
  expect_equal(Nestimate:::.te_extract_counts(hot), counts$weights)
})

test_that("a resample that loses a state counts as below the threshold", {
  scores <- data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
                       C = rep(c(0, 0, 1, 1), 6L))
  net <- build_network(scores, method = "cor")
  # some subsets leave `B` constant, which the estimator drops; stopped with
  # `subscript out of bounds` before
  stability <- suppressWarnings(suppressMessages(
    centrality_stability(net, iter = 3L, seed = 1)))
  expect_s3_class(stability, "net_stability")
  expect_named(stability$cs, c("InStrength", "OutStrength", "Betweenness"))
})


# ---- Round six: declared 0/1 positions, nested one-hot, association networks ----

.round6_onehot <- function() {
  data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
             C = rep(c(0, 0, 1, 1), 6L), actor = rep(paste0("p", 1:6), each = 4L))
}

test_that("co_occurrence reads declared 0/1 positions as the states 0 and 1", {
  u <- data.frame(T1 = c(1, 0), T2 = c(0, 0), T3 = c(1, 1))
  # rows (1,0,1) and (0,0,1): unordered position pairs 1-0, 1-1, 0-1 and 0-0,
  # 0-1, 0-1, so four mixed pairs, one 0-0 and one 1-1
  expected <- matrix(c(1, 4, 4, 1), 2L, dimnames = list(c("0", "1"), c("0", "1")))
  expect_equal(build_network(u, method = "co_occurrence",
                             state_cols = names(u))$weights, expected)
  expect_equal(build_network(u, method = "co_occurrence",
                             params = list(cols = names(u)))$weights, expected)
})

test_that("an association network carries no guessed one-hot declaration", {
  items <- data.frame(x = c(1, 1, 0), y = c(1, 0, 1), z = c(0, 1, 1))
  net <- build_network(items, method = "cor")
  # each pair: deviations (1, 1, -2)/3 against a permutation of them, covariance
  # sum -3/9 over variance sums 6/9, so r = -0.5
  expect_equal(unname(net$weights[upper.tri(net$weights)]), rep(-0.5, 3L))
  expect_null(net$params$codes)
  expect_null(net$params$window_size)
  # its rows are variables measured together, not transactions or sequences
  expect_error(association_rules(net, min_support = 0.1, min_confidence = 0.1),
               class = "nestimate_bad_input")
  expect_error(extract_sequences(net), class = "nestimate_bad_input")
  rules <- association_rules(items, format = "binary", min_support = 0.1,
                             min_confidence = 0.1, min_lift = 0)
  # three transactions, each pair together in one of them
  expect_equal(unique(summary(rules)$support), 1 / 3)
  expect_equal(unique(summary(rules)$confidence), 1 / 2)
})

test_that("a binary MGM network resamples its rows, not one sequence", {
  a <- rep(c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0), 6L)
  b <- a
  flip <- c(3L, 14L, 25L, 36L, 47L, 58L)
  b[flip] <- 1 - b[flip]
  scores <- data.frame(A = a, B = b, C = rep(c(1, 1, 0, 0, 1, 0, 1, 0, 0, 1), 6L))
  fit <- function(rows) {
    build_network(scores[rows, ], method = "mgm", type = rep("c", 3L),
                  level = rep(2L, 3L))$weights
  }
  net <- build_network(scores, method = "mgm", type = rep("c", 3L),
                       level = rep(2L, 3L))
  boot <- bootstrap_network(net, iter = 5L, seed = 1)
  # the bootstrap mean is the mean of the fits on rows drawn with replacement
  set.seed(1)
  refits <- lapply(seq_len(5L), \(i) fit(sample(nrow(scores), replace = TRUE)))
  expect_equal(boot$mean, Reduce(`+`, refits) / 5, ignore_attr = TRUE)
  rel <- network_reliability(net, iter = 3L, seed = 1)
  expect_s3_class(rel, "net_reliability")
  # the summary is the mean over the split-half iterations
  expect_equal(subset(rel$summary, metric == "mean_dev")$mean,
               mean(rel$iterations$mean_dev))
  expect_equal(subset(rel$summary, metric == "cor")$mean,
               mean(rel$iterations$cor))
  expect_true(all(rel$iterations$cor >= -1 & rel$iterations$cor <= 1))
  # each iteration compares the fits on a seeded half and on the rest
  set.seed(1)
  halves <- lapply(seq_len(3L), \(i) {
    first <- sample.int(nrow(scores), 30L)
    list(fit(first), fit(setdiff(seq_len(nrow(scores)), first)))
  })
  expect_equal(rel$iterations$mean_dev,
               vapply(halves, \(h) mean(abs(h[[1]] - h[[2]])), numeric(1)))
  expect_equal(rel$iterations$cor,
               vapply(halves, \(h) stats::cor(c(h[[1]]), c(h[[2]])), numeric(1)))
  expect_equal(rel$iterations$median_dev,
               vapply(halves, \(h) stats::median(abs(h[[1]] - h[[2]])), numeric(1)))
  expect_equal(rel$iterations$max_dev,
               vapply(halves, \(h) max(abs(h[[1]] - h[[2]])), numeric(1)))
  expect_equal(subset(rel$summary, metric == "max_dev")$mean,
               mean(rel$iterations$max_dev))
})

test_that("sequence readers refuse the one-hot networks of a group", {
  q <- .round6_onehot()
  q$cohort <- rep(c("g1", "g2"), each = 12L)
  group <- build_network(q, method = "frequency", format = "onehot", actor = "actor",
                         group = "cohort", window_size = 2L)
  expect_error(sequence_plot(group), class = "nestimate_bad_input")
  expect_error(distribution_plot(group), class = "nestimate_bad_input")
  expect_error(sequence_compare(group, sub = 2L, min_freq = 1L, iter = 2L,
                                test = "none"), class = "nestimate_bad_input")
  expect_error(build_hon(group, max_order = 1L), class = "nestimate_bad_input")
  expect_error(build_mogen(group, max_order = 1L), class = "nestimate_bad_input")
  expect_error(build_hypa(group, order = 1L, min_count = 1L),
               class = "nestimate_bad_input")
  expect_error(state_frequencies(group), class = "nestimate_bad_input")
  expect_error(path_counts(group), class = "nestimate_bad_input")
})

test_that("a bare cograph_network keeps its one-hot declaration", {
  hot <- build_network(.round6_onehot(), method = "frequency", format = "onehot",
                       actor = "actor", window_size = 2L)
  cg <- hot
  class(cg) <- "cograph_network"
  expect_error(extract_sequences(cg), class = "nestimate_bad_input")
  expect_error(build_hon(cg, max_order = 1L), class = "nestimate_bad_input")
  expect_identical(.as_netobject(cg)$params$codes, c("A", "B", "C"))
  rules <- association_rules(cg, min_support = 0.1, min_confidence = 0.1)
  expect_identical(rules$items, c("A", "B", "C"))
  # windows of two rows per actor are {A, B} and {A, C}: 12 windows over six
  # actors, A in all of them, B and C in half, AB and AC each in half
  expect_identical(rules$n_transactions, 12L)
  tab <- summary(rules)
  expect_equal(tab$support, rep(0.5, 4L))
  expect_equal(subset(tab, antecedent == "B" & consequent == "A")$confidence, 1)
  expect_equal(subset(tab, antecedent == "A" & consequent == "B")$confidence, 0.5)
})

test_that("a sequence reader's advice names only routes it reads", {
  ua <- data.frame(T1 = rep(c(1, 0, 1, 1, 0, 0), 5L), T2 = rep(c(0, 0, 1, 0, 1, 1), 5L),
                   T3 = rep(c(1, 1, 0, 1, 0, 1), 5L))
  err <- tryCatch(build_hon(ua, max_order = 1L), error = identity)
  expect_s3_class(err, "nestimate_bad_input")
  expect_false(grepl("format = \"onehot\"", conditionMessage(err), fixed = TRUE))
  expect_s3_class(build_hon(build_network(ua, method = "frequency",
                                          state_cols = names(ua)), max_order = 1L),
                  "net_hon")
  err <- tryCatch(build_network(ua, method = "attention"), error = identity)
  expect_false(grepl("format = \"onehot\"", conditionMessage(err), fixed = TRUE))
  expect_s3_class(build_network(ua, method = "attention", state_cols = names(ua)),
                  "netobject")
})

test_that("cooccurrence() reads logical items under format = 'binary' and refuses networks", {
  items <- data.frame(x = c(1, 1, 0), y = c(1, 0, 1), z = c(0, 1, 1))
  expected <- cooccurrence(items, format = "binary")$weights
  # each item is in two rows, each pair together in one
  expect_equal(unname(diag(expected)), c(2, 2, 2))
  expect_equal(unname(expected[upper.tri(expected)]), c(1, 1, 1))
  expect_equal(cooccurrence(items > 0, format = "binary")$weights, expected)
  expect_equal(cooccurrence(as.data.frame(items > 0), format = "binary")$weights,
               expected)
  expect_error(cooccurrence(items > 0), class = "nestimate_bad_input")
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"))
  expect_error(cooccurrence(build_network(seqs, method = "frequency")),
               class = "nestimate_bad_input")
})

test_that("duplicated column names are refused, so no added column is trusted", {
  d <- data.frame(week_a = c("A", "B"), week_b = c("B", "A"), week_c = c("C", "B"))
  x <- extract_sequences(d, cols = names(d))
  x[[4]] <- c("NOVEL", "NOVEL")
  names(x)[4] <- "week_a"
  expect_null(.marked_positions(x))
  expect_error(build_hon(x, max_order = 1L), class = "nestimate_bad_input")
  expect_error(frequencies(x), class = "nestimate_bad_input")
  expect_error(build_network(x, method = "frequency"), class = "nestimate_bad_input")
})

test_that("one 0/1 column, a 0/1 matrix and a 0/1 run for build_mcml are not guessed", {
  u <- data.frame(T1 = c(1, 0), T2 = c(0, 0), T3 = c(1, 1))
  expect_error(build_network(u["T1"], method = "frequency"),
               class = "nestimate_bad_input")
  expect_error(state_frequencies(u["T1"]), class = "nestimate_bad_input")
  expect_error(sequence_metrics(as.matrix(u)), class = "nestimate_bad_input")
  expect_error(sequence_plot(as.matrix(u)), class = "nestimate_bad_input")
  expect_error(path_counts(as.matrix(u)), class = "nestimate_bad_input")
  expect_error(build_hon(as.matrix(u) > 0, max_order = 1L),
               class = "nestimate_bad_input")
  expect_error(build_mcml(u, list(G = "0", H = "1"), type = "raw"),
               class = "nestimate_bad_input")
  macro <- build_mcml(u, list(G = "0", H = "1"), type = "raw", cols = names(u))
  # transitions 1->0, 0->1 in row one and 0->0, 0->1 in row two
  expect_equal(unname(macro$macro$weights), matrix(c(1, 1, 2, 0), 2L))
})

test_that("case-dropping and split halves survive a resample that loses a state", {
  scores <- data.frame(A = rep(c(1, 0, 1, 0, 1), 12L), B = rep(c(0, 1, 1, 0, 0), 12L),
                       C = c(1, rep(0, 59L)))
  net <- suppressMessages(build_network(scores, method = "cor"))
  # most subsets lose the one row where `C` is 1, and C is then constant
  expect_warning(
    cd <- suppressMessages(casedrop_reliability(net, iter = 3L, seed = 1)),
    class = "nestimate_resample_failed")
  expect_s3_class(cd, "net_casedrop_reliability")
  # the CS-coefficient is one of the drop proportions, or 0 when none qualifies
  expect_true(cd$cs %in% c(0, cd$drop_prop))
  expect_true(all(cd$correlations >= -1 & cd$correlations <= 1, na.rm = TRUE))
  # each correlation is the Spearman correlation of the off-diagonal edges with
  # those of a refit on the seeded subset; a subset that loses C is NA
  off <- !diag(3L)
  refit <- function(rows) {
    sub <- suppressMessages(build_network(scores[rows, ], method = "cor"))
    if (!all(c("A", "B", "C") %in% rownames(sub$weights))) return(NA_real_)
    stats::cor(net$weights[off],
               sub$weights[c("A", "B", "C"), c("A", "B", "C")][off],
               method = "spearman")
  }
  set.seed(1)
  recomputed <- vapply(cd$drop_prop, \(p) vapply(seq_len(3L), \(i)
    refit(sample(60L, 60L - floor(60L * p))), numeric(1)), numeric(3))
  expect_equal(unname(cd$correlations), recomputed)
  # CS is the largest proportion whose draws all reach 0.7 (NA counts below)
  reached <- colMeans(!is.na(recomputed) & recomputed >= 0.7) >= 0.95
  expect_equal(cd$cs, max(cd$drop_prop[reached]))
  expect_equal(cd$cs, 0.5)
  expect_warning(
    rel <- suppressMessages(network_reliability(net, iter = 3L, seed = 1)),
    class = "nestimate_resample_failed")
  expect_equal(subset(rel$summary, metric == "mean_dev")$mean,
               mean(rel$iterations$mean_dev, na.rm = TRUE))
})

# ---- Round seven ----

test_that("wtna refuses positions declared beside its one-hot reading", {
  u <- data.frame(T1 = c(1, 0), T2 = c(0, 0), T3 = c(1, 1))
  expect_error(build_network(u, method = "wtna", state_cols = names(u)),
               class = "nestimate_bad_input")
  expect_error(build_network(u, method = "wtna_cooccurrence",
                             params = list(cols = names(u))),
               class = "nestimate_bad_input")
})

test_that("a fitted result is not read as sequences or transactions", {
  x <- stats::t.test(1:6, 2:7)
  expect_error(build_hon(x, max_order = 1L), class = "nestimate_bad_input")
  expect_error(path_counts(x), class = "nestimate_bad_input")
  expect_error(state_frequencies(x), class = "nestimate_bad_input")
  expect_error(cooccurrence(x), class = "nestimate_bad_input")
  # a plain list of state vectors still reads: A->B once, B->A twice
  hon <- build_hon(list(c("A", "B", "A"), c("B", "A")), max_order = 1L)
  expect_s3_class(hon, "net_hon")
  expect_identical(hon$ho_edges$path, c("A -> B", "B -> A"))
  expect_equal(hon$ho_edges$count, c(1, 2))
})

test_that("build_mmm() refuses an association network", {
  scores <- data.frame(A = rep(c(1, 0, 1, 0, 1), 12L), B = rep(c(0, 1, 1, 0, 0), 12L),
                       C = rep(c(1, 1, 0, 0, 1), 12L))
  net <- build_network(scores, method = "mgm", type = rep("c", 3L),
                       level = rep(2L, 3L))
  expect_error(build_mmm(net, k = 2L, n_starts = 1L, max_iter = 5L, seed = 1),
               class = "nestimate_bad_input")
  expect_error(cluster_mmm(net, k = 2L, n_starts = 1L, max_iter = 5L, seed = 1),
               class = "nestimate_bad_input")
})

test_that("association_rules() advice on counts can be followed", {
  counts <- matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 1), 3L, byrow = TRUE,
                   dimnames = list(NULL, c("x", "y", "z")))
  err <- tryCatch(association_rules(counts), error = identity)
  expect_s3_class(err, "nestimate_bad_input")
  expect_true(grepl("format = \"binary\"", conditionMessage(err), fixed = TRUE))
  rules <- association_rules(counts > 0, format = "binary", min_support = 0.1,
                             min_confidence = 0.1, min_lift = 0)
  expect_identical(rules$n_transactions, 3L)
  # transactions {x, y}, {x, y, z}, {y, z}: x is in two, both with y
  r <- summary(rules)
  expect_equal(r$support[r$antecedent == "x" & r$consequent == "y"], 2 / 3)
  expect_equal(r$confidence[r$antecedent == "x" & r$consequent == "y"], 1)
  expect_equal(r$confidence[r$antecedent == "y" & r$consequent == "x"], 2 / 3)
})

# ---- Round eight ----

test_that("wtna reads logical one-hot columns by its name, as numeric ones", {
  hot <- data.frame(A = c(TRUE, FALSE, FALSE, TRUE), B = c(FALSE, TRUE, FALSE, TRUE),
                    C = c(FALSE, FALSE, TRUE, FALSE))
  # windows of one row: A -> B, B -> C, then C -> A and C -> B
  expect_warning(net <- build_network(hot, method = "wtna", window_size = 1L))
  expect_equal(unname(net$weights),
               matrix(c(0, 0, 1, 1, 0, 1, 0, 1, 0), 3L))
  expect_identical(net$params$codes, c("A", "B", "C"))
  one <- data.frame(A = c(TRUE, FALSE, FALSE, TRUE), B = c(FALSE, TRUE, FALSE, FALSE),
                    C = c(FALSE, FALSE, TRUE, FALSE))
  # one code per window of one row: nothing co-occurs with anything else
  expect_warning(co <- build_network(one, method = "wtna_cooccurrence",
                                     window_size = 1L))
  expect_equal(unname(co$weights), diag(c(2, 1, 1)))
  expect_identical(co$params$window_size, 1L)
  expect_error(build_hon(co, max_order = 1L), class = "nestimate_bad_input")
})

test_that("association_rules() refuses a classed list of character fields", {
  expect_error(association_rules(utils::packageDescription("base")),
               class = "nestimate_bad_input")
  expect_error(association_rules(list(c("a", "b"), 1:2)),
               class = "nestimate_bad_input")
  rules <- association_rules(list(c("a", "b"), c("a", "c"), c("a", "b")),
                             min_support = 0.1, min_confidence = 0.1)
  expect_identical(rules$n_transactions, 3L)
})

test_that("build_mmm() refuses a bare cograph_network of an association model", {
  scores <- data.frame(A = rep(c(1, 0, 1, 0, 1), 12L), B = rep(c(0, 1, 1, 0, 0), 12L),
                       C = rep(c(1, 1, 0, 0, 1), 12L))
  net <- build_network(scores, method = "mgm", type = rep("c", 3L),
                       level = rep(2L, 3L))
  class(net) <- "cograph_network"
  expect_error(build_mmm(net, k = 2L, n_starts = 1L, max_iter = 5L, seed = 1),
               class = "nestimate_bad_input")
  expect_error(cluster_mmm(net, k = 2L, n_starts = 1L, max_iter = 5L, seed = 1),
               class = "nestimate_bad_input")
  seqs <- data.frame(T1 = c("a", "b", "a", "c", "a", "b"),
                     T2 = c("b", "a", "c", "a", "b", "c"),
                     T3 = c("c", "c", "b", "b", "a", "a"))
  freq <- build_network(seqs, method = "frequency")
  class(freq) <- "cograph_network"
  fit <- suppressWarnings(build_mmm(freq, k = 2L, n_starts = 1L, max_iter = 5L,
                                    seed = 1))
  # the six sequences are shared out over the two components
  expect_s3_class(fit, "net_mmm")
  expect_identical(fit$n_sequences, 6L)
  expect_equal(sum(vapply(fit$models, \(m) nrow(m$data), integer(1))), 6L)
})

# ---- Round nine ----

.round9_units <- function() {
  # u1 and u2 run A -> B -> C three times, u3 and u4 run C -> B -> A
  hot <- function(order) diag(3L)[rep(order, 3L), , drop = FALSE]
  rows <- rbind(hot(1:3), hot(1:3), hot(3:1), hot(3:1))
  data.frame(id = rep(paste0("u", 1:4), each = 9L),
             A = rows[, 1], B = rows[, 2], C = rows[, 3])
}

test_that("build_mmm() counts a one-hot network within its windows, per unit", {
  x <- .round9_units()
  net <- suppressWarnings(build_network(x, method = "frequency", format = "onehot",
                                        codes = c("A", "B", "C"), actor = "id",
                                        window_size = 1L))
  # per unit A -> B 3, B -> C 3, C -> A 2, or the reverse
  expect_equal(unname(net$weights), matrix(c(0, 6, 4, 6, 0, 6, 4, 6, 0), 3L))
  fit <- build_mmm(net, k = 2L, n_starts = 5L, seed = 1)
  expect_identical(fit$n_sequences, 4L)
  expect_identical(as.vector(table(fit$assignments)), c(2L, 2L))
  expect_identical(fit$assignments[1], fit$assignments[2])
  expect_identical(fit$assignments[3], fit$assignments[4])
  # the cluster networks add up to the network
  nets <- build_network(fit)
  expect_equal(Reduce(`+`, lapply(nets, \(n) n$weights)), net$weights)
  # the likelihood of the fitted parameters on the hand counts: u1, u2 start
  # in A with A -> B 3, B -> C 3, C -> A 2; u3, u4 start in C, reversed
  ahead <- matrix(0, 3L, 3L, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  ahead["A", "B"] <- 3
  ahead["B", "C"] <- 3
  ahead["C", "A"] <- 2
  units <- list(list(n = ahead, first = "A"), list(n = ahead, first = "A"),
                list(n = t(ahead), first = "C"), list(n = t(ahead), first = "C"))
  ll <- sum(vapply(units, \(u) log(sum(vapply(seq_len(2L), \(k) {
    m <- fit$models[[k]]
    fit$mixing[k] * m$initial[[u$first]] * prod(m$weights^u$n)
  }, numeric(1)))), numeric(1)))
  expect_equal(fit$log_likelihood, ll)
  # 2 x 3 x 2 transition, 2 x 2 initial and 1 mixing parameters
  expect_identical(fit$n_params, 17L)
  expect_equal(fit$BIC, -2 * ll + 17 * log(4))
  bare <- net
  class(bare) <- "cograph_network"
  expect_identical(build_mmm(bare, k = 2L, n_starts = 5L, seed = 1)$assignments,
                   fit$assignments)
  # without an actor a bare one-hot network is one unit
  alone <- suppressWarnings(build_network(x[c("A", "B", "C")], method = "frequency",
                                          format = "onehot", window_size = 1L))
  class(alone) <- "cograph_network"
  expect_error(build_mmm(alone, k = 2L, seed = 1), class = "nestimate_bad_input")
})

test_that("rebuilding an empty mixture component is refused by name", {
  seqs <- data.frame(T1 = rep("A", 4L), T2 = rep("B", 4L), T3 = rep("C", 4L))
  fit <- suppressWarnings(build_mmm(build_network(seqs, method = "frequency"),
                                    k = 2L, n_starts = 1L, max_iter = 5L, seed = 1))
  expect_identical(fit$n_sequences, 4L)
  expect_error(build_network(fit), class = "nestimate_empty_cluster")
})

test_that("wtna keeps declared metadata out of its codes", {
  x <- data.frame(A = c(TRUE, FALSE, FALSE, TRUE), B = c(FALSE, TRUE, FALSE, TRUE),
                  C = c(FALSE, FALSE, TRUE, FALSE), age = rep(20, 4L))
  net <- suppressWarnings(build_network(x, method = "wtna", metadata_cols = "age",
                                        window_size = 1L))
  # windows of one row: A -> B, B -> C, then C -> A and C -> B
  expect_equal(unname(net$weights), matrix(c(0, 0, 1, 1, 0, 1, 0, 1, 0), 3L))
  expect_identical(net$params$window_size, 1L)
  co <- suppressWarnings(build_network(x, method = "wtna_cooccurrence",
                                       metadata_cols = "age", window_size = 1L))
  # A and B share row 4; each code occurs in its own rows
  expect_equal(unname(co$weights), matrix(c(2, 1, 0, 1, 2, 0, 0, 0, 1), 3L))
})

test_that("entropy_bayes() reads counts only from a network that holds them", {
  same <- data.frame(A = 1:6, B = 1:6, C = 1:6)
  expect_warning(cor_net <- build_network(same, method = "cor"),
                 class = "nestimate_singular_neighbours")
  expect_error(entropy_bayes(cor_net, draws = 100L, seed = 1),
               class = "nestimate_bad_input")
  class(cor_net) <- "cograph_network"
  expect_error(entropy_bayes(cor_net, draws = 100L, seed = 1),
               class = "nestimate_bad_input")
  y <- data.frame(id = rep(paste0("u", 1:4), each = 3L), A = rep(c(1, 0, 0), 4L),
                  B = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1),
                  C = c(0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0))
  rel <- suppressWarnings(build_network(y, method = "relative", format = "onehot",
                                        actor = "id", window_size = 1L))
  freq <- suppressWarnings(build_network(y, method = "frequency", format = "onehot",
                                         actor = "id", window_size = 1L))
  # ABC, ACB, ABC, ACB: eight transitions, recounted per unit from either
  hand <- matrix(c(0, 0, 0, 2, 0, 2, 2, 2, 0), 3L)
  expect_equal(unname(Nestimate:::.te_extract_counts(freq)), hand)
  expect_equal(unname(Nestimate:::.te_extract_counts(rel)), hand)
  class(rel) <- "cograph_network"
  expect_equal(unname(Nestimate:::.te_extract_counts(rel)), hand)
  # a deterministic chain has whole-number probabilities, not counts of one
  seqs <- data.frame(T1 = rep("A", 4L), T2 = rep("B", 4L), T3 = rep("C", 4L))
  chain <- build_network(seqs, method = "relative")
  four <- matrix(c(0, 0, 0, 4, 0, 0, 0, 4, 0), 3L)
  expect_equal(unname(Nestimate:::.te_extract_counts(chain)), four)
  # scaled counts are not counts either
  scaled <- build_network(seqs, method = "frequency", scaling = "max")
  expect_equal(unname(Nestimate:::.te_extract_counts(scaled)), four)
  # the network's own start and end states and concatenation are recounted
  one <- data.frame(T1 = "A", T2 = "B", T3 = "C")
  ends <- Nestimate:::.te_extract_counts(
    build_network(one, method = "relative", start = TRUE, end = TRUE))
  expect_identical(rownames(ends), c("Start", "A", "B", "C", "End"))
  expect_equal(unname(ends), rbind(cbind(0, diag(4L)), 0))
  two <- data.frame(T1 = c("A", "C"), T2 = c("B", "D"))
  joined <- Nestimate:::.te_extract_counts(
    build_network(two, method = "relative", params = list(concat = 2L)))
  # A B C D read as one sequence: the bridge B -> C counts
  expect_equal(unname(joined), rbind(cbind(0, diag(3L)), 0))
  # a network with weights that are not counts and no counting rule is refused
  expect_error(Nestimate:::.te_extract_counts(build_network(seqs, method = "attention")),
               class = "nestimate_bad_input")
  # relative co-occurrence: A and B each present in two windows of one row
  pairs <- data.frame(A = c(1, 1, 0, 0), B = c(0, 0, 1, 1))
  co <- suppressWarnings(build_network(pairs, method = "co_occurrence",
                                       format = "onehot", window_size = 1L,
                                       type = "relative"))
  expect_equal(unname(Nestimate:::.te_extract_counts(co)), diag(c(2, 2)))
  # weighted counts divide by sequence length; the posterior needs the raw ones
  weighted <- build_network(seqs, method = "frequency",
                            params = list(weighted = TRUE))
  expect_equal(unname(Nestimate:::.te_extract_counts(weighted)), four)
  skip_if_not_installed("tna")
  # a relative tna model declares its type; its counts come from its sequences
  expect_equal(unname(Nestimate:::.te_extract_counts(tna::tna(seqs))), four)
  # counting types hold their counts: reverse B -> A, C -> B; co-occurrence
  # every pair of A, B, C in each of the four rows
  reverse <- Nestimate:::.te_extract_counts(tna::build_model(seqs, type = "reverse"))
  expect_equal(unname(reverse), t(four))
  co_tna <- Nestimate:::.te_extract_counts(
    tna::build_model(seqs, type = "co-occurrence"))
  expect_equal(unname(co_tna), 4 * (1 - diag(3L)))
  # scaled and weighted models are recounted from their sequences
  expect_equal(unname(Nestimate:::.te_extract_counts(
    tna::ftna(seqs, scaling = "max"))), four)
  expect_equal(unname(Nestimate:::.te_extract_counts(
    tna::ftna(seqs, params = list(weighted = TRUE)))), four)
  # a weighted model of any type is not read as counts: reverse and
  # co-occurrence cannot be replayed, so they are refused
  expect_error(Nestimate:::.te_extract_counts(tna::build_model(
    seqs, type = "reverse", params = list(weighted = TRUE))),
    class = "nestimate_bad_input")
  expect_error(Nestimate:::.te_extract_counts(tna::build_model(
    seqs, type = "co-occurrence", params = list(weighted = TRUE))),
    class = "nestimate_bad_input")
  # 3-grams of A B C D: ABC and BCD, so B -> C counts twice per row
  abcd <- data.frame(T1 = rep("A", 4L), T2 = rep("B", 4L), T3 = rep("C", 4L),
                     T4 = rep("D", 4L))
  grams <- Nestimate:::.te_extract_counts(
    tna::build_model(abcd, type = "n-gram", params = list(n_gram = 3)))
  expect_equal(unname(grams), rbind(cbind(0, diag(c(4, 8, 4))), 0))
  # decay and gap weights are not counts, and windows are not replayed
  expect_error(Nestimate:::.te_extract_counts(
    tna::build_model(seqs, type = "attention")), class = "nestimate_bad_input")
  expect_error(Nestimate:::.te_extract_counts(
    tna::build_model(seqs, type = "gap")), class = "nestimate_bad_input")
  units <- data.frame(id = rep(c("u", "v"), each = 4L), A = rep(c(1, 0, 0, 1), 2L),
                      B = rep(c(0, 1, 0, 0), 2L), C = rep(c(0, 0, 1, 0), 2L))
  windows <- tna::import_onehot(units, cols = c("A", "B", "C"), actor = "id",
                                window_size = 2L)
  expect_error(Nestimate:::.te_extract_counts(tna::tna(windows)),
               class = "nestimate_bad_input")
  # windows {A, B} then {C, A} in each unit: A, B -> C, A twice each
  counted <- Nestimate:::.te_extract_counts(tna::ftna(windows))
  expect_equal(unname(counted), matrix(c(2, 2, 0, 0, 0, 0, 2, 2, 0), 3L))
})

test_that("an MMM fit of a one-hot co-occurrence network rebuilds its clusters", {
  x <- .round9_units()
  net <- suppressWarnings(build_network(x, method = "wtna_cooccurrence",
                                        codes = c("A", "B", "C"), actor = "id",
                                        window_size = 1L))
  # one code per row: each code co-occurs with itself in its 12 rows
  expect_equal(unname(net$weights), diag(12, 3L))
  fit <- build_mmm(net, k = 2L, n_starts = 5L, seed = 1)
  nets <- build_network(fit)
  expect_equal(Reduce(`+`, lapply(nets, \(n) n$weights)), net$weights)
  routed <- cluster_network(net, k = 2L, cluster_by = "mmm", n_starts = 5L,
                            seed = 1)
  expect_equal(Reduce(`+`, lapply(routed, \(n) n$weights)), net$weights)
  # a bare network keeps its declared method for the rebuild
  bare <- net
  class(bare) <- "cograph_network"
  bare_fit <- build_mmm(bare, k = 2L, n_starts = 5L, seed = 1)
  expect_identical(bare_fit$network_method, "co_occurrence")
  expect_equal(Reduce(`+`, lapply(build_network(bare_fit), \(n) n$weights)),
               net$weights)
  # the declaration in meta$tna alone is enough
  meta_only <- bare
  meta_only$method <- NULL
  meta_only$meta$tna$method <- "wtna_cooccurrence"
  meta_fit <- build_mmm(meta_only, k = 2L, n_starts = 5L, seed = 1)
  expect_identical(meta_fit$network_method, "co_occurrence")
  skip_if_not_installed("tna")
  # a windowed tna model's windows cannot be replayed: refused, not miscounted
  windows <- tna::import_onehot(x, cols = c("A", "B", "C"), actor = "id",
                                window_size = 1L)
  expect_error(build_mmm(tna::ftna(windows), k = 2L, seed = 1),
               class = "nestimate_bad_input")
})

test_that("a one-hot unit whose first window holds two codes has no first state", {
  x <- .round9_units()
  # u1 and u2 open with A and B together
  x$B[x$id %in% c("u1", "u2") & !duplicated(x$id)] <- 1
  net <- suppressWarnings(build_network(x, method = "frequency", format = "onehot",
                                        codes = c("A", "B", "C"), actor = "id",
                                        window_size = 1L))
  fit <- build_mmm(net, k = 2L, n_starts = 5L, seed = 1)
  counts <- Nestimate:::.precompute_per_sequence_windowed(fit$data, "relative",
                                                          list(), fit$states)
  first <- c(NA, NA, "C", "C")
  ll <- sum(vapply(seq_len(4L), \(u) log(sum(vapply(seq_len(2L), \(k) {
    m <- fit$models[[k]]
    start <- if (is.na(first[u])) 1 else m$initial[[first[u]]]
    fit$mixing[k] * start * prod(t(m$weights)^counts[u, ])
  }, numeric(1)))), numeric(1)))
  expect_equal(fit$log_likelihood, ll)
})

test_that("build_clusters() refuses a bare cograph_network of an association model", {
  scores <- data.frame(A = rep(c(1, 0, 1, 0, 1), 12L), B = rep(c(0, 1, 1, 0, 0), 12L),
                       C = rep(c(1, 1, 0, 0, 1), 12L))
  net <- build_network(scores, method = "mgm", type = rep("c", 3L),
                       level = rep(2L, 3L))
  class(net) <- "cograph_network"
  expect_error(build_clusters(net, k = 2L), class = "nestimate_bad_input")
  x <- .round9_units()
  hot <- suppressWarnings(build_network(x, method = "frequency", format = "onehot",
                                        codes = c("A", "B", "C"), actor = "id",
                                        window_size = 1L))
  class(hot) <- "cograph_network"
  # four units, not 36 rows
  cl <- build_clusters(hot, k = 2L)
  expect_length(cl$assignments, 4L)
  expect_identical(cl$assignments[1], cl$assignments[2])
  expect_false(identical(cl$assignments[1], cl$assignments[3]))
})

# ---- Round thirteen ----

test_that("build_mmm() refuses inputs it would read as false sequences", {
  seqs <- data.frame(T1 = c("A", "A", "C", "C"), T2 = rep("B", 4L),
                     T3 = c("C", "C", "A", "A"), Group = c("g1", "g2", "g1", "g2"))
  groups <- build_network(seqs, method = "frequency", group = "Group")
  # its edge tables would become sequences with group labels as states
  expect_error(build_mmm(groups, k = 2L, seed = 1), class = "nestimate_bad_input")
  expect_error(cluster_network(groups, k = 2L, cluster_by = "mmm", seed = 1),
               class = "nestimate_bad_input")
  scores <- data.frame(A = 1:8, B = c(2, 1, 4, 3, 6, 5, 8, 7),
                       C = c(3, 5, 1, 7, 2, 8, 4, 6))
  cor_net <- build_network(scores, method = "cor")
  class(cor_net) <- "cograph_network"
  cor_net$meta$tna <- NULL
  # only $method declares it; the symmetric matrix would infer co_occurrence
  expect_error(build_mmm(cor_net, k = 2L, seed = 1), class = "nestimate_bad_input")
})

test_that("an MMM rebuild keeps the source network's settings", {
  seqs <- data.frame(T1 = c("A", "A", "D", "D"), T2 = c("B", "B", "C", "C"),
                     T3 = c("C", "C", "B", "B"), T4 = c("D", "D", "A", "A"))
  grams <- build_network(seqs, method = "ngram", params = list(n_gram = 3L))
  # ABC, BCD per forward row and DCB, CBA per reverse row: 16 in all
  expect_equal(sum(grams$weights), 16)
  fit <- build_mmm(grams, k = 2L, n_starts = 5L, seed = 1)
  expect_equal(Reduce(`+`, lapply(build_network(fit), \(n) n$weights)),
               grams$weights)
  routed <- cluster_network(grams, k = 2L, cluster_by = "mmm", n_starts = 5L,
                            seed = 1)
  expect_equal(Reduce(`+`, lapply(routed, \(n) n$weights)), grams$weights)
})

test_that("a network converted from a tna model counts only when its data agree", {
  skip_if_not_installed("tna")
  skip_if_not_installed("cograph")
  seqs <- data.frame(T1 = rep("A", 3L), T2 = rep("B", 3L), T3 = rep("C", 3L))
  converted <- cograph::as_cograph(tna::tna(seqs))
  # its probabilities of 1 are not counts; its sequences reproduce them, and
  # hold A -> B and B -> C three times
  expect_equal(unname(Nestimate:::.te_extract_counts(converted)),
               matrix(c(0, 0, 0, 3, 0, 0, 0, 3, 0), 3L))
  # a windowed model's conversion keeps no windows: its slots recount to other
  # weights, so it is refused rather than miscounted
  x <- .round9_units()
  windows <- tna::import_onehot(x, cols = c("A", "B", "C"), actor = "id",
                                window_size = 1L)
  lost <- cograph::as_cograph(tna::ftna(windows))
  expect_error(entropy_bayes(lost, draws = 100L, seed = 1),
               class = "nestimate_bad_input")
  expect_error(build_mmm(lost, k = 2L, seed = 1), class = "nestimate_bad_input")
})
