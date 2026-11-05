# One-hot data through association rules, clustering and resampling.
# Expected values are worked out by hand in the comments.

.hot_actors <- function() {
  # Three actors, rows in order: a = A B A B, b = C C B C, c = A C A A.
  data.frame(actor = rep(c("a", "b", "c"), each = 4L),
             A = c(1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1),
             B = c(0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0),
             C = c(0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0))
}

test_that("association rules on a one-hot network use codes and windows", {
  d <- data.frame(A = c(1, 0, 0, 0), B = c(0, 1, 0, 0), C = c(0, 0, 1, 0))
  tumbling <- suppressWarnings(build_network(d, method = "cna",
                                             codes = c("A", "B", "C"),
                                             window_size = 2))
  # Windows {A, B} and {C}, then an empty one has no rows left to form.
  parsed <- .ar_parse_input(tumbling)
  expect_identical(parsed$items, c("A", "B", "C"))
  expect_identical(nrow(parsed$matrix), 2L)
  expect_identical(unname(parsed$matrix[1L, ]), c(TRUE, TRUE, FALSE))
  # Overlapping: {A, B}, {B, C}, {C}.
  sliding <- suppressWarnings(build_network(d, method = "cna",
                                            codes = c("A", "B", "C"),
                                            window_size = 2,
                                            mode = "overlapping"))
  parsed <- .ar_parse_input(sliding)
  expect_identical(unname(colSums(parsed$matrix)), c(1, 2, 2))
  # Windows stop at the end of each actor: three actors of 4 rows, 2 windows each.
  by_actor <- build_network(.hot_actors(), method = "cna", actor = "actor",
                            codes = c("A", "B", "C"), window_size = 2)
  expect_identical(nrow(.ar_parse_input(by_actor)$matrix), 6L)
  # wtna() networks read the same way.
  expect_identical(.ar_parse_input(wtna(.hot_actors(), actor = "actor",
                                        method = "cooccurrence",
                                        codes = c("A", "B", "C"),
                                        window_size = 2))$items,
                   c("A", "B", "C"))
})

test_that("windows of .prepare_onehot() of different lengths are one run", {
  # Rows, split into windows of 2: id 1 {A}, {B} | {A}; id 2 {B}, {C} | {B};
  # id 3 {A}, {A} | {C}.
  d <- data.frame(id = rep(1:3, each = 3L),
                  A = c(1, 0, 1, 0, 0, 0, 1, 1, 0),
                  B = c(0, 1, 0, 1, 0, 1, 0, 0, 0),
                  C = c(0, 0, 0, 0, 1, 0, 0, 0, 1))
  wide <- .prepare_onehot(d, cols = c("A", "B", "C"), actor = "id",
                         window_size = 2)
  # Two rows of three codes, then one row: W0 has 6 columns and W1 has 3.
  expect_identical(ncol(wide), 9L)
  expect_null(.build_sequence_cols(wide))
  clusters <- build_clusters(wide, k = 2, method = "ward.D2")
  expect_identical(ncol(clusters$data), 9L)
  expect_true(isTRUE(attr(clusters$data, "windowed")))
  # Between the two windows: id 1 A->A, B->A; id 2 B->B, C->B; id 3 A->C twice.
  net <- build_network(wide, method = "frequency")
  expect_identical(ncol(net$data), 9L)
  expected <- matrix(c(1, 0, 2,
                       1, 1, 0,
                       0, 1, 0), 3L, byrow = TRUE,
                     dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_equal(net$weights[c("A", "B", "C"), c("A", "B", "C")], expected)
})

test_that("an overlapping window covers exactly window_size rows", {
  one <- data.frame(A = c(0, 0, 0, 1, 0, 0, 0, 0), B = 0)
  wide <- .prepare_onehot(one, cols = c("A", "B"), window_size = 3,
                         window_type = "overlapping")
  first_code <- grep("_T1$", names(wide), value = TRUE)
  # Eight rows give six windows, ending at rows 3 to 8; three hold row 4.
  expect_length(first_code, 6L)
  expect_identical(unname(unlist(wide[1L, first_code])),
                   c(NA, "A", "A", "A", NA, NA))
  # A group shorter than a window has no complete window.
  short <- .prepare_onehot(data.frame(id = c(1, 2, 2, 2), A = 1, B = 0),
                          cols = c("A", "B"), actor = "id", window_size = 3,
                          window_type = "overlapping")
  expect_identical(nrow(short), 1L)
})

test_that("resampling counts the windows the network counts", {
  d <- data.frame(id = rep(1:3, each = 3L),
                  A = c(1, 0, 1, 0, 0, 1, 1, 1, 0),
                  B = c(0, 1, 0, 1, 0, 0, 0, 1, 1),
                  C = c(0, 0, 0, 0, 1, 0, 1, 0, 0))
  wide <- .prepare_onehot(d, cols = c("A", "B", "C"), actor = "id",
                         window_size = 2)
  invisible(lapply(c("frequency", "co_occurrence"), function(method) {
    net <- build_network(wide, method = method)
    states <- net$nodes$label
    counts <- .precompute_per_sequence(net$data, method, net$params, states)
    expect_identical(nrow(counts), 3L)
    expect_equal(matrix(colSums(counts), 3L, 3L, byrow = TRUE),
                 unname(net$weights[states, states]), info = method)
  }))
  # With identical units every resample is the network itself.
  same <- data.frame(id = rep(1:2, each = 3L), A = c(1, 0, 1, 1, 0, 1),
                     B = c(0, 1, 0, 0, 1, 0))
  net <- build_network(.prepare_onehot(same, cols = c("A", "B"), actor = "id",
                                      window_size = 2),
                       method = "frequency")
  boot <- bootstrap_network(net, iter = 5, seed = 1)
  expect_equal(boot$mean, net$weights)
})

test_that("a one-hot network clusters and rebuilds by unit", {
  net <- build_network(.hot_actors(), method = "cna", actor = "actor",
                       codes = c("A", "B", "C"), window_size = 2)
  clusters <- build_clusters(net, k = 2, method = "ward.D2")
  # One row per actor holding states, not one row per time point holding 0 or 1.
  expect_identical(nrow(clusters$data), 3L)
  expect_false(any(c("0", "1") %in% unlist(clusters$data)))
  # The clusters partition the actors, so their networks add up to the whole.
  nets <- build_network(clusters, method = "co_occurrence")
  codes <- c("A", "B", "C")
  total <- matrix(0, 3L, 3L, dimnames = list(codes, codes))
  parts <- lapply(nets, function(n) {
    states <- rownames(n$weights)
    total[states, states] <<- total[states, states] + n$weights
  })
  expect_equal(total, net$weights[codes, codes])
  # Without a unit there is nothing to cluster.
  alone <- suppressWarnings(build_network(.hot_actors()[c("A", "B", "C")],
                                          method = "cna", window_size = 2))
  expect_error(build_clusters(alone, k = 2), class = "nestimate_bad_input")
  mmm <- suppressWarnings(build_mmm(net, k = 2, n_starts = 1, max_iter = 5,
                                    seed = 1))
  expect_identical(nrow(mmm$data), 3L)
})

test_that("resampling verbs run on a one-hot network with an actor", {
  d <- data.frame(actor = rep(paste0("p", 1:6), each = 6L),
                  A = rep(c(1, 0, 0), 12L), B = rep(c(0, 1, 0), 12L),
                  C = rep(c(0, 0, 1), 12L))
  d$A[c(4, 11, 20)] <- 0
  net <- build_network(d, method = "cna", actor = "actor",
                       codes = c("A", "B", "C"), window_size = 2)
  expect_no_condition(stability <- suppressWarnings(
    centrality_stability(net, iter = 3, seed = 1)))
  expect_s3_class(stability, "net_stability")
  expect_s3_class(network_reliability(net, iter = 3, seed = 1), "net_reliability")
  expect_s3_class(suppressWarnings(casedrop_reliability(net, iter = 3, seed = 1)),
                  "net_casedrop_reliability")
  # The units stay the actors.
  counts <- .precompute_per_sequence(
    .resampling_transition_data(net$data, net$metadata, net$params),
    "co_occurrence", net$params, net$nodes$label)
  expect_identical(nrow(counts), 6L)
})

test_that("one-hot data with no actor is one unit and cannot be resampled", {
  d <- data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
                  C = rep(c(0, 0, 1, 1), 6L))
  codes <- c("A", "B", "C")
  nets <- list(
    wtna(d),
    wtna(d, method = "cooccurrence"),
    build_network(d, method = "wtna", codes = codes),
    build_network(d, method = "cna", codes = codes, window_size = 2),
    build_network(d, method = "frequency", codes = codes),
    build_network(.prepare_onehot(d, codes, window_size = 2), method = "frequency")
  )
  invisible(lapply(nets, function(net) {
    expect_error(bootstrap_network(net, iter = 3, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(permutation(net, net, iter = 3, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(network_reliability(net, iter = 3, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(centrality_stability(net, iter = 3, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(casedrop_reliability(net, iter = 3, seed = 1),
                 class = "nestimate_single_sequence")
  }))
  mixed <- wtna(d, method = "both")
  expect_error(bootstrap_network(mixed, iter = 3, seed = 1),
               class = "nestimate_single_sequence")
  expect_error(permutation(mixed, mixed, iter = 3, seed = 1),
               class = "nestimate_single_sequence")

  # An actor column with one value is still one unit.
  d$actor <- 1L
  lone <- build_network(d, method = "cna", actor = "actor", codes = codes,
                        window_size = 2)
  expect_error(bootstrap_network(lone, iter = 3, seed = 1),
               class = "nestimate_single_sequence")

  # With units it runs.
  d$actor <- rep(1:3, each = 8L)
  grouped <- build_network(d, method = "cna", actor = "actor", codes = codes,
                           window_size = 2)
  expect_s3_class(bootstrap_network(grouped, iter = 3, seed = 1), "net_bootstrap")
})

test_that("one sequence cannot be resampled, on any verb or method", {
  long <- data.frame(action = c("A", "B", "C", "A", "B"))
  nets <- list(
    suppressMessages(build_tna(long, action = "action")),
    suppressMessages(build_ftna(long, action = "action")),
    suppressMessages(build_cna(long, action = "action")),
    suppressMessages(build_atna(long, action = "action")),
    suppressMessages(build_network(long, method = "ngram", action = "action")),
    suppressMessages(build_network(long, method = "gap", action = "action")),
    suppressMessages(build_network(long, method = "reverse", action = "action")),
    build_tna(data.frame(T1 = "A", T2 = "B", T3 = "C", T4 = "A"))
  )
  invisible(lapply(nets, function(net) {
    expect_error(bootstrap_network(net, iter = 2L, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(permutation(net, net, iter = 2L, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(network_reliability(net, iter = 2L, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(centrality_stability(net, iter = 2L, seed = 1),
                 class = "nestimate_single_sequence")
    expect_error(casedrop_reliability(net, iter = 2L, seed = 1),
                 class = "nestimate_single_sequence")
  }))
  # Building one stays legitimate, with a note.
  expect_message(build_tna(long, action = "action"),
                 class = "nestimate_single_sequence")

  # A group of one sequence is named.
  wide <- data.frame(T1 = c("A", "B", "A"), T2 = c("B", "A", "B"),
                     T3 = c("A", "A", "B"), g = c("x", "x", "y"))
  grouped <- build_network(wide, method = "relative", group = "g")
  expect_error(bootstrap_network(grouped, iter = 2L, seed = 1), "^y: ",
               class = "nestimate_single_sequence")
})

test_that("a column's name does not overrule the declared states", {
  # Eight rows of five positions: 8 x 4 = 32 transitions, whatever else is there.
  set.seed(2)
  wide <- as.data.frame(matrix(sample(c("x", "y", "z"), 40, TRUE), 8, 5,
                               dimnames = list(NULL, paste0("T", 1:5))))
  covariate <- rep(c("p", "q"), 4L)
  frame <- cbind(wide, Action = covariate)
  positions <- paste0("T", 1:5)
  reference <- build_network(wide, method = "frequency")$weights

  # Declared by state_cols, format or metadata_cols: the name is not read.
  by_cols <- build_network(frame, method = "frequency", state_cols = positions)
  expect_equal(by_cols$weights, reference)
  by_meta <- build_network(frame, method = "frequency", metadata_cols = "Action")
  expect_equal(by_meta$weights, reference)
  by_format <- suppressMessages(
    build_network(frame, method = "frequency", format = "wide"))
  expect_equal(by_format$weights, reference)
  expect_message(build_network(frame, method = "frequency", format = "wide"),
                 "aside", class = "nestimate_column_role")
  expect_no_message(build_network(frame, method = "frequency",
                                  state_cols = positions),
                    class = "nestimate_column_role")

  # Undeclared, it could be either: the caller says which.
  expect_error(build_network(frame, method = "frequency"),
               class = "nestimate_bad_input")
  # Stating the event column reads it as an event log.
  frame$Actor <- rep(1:2, each = 4L)
  explicit <- suppressMessages(build_network(
    frame, method = "frequency", action = "Action", actor = "Actor"))
  # p q p q within each actor: 3 transitions each.
  expect_equal(sum(explicit$weights), 6)
})

test_that("no role comes from a column's name", {
  log <- data.frame(Actor = rep(1:2, each = 4L), Action = rep(c("a", "b"), 4L))
  # Named Action and undeclared: the call asks, in any spelling.
  expect_error(build_network(log, method = "frequency", actor = "Actor"),
               class = "nestimate_bad_input")
  names(log)[2L] <- "ACTION"
  expect_error(build_network(log, method = "frequency", actor = "Actor"),
               class = "nestimate_bad_input")
  names(log)[2L] <- "Action"
  net <- build_network(log, method = "frequency", actor = "Actor",
                       action = "Action")
  expect_equal(sum(net$weights), 6)

  # Time is not used unless set, and the call says so. Rows reversed in each
  # actor: b a b a read in row order gives b->a 2 + a->b 1 per actor.
  log <- data.frame(Actor = rep(1:2, each = 4L), Time = rep(4:1, 2L),
                    Action = rep(c("b", "a", "b", "a"), 2L))
  expect_message(by_row <- build_network(log, method = "frequency",
                                         actor = "Actor", action = "Action"),
                 "`Time`", class = "nestimate_column_role")
  expect_equal(by_row$weights["b", "a"], 4)
  expect_equal(by_row$weights["a", "b"], 2)
  by_time <- build_network(log, method = "frequency", actor = "Actor",
                           action = "Action", time = "Time")
  expect_equal(by_time$weights["a", "b"], 4)
  expect_no_message(build_network(log, method = "frequency", actor = "Actor",
                                  action = "Action", metadata_cols = "Time"),
                    class = "nestimate_column_role")

  # frequencies() and magnitude_difference() read no default names either.
  log <- data.frame(T1 = c("a", "b"), T2 = c("b", "a"), Action = c("p", "q"))
  expect_equal(sum(frequencies(log, cols = c("T1", "T2"))), 2)
})

test_that("an undeclared identifier beside 0/1 columns stops the call", {
  set.seed(1)
  d <- data.frame(A = rbinom(30, 1, .4), B = rbinom(30, 1, .4),
                  C = rbinom(30, 1, .4), student = rep(1:3, each = 10L))
  # Read as wide, each of its values would be a node.
  expect_error(build_network(d, method = "relative"),
               class = "nestimate_bad_input")
  expect_error(build_network(d, method = "wtna"), class = "nestimate_bad_input")
  net <- build_network(d, method = "frequency", actor = "student")
  expect_identical(rownames(net$weights), c("A", "B", "C"))
  net <- build_network(d, method = "frequency", session = "student")
  expect_identical(rownames(net$weights), c("A", "B", "C"))
  # Other methods keep such variables, and numbered numeric positions are states.
  expect_true("student" %in% rownames(build_network(d, method = "cor")$weights))
  wide <- data.frame(T1 = c(0, 1, 0), T2 = c(1, 0, 1), T3 = c(2, 1, 0))
  expect_identical(rownames(build_network(wide, method = "frequency")$weights),
                   c("0", "1", "2"))
})

test_that("a character column left out of the sequence is announced once", {
  set.seed(2)
  wide <- as.data.frame(matrix(sample(c("x", "y", "z"), 40, TRUE), 8, 5,
                               dimnames = list(NULL, paste0("T", 1:5))))
  wide$final_state <- rep(c("x", "y"), 4L)
  seen <- character()
  net <- withCallingHandlers(
    build_network(wide, method = "frequency"),
    nestimate_column_role = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_length(seen, 1L)
  expect_match(seen, "`final_state`")
  # Declaring it counts it: 8 rows x 5 transitions.
  all_states <- build_network(wide, method = "frequency",
                              state_cols = names(wide))
  expect_equal(sum(all_states$weights), 40)
  # A declared covariate, and a numeric one, say nothing.
  expect_no_message(build_network(wide, method = "frequency",
                                  metadata_cols = "final_state"),
                    class = "nestimate_column_role")
  wide$final_state <- 1:8
  expect_no_message(build_network(wide, method = "frequency"),
                    class = "nestimate_column_role")
})

test_that("a resample is scaled the way the network was", {
  # Identical rows: every resample holds the same counts, so its scaled matrix
  # must equal the network's, zeros and all, on every path.
  d <- data.frame(T1 = rep("A", 4L), T2 = rep("B", 4L), T3 = rep("A", 4L),
                  T4 = rep("C", 4L), T5 = rep("A", 4L), T6 = rep("B", 4L))
  invisible(lapply(c("frequency", "relative", "co_occurrence", "attention"),
                   function(method) {
    invisible(lapply(c("minmax", "max", "rank"), function(scaling) {
      net <- build_network(d, method = method, scaling = scaling)
      boot <- bootstrap_network(net, iter = 3, seed = 1)
      expect_equal(boot$mean, net$weights,
                   info = paste(method, scaling))
    }))
  }))
  # Counts 0, 2, 4, 0 with zeros in the scale: 0, .5, 1, 0.
  expect_equal(as.vector(t(.postprocess_counts(c(0, 2, 4, 0), 2L, FALSE,
                                               "minmax", 0,
                                               include_zeros = TRUE))),
               c(0, 0.5, 1, 0))
  expect_true(.scaling_includes_zeros("attention"))
  expect_false(.scaling_includes_zeros("glasso"))
})

test_that("a 0/1 column stays a code, whatever it is called", {
  set.seed(1)
  binary <- function(names, n = 30L) {
    as.data.frame(stats::setNames(lapply(names, \(k) stats::rbinom(n, 1, 0.4)),
                                  names))
  }
  action <- suppressWarnings(build_network(binary(c("Action", "B", "C")),
                                           method = "cna", window_size = 2))
  expect_identical(rownames(action$weights), c("Action", "B", "C"))
  session <- suppressWarnings(build_network(binary(c("A", "B", "Session")),
                                            method = "cna", window_size = 2))
  expect_identical(rownames(session$weights), c("A", "B", "Session"))
  expect_null(session$params$actor)
  # A grouping column is used only when set.
  grouped <- binary(c("A", "B", "C"))
  grouped$session_id <- rep(1:3, each = 10L)
  expect_error(build_network(grouped, method = "cna", window_size = 2),
               class = "nestimate_bad_input")
  net <- build_network(grouped, method = "cna", window_size = 2,
                       session = "session_id")
  expect_identical(net$params$actor, "session_id")
})

test_that("one-hot values and identifiers are checked", {
  d <- data.frame(A = c(1, 0, 1, 0), B = c(0, 1, 0, 1), C = c(1, 1, 0, 0))
  codes <- c("A", "B", "C")
  counted <- d
  counted$A[1L] <- 2
  expect_error(build_network(counted, method = "cna", codes = codes,
                             window_size = 2),
               class = "nestimate_bad_input")
  labelled <- d
  labelled$C <- ifelse(labelled$C == 1, "Yes", "No")
  expect_error(build_network(labelled, method = "cna", codes = codes,
                             window_size = 2),
               class = "nestimate_bad_input")
  expect_error(wtna(counted, codes = codes), class = "nestimate_bad_input")
  # Missing identifiers are refused, not dropped.
  orphan <- d
  orphan$actor <- c("x", NA, "y", "y")
  expect_error(build_network(orphan, method = "cna", actor = "actor",
                             codes = codes, window_size = 2),
               class = "nestimate_bad_input")
  expect_error(wtna(orphan, actor = "actor", codes = codes),
               class = "nestimate_bad_input")
  # 0, 1, TRUE, FALSE and a missing cell stay valid; a missing cell is inactive.
  logical_d <- as.data.frame(lapply(d, as.logical))
  logical_d$A[2L] <- NA
  net <- suppressWarnings(build_network(logical_d, method = "cna", codes = codes,
                                        window_size = 1))
  expect_identical(unname(net$weights["A", "A"]), 2)
})

test_that("wtna networks resample their units and keep the time order", {
  codes <- c("A", "B", "C")
  # A B C A B C A B C in each of six units: 3 A->B, 3 B->C, 2 C->A per unit,
  # so every resample of six units gives the network exactly.
  d <- data.frame(A = rep(c(1, 0, 0), 18L), B = rep(c(0, 1, 0), 18L),
                  C = rep(c(0, 0, 1), 18L), unit = rep(1:6, each = 9L))
  frequency <- wtna(d, method = "transition", actor = "unit", window_size = 1L)
  expect_equal(unname(frequency$weights),
               rbind(c(0, 18, 0), c(0, 0, 18), c(12, 0, 0)))
  boot <- bootstrap_network(frequency, iter = 20, seed = 1)
  # Rows shuffled one by one gave a mean of about 6 for A->A here.
  expect_equal(unname(boot$mean), unname(frequency$weights))
  expect_equal(max(boot$sd), 0)

  nets <- list(
    wtna(d, method = "transition", type = "relative", actor = "unit",
         window_size = 1L),
    wtna(d, method = "transition", actor = "unit", window_size = 2L),
    wtna(d, method = "cooccurrence", actor = "unit", window_size = 2L),
    wtna(d, method = "cooccurrence", type = "relative", actor = "unit",
         window_size = 2L),
    wtna(d, method = "transition", actor = "unit", window_size = 2L,
         mode = "overlapping"),
    build_network(d, method = "wtna", codes = codes, actor = "unit",
                  window_size = 2L),
    build_network(d, method = "wtna_cooccurrence", codes = codes, actor = "unit",
                  window_size = 2L)
  )
  invisible(lapply(nets, function(net) {
    # Identical units: a resample on another scale or a broken window shows here.
    boot <- bootstrap_network(net, iter = 5, seed = 1)
    expect_equal(unname(boot$mean), unname(net$weights))
    expect_equal(max(boot$sd), 0)
    perm <- permutation(net, net, iter = 5, seed = 1)
    expect_equal(max(abs(perm$diff)), 0)
    expect_s3_class(network_reliability(net, iter = 3, seed = 1), "net_reliability")
    # All three codes have the same centrality here, and it warns so.
    expect_s3_class(suppressWarnings(suppressMessages(
      centrality_stability(net, iter = 3, seed = 1))), "net_stability")
    expect_s3_class(suppressWarnings(
      casedrop_reliability(net, iter = 3, seed = 1)), "net_casedrop_reliability")
  }))
})

test_that("overlapping means one count on every route", {
  # A B C A with window 2 gives the windows {A,B}, {B,C}, {C,A}. Windows 1-2 form
  # a block (A 1, B 2, C 1) and window 3 the next (A 1, C 1).
  d <- data.frame(A = c(1, 0, 0, 1), B = c(0, 1, 0, 0), C = c(0, 0, 1, 0))
  transitions <- wtna(d, window_size = 2L, mode = "overlapping")
  expect_equal(unname(transitions$weights),
               rbind(c(1, 0, 1), c(2, 0, 2), c(1, 0, 1)))
  # Within a block: (1,2,1) x (1,2,1) plus (1,0,1) x (1,0,1).
  cooccurrence <- wtna(d, method = "cooccurrence", window_size = 2L,
                       mode = "overlapping")
  expect_equal(unname(cooccurrence$weights),
               rbind(c(2, 2, 2), c(2, 4, 2), c(2, 2, 2)))

  set.seed(4)
  codes <- c("A", "B", "C")
  d <- as.data.frame(matrix(rbinom(3 * 60, 1, .4), ncol = 3,
                            dimnames = list(NULL, codes)))
  d$unit <- rep(1:6, each = 10L)
  invisible(lapply(c(2L, 3L, 5L), function(w) {
    prepared <- .prepare_onehot(d, codes, actor = "unit", window_size = w,
                                window_type = "overlapping")
    expect_equal(
      unname(wtna(d, codes = codes, actor = "unit", window_size = w,
                  mode = "overlapping")$weights),
      unname(build_network(prepared, method = "frequency")$weights[codes, codes]))
    expect_equal(
      unname(build_network(d, method = "cna", codes = codes, actor = "unit",
                           window_size = w, mode = "overlapping")$weights),
      unname(build_network(prepared, method = "cna")$weights[codes, codes]))
    # Per-unit counts add up to the network.
    net <- build_network(d, method = "frequency", codes = codes, actor = "unit",
                         window_size = w, mode = "overlapping")
    counts <- .precompute_per_sequence(
      .resampling_transition_data(net$data, net$metadata, net$params),
      "frequency", net$params, codes)
    expect_equal(sum(counts), sum(net$weights))
  }))

  # A group shorter than the window gives no windows.
  short <- wtna(d[1:3, codes], window_size = 5L, mode = "overlapping")
  expect_equal(sum(short$weights), 0)
})
