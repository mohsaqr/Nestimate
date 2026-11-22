# Refusal parity (layer 7): one degenerate or unsupported input passed through every verb that accepts that kind
# of object. Each cell is one call; its outcome is the condition class the call raised (the expected class when it
# is among the raised classes, else the first class) or the class of the result it returned. The outcome is compared
# as text with the class the verb's siblings raise on the same input, or with the result class when the input is
# legitimately supported. A wrong cell is recorded with `defect =`, which keeps the job green and fails once the
# cell is fixed (the marker is then stale). Nothing here is a reference computation: the authority is internal.

rf_reference <- "internal: refusal parity"

# The outcome of one call, as one string. `allow_warning` names warning classes the cell is known to raise (each is
# muffled, every other warning reaches testthat); `allow_message` the same for messages.
rf_outcome <- function(expr, expected, allow_warning = character(), allow_message = character()) {
  muffle_listed <- function(listed, restart) {
    function(condition) if (any(class(condition) %in% listed)) invokeRestart(restart)
  }
  caught <- tryCatch(
    withCallingHandlers(list(value = expr),
                        warning = muffle_listed(allow_warning, "muffleWarning"),
                        message = muffle_listed(allow_message, "muffleMessage")),
    error = \(condition) list(error = condition))
  if (is.null(caught$error)) return(class(caught$value)[[1L]])
  raised <- class(caught$error)
  hit <- intersect(expected, raised)
  if (length(hit)) hit[[1L]] else raised[[1L]]
}

# Assert and record one cell. `expected` is a condition class (`nestimate_...`) or the class of a legitimate result;
# `defect` describes a cell that does not give it (the disagreement is recorded, the test stays green).
rf_cell <- function(expr, expected, covers, input, verb = covers, defect = NULL, note = NULL,
                    allow_warning = character(), allow_message = character()) {
  got <- rf_outcome(expr, expected, allow_warning = allow_warning, allow_message = allow_message)
  expect_reference(got, expected, reference = rf_reference, authority = "internal", covers = covers,
                   quantity = sprintf("%s on %s", verb, input), tolerance = 0, defect = defect, note = note)
  invisible(got)
}

# Run a table of cells: one row per cell with columns `input`, `verb`, `expected`, and optionally `defect`, `note`,
# `allow_warning` (NA when none). `inputs` and `verbs` are named lists; each verb is a one-argument function and
# `covers` maps a verb name to the exported function(s) it calls.
rf_run <- function(cells, inputs, verbs, covers) {
  column <- function(name) if (name %in% names(cells)) cells[[name]] else rep(NA_character_, nrow(cells))
  missing_to_null <- function(value) if (is.na(value)) NULL else value
  stopifnot(
    "`cells` needs columns input, verb, expected" = all(c("input", "verb", "expected") %in% names(cells)),
    "every input of `cells` must be in `inputs`" = all(cells$input %in% names(inputs)),
    "every verb of `cells` must be in `verbs` and `covers`" = all(cells$verb %in% names(verbs)) &&
      all(cells$verb %in% names(covers))
  )
  invisible(Map(\(input, verb, expected, defect, note, allow_warning) {
    rf_cell(verbs[[verb]](inputs[[input]]), expected = expected, covers = covers[[verb]], input = input,
            verb = verb, defect = missing_to_null(defect), note = missing_to_null(note),
            allow_warning = if (is.na(allow_warning)) character() else strsplit(allow_warning, ",", fixed = TRUE)[[1L]])
  }, cells$input, cells$verb, cells$expected, column("defect"), column("note"), column("allow_warning")))
}

# A cell table from one verb and a named character vector input = expected; a later entry for an input replaces an
# earlier one, so a shared default vector can be specialised.
rf_cells <- function(verb, expected, defect = NULL, note = NULL, allow_warning = NULL) {
  expected <- expected[!duplicated(names(expected), fromLast = TRUE)]
  out <- data.frame(input = names(expected), verb = verb, expected = unname(expected), stringsAsFactors = FALSE)
  out$defect <- if (is.null(defect)) NA_character_ else unname(defect[out$input]) %||% NA_character_
  out$note <- if (is.null(note)) NA_character_ else unname(note[out$input])
  out$allow_warning <- if (is.null(allow_warning)) NA_character_ else unname(allow_warning[out$input])
  out
}

# ---- Inputs -----------------------------------------------------------------------------------------------------

# Wide sequences over four states, columns T1..Tn; `group` adds a column g splitting the rows in two halves.
rf_wide <- function(n = 30L, length = 10L, seed = 11L, states = c("A", "B", "C", "D"), group = FALSE) {
  set.seed(seed)
  wide <- as.data.frame(matrix(sample(states, n * length, replace = TRUE), n, length), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(length))
  if (group) wide$g <- rep(c("x", "y"), each = n %/% 2L, length.out = n)
  wide
}

# An event log of 12 chats of 25 messages with times, two teams: the input of method = "latency".
rf_event_log <- function(seed = 7L) {
  set.seed(seed)
  t0 <- as.POSIXct("2024-03-01", tz = "UTC")
  do.call(rbind, lapply(seq_len(12L), \(i) {
    data.frame(team = if (i <= 6L) "red" else "blue", chat = paste0("chat", i),
               sender = sample(c("A", "B", "C", "D"), 25L, replace = TRUE),
               at = t0 + cumsum(stats::rexp(25L, 1 / 60)), stringsAsFactors = FALSE)
  }))
}

# One-hot rows (time points) over codes A, B, C with no actor and no session: one unit.
rf_onehot <- function() {
  data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L), C = rep(c(0, 0, 1, 1), 6L))
}

# Numeric variables (V2 correlated with V1) and their 0/1 split: the inputs of the association methods.
rf_numeric <- function(seed = 3L) {
  set.seed(seed)
  num <- as.data.frame(matrix(stats::rnorm(200L * 4L), 200L, 4L))
  names(num) <- paste0("V", seq_len(4L))
  num$V2 <- num$V2 + num$V1
  num
}

# Nine states in three clusters with unequal state probabilities, so the within-cluster layers vary.
rf_states9 <- c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3")
rf_clusters9 <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"), C = c("c1", "c2", "c3"))
rf_wide9 <- function(seed = 21L) {
  set.seed(seed)
  probability <- sort(stats::runif(9L))
  wide <- as.data.frame(matrix(sample(rf_states9, 40L * 12L, replace = TRUE, prob = probability / sum(probability)),
                               40L, 12L), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(12L))
  wide
}

# A weighted edge list over the nine states: a cycle inside each cluster, edges between the clusters, one heavy edge.
rf_edges9 <- function() {
  data.frame(from = c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3", "a1", "b1", "c1", "a2"),
             to = c("a2", "a3", "a1", "b2", "b3", "b1", "c2", "c3", "c1", "b1", "c1", "a1", "a3"),
             w = c(rep(1, 12L), 50))
}

# A bare cograph_network: a weight matrix with nodes, no $data, no method.
rf_bare <- function(weights = matrix(c(0, .5, .3, .2, .4, 0, .4, .2, .1, .6, 0, .3, .3, .3, .4, 0), 4L, 4L,
                                     byrow = TRUE, dimnames = list(LETTERS[1:4], LETTERS[1:4]))) {
  structure(list(weights = weights,
                 nodes = data.frame(id = seq_len(nrow(weights)), label = rownames(weights), name = rownames(weights)),
                 directed = TRUE), class = "cograph_network")
}

# A three-level panel (20 people, 3 days of 10 beeps, 3 variables, person means apart): the input of build_mlvar().
rf_panel <- function(seed = 5L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(20L), \(person) {
    # a person mean per variable, so the random intercepts have variance and the between network is estimable
    series <- sweep(matrix(stats::rnorm(30L * 3L), 30L, 3L), 2L, stats::rnorm(3L, sd = 1.5), "+")
    colnames(series) <- paste0("V", seq_len(3L))
    data.frame(id = person, day = rep(seq_len(3L), each = 10L), beep = rep(seq_len(10L), 3L), series)
  }))
}

# The one-sequence note of a build is documented (class nestimate_single_sequence); only it is muffled.
rf_quiet_build <- function(expr) suppressMessages(expr, classes = "nestimate_single_sequence")
