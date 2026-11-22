# ---- tna-parity sequence estimators: ngram, gap, reverse ----
#
# These mirror `tna::build_model()` types "n-gram", "gap" and "reverse".
# All three operate on a wide sequence matrix (rows = sequences, columns =
# time points) and count state pairs at chosen column offsets. The shared
# kernel `.count_column_pairs()` does the tabulation; each estimator only
# decides which column pairs to count and with what weight.

#' Coerce sequence input (wide or long) to a wide character matrix
#'
#' @return List with `matrix` (character, NA-padded) and `states`.
#' @noRd
.sequence_matrix <- function(data, format = "auto", action = NULL,
                             id = NULL, time = NULL, cols = NULL,
                             alphabet = NULL, concat = 1L,
                             begin_state = NULL, end_state = NULL,
                             end_column = "always") {
  .assert("`data` must be a data frame" = is.data.frame(data))
  if (identical(format, "auto")) {
    format <- .format_from_action(action, data)
  }
  if (identical(format, "long")) {
    id_col <- if (is.null(id)) ".sequence" else id[1L]
    long <- as.data.frame(data)
    if (is.null(id)) long$.sequence <- 1L
    if (length(id) > 1L) {
      long[[".sequence"]] <- .group_label(long[, id, drop = FALSE])
      id_col <- ".sequence"
    }
    time_col <- intersect(time, names(long))
    if (length(time_col) == 0L) time_col <- NULL
    wide <- long_to_wide(long, id_col = id_col, time_col = time_col,
                         action_col = action, fill_na = TRUE)
    data <- wide
    id <- id_col
    cols <- NULL
  }
  .prepare_sequence_wide(data, id = id, cols = cols, alphabet = alphabet,
                         concat = concat, begin_state = begin_state,
                         end_state = end_state, end_column = end_column)
}

#' Count state pairs across a set of column offsets
#'
#' @param mat Character matrix (rows = sequences).
#' @param states Character vector of states (matrix dimnames).
#' @param from_cols,to_cols Integer vectors of equal length: column pairs.
#' @param pair_weight Numeric vector, one weight per column pair.
#' @param row_weight Optional numeric vector, one weight per row of `mat`.
#' @return Numeric `states x states` matrix of weighted pair counts.
#' @noRd
.count_column_pairs <- function(mat, states, from_cols, to_cols, pair_weight,
                                row_weight = NULL) {
  n_states <- length(states)
  nbins <- n_states * n_states
  empty <- matrix(0, nrow = n_states, ncol = n_states,
                  dimnames = list(states, states))
  if (n_states == 0L || length(from_cols) == 0L || nrow(mat) == 0L) {
    return(empty)
  }
  row_weight <- row_weight %||% rep(1, nrow(mat))
  contributions <- Map(function(fc, tc, w) {
    from <- match(mat[, fc], states)
    to <- match(mat[, tc], states)
    ok <- !is.na(from) & !is.na(to)
    if (!any(ok)) return(numeric(nbins))
    idx <- (from[ok] - 1L) * n_states + to[ok]
    sums <- rowsum(w * row_weight[ok], idx, reorder = FALSE)
    out <- numeric(nbins)
    out[as.integer(rownames(sums))] <- sums[, 1L]
    out
  }, from_cols, to_cols, pair_weight)
  counts <- Reduce(`+`, contributions, numeric(nbins))
  matrix(counts, nrow = n_states, ncol = n_states, byrow = TRUE,
         dimnames = list(states, states))
}

#' Assemble the estimator return list shared by the three estimators
#' @noRd
.sequence_estimator_result <- function(mat_counts, data, states, format,
                                       action, id, time, cols, concat,
                                       begin_state, end_state) {
  resolved_format <- if (identical(format, "auto")) {
    .format_from_action(action, data)
  } else format
  initial <- .compute_initial_probs(
    data, states, format = resolved_format,
    action = action, id = id, time = time, cols = cols, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  list(
    matrix = mat_counts,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = mat_counts,
    initial = initial
  )
}

#' n-gram estimator (tna type "n-gram")
#'
#' Every window of `n_gram` consecutive positions contributes each of its
#' adjacent pairs once, so interior pairs are counted up to `n_gram - 1`
#' times. Matches `tna::build_model(type = "n-gram")`.
#' @noRd
.estimator_ngram <- function(data, format = "auto", action = NULL,
                             id = NULL, time = NULL, cols = NULL,
                             alphabet = NULL, concat = 1L,
                             begin_state = NULL, end_state = NULL,
                             n_gram = 2L, ...) {
  .assert(
    "`n_gram` must be a single integer >= 2" =
      is.numeric(n_gram) && length(n_gram) == 1L && !is.na(n_gram) &&
        n_gram >= 2
  )
  n_gram <- as.integer(n_gram)
  # The windows run over the columns, so the end state gets a column of its
  # own only when a sequence has no room for it (tna's rule).
  prep <- .sequence_matrix(data, format, action, id, time, cols, alphabet,
                           concat, begin_state, end_state,
                           end_column = "needed")
  mat <- prep$matrix
  p <- ncol(mat)
  n_windows <- p - n_gram + 1L
  if (n_windows >= 1L) {
    # Pair (j, j+1) is counted once for each window that contains it.
    window_starts <- seq_len(n_windows)
    pair_starts <- unlist(lapply(window_starts, \(i) seq(i, i + n_gram - 2L)))
    multiplicity <- tabulate(pair_starts, nbins = p - 1L)
    from_cols <- which(multiplicity > 0L)
    counts <- .count_column_pairs(mat, prep$states, from_cols, from_cols + 1L,
                                  multiplicity[from_cols])
  } else {
    counts <- .count_column_pairs(mat, prep$states, integer(0), integer(0),
                                  numeric(0))
  }
  .sequence_estimator_result(counts, data, prep$states, format, action, id,
                             time, cols, concat, begin_state, end_state)
}

#' Gap estimator (tna type "gap")
#'
#' Pairs `(i, j)` with `1 <= j - i <= max_gap + 1` are counted with weight
#' `1 / (j - i)`. Matches `tna::build_model(type = "gap")`, where
#' `max_gap = 1` (default) allows one intervening state.
#' @noRd
.estimator_gap <- function(data, format = "auto", action = NULL,
                           id = NULL, time = NULL, cols = NULL,
                           alphabet = NULL, concat = 1L,
                           begin_state = NULL, end_state = NULL,
                           max_gap = 1L, ...) {
  .assert(
    "`max_gap` must be a single non-negative integer" =
      is.numeric(max_gap) && length(max_gap) == 1L && !is.na(max_gap) &&
        max_gap >= 0
  )
  max_gap <- as.integer(max_gap)
  prep <- .sequence_matrix(data, format, action, id, time, cols, alphabet,
                           concat, begin_state, end_state)
  mat <- prep$matrix
  p <- ncol(mat)
  distances <- seq_len(max_gap + 1L)
  pairs <- expand.grid(from = seq_len(max(p - 1L, 0L)), dist = distances)
  pairs$to <- pairs$from + pairs$dist
  pairs <- pairs[pairs$to <= p, , drop = FALSE]
  counts <- .count_column_pairs(mat, prep$states, pairs$from, pairs$to,
                                1 / pairs$dist)
  .sequence_estimator_result(counts, data, prep$states, format, action, id,
                             time, cols, concat, begin_state, end_state)
}

#' Reverse (reply-network) estimator (tna type "reverse")
#'
#' Each adjacent pair is counted from the later state to the earlier one:
#' the transpose of the frequency matrix. With `weighted = TRUE` each pair is
#' weighted by the inverse of its sequence length, as in tna.
#' @noRd
.estimator_reverse <- function(data, format = "auto", action = NULL,
                               id = NULL, time = NULL, cols = NULL,
                               alphabet = NULL, concat = 1L,
                               begin_state = NULL, end_state = NULL,
                               weighted = FALSE, ...) {
  .assert("`weighted` must be TRUE or FALSE" = isTRUE(weighted) ||
              isFALSE(weighted))
  prep <- .sequence_matrix(data, format, action, id, time, cols, alphabet,
                           concat, begin_state, end_state)
  mat <- prep$matrix
  p <- ncol(mat)
  from_cols <- seq_len(max(p - 1L, 0L)) + 1L
  row_weight <- if (weighted) 1 / rowSums(!is.na(mat)) else NULL
  counts <- .count_column_pairs(mat, prep$states, from_cols, from_cols - 1L,
                                rep(1, length(from_cols)), row_weight)
  .sequence_estimator_result(counts, data, prep$states, format, action, id,
                             time, cols, concat, begin_state, end_state)
}


# ---- Response-latency reply network ----

# Prefix of the time columns the latency estimator reads beside the states
# (build_network() adds them to `data` from the prepared event log).
.LATENCY_TIME_PREFIX <- ".latency_time_"
# Column of the per-sequence cap unit (actor or session) for latency_scope
# "actor" and "session", beside the time columns.
.LATENCY_UNIT_COL <- ".latency_unit"
.LATENCY_SCOPES <- c("global", "actor", "session", "responder")
# A unit with fewer sender-change latencies than this takes the global cap.
.LATENCY_MIN_UNIT_N <- 5L

#' Latency settings: validate floor, cap and aggregate
#'
#' @return List with `floor` (number), `cap` (a percentile in (0, 1] or
#'   NULL), `cap_seconds` (seconds or NULL) and `aggregate`.
#' @noRd
.latency_settings <- function(latency_floor = 0.5, latency_cap = 0.95,
                              latency_aggregate = "mean",
                              latency_scope = "global") {
  bad <- function(msg) {
    stop(errorCondition(msg, class = "nestimate_bad_input", call = NULL))
  }
  if (!is.numeric(latency_floor) || length(latency_floor) != 1L ||
      !is.finite(latency_floor) || latency_floor < 0 || latency_floor > 1) {
    bad("`latency_floor` must be a single number between 0 and 1.")
  }
  if (!is.character(latency_aggregate) || length(latency_aggregate) != 1L ||
      !latency_aggregate %in% c("sum", "mean")) {
    bad("`latency_aggregate` must be \"sum\" or \"mean\".")
  }
  if (!is.character(latency_scope) || length(latency_scope) != 1L ||
      !latency_scope %in% .LATENCY_SCOPES) {
    bad(paste0("`latency_scope` must be one of ",
               paste0("\"", .LATENCY_SCOPES, "\"", collapse = ", "), "."))
  }
  cap_msg <- paste0("`latency_cap` must be a percentile in (0, 1] (e.g. ",
                    "0.95) or a positive number of seconds written as a ",
                    "string (e.g. \"120s\").")
  if (is.character(latency_cap)) {
    if (length(latency_cap) != 1L || is.na(latency_cap) ||
        !grepl("^\\s*[0-9]*\\.?[0-9]+\\s*s\\s*$", latency_cap)) {
      bad(cap_msg)
    }
    seconds <- as.numeric(sub("^\\s*([0-9]*\\.?[0-9]+)\\s*s\\s*$", "\\1",
                              latency_cap))
    if (!is.finite(seconds) || seconds <= 0) bad(cap_msg)
    return(list(floor = latency_floor, cap = NULL, cap_seconds = seconds,
                aggregate = latency_aggregate, scope = latency_scope))
  }
  if (!is.numeric(latency_cap) || length(latency_cap) != 1L ||
      !is.finite(latency_cap) || latency_cap <= 0 || latency_cap > 1) {
    bad(cap_msg)
  }
  list(floor = latency_floor, cap = latency_cap, cap_seconds = NULL,
       aggregate = latency_aggregate, scope = latency_scope)
}

#' Split a latency frame into its state and time parts
#'
#' @return List with `states` (data frame), `times` (numeric matrix of
#'   seconds, same rows and as many columns as the state positions) and
#'   `units` (the cap unit of each row, or NULL).
#' @noRd
.latency_split <- function(data, cols = NULL) {
  time_cols <- names(data)[startsWith(names(data), .LATENCY_TIME_PREFIX)]
  if (length(time_cols) == 0L) {
    stop(errorCondition(paste0(
      "method = \"latency\" needs the event times: build it from an event ",
      "log with `action` and `time` set."),
      class = "nestimate_bad_input", call = NULL))
  }
  state_cols <- cols %||% setdiff(names(data), c(time_cols, .LATENCY_UNIT_COL))
  units <- if (.LATENCY_UNIT_COL %in% names(data)) {
    as.character(data[[.LATENCY_UNIT_COL]])
  }
  times <- as.matrix(data[, time_cols, drop = FALSE])
  storage.mode(times) <- "double"
  states <- data[, state_cols, drop = FALSE]
  if (ncol(times) != ncol(states)) {
    stop(errorCondition(sprintf(paste0(
      "method = \"latency\": %d state position(s) but %d time column(s)."),
      ncol(states), ncol(times)),
      class = "nestimate_bad_input", call = NULL))
  }
  list(states = states, times = times, units = units)
}

#' Lag-1 reply pairs of a sequence matrix with their latencies
#'
#' The pairs of `method = "reverse"`: position j+1 (the responder) to
#' position j (the previous sender), within each row.
#' @return Data frame with integer `from`, `to` (state indices), numeric
#'   `latency` (seconds) and the sequence `row`, one row per pair.
#' @noRd
.latency_pairs <- function(mat, times, states) {
  p <- ncol(mat)
  if (p < 2L || nrow(mat) == 0L) {
    return(data.frame(from = integer(0), to = integer(0),
                      latency = numeric(0), row = integer(0)))
  }
  later <- seq_len(p - 1L) + 1L
  responder <- match(as.vector(mat[, later, drop = FALSE]), states)
  previous <- match(as.vector(mat[, later - 1L, drop = FALSE]), states)
  latency <- as.vector(times[, later, drop = FALSE]) -
    as.vector(times[, later - 1L, drop = FALSE])
  ok <- !is.na(responder) & !is.na(previous)
  if (anyNA(latency[ok])) {
    stop(errorCondition(paste0(
      "method = \"latency\": an event has no time. Every event of the log ",
      "needs a timestamp."),
      class = "nestimate_bad_input", call = NULL))
  }
  row <- rep(seq_len(nrow(mat)), times = p - 1L)
  data.frame(from = responder[ok], to = previous[ok], latency = latency[ok],
             row = row[ok])
}

#' The reply pairs a latency network is built from
#'
#' Reads the states as `.estimator_latency()` does (cleaned, with
#' `alphabet`, when given, as the states; a state outside it forms no
#' pair), so the cap prepass of a grouped build and each group's estimator
#' select the same replies.
#' @param state_frame Data frame of the state positions.
#' @param times Numeric matrix of seconds, laid out as `state_frame`.
#' @return List with `states` (the alphabet) and `pairs` (from
#'   `.latency_pairs()`).
#' @noRd
.latency_reply_pairs <- function(state_frame, times, alphabet = NULL,
                                 id = NULL) {
  prep <- .sequence_matrix(state_frame, format = "wide", id = id,
                           alphabet = alphabet)
  list(states = prep$states,
       pairs = .latency_pairs(prep$matrix, times, prep$states))
}

#' Resolve a latency cap to seconds
#'
#' A percentile is taken over the latencies of the sender-change pairs
#' (responder differs from the previous sender), with quantile type 7.
#' @noRd
.latency_cap_seconds <- function(settings, latencies_changed) {
  if (!is.null(settings$cap_seconds)) return(settings$cap_seconds)
  if (length(latencies_changed) == 0L) {
    stop(errorCondition(paste0(
      "method = \"latency\": the data hold no reply from one sender to ",
      "another, so a percentile `latency_cap` is undefined. Give the cap in ",
      "seconds (e.g. latency_cap = \"120s\")."),
      class = "nestimate_bad_input", call = NULL))
  }
  cap <- stats::quantile(latencies_changed, probs = settings$cap,
                         names = FALSE, type = 7)
  if (!(cap > 0)) {
    stop(errorCondition(sprintf(paste0(
      "method = \"latency\": the %g percentile of the reply latencies is 0 ",
      "seconds. Give the cap in seconds (e.g. latency_cap = \"120s\")."),
      settings$cap),
      class = "nestimate_bad_input", call = NULL))
  }
  cap
}

#' The cap unit of each reply pair under a latency scope
#' @noRd
.latency_pair_units <- function(pairs, scope, units, states) {
  switch(scope,
    global = NULL,
    responder = states[pairs$from],
    {
      if (is.null(units)) {
        stop(errorCondition(sprintf(paste0(
          "latency_scope = \"%s\" needs the %s of each sequence; rebuild ",
          "the network with build_network(method = \"latency\")."),
          scope, scope), class = "nestimate_bad_input", call = NULL))
      }
      units[pairs$row]
    })
}

#' Per-unit latency caps
#'
#' The percentile of each unit's sender-change latencies. A unit with fewer
#' than `.LATENCY_MIN_UNIT_N` of them (or a percentile of 0 s) is left out,
#' and takes the global cap; their number is said in a message of class
#' `nestimate_latency_fallback`.
#' @param latency,changed,unit One value per reply pair.
#' @return Named numeric vector of caps in seconds.
#' @noRd
.latency_unit_caps <- function(settings, latency, changed, unit, scope) {
  all_units <- unique(unit)
  by_unit <- split(latency[changed], factor(unit[changed], levels = all_units))
  caps <- vapply(by_unit, function(v) {
    if (length(v) < .LATENCY_MIN_UNIT_N) return(NA_real_)
    cap <- stats::quantile(v, probs = settings$cap, names = FALSE, type = 7)
    if (cap > 0) cap else NA_real_
  }, numeric(1))
  n_fallback <- sum(is.na(caps))
  if (n_fallback > 0L) {
    cond <- simpleMessage(sprintf(paste0(
      "latency_scope = \"%s\": %d of %d unit(s) have fewer than %d replies ",
      "to another sender and use the global cap.\n"),
      scope, n_fallback, length(caps), .LATENCY_MIN_UNIT_N))
    class(cond) <- c("nestimate_latency_fallback", class(cond))
    message(cond)
  }
  caps <- caps[!is.na(caps)]
  # Sorted by unit, so the caps do not depend on the order units appear in
  caps[order(names(caps))]
}

#' Resolve the global cap and the per-unit caps of a set of reply pairs
#' @param pairs Data frame from `.latency_pairs()` plus a `unit` column
#'   (NULL for the global scope).
#' @return List with `global` (seconds) and `units` (named seconds or NULL).
#' @noRd
.latency_resolve_caps <- function(settings, latency, changed, unit) {
  global <- .latency_cap_seconds(settings, latency[changed])
  units <- if (!identical(settings$scope, "global") &&
               is.null(settings$cap_seconds)) {
    .latency_unit_caps(settings, latency, changed, unit, settings$scope)
  }
  list(global = global, units = units)
}

#' Response-latency estimator
#'
#' The reply pairs of `"reverse"` (responder -> previous sender, lag 1 within
#' each sequence), each weighted by how fast the reply came:
#' `w = 1 - (1 - latency_floor) * min(t / cap, 1)`, with `t` the seconds
#' between the two events. `latency_aggregate = "mean"` (default) averages
#' the weights of a pair, `"sum"` adds them. The diagonal (a sender continuing) is
#' weighted by the same rule. `latency_cap_seconds` (and, for a non-global
#' `latency_scope`, `latency_unit_caps`), when given, are the caps already
#' resolved on the full data and win over `latency_cap`.
#' @noRd
.estimator_latency <- function(data, format = "auto", action = NULL,
                               id = NULL, time = NULL, cols = NULL,
                               alphabet = NULL, concat = 1L,
                               begin_state = NULL, end_state = NULL,
                               weighted = FALSE,
                               latency_floor = 0.5, latency_cap = 0.95,
                               latency_aggregate = "mean",
                               latency_scope = "global",
                               latency_cap_seconds = NULL,
                               latency_unit_caps = NULL, ...) {
  .assert("`data` must be a data frame" = is.data.frame(data))
  if (!isFALSE(weighted)) {
    stop(errorCondition(paste0(
      "`weighted` (inverse sequence length) does not apply to method = ",
      "\"latency\", whose replies are weighted by their speed."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!identical(as.integer(concat %||% 1L), 1L) || !is.null(begin_state) ||
      !is.null(end_state)) {
    stop(errorCondition(paste0(
      "method = \"latency\" reads the times of consecutive events; `concat`, ",
      "`start` and `end` have no times and do not apply."),
      class = "nestimate_bad_input", call = NULL))
  }
  settings <- .latency_settings(latency_floor, latency_cap, latency_aggregate,
                                latency_scope)
  if (!is.null(latency_cap_seconds)) {
    .assert("`latency_cap_seconds` must be a single positive number" =
              is.numeric(latency_cap_seconds) &&
              length(latency_cap_seconds) == 1L &&
              is.finite(latency_cap_seconds) && latency_cap_seconds > 0)
  }
  parts <- .latency_split(data, cols)
  reply <- .latency_reply_pairs(parts$states, parts$times,
                                alphabet = alphabet, id = id)
  states <- reply$states
  pairs <- reply$pairs
  # A cap in seconds applies to every reply, whatever the scope
  scope <- if (is.null(settings$cap_seconds)) settings$scope else "global"
  unit <- .latency_pair_units(pairs, scope, parts$units, states)
  caps <- if (!is.null(latency_cap_seconds)) {
    # resolved when the network was built: held fixed under resampling
    list(global = latency_cap_seconds, units = latency_unit_caps)
  } else {
    settings$scope <- scope
    .latency_resolve_caps(settings, pairs$latency, pairs$from != pairs$to,
                          unit)
  }
  pair_cap <- rep(caps$global, nrow(pairs))
  if (!identical(scope, "global") && length(caps$units) > 0L) {
    unit_cap <- unname(caps$units[unit])
    pair_cap <- ifelse(is.na(unit_cap), pair_cap, unit_cap)
  }
  pairs$weight <- 1 - (1 - settings$floor) *
    pmin(pairs$latency / pair_cap, 1)

  n_states <- length(states)
  cell <- (pairs$from - 1L) * n_states + pairs$to
  to_matrix <- function(values) {
    out <- numeric(n_states * n_states)
    out[as.integer(names(values))] <- values
    matrix(out, n_states, n_states, byrow = TRUE,
           dimnames = list(states, states))
  }
  cell_f <- factor(cell, levels = sort(unique(cell)))
  n_replies <- to_matrix(vapply(split(pairs$weight, cell_f), length,
                                numeric(1)))
  weight_sum <- to_matrix(vapply(split(pairs$weight, cell_f), sum,
                                 numeric(1)))
  median_latency <- to_matrix(vapply(split(pairs$latency, cell_f),
                                     stats::median, numeric(1)))
  median_latency[n_replies == 0] <- NA_real_
  weights <- if (identical(settings$aggregate, "mean")) {
    ifelse(n_replies > 0, weight_sum / pmax(n_replies, 1), 0)
  } else {
    weight_sum
  }
  dimnames(weights) <- list(states, states)

  result <- .sequence_estimator_result(weights, parts$states, states, "wide",
                                       NULL, id, NULL, cols, 1L, NULL, NULL)
  # The weights are not counts: verbs that read transition counts
  # (bayes_compare(), entropy_bayes()) must not take them for counts.
  result$frequency_matrix <- NULL
  result$n_replies <- n_replies
  result$median_latency <- median_latency
  result$latency_cap_seconds <- caps$global
  result$latency_unit_caps <- caps$units
  result$latency_time <- as.data.frame(parts$times)
  if (!is.null(parts$units)) result$latency_time[[.LATENCY_UNIT_COL]] <- parts$units
  result
}

#' Event times of a prepared log as latency time columns (seconds)
#'
#' `prepare()` returns `time_data` aligned with `sequence_data`, one
#' `time_T<j>` column per position. Every event with a state needs a time:
#' `prepare()` sorts an event without one last and starts a session there,
#' which would drop its reply, so such events are refused here, by row.
#' @noRd
.latency_time_frame <- function(prepared) {
  .latency_stop_missing_times(prepared)
  td <- prepared$time_data
  n_pos <- ncol(prepared$sequence_data)
  if (!is.data.frame(td) || ncol(td) != n_pos ||
      nrow(td) != nrow(prepared$sequence_data)) {
    stop(errorCondition(paste0(
      "method = \"latency\": the event times could not be aligned with the ",
      "sequences. Check that `time` holds a timestamp for every event."),
      class = "nestimate_bad_input", call = NULL))
  }
  out <- as.data.frame(lapply(td, as.numeric))
  names(out) <- paste0(.LATENCY_TIME_PREFIX, seq_len(n_pos))
  out
}

#' Refuse a latency event log with an event that has no time
#'
#' An event counts when it has a state (an event without an action forms
#' no reply). Its row is read from `prepare()`'s `long_data`, whose row
#' names are those of the data given.
#' @noRd
.latency_stop_missing_times <- function(prepared) {
  long <- prepared$long_data
  if (!is.data.frame(long) || !".parsed_time" %in% names(long)) {
    return(invisible(NULL))
  }
  no_time <- is.na(long$.parsed_time)
  if (!any(no_time)) return(invisible(NULL))
  seq_row <- match(long$.session_id, prepared$meta_data$.session_id)
  state <- .clean_states(as.character(
    as.matrix(prepared$sequence_data)[cbind(seq_row, long$.sequence)]))
  missing <- no_time & !is.na(state)
  if (!any(missing)) return(invisible(NULL))
  rows <- rownames(long)[missing]
  shown <- paste0(paste(utils::head(rows, 5L), collapse = ", "),
                  if (length(rows) > 5L) ", ..." else "")
  stop(errorCondition(sprintf(paste0(
    "method = \"latency\" needs the time of every event: %d event(s) have ",
    "no timestamp (row(s) %s). Give them a time or remove them."),
    length(rows), shown),
    class = "nestimate_bad_input", call = NULL))
}

#' The latency caps of a grouped build
#'
#' Prepares each group's rows as its own network will, and resolves the
#' caps once over the reply pairs of all groups: the global cap pools every
#' group, a per-actor or per-responder cap pools that unit over the groups,
#' and a session is a unit of one group.
#' @param group_rows List of data frames, one per group.
#' @return A list with one element per group, each a list with `global`
#'   (seconds) and `units` (named seconds, or NULL).
#' @noRd
.latency_group_caps <- function(group_rows, params, action, actor, time,
                                session, order, time_threshold, timezone) {
  settings <- .latency_settings(params$latency_floor, params$latency_cap,
                                params$latency_aggregate,
                                params$latency_scope %||% "global")
  scope <- if (is.null(settings$cap_seconds)) settings$scope else "global"
  per_group <- lapply(seq_along(group_rows), function(i) {
    prepared <- .prepare_event_log(group_rows[[i]], params, action = action,
                                   actor = actor, time = time,
                                   session = session, order = order,
                                   time_threshold = time_threshold,
                                   timezone = timezone)
    times <- as.matrix(.latency_time_frame(prepared))
    # the replies each group's estimator reads: the same alphabet (the
    # caller's, or the one pooled over the groups) drops the same pairs
    reply <- .latency_reply_pairs(prepared$sequence_data, times,
                                  alphabet = .param_get(params, "alphabet"),
                                  id = .param_get(params, "id"))
    pairs <- reply$pairs
    units <- .latency_row_units(prepared, scope, actor)
    unit <- .latency_pair_units(pairs, scope, units, reply$states)
    data.frame(group = rep(i, nrow(pairs)), latency = pairs$latency,
               changed = pairs$from != pairs$to,
               unit = if (is.null(unit)) rep(NA_character_, nrow(pairs)) else
                 unit, stringsAsFactors = FALSE)
  })
  pooled <- do.call(rbind, per_group)
  global <- .latency_cap_seconds(settings, pooled$latency[pooled$changed])
  if (identical(scope, "global")) {
    return(rep(list(list(global = global, units = NULL)), length(group_rows)))
  }
  # a session belongs to one group: key it by group so labels cannot meet
  key <- if (identical(scope, "session")) {
    paste(pooled$group, pooled$unit, sep = "\r")
  } else {
    pooled$unit
  }
  unit_caps <- .latency_unit_caps(settings, pooled$latency, pooled$changed,
                                  key, scope)
  lapply(seq_along(group_rows), function(i) {
    caps <- if (identical(scope, "session")) {
      prefix <- paste0(i, "\r")
      own <- unit_caps[startsWith(names(unit_caps), prefix)]
      names(own) <- substring(names(own), nchar(prefix) + 1L)
      own
    } else {
      unit_caps
    }
    list(global = global, units = caps)
  })
}

#' The cap unit of each sequence of a prepared log
#'
#' The actor or the session label, one per row of `sequence_data`; NULL for
#' the global and responder scopes. Several actor columns are joined by
#' "-", with a "-" or a backslash inside a value escaped by a backslash, so
#' two different combinations never share a unit: ("a-b", "c") is a\-b-c
#' and ("a", "b-c") is a-b\-c. The label depends on the row alone, so the
#' groups of a grouped build name a unit alike.
#' @noRd
.latency_row_units <- function(prepared, scope, actor) {
  if (!scope %in% c("actor", "session")) return(NULL)
  md <- prepared$meta_data
  if (identical(scope, "actor")) {
    if (is.null(actor)) {
      stop(errorCondition(paste0(
        "latency_scope = \"actor\" needs `actor`, the column whose values ",
        "get their own cap."), class = "nestimate_bad_input", call = NULL))
    }
    if (length(actor) == 1L) return(as.character(md[[actor]]))
    escaped <- lapply(md[actor], function(v) {
      gsub("-", "\\-", gsub("\\", "\\\\", as.character(v), fixed = TRUE),
           fixed = TRUE)
    })
    return(do.call(paste, c(unname(escaped), sep = "-")))
  }
  as.character(md[[".session_label"]])
}

#' Refuse to resample a latency network
#'
#' Only bootstrap_network() resamples latency networks (validated by
#' simulation); the verbs that shuffle, drop or otherwise resample sequences
#' refuse a network built with `method = "latency"`, and a group that holds
#' one.
#' @param nets List of the verb's network inputs; anything that is not a
#'   `netobject` or `netobject_group` passes.
#' @param verb The verb's name, for the message.
#' @noRd
.stop_latency_resampling <- function(nets, verb) {
  is_latency <- function(obj) {
    if (inherits(obj, "netobject")) return(identical(obj$method, "latency"))
    if (inherits(obj, "netobject_group")) {
      return(any(vapply(obj, is_latency, logical(1))))
    }
    FALSE
  }
  if (!any(vapply(nets, is_latency, logical(1)))) return(invisible(NULL))
  stop(errorCondition(sprintf(paste0(
    "%s does not take a network built with method = \"latency\"; ",
    "bootstrap_network() resamples latency networks."), verb),
    class = c("nestimate_latency_resampling", "nestimate_bad_input"),
    call = NULL))
}
