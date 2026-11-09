# Sequence data conversion. This file holds frequencies() (long/wide ->
# transition count matrix) and convert_sequence_format() (long/wide ->
# frequency / one-hot / edgelist / follows), each with its own Rd topic.
#
# NOTE: do not reintroduce a file-header roxygen block with
# `@name frequencies`. roxygen merges it into frequencies.Rd, overwriting
# the function's own title with the header's and duplicating the
# description.

#' Build a Transition Frequency Matrix
#'
#' @description
#' `frequencies()` counts how often each state is followed by each other
#' state across all sequences and returns the counts as a square integer
#' matrix. It reads an event log (one row per event, `action` set) or a wide
#' table (one row per sequence, one column per position).
#'
#' @details
#' In an \strong{event log}, the sequences are defined by the `id`
#' column(s), and the events of each sequence are ordered by `time` (events
#' sharing a time keep their row order) or, without `time`, by row order.
#' Each pair of consecutive events is one transition.
#'
#' In a \strong{wide table}, each row is a sequence and the state columns are
#' consecutive positions; each pair of adjacent positions is one transition.
#'
#' In both, a missing state (`NA` or an empty string) breaks the sequence: no
#' transition is counted into or out of it, and the states on either side
#' are not joined. [convert_sequence_format()] closes such gaps instead.
#'
#' @param data A data frame of sequences: an event log or a wide table.
#' @param action The column holding each event's state; setting it declares
#'   an event log. Default `NULL`. Column names give no role: a column called
#'   `Action` with nothing declared raises an error of class
#'   `nestimate_bad_input`, and so does an `action` not in `data`.
#' @param id The column(s) identifying sequences. Default `NULL`. In an
#'   event log, each unique combination of their values is one sequence, and
#'   without `id` all rows form a single sequence. In a wide table, they are
#'   left out of the state columns. A column not in `data` raises an error of
#'   class `nestimate_bad_input`.
#' @param time The column that orders the events of each sequence in an
#'   event log. Default `NULL`: events keep their row order. A wide table
#'   does not use it and leaves it out of the positions. A column not in
#'   `data`, or one also listed in `cols`, raises an error of class
#'   `nestimate_bad_input`.
#' @param cols The state columns of a wide table; setting them declares a
#'   wide table. Default `NULL`. With neither `action` nor `cols`, the frame
#'   must be one column or hold position columns named as a numbered run
#'   (`T1`, `T2`, ..., or `1`, `2`, ...), with nothing beside them but the
#'   `id` (and `time`) columns; anything else raises an error of class
#'   `nestimate_bad_input`, as does a listed column not in `data`. A single
#'   column is read as one position per row, so it counts no transitions;
#'   declare an event log with `action`.
#' @param format How to read `data`. `"auto"` (default) reads an event log
#'   when `action` is set and a wide table otherwise. `"long"` forces an
#'   event log and needs `action`. `"wide"` forces a wide table and reads
#'   every column outside `id` and `time` (or the columns in `cols`) as a
#'   position, without the check on undeclared columns; it ignores `action`,
#'   so that column is counted as a position unless `cols` leaves it out.
#'   `"long"` ignores `cols`.
#'
#' @return A square integer matrix of transition counts, of class
#'   `c("nest_transition_counts", "matrix", "array")`, where `mat[i, j]` is
#'   the number of times state i was followed by state j; the row and column
#'   names are the sorted unique states. It behaves as an ordinary matrix and
#'   can be passed directly to `tna::tna()`. Its `summary()` returns the same
#'   counts as a data frame with one row per non-zero transition and columns
#'   `from`, `to` and `count`.
#'
#' @examples
#' # Wide format
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' freq <- frequencies(seqs, format = "wide")
#'
#' # Long format
#' long <- data.frame(
#'   Actor = rep(1:2, each = 3), Time = rep(1:3, 2),
#'   Action = c("A","B","C","B","A","C")
#' )
#' freq <- frequencies(long, action = "Action", id = "Actor")
#'
#' @seealso \code{\link{convert_sequence_format}} for converting to other
#'   representations (frequency counts, one-hot, edge lists).
#'
#' @export
frequencies <- function(data,
                        action = NULL,
                        id = NULL,
                        time = NULL,
                        cols = NULL,
                        format = c("auto", "long", "wide")) {
  .assert("`data` must be a data frame" = is.data.frame(data))
  .assert("`action` must be NULL or a single column name" =
            is.null(action) || (is.character(action) && length(action) == 1))
  .assert("`time` must be NULL or a single column name" =
            is.null(time) || (is.character(time) && length(time) == 1))
  .assert("`id` must be NULL or a character vector of column names" =
            is.null(id) || is.character(id))
  .assert("`cols` must be NULL or a character vector of column names" =
            is.null(cols) || is.character(cols))
  format <- match.arg(format)
  # without this the event-log reader stopped with "Action column 'NULL' not
  # found in data."
  if (identical(format, "long") && is.null(action)) {
    stop(errorCondition(paste0(
      "frequencies(): `format = \"long\"` reads an event log, so `action` ",
      "must name the column holding the states."),
      class = "nestimate_bad_input", call = NULL))
  }
  .check_declared_columns(data, "frequencies()", action = action, id = id,
                          time = time, cols = cols)
  if (is.null(action) && is.null(cols) && identical(format, "auto")) {
    .refuse_undeclared_roles(data, "frequencies()", taken = c(id, time),
                             whole_frame = TRUE)
  }
  # A wide table never reads a declared id or time column as a position.
  reads_wide <- identical(format, "wide") ||
    (identical(format, "auto") && is.null(action))
  if (reads_wide) {
    .refuse_role_positions(cols, c(id, time), "cols", "`id`/`time`")
    if (!is.null(time)) data <- data[setdiff(names(data), time)]
  }

  mat <- .count_transitions(
    data, format = format, action = action, id = id, time = time, cols = cols
  )
  if (is.matrix(mat)) {
    class(mat) <- c("nest_transition_counts", class(mat))
  }
  mat
}


#' @rdname frequencies
#' @param object For the `summary()` method: a `nest_transition_counts`
#'   matrix from `frequencies()`.
#' @param ... For the `summary()` method: ignored.
#' @export
summary.nest_transition_counts <- function(object, ...) {
  .matrix_to_long_df(unclass(object), value_col = "count",
                     include = "positive")
}


#' Convert Sequence Data to Different Formats
#'
#' @description
#' `convert_sequence_format()` turns wide or long sequence data into one of
#' four data frames: per-sequence state counts, per-sequence state presence
#' (one-hot), the list of consecutive transitions, or each event paired with
#' the one before it. It reads an event log (`action` set) or a wide table
#' (one row per sequence, one column per position). The counts can serve as
#' per-sequence behaviour profiles for the association estimators of
#' [build_network()] and for [build_mcml_pc()]. The id columns (those that
#' `id_col` names, and `.rid` for a wide table without `id_col`) are not
#' variables and must be set aside there (`id_col` in [build_mcml_pc()],
#' `actor` in [build_network()]); the `rid` column is
#' dropped automatically.
#'
#' @details
#' A missing state (`NA` or an empty string) is dropped and the states on
#' either side are joined, so a sequence `A, NA, B` gives the transition
#' `A -> B`. [frequencies()] counts no transition across a missing state
#' instead.
#'
#' @param data A data frame of sequences: an event log or a wide table.
#' @param seq_cols The state columns of a wide table; setting them declares
#'   a wide table. Default `NULL`: every column except `id_col` and `time`.
#'   With neither `action` nor `seq_cols`, the frame must be one column or
#'   hold position columns named as a numbered run (`T1`, `T2`, ..., or `1`,
#'   `2`, ...), with nothing beside them but the `id_col` (and `time`)
#'   columns; anything else raises an error of class `nestimate_bad_input`,
#'   as does a column listed both here and in `id_col` or `time`.
#' @param id_col The identifier column(s), copied to every output row. In an
#'   event log, each unique combination of their values is one sequence,
#'   and it is required: an event log without `id_col` raises an error of
#'   class `nestimate_bad_input` (no column is taken by position). In a wide
#'   table, each row is one sequence and no id is inferred from the data:
#'   without `id_col`, a row-index column `.rid` is added. Default `NULL`.
#' @param action The column holding each event's state; setting it declares
#'   an event log. Default `NULL`.
#' @param time The column that orders the events within each sequence of an
#'   event log; without it, events keep their row order. A wide table leaves
#'   it out of the state columns. Default `NULL`.
#' @param format The output format: `"frequency"` (default), the count of
#'   each state per sequence; `"onehot"`, whether each state occurs in the
#'   sequence (0/1); `"edgelist"`, each pair of consecutive states (`from`,
#'   `to`); `"follows"`, each event (`act`) paired with the event before it
#'   (`follows`).
#'
#' @return A data frame with the id column(s) first. For `"frequency"` and
#'   `"onehot"`: one row per sequence, a `rid` column identifying the
#'   sequence, and one integer column per state, holding counts or 0/1. For
#'   `"edgelist"`: one row per transition, with columns `from` and `to`. For
#'   `"follows"`: one row per event after the first of its sequence, with
#'   columns `act` and `follows`. An event log's rows come out sorted by
#'   `id_col`.
#'
#' @examples
#' # Wide format input
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' convert_sequence_format(seqs, format = "frequency")
#' convert_sequence_format(seqs, format = "edgelist")
#'
#' @seealso \code{\link{frequencies}} for building transition frequency matrices.
#'
#' @export
convert_sequence_format <- function(data,
                                    seq_cols = NULL,
                                    id_col = NULL,
                                    action = NULL,
                                    time = NULL,
                                    format = c("frequency", "onehot",
                                               "edgelist", "follows")) {
  .assert("`data` must be a data frame" = is.data.frame(data))
  format <- match.arg(format)

  .check_declared_columns(data, "convert_sequence_format()", action = action,
                          id_col = id_col, time = time, seq_cols = seq_cols)
  is_long <- !is.null(action)
  if (is.null(action) && is.null(seq_cols)) {
    .refuse_undeclared_roles(data, "convert_sequence_format()",
                             taken = c(id_col, time), wide = "`seq_cols`",
                             whole_frame = TRUE)
  }

  if (is_long && is.null(id_col)) {
    stop(errorCondition(
      paste("an event log needs its sequences declared: set `id_col` to the",
            "column(s) identifying them; no column is taken by position."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is_long) {
    .refuse_role_positions(seq_cols, c(id_col, time), "seq_cols",
                           "`id_col`/`time`")
  }

  if (is_long) {
    long <- .standardize_long(data, id_col, action, time)
    id_col <- long$id_col
    long <- long$data
  } else {
    if (is.null(id_col)) {
      # No id column declared. Never infer one from the data: value-based id
      # detection silently drops the first time point whenever a start-only
      # state does not recur later (audit A11-F01; Codex review). Synthesize a
      # row-index id (the same safe pattern as wide_to_long()) and treat every
      # column as a state, so no event is ever dropped. Callers with a genuine
      # leading id column must pass `id_col` explicitly.
      data[[".rid"]] <- seq_len(nrow(data))
      id_col <- ".rid"
    }
    if (is.null(seq_cols)) seq_cols <- setdiff(names(data), c(id_col, time))
    .validate_cols(data, c(id_col, seq_cols))
    long <- .pivot_wide_to_long(data, id_col, seq_cols)
  }

  switch(format,
    frequency = .fmt_frequency(long, id_col),
    onehot    = .fmt_onehot(long, id_col),
    edgelist  = .fmt_edgelist(long, id_col),
    follows   = .fmt_follows(long, id_col)
  )
}


# ---- Internal helpers for frequencies() ----

#' @noRd
.build_pairs_long <- function(data, action, id, time) {
  if (is.null(action) || !action %in% names(data)) {
    stop(errorCondition(paste0("Action column '", action %||% "NULL",
                               "' not found in data."),
                               class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(id)) {
    missing_ids <- setdiff(id, names(data))
    if (length(missing_ids) > 0) {
      stop(errorCondition(paste0("ID column(s) not found: ",
                                 paste(missing_ids, collapse = ", ")),
                          class = "nestimate_bad_input", call = NULL))
    }
  }

  if (is.null(id)) {
    groups <- list(data)
  } else if (length(id) == 1L) {
    groups <- split(data, data[[id]], drop = TRUE)
  } else {
    group_key <- interaction(data[, id, drop = FALSE], drop = TRUE)
    groups <- split(data, group_key, drop = TRUE)
  }

  pairs <- lapply(groups, function(g) {
    if (!is.null(time) && time %in% names(g)) {
      g <- g[order(g[[time]]), ]
    }
    actions <- g[[action]]
    n <- length(actions)
    if (n < 2L) {
      return(data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE))
    }
    from <- actions[-n]
    to <- actions[-1L]
    valid <- !is.na(from) & !is.na(to)
    data.frame(from = from[valid], to = to[valid], stringsAsFactors = FALSE)
  })

  do.call(rbind, pairs)
}

#' @noRd
.build_pairs_wide <- function(data, id, cols) {
  if (!is.null(cols)) {
    state_cols <- cols
  } else if (!is.null(id)) {
    state_cols <- setdiff(names(data), id)
  } else {
    state_cols <- names(data)
  }

  missing_cols <- setdiff(state_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(errorCondition(paste0("Columns not found: ",
                               paste(missing_cols, collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }
  if (length(state_cols) < 2L) {
    stop(errorCondition(
      "At least 2 state columns are required for wide format.",
      class = "nestimate_bad_input", call = NULL))
  }

  state_data <- data[, state_cols, drop = FALSE]

  pairs <- lapply(seq_len(nrow(state_data)), function(i) {
    row_vals <- unlist(state_data[i, ], use.names = FALSE)
    row_vals <- row_vals[!is.na(row_vals)]
    n <- length(row_vals)
    if (n < 2L) {
      return(data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE))
    }
    data.frame(from = row_vals[-n], to = row_vals[-1L],
               stringsAsFactors = FALSE)
  })

  do.call(rbind, pairs)
}


# ---- Internal helpers for convert_sequence_format() ----

#' Refuse a column declared both as a position and as another role
#' @noRd
.refuse_role_positions <- function(positions, roles, arg, role_args) {
  clash <- intersect(positions, roles)
  if (length(clash) > 0L) {
    stop(errorCondition(
      sprintf(paste("column(s) %s are declared both as positions (`%s`)",
                    "and as %s; a column has one role."),
              paste(clash, collapse = ", "), arg, role_args),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(NULL)
}

#' @noRd
.validate_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(errorCondition(paste0("Columns not found: ",
                               paste(missing, collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }
}

#' Standardize long format input for convert_sequence_format
#' @noRd
.standardize_long <- function(data, id_col, action, time) {
  .validate_cols(data, c(id_col, action))

  # Order by id columns + time
  order_cols <- id_col
  if (!is.null(time) && time %in% names(data)) {
    order_cols <- c(order_cols, time)
  }
  data <- data[do.call(order, data[, order_cols, drop = FALSE]), ]

  # Create sequence identifier
  if (length(id_col) == 1L) {
    rid_vals <- data[[id_col]]
  } else {
    rid_vals <- interaction(data[, id_col, drop = FALSE], drop = TRUE)
  }

  long <- data[, id_col, drop = FALSE]
  long$rid <- rid_vals
  long$act <- data[[action]]
  long <- long[!is.na(long$act) & long$act != "", ]

  list(data = long, id_col = id_col)
}

#' Pivot wide data to standardized long format
#' @noRd
.pivot_wide_to_long <- function(data, id_col, seq_cols) {
  rid_vals <- seq_len(nrow(data))
  parts <- lapply(seq_cols, function(col) {
    df <- data[, id_col, drop = FALSE]
    df$rid <- rid_vals
    df$act <- data[[col]]
    df
  })
  long <- do.call(rbind, parts)
  long <- long[!is.na(long$act) & long$act != "", ]
  as.data.frame(long)
}

#' @noRd
.fmt_frequency <- function(long, id_col) {
  tab <- table(long$rid, long$act)
  freq_df <- as.data.frame.matrix(tab)
  # `rid` is whatever identifies a sequence: an integer row index for wide
  # data, but the raw id value for long data -- often character, sometimes a
  # factor from a multi-column interaction. table() hands its levels back as
  # character rownames, so align on the character form. Coercing them with
  # as.integer() turned every character id into NA and the merge then
  # returned zero rows, silently.
  id_map <- unique(long[, c(id_col, "rid"), drop = FALSE])
  freq_df <- freq_df[match(as.character(id_map$rid), rownames(freq_df)), ,
                     drop = FALSE]
  result <- cbind(id_map, freq_df)
  rownames(result) <- NULL
  as.data.frame(result)
}

#' @noRd
.fmt_onehot <- function(long, id_col) {
  uniq <- unique(long[, c(id_col, "rid", "act"), drop = FALSE])
  all_acts <- sort(unique(uniq$act))
  id_map <- unique(uniq[, c(id_col, "rid"), drop = FALSE])
  mat <- matrix(0L, nrow = nrow(id_map), ncol = length(all_acts),
                dimnames = list(NULL, all_acts))
  # Loop kept: one indexed assignment per sequence, by action name; a
  # two-column index matrix would need its own rule for repeated row ids.
  for (i in seq_len(nrow(id_map))) {
    acts_i <- uniq$act[uniq$rid == id_map$rid[i]]
    mat[i, acts_i] <- 1L
  }
  result <- cbind(id_map, as.data.frame(mat))
  as.data.frame(result)
}

#' @noRd
.fmt_edgelist <- function(long, id_col) {
  groups <- split(long, long$rid)
  parts <- lapply(groups, function(g) {
    n <- nrow(g)
    if (n < 2L) {
      return(data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE))
    }
    edge <- g[seq_len(n - 1L), id_col, drop = FALSE]
    edge$from <- g$act[-n]
    edge$to <- g$act[-1L]
    edge
  })
  result <- do.call(rbind, parts)
  result <- result[, c(id_col, "from", "to"), drop = FALSE]
  rownames(result) <- NULL
  as.data.frame(result)
}

#' @noRd
.fmt_follows <- function(long, id_col) {
  groups <- split(long, long$rid)
  parts <- lapply(groups, function(g) {
    n <- nrow(g)
    if (n < 2L) {
      return(data.frame(act = character(0), follows = character(0),
                        stringsAsFactors = FALSE))
    }
    fol <- g[2:n, c(id_col), drop = FALSE]
    fol$act <- g$act[-1L]
    fol$follows <- g$act[-n]
    fol
  })
  result <- do.call(rbind, parts)
  result <- result[, c(id_col, "act", "follows"), drop = FALSE]
  rownames(result) <- NULL
  as.data.frame(result)
}
