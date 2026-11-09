#' Prepare Event Log Data for Network Estimation
#'
#' @description
#' `prepare()` turns an event log (one row per event: the action, and
#' optionally who did it and when) into sequences, one row per session and
#' one column per position (`T1`, `T2`, ...). It parses the timestamps,
#' orders the events, splits each actor's events into sessions by explicit
#' session columns and by time gaps, and keeps every other column as one
#' row of metadata per session. [build_network()] runs the same preparation
#' itself when given `action`, with the same arguments, so `prepare()` is
#' for inspecting the sequences or passing them to the sequence verbs.
#'
#' @param data The event log: a data frame with one row per event.
#' @param actor The column, or columns, identifying who performed each event
#'   (e.g. `"student"` or `c("student", "group")`); several columns are
#'   combined, each observed combination being one actor. Missing or `NULL`
#'   treats all events as one actor.
#' @param action The column holding each event's action (state or code).
#' @param time The timestamp column. Default `NULL`, in which case the row
#'   order (or `order`) gives the order of the events and no session is split
#'   by time. Numeric values are read as Unix timestamps, and POSIXct and
#'   POSIXlt values are used as they are, in their own time zone. Text may be
#'   ISO 8601 (with `Z` or an offset), `YYYY-MM-DD HH:MM[:SS]` with `-`, `/`
#'   or `.`, day-first or month-first dates, dates with month names, or
#'   compact `YYYYMMDDHHMMSS`, more than 40 layouts in all; text that matches
#'   none but reads as a number is taken as a Unix timestamp. A value must
#'   match a layout completely, with a four-digit year, so a two-digit year
#'   raises an error of class `nestimate_bad_input` instead of being misread,
#'   as does any value no layout reads. Day-first or month-first is decided
#'   once for the whole column: any first field above 12 means day-first, any
#'   second field above 12 means month-first. When the column cannot decide,
#'   values are read day-first with a warning of class
#'   `nestimate_ambiguous_dates`; `custom_format` settles it. An empty or
#'   missing timestamp is not an error: that event is sorted after its
#'   session's timed events and starts a session of its own.
#' @param order The column that orders the events: with `time`, it breaks
#'   ties between events with the same timestamp; without `time`, it gives
#'   the order of each session's events. Default `NULL` (the row order).
#' @param session The column, or columns, identifying explicit sessions
#'   (e.g. `"course"` or `c("course", "semester")`). Each actor's events in
#'   one session form a sequence, further split by time gaps when `time` is
#'   set. Default `NULL`.
#' @param time_threshold The longest gap, in seconds, between consecutive
#'   events of one session; a longer gap starts a new session, numbered
#'   `s1`, `s2`, ... within its actor and session. Used only with `time`.
#'   `FALSE` turns the split off, so each actor (or actor-session) forms one
#'   sequence however long the gaps. Default `900` (15 minutes).
#' @param custom_format A [strptime()] format tried before the built-in
#'   layouts; values it does not match fall through to them. Default `NULL`.
#' @param is_unix_time Whether `time` holds Unix timestamps. `TRUE` reads a
#'   text column as numbers and never as a date layout (a 12- or 14-digit
#'   value would otherwise match the compact `YYYYMMDDHHMM[SS]` layouts); a
#'   value that is not a number then raises an error of class
#'   `nestimate_bad_input`. With `FALSE` (default), text is matched to the
#'   date layouts first and read as a Unix timestamp only when it matches
#'   none. A numeric `time` column is always read as Unix timestamps.
#' @param unix_time_unit The unit of Unix timestamps: `"seconds"` (default),
#'   `"milliseconds"` or `"microseconds"`. A numeric `time` column that gives
#'   dates after the year 2500 raises a warning of class
#'   `nestimate_time_unit`, the usual sign of the wrong unit.
#' @param timezone The Olson time zone (see [OlsonNames()]) in which
#'   timestamps without zone information are read and Unix timestamps are
#'   expressed. Timestamps that end in `Z`, `UTC`, `GMT` or a numeric offset
#'   such as `+02:00` are converted from that offset, so parsing does not
#'   depend on the machine's time zone. Default `"UTC"`.
#'
#' @details
#' Within each actor and session, events are sorted by time and then by
#' `order`. Sessions are identified by the observed combinations of the
#' `actor` and `session` values, so identifiers containing separator
#' characters stay distinct and high-cardinality identifiers cannot
#' overflow. A missing value in any of these columns raises an error of
#' class `nestimate_bad_input`: drop or relabel those events first.
#'
#' Every other column is aggregated into the session metadata: numeric
#' columns by their mean, other columns by their most frequent value. A tie
#' goes to the tied value that occurs first in the session, and a message
#' names the columns where that happened.
#'
#' @return A `nestimate_data` object. `print()` gives the numbers of
#'   sessions and actions, the longest sequence, the number of actors when
#'   `actor` was given and whether times were parsed, and returns its input
#'   invisibly. [extract_sequences()], [sequence_patterns()] and
#'   [sequence_metrics()] read it directly. It holds:
#' \describe{
#'   \item{sequence_data}{The sequences as a data frame: one row per session,
#'     one column per position (`T1`, `T2`, ...), padded with `NA`.}
#'   \item{long_data}{The events in their sorted order, with the parsed time
#'     and the working columns that identify each event's session and
#'     position (`.session_id`, `.session_label`, `.sequence`, ...).}
#'   \item{meta_data}{One row per session, in the row order of
#'     `sequence_data`: `.session_id`, `.session_label`, the actor column
#'     (with several, their combined label `.actor`, joined by `"-"`, or by
#'     the first of `"|"`, `" | "`, `" / "`, `" :: "` that keeps distinct
#'     actors apart when values contain `"-"`), the session column(s)
#'     under their own names, and every other column aggregated per
#'     session.}
#'   \item{time_data}{The parsed times in the layout of `sequence_data`
#'     (columns `time_T1`, `time_T2`, ...), in `timezone` (or the zone of a
#'     date-time `time` column), or `NULL` without `time`.}
#'   \item{statistics}{A list with `total_sessions`, `total_actions` and
#'     `max_sequence_length`, plus `unique_actors` (distinct combinations of
#'     the `actor` columns) only when `actor` was given (without it, every
#'     event belongs to one actor).}
#' }
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   student = rep(1:3, each = 5),
#'   code = sample(c("read", "write", "test"), 15, replace = TRUE),
#'   timestamp = seq.POSIXt(as.POSIXct("2024-01-01"), by = "min", length.out = 15)
#' )
#' prepared <- prepare(df, actor = "student", action = "code",
#'                     time = "timestamp")
#' prepared
#'
#' # build_network() runs the same preparation itself
#' net <- build_network(df, actor = "student", action = "code",
#'                      time = "timestamp", method = "relative")
#' net
#'
#' @seealso \code{\link{build_network}}
#'
#' @export
prepare <- function(data,
                         actor,
                         action,
                         time = NULL,
                         order = NULL,
                         session = NULL,
                         time_threshold = 900,
                         custom_format = NULL,
                         is_unix_time = FALSE,
                         unix_time_unit = c("seconds", "milliseconds",
                                            "microseconds"),
                         timezone = "UTC") {
  .assert(is.data.frame(data))
  .check_timezone(timezone)
  .assert(is.character(action), length(action) == 1, action %in% names(data))
  # FALSE disables gap-based splitting. Inf makes the `gaps > time_threshold`
  # test below never fire, so no separate code path is needed.
  if (isFALSE(time_threshold)) {
    time_threshold <- Inf
  }
  .assert(is.numeric(time_threshold), length(time_threshold) == 1,
          !is.na(time_threshold), time_threshold > 0)
  unix_time_unit <- match.arg(unix_time_unit)

  df <- as.data.frame(data)
  n <- nrow(df)

  # ---- Actor ----
  default_actor <- FALSE
  if (missing(actor) || is.null(actor)) {
    df$.actor <- "all"
    actor_col <- ".actor"
    default_actor <- TRUE
  } else {
    .assert(is.character(actor), all(actor %in% names(df)))
    if (length(actor) > 1L) {
      # Display label only. Identity comes from the underlying columns below.
      df$.actor <- .distinct_label(df[, actor, drop = FALSE])
      actor_col <- ".actor"
    } else {
      actor_col <- actor
    }
  }
  actor_key_cols <- if (default_actor) ".actor" else actor

  # ---- Session (explicit grouping) ----
  if (!is.null(session)) {
    .assert(is.character(session), all(session %in% names(df)))
    if (length(session) > 1L) {
      df$.session_explicit <- .distinct_label(df[, session, drop = FALSE])
    } else {
      df$.session_explicit <- df[[session]]
    }
  }

  # ---- Order (tiebreaker) ----
  if (is.null(order)) {
    df$.order <- seq_len(n)
    order_col <- ".order"
  } else {
    .assert(is.character(order), length(order) == 1, order %in% names(df))
    order_col <- order
  }

  # ---- Build base grouping key (actor + session) ----
  # Identity is an integer over observed combinations of the *original*
  # columns, so it cannot overflow and no separator can collide. The label is
  # kept separately for display.
  key_cols <- df[, c(actor_key_cols, session), drop = FALSE]
  df$.base_group <- .observed_group_id(key_cols, context = "actor/session")
  df$.base_label <- if (is.null(session)) {
    as.character(df[[actor_col]])
  } else {
    .distinct_label(list(actor = df[[actor_col]],
                         session = df$.session_explicit),
                    seps = c(" | ", " || ", " / ", " :: "))
  }

  # ---- Time parsing + inferred session detection ----
  if (!is.null(time)) {
    .assert(is.character(time), length(time) == 1, time %in% names(df))

    # Auto-detect numeric as unix
    if (is.numeric(df[[time]])) {
      is_unix_time <- TRUE
    }

    if (isTRUE(is_unix_time) && !is.numeric(df[[time]])) {
      # Declared Unix time: read the text as numbers, never as a date layout
      # (a 12- or 14-digit value also matches the compact YYYYMMDDHHMM[SS]).
      as_num <- suppressWarnings(as.numeric(trimws(df[[time]])))
      bad <- !is.na(df[[time]]) & nzchar(trimws(df[[time]])) & is.na(as_num)
      if (any(bad)) {
        stop(errorCondition(sprintf(paste(
          "`is_unix_time = TRUE`, but `%s` holds values that are not",
          "numbers, e.g. %s."), time,
          paste(utils::head(df[[time]][bad], 3), collapse = ", ")),
          class = "nestimate_bad_input", call = NULL))
      }
      df[[time]] <- as_num
    }
    parsed <- .parse_time(df[[time]], custom_format = custom_format,
                          is_unix_time = is_unix_time,
                          unix_time_unit = unix_time_unit,
                          timezone = timezone)
    df$.parsed_time <- parsed

    # Numeric times in the wrong unit parse "successfully" to dates millennia
    # away and turn every event into its own session. Say so.
    if (is.numeric(df[[time]]) && any(!is.na(parsed)) &&
        max(parsed, na.rm = TRUE) > as.POSIXct("2500-01-01", tz = "UTC")) {
      warning(warningCondition(sprintf(paste0(
        "`%s` read as Unix %s gives dates after the year 2500 (latest %s). ",
        "Set `unix_time_unit` (\"milliseconds\" or \"microseconds\")."),
        time, unix_time_unit, format(max(parsed, na.rm = TRUE), "%Y")),
        class = "nestimate_time_unit", call = NULL))
    }

    # Sort by base_group + time + order
    df <- df[base::order(df$.base_group, df$.parsed_time, df[[order_col]]), ]

    # Infer sub-sessions from time gaps within each base group
    df$.inferred_nr <- ave(
      as.numeric(df$.parsed_time),
      df$.base_group,
      FUN = function(t) {
        gaps <- c(NA_real_, diff(t))
        new_session <- is.na(gaps) | gaps > time_threshold
        cumsum(new_session)
      }
    )

    # Both parts are integers, so this key is unambiguous.
    df$.session_id <- paste0(df$.base_group, " s", df$.inferred_nr)
    df$.session_label <- paste0(df$.base_label, " s", df$.inferred_nr)

  } else {
    # No time: sort by base_group + order
    df <- df[base::order(df$.base_group, df[[order_col]]), ]

    # Each base group = one session
    df$.session_id <- as.character(df$.base_group)
    df$.session_label <- df$.base_label
  }

  # ---- Sequence numbering within sessions ----
  df$.sequence <- ave(seq_len(nrow(df)), df$.session_id, FUN = seq_along)

  # ---- Pivot to wide ----
  sessions <- unique(df$.session_id)
  max_len <- max(df$.sequence)

  # Build wide sequence data via matrix indexing (fast)
  seq_mat <- matrix(NA_character_, nrow = length(sessions), ncol = max_len)
  session_idx <- match(df$.session_id, sessions)
  seq_mat[cbind(session_idx, df$.sequence)] <- as.character(df[[action]])

  sequence_data <- as.data.frame(seq_mat, stringsAsFactors = FALSE)
  names(sequence_data) <- paste0("T", seq_len(max_len))

  # Build wide time data (if time was provided)
  if (!is.null(time)) {
    time_mat <- matrix(NA_real_, nrow = length(sessions), ncol = max_len)
    time_mat[cbind(session_idx, df$.sequence)] <- as.numeric(df$.parsed_time)
    # Convert every column first, then build the frame once. Assigning into a
    # data.frame column-by-column copies the whole frame on each iteration,
    # which is quadratic in max_len -- 27k+ columns for a single long sequence.
    # Keep the zone of the parsed times (`timezone`, or a POSIXct input's
    # own), so the display does not depend on the machine's local zone.
    time_tz <- attr(df$.parsed_time, "tzone")
    if (is.null(time_tz) || !nzchar(time_tz[1L])) time_tz <- timezone
    time_cols <- lapply(
      seq_len(max_len),
      function(j) as.POSIXct(time_mat[, j], origin = "1970-01-01",
                             tz = time_tz[1L])
    )
    names(time_cols) <- paste0("time_T", seq_len(max_len))
    time_data <- as.data.frame(time_cols, stringsAsFactors = FALSE)
  } else {
    time_data <- NULL
  }

  # Meta data. Keyed on the unambiguous .session_id; .session_label carries the
  # readable actor/session name for reporting.
  meta_data <- data.frame(.session_id = sessions, stringsAsFactors = FALSE)
  label_map <- df$.session_label[match(sessions, df$.session_id)]
  meta_data$.session_label <- label_map
  if (!default_actor) {
    actor_map <- df[!duplicated(df$.session_id),
                    c(".session_id", actor_col), drop = FALSE]
    meta_data <- merge(meta_data, actor_map, by = ".session_id", sort = FALSE)
  }
  # The session column(s) keep their own names, so a sequence can be traced
  # back to its session without parsing .session_label. Every row of a
  # session carries the same value, so the first row is the session's value.
  if (!is.null(session)) {
    session_map <- df[!duplicated(df$.session_id),
                      c(".session_id", session), drop = FALSE]
    meta_data <- merge(meta_data, session_map, by = ".session_id", sort = FALSE)
  }

  # Aggregate extra columns per session
  special_cols <- c(action, actor_col, time, order_col, session,
                    grep("^\\.", names(df), value = TRUE))
  extra_cols <- setdiff(names(df), special_cols)
  if (length(extra_cols) > 0L) {
    agg <- .aggregate_metadata(df, session_col = ".session_id",
                               extra_cols = extra_cols)
    meta_data <- merge(meta_data, agg, by = ".session_id", sort = FALSE)
  }
  # merge() leaves the row order unspecified; row i must describe sequence i.
  meta_data <- meta_data[match(sessions, meta_data$.session_id), , drop = FALSE]
  rownames(meta_data) <- NULL

  # Statistics
  stats <- list(
    total_sessions = length(sessions),
    total_actions = nrow(df),
    max_sequence_length = max_len
  )
  if (!default_actor) {
    stats$unique_actors <- nrow(unique(df[, actor, drop = FALSE]))
  }

  rownames(sequence_data) <- NULL

  structure(
    list(
      sequence_data = sequence_data,
      long_data = df,
      meta_data = meta_data,
      time_data = time_data,
      statistics = stats
    ),
    class = "nestimate_data"
  )
}


#' @rdname prepare
#' @param x For the `print()` method: a `nestimate_data` object.
#' @param ... Not used.
#' @export
print.nestimate_data <- function(x, ...) {
  cat("Prepared Data for Network Estimation\n")
  s <- x$statistics
  cat(sprintf("  Sessions: %d  |  Actions: %d  |  Max length: %d\n",
              s$total_sessions, s$total_actions, s$max_sequence_length))
  if (!is.null(s$unique_actors)) {
    cat(sprintf("  Actors: %d\n", s$unique_actors))
  }
  if (!is.null(x$time_data)) {
    cat("  Time data: available\n")
  }
  invisible(x)
}


# ---- Time parsing ----

#' Validate an Olson time zone name
#' @noRd
.check_timezone <- function(timezone) {
  .assert(
    "`timezone` must be a single Olson time zone name (see OlsonNames())" =
      is.character(timezone) && length(timezone) == 1L &&
        !is.na(timezone) && timezone %in% OlsonNames()
  )
  invisible(timezone)
}

#' Parse character timestamps element-wise through a list of formats
#'
#' Each format is tried only on the elements still unparsed, so a column
#' mixing several layouts is resolved completely instead of stopping at the
#' first format that matches any element.
#'
#' @param time Character vector.
#' @param formats Character vector of \code{strptime} formats, tried in order.
#' @param timezone Olson time zone name.
#' @param strict Require each value to match the whole of a format with a
#'   four-digit year. Off for a user's own `custom_format`.
#' @return POSIXct vector (NA where no format matched).
#' @noRd
.parse_time_formats <- function(time, formats, timezone, strict = FALSE) {
  out <- rep(NA_real_, length(time))
  # Sequential fill: each format only sees what earlier formats left NA, so
  # the loop is over formats (short), not elements.
  for (fmt in formats) {
    idx <- which(is.na(out) & !is.na(time))
    if (length(idx) == 0L) break
    # strptime ignores trailing characters and %Y takes 1 to 4 digits, so
    # "%Y/%m/%d" alone reads 15/01/2024 as the year 15.
    if (strict) idx <- idx[grepl(.strptime_regex(fmt), time[idx])]
    if (length(idx) == 0L) next
    parsed <- as.POSIXct(strptime(time[idx], format = fmt, tz = timezone))
    ok <- !is.na(parsed)
    out[idx[ok]] <- as.numeric(parsed[ok])
  }
  as.POSIXct(out, origin = "1970-01-01", tz = timezone)
}

#' Regular expression for the complete strings a strptime format may accept
#'
#' Fields run together ("%Y%m%d%H%M%S") have a fixed width, otherwise a Unix
#' time written as a string could pass for a compact date.
#' @param fmt A single strptime format.
#' @return A single anchored regular expression.
#' @noRd
.strptime_regex <- function(fmt) {
  field <- if (grepl("[^%A-Za-z]", fmt)) "[0-9]{1,2}" else "[0-9]{2}"
  seconds <- paste0(field, "(\\.[0-9]+)?")
  tokens <- c("%Y" = "[0-9]{4}", "%m" = field, "%d" = field, "%H" = field,
              "%M" = field, "%OS" = seconds, "%S" = seconds,
              "%b" = "[[:alpha:]]+", "%B" = "[[:alpha:]]+",
              "%z" = "\\s*[+-][0-9]{4}")
  literal <- gsub(" ", "\\s+", gsub(".", "\\.", fmt, fixed = TRUE),
                  fixed = TRUE)
  body <- Reduce(\(pattern, token) gsub(token, tokens[[token]], pattern,
                                        fixed = TRUE),
                 names(tokens), init = literal)
  paste0("^", body, "$")
}

#' Decide once per column whether numeric dates are day-first or month-first
#'
#' Trying each value separately would read a month-first column day-first
#' wherever the day is 12 or less, giving one column two layouts.
#' @param time Character vector of naive time strings.
#' @param formats The candidate strptime formats, day-first before month-first.
#' @return `formats`, without the day-first numeric ones for a month-first
#'   column. Warns when the column cannot decide.
#' @noRd
.numeric_date_formats <- function(time, formats) {
  pattern <- "^([0-9]{1,2})[-/.]([0-9]{1,2})[-/.][0-9]{4}([ T].*)?$"
  dated <- time[grepl(pattern, time)]
  if (length(dated) == 0L) return(formats)
  first <- as.integer(sub(pattern, "\\1", dated))
  second <- as.integer(sub(pattern, "\\2", dated))
  day_first <- any(first > 12L)
  month_first <- any(second > 12L)
  if (month_first && !day_first) {
    return(formats[!grepl("^%d[-/.]%m[-/.]%Y", formats)])
  }
  undecided <- first <= 12L & second <= 12L & first != second
  if ((day_first && month_first && any(undecided)) ||
      (!day_first && !month_first && any(undecided))) {
    warning(warningCondition(
      sprintf(paste0("%d time value(s) could be day-first or month-first ",
                     "(e.g. %s); read as day-first. Set `custom_format` to ",
                     "settle it, e.g. \"%%m/%%d/%%Y\"."),
              sum(undecided), dated[which(undecided)[1L]]),
      class = "nestimate_ambiguous_dates", call = NULL
    ))
  }
  formats
}

#' Parse time values into POSIXct, independent of the machine time zone
#'
#' Ported from \code{tna:::parse_time()}. Handles POSIXct input, numeric Unix
#' timestamps, ISO-8601 strings with \code{Z}/\code{UTC}/\code{GMT} markers or
#' numeric offsets, a wide set of naive layouts, and a numeric-string Unix
#' fallback. Naive strings are interpreted in \code{timezone}; offset strings
#' are converted from their offset.
#'
#' @param time Vector of time values.
#' @param custom_format Character or NULL. Custom strptime format tried first.
#' @param is_unix_time Logical. Treat numeric as Unix timestamp.
#' @param unix_time_unit Character. "seconds", "milliseconds", "microseconds".
#' @param timezone Olson time zone name.
#' @return POSIXct vector in \code{timezone}.
#' @noRd
.parse_time <- function(time, custom_format = NULL, is_unix_time = FALSE,
                        unix_time_unit = "seconds", timezone = "UTC") {
  .check_timezone(timezone)
  unix_divisor <- switch(unix_time_unit,
    seconds = 1, milliseconds = 1000, microseconds = 1e6,
    stop("'unix_time_unit' must be seconds, milliseconds or microseconds.",
         call. = FALSE)
  )
  as_unix <- function(x) {
    as.POSIXct(x / unix_divisor, origin = "1970-01-01", tz = timezone)
  }

  if (inherits(time, c("POSIXct", "POSIXlt"))) {
    return(as.POSIXct(time))
  }
  if (is.numeric(time) && is_unix_time) {
    return(as_unix(time))
  }

  time_original <- time
  time <- trimws(as.character(time))
  time_empty <- is.na(time) | !nzchar(time)
  time[time_empty] <- NA_character_
  parsed <- as.POSIXct(rep(NA_real_, length(time)), origin = "1970-01-01",
                       tz = timezone)

  # 1. Custom format first; unmatched values fall through.
  if (!is.null(custom_format)) {
    custom_parsed <- .parse_time_formats(time, custom_format, timezone)
    ok <- !is.na(custom_parsed)
    parsed[ok] <- custom_parsed[ok]
  }

  # 2. Offset-bearing strings: normalise Z/UTC/GMT and "+02:00" to "+0200"
  #    so %z parses them consistently, then convert from that offset.
  time_offset <- sub("(?:Z|UTC|GMT)$", "+0000", time,
                     ignore.case = TRUE, perl = TRUE)
  time_offset <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", time_offset,
                     perl = TRUE)
  has_timezone <- grepl(
    # A numeric offset follows a time; on its own, "-2024" ends DD-MM-YYYY.
    "(?:Z|UTC|GMT|:[0-9]{2}(?:\\.[0-9]+)?\\s*[+-][0-9]{2}:?[0-9]{2}|[A-Za-z]{2,})$",
    time, ignore.case = TRUE, perl = TRUE
  )
  offset_formats <- c(
    "%Y-%m-%dT%H:%M:%OS%z", "%Y-%m-%d %H:%M:%OS%z",
    "%Y-%m-%dT%H:%M%z", "%Y-%m-%d %H:%M%z"
  )
  offset_idx <- which(is.na(parsed) & has_timezone & !time_empty)
  if (length(offset_idx) > 0L) {
    offset_parsed <- .parse_time_formats(time_offset[offset_idx],
                                         offset_formats, timezone,
                                         strict = TRUE)
    ok <- !is.na(offset_parsed)
    parsed[offset_idx[ok]] <- offset_parsed[ok]
  }

  # 3. Naive strings, interpreted in `timezone`.
  # %OS, not %S: %S drops a fraction of a second, which ties or reorders events
  # within the same second. Compact "%Y%m%d%H%M%S" has no separator for one.
  formats <- c(
    "%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M",
    "%Y/%m/%d %H:%M:%OS", "%Y/%m/%d %H:%M",
    "%Y.%m.%d %H:%M:%OS", "%Y.%m.%d %H:%M",
    "%Y-%m-%dT%H:%M:%OS", "%Y-%m-%dT%H:%M",
    "%Y%m%d%H%M%S", "%Y%m%d%H%M",
    "%d-%m-%Y %H:%M:%OS", "%d-%m-%Y %H:%M",
    "%d/%m/%Y %H:%M:%OS", "%d/%m/%Y %H:%M",
    "%d.%m.%Y %H:%M:%OS", "%d.%m.%Y %H:%M",
    "%d-%m-%YT%H:%M:%OS", "%d-%m-%YT%H:%M",
    "%m-%d-%Y %H:%M:%OS", "%m-%d-%Y %H:%M",
    "%m/%d/%Y %H:%M:%OS", "%m/%d/%Y %H:%M",
    "%m.%d.%Y %H:%M:%OS", "%m.%d.%Y %H:%M",
    "%m-%d-%YT%H:%M:%OS", "%m-%d-%YT%H:%M",
    "%d %b %Y %H:%M:%OS", "%d %b %Y %H:%M",
    "%d %B %Y %H:%M:%OS", "%d %B %Y %H:%M",
    "%b %d %Y %H:%M:%OS", "%b %d %Y %H:%M",
    "%B %d %Y %H:%M:%OS", "%B %d %Y %H:%M",
    "%Y-%m-%d", "%Y/%m/%d", "%Y.%m.%d",
    "%d-%m-%Y", "%d/%m/%Y", "%d.%m.%Y",
    "%m-%d-%Y", "%m/%d/%Y", "%m.%d.%Y",
    "%d %b %Y", "%d %B %Y", "%b %d %Y", "%B %d %Y"
  )
  naive_idx <- which(is.na(parsed) & !has_timezone & !time_empty)
  if (length(naive_idx) > 0L) {
    naive_formats <- .numeric_date_formats(time[naive_idx], formats)
    naive_parsed <- .parse_time_formats(time[naive_idx], naive_formats,
                                        timezone, strict = TRUE)
    ok <- !is.na(naive_parsed)
    parsed[naive_idx[ok]] <- naive_parsed[ok]
  }

  # 4. Anything left: numeric strings as Unix time.
  unresolved <- which(is.na(parsed) & !time_empty)
  if (length(unresolved) > 0L) {
    numeric_time <- suppressWarnings(as.numeric(time[unresolved]))
    ok <- !is.na(numeric_time)
    parsed[unresolved[ok]] <- as_unix(numeric_time[ok])
  }

  invalid <- which(is.na(parsed) & !time_empty)
  if (length(invalid) > 0L) {
    stop(errorCondition(
      paste0("Could not parse time values. Sample: ",
             paste(utils::head(time_original[invalid], 3), collapse = ", "),
             ". Supported: ISO-8601 (with Z/offset), YYYY-MM-DD HH:MM[:SS], ",
             "DD-MM-YYYY, MM-DD-YYYY (four-digit year), month names, ",
             "compact YYYYMMDDHHMMSS, ",
             "Unix timestamps. Use 'custom_format' otherwise."),
      class = "nestimate_bad_input", call = NULL))
  }
  parsed
}


#' Aggregate extra columns per session
#'
#' Numeric columns: mean. Character/factor columns: mode (most frequent value);
#' if tied, uses the first occurring value and emits a message.
#' @noRd
.aggregate_metadata <- function(df, session_col, extra_cols) {
  sessions <- unique(df[[session_col]])
  tie_counts <- integer(0)

  agg_list <- lapply(extra_cols, function(col) {
    vals <- df[[col]]
    if (is.numeric(vals)) {
      agg <- tapply(vals, df[[session_col]], mean, na.rm = TRUE)[sessions]
      agg[is.nan(agg)] <- NA_real_
      agg
    } else {
      n_ties <- 0L
      result <- tapply(vals, df[[session_col]], function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0L) return(NA_character_)
        tab <- table(v)
        max_count <- max(tab)
        modes <- names(tab)[tab == max_count]
        if (length(modes) > 1L) {
          n_ties <<- n_ties + 1L
          v[v %in% modes][1L]
        } else {
          modes
        }
      })[sessions]
      if (n_ties > 0L) tie_counts[[col]] <<- n_ties
      result
    }
  })
  names(agg_list) <- extra_cols

  if (length(tie_counts) > 0L) {
    parts <- vapply(names(tie_counts), function(col) {
      sprintf("'%s' (%d sessions)", col, tie_counts[[col]])
    }, character(1))
    message("Metadata aggregated per session: ties resolved by first ",
            "occurrence in ", paste(parts, collapse = ", "))
  }

  result <- data.frame(
    .session_id = sessions,
    stringsAsFactors = FALSE
  )
  for (col in extra_cols) {
    result[[col]] <- unname(agg_list[[col]])
  }
  result
}

#' Display label for several grouping columns that never merges two groups
#'
#' Pasting with a separator that occurs in a value can give two different
#' combinations the same label (`"1-a" + "b"` and `"1" + "a-b"` both give
#' `"1-a-b"`). The first separator whose labels map one-to-one onto the
#' combinations is used (the first, `"-"`, whenever nothing collides); when
#' every candidate collides, colliding labels get the group's number
#' appended.
#' @noRd
.distinct_label <- function(cols, seps = c("-", "|", " | ", " / ", " :: "),
                            context = "actor/session") {
  cols <- lapply(as.list(cols), as.character)
  id <- .observed_group_id(cols, context = context)
  cols <- unname(cols)
  one_to_one <- function(label) {
    all(tapply(id, label, function(i) length(unique(i)) == 1L))
  }
  labels <- lapply(seps, function(s) do.call(paste, c(cols, list(sep = s))))
  ok <- vapply(labels, one_to_one, logical(1))
  if (any(ok)) return(labels[[which(ok)[1L]]])
  label <- labels[[1L]]
  clash <- tapply(id, label, function(i) length(unique(i)) > 1L)
  label_clash <- label %in% names(clash)[clash]
  label[label_clash] <- paste0(label[label_clash], " #", id[label_clash])
  label
}
