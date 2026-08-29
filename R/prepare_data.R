#' Prepare Event Log Data for Network Estimation
#'
#' @description
#' Converts event log data (actor, action, time) into wide sequence format
#' suitable for \code{\link{build_network}}. Automatically parses timestamps,
#' detects sessions from time gaps, and handles tie-breaking.
#'
#' @param data Data frame with event log columns.
#' @param actor Character or character vector. Column name(s) identifying who
#'   performed the action (e.g. \code{"student"} or
#'   \code{c("student", "group")}). If missing, all data is treated as one
#'   actor.
#' @param action Character. Column name containing the action/state/code.
#' @details
#' Sessions are identified by the observed combinations of the \code{actor} and
#' \code{session} columns, so identifiers containing separator characters stay
#' distinct and high-cardinality identifiers cannot overflow. Missing values in
#' any grouping column raise an error: drop or relabel those events first.
#' @param time Character or NULL. Column name containing timestamps.
#'   Supports ISO8601, Unix timestamps (numeric), and 40+ date/time formats.
#'   If NULL, row order defines the sequence. Default: NULL.
#' @param order Character or NULL. Column name for tie-breaking when
#'   timestamps are identical. If NULL, original row order is used.
#'   Default: NULL.
#' @param session Character, character vector, or NULL. Column name(s) for
#'   explicit session grouping (e.g. \code{"course"} or
#'   \code{c("course", "semester")}). When combined with \code{time},
#'   sessions are further split by time gaps. Default: NULL.
#' @param time_threshold Numeric or FALSE. Maximum gap in seconds between
#'   consecutive events before a new session starts. Only used when
#'   \code{time} is provided. Set to \code{FALSE} to switch session-interval
#'   splitting off, so each actor (or actor-session) forms a single sequence
#'   however long the gaps are. Default: 900 (15 minutes).
#' @param custom_format Character or NULL. Custom \code{strptime} format
#'   string for parsing timestamps. Default: NULL (auto-detect).
#' @param is_unix_time Logical. If TRUE, treat numeric time values as Unix
#'   timestamps. Default: FALSE (auto-detected for numeric columns).
#' @param unix_time_unit Character. Unit for Unix timestamps:
#'   \code{"seconds"}, \code{"milliseconds"}, or \code{"microseconds"}.
#'   Default: \code{"seconds"}.
#' @param timezone Character. An Olson time zone (see
#'   \code{\link[base]{OlsonNames}}) used to interpret timestamps that carry
#'   no zone information, and in which Unix timestamps are expressed.
#'   Timestamps that end in \code{Z}, \code{UTC}, \code{GMT} or a numeric
#'   offset such as \code{+02:00} are converted from that offset. Parsing is
#'   therefore independent of the machine's local time zone.
#'   Default: \code{"UTC"}.
#'
#' @return A list with class \code{"nestimate_data"} containing:
#' \describe{
#'   \item{sequence_data}{Data frame in wide format (one row per session,
#'     columns T1, T2, ...).}
#'   \item{long_data}{The processed long-format data with session IDs.}
#'   \item{meta_data}{Session-level metadata (session ID, actor).}
#'   \item{time_data}{Parsed time values in wide format (if time provided).}
#'   \item{statistics}{List with total_sessions, total_actions,
#'     max_sequence_length, unique_actors, etc.}
#' }
#'
#' @examples
#' df <- data.frame(
#'   student = rep(1:3, each = 5),
#'   code = sample(c("read", "write", "test"), 15, replace = TRUE),
#'   timestamp = seq.POSIXt(as.POSIXct("2024-01-01"), by = "min", length.out = 15)
#' )
#' prepared <- prepare(df, actor = "student", action = "code",
#'                     time = "timestamp")
#' net <- build_network(prepared$sequence_data, method = "relative")
#'
#' @seealso \code{\link{build_network}}, \code{\link{prepare_onehot}}
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
  stopifnot(is.data.frame(data))
  .check_timezone(timezone)
  stopifnot(is.character(action), length(action) == 1, action %in% names(data))
  # FALSE disables gap-based splitting. Inf makes the `gaps > time_threshold`
  # test below never fire, so no separate code path is needed.
  if (isFALSE(time_threshold)) {
    time_threshold <- Inf
  }
  stopifnot(is.numeric(time_threshold), length(time_threshold) == 1,
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
    stopifnot(is.character(actor), all(actor %in% names(df)))
    if (length(actor) > 1L) {
      # Display label only. Identity comes from the underlying columns below.
      df$.actor <- .group_label(df[, actor, drop = FALSE], sep = "-")
      actor_col <- ".actor"
    } else {
      actor_col <- actor
    }
  }
  actor_key_cols <- if (default_actor) ".actor" else actor

  # ---- Session (explicit grouping) ----
  if (!is.null(session)) {
    stopifnot(is.character(session), all(session %in% names(df)))
    if (length(session) > 1L) {
      df$.session_explicit <- .group_label(df[, session, drop = FALSE],
                                           sep = "-")
    } else {
      df$.session_explicit <- df[[session]]
    }
  }

  # ---- Order (tiebreaker) ----
  if (is.null(order)) {
    df$.order <- seq_len(n)
    order_col <- ".order"
  } else {
    stopifnot(is.character(order), length(order) == 1, order %in% names(df))
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
    .group_label(list(df[[actor_col]], df$.session_explicit))
  }

  # ---- Time parsing + inferred session detection ----
  if (!is.null(time)) {
    stopifnot(is.character(time), length(time) == 1, time %in% names(df))

    # Auto-detect numeric as unix
    if (is.numeric(df[[time]])) {
      is_unix_time <- TRUE
    }

    parsed <- .parse_time(df[[time]], custom_format = custom_format,
                          is_unix_time = is_unix_time,
                          unix_time_unit = unix_time_unit,
                          timezone = timezone)
    df$.parsed_time <- parsed

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
    time_cols <- lapply(
      seq_len(max_len),
      function(j) as.POSIXct(time_mat[, j], origin = "1970-01-01")
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

  # Aggregate extra columns per session
  special_cols <- c(action, actor_col, time, order_col, session,
                    grep("^\\.", names(df), value = TRUE))
  extra_cols <- setdiff(names(df), special_cols)
  if (length(extra_cols) > 0L) {
    agg <- .aggregate_metadata(df, session_col = ".session_id",
                               extra_cols = extra_cols)
    meta_data <- merge(meta_data, agg, by = ".session_id", sort = FALSE)
  }

  # Statistics
  stats <- list(
    total_sessions = length(sessions),
    total_actions = nrow(df),
    max_sequence_length = max_len
  )
  if (!default_actor) {
    stats$unique_actors <- length(unique(df[[actor_col]]))
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


#' Print Method for nestimate_data
#' @param x A \code{nestimate_data} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @examples
#' events <- data.frame(
#'   actor  = c("u1","u1","u1","u2","u2","u2"),
#'   action = c("A","B","C","B","A","C"),
#'   time   = c(1,2,3,1,2,3)
#' )
#' nd <- prepare(events, action = "action",
#'               actor = "actor", time = "time")
#' print(nd)
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
  stopifnot(
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
#' @return POSIXct vector (NA where no format matched).
#' @noRd
.parse_time_formats <- function(time, formats, timezone) {
  out <- rep(NA_real_, length(time))
  # Sequential fill: each format only sees what earlier formats left NA, so
  # the loop is over formats (short), not elements.
  for (fmt in formats) {
    idx <- which(is.na(out) & !is.na(time))
    if (length(idx) == 0L) break
    parsed <- as.POSIXct(strptime(time[idx], format = fmt, tz = timezone))
    ok <- !is.na(parsed)
    out[idx[ok]] <- as.numeric(parsed[ok])
  }
  as.POSIXct(out, origin = "1970-01-01", tz = timezone)
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
    "(?:Z|UTC|GMT|[+-][0-9]{2}:?[0-9]{2}|[A-Za-z]{2,})$",
    time, ignore.case = TRUE, perl = TRUE
  )
  offset_formats <- c(
    "%Y-%m-%dT%H:%M:%OS%z", "%Y-%m-%d %H:%M:%OS%z",
    "%Y-%m-%dT%H:%M%z", "%Y-%m-%d %H:%M%z"
  )
  offset_idx <- which(is.na(parsed) & has_timezone & !time_empty)
  if (length(offset_idx) > 0L) {
    offset_parsed <- .parse_time_formats(time_offset[offset_idx],
                                         offset_formats, timezone)
    ok <- !is.na(offset_parsed)
    parsed[offset_idx[ok]] <- offset_parsed[ok]
  }

  # 3. Naive strings, interpreted in `timezone`.
  formats <- c(
    "%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
    "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M",
    "%Y.%m.%d %H:%M:%S", "%Y.%m.%d %H:%M",
    "%Y-%m-%dT%H:%M:%OS", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%dT%H:%M",
    "%Y%m%d%H%M%S", "%Y%m%d%H%M",
    "%d-%m-%Y %H:%M:%S", "%d-%m-%Y %H:%M",
    "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M",
    "%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M",
    "%d-%m-%YT%H:%M:%S", "%d-%m-%YT%H:%M",
    "%m-%d-%Y %H:%M:%S", "%m-%d-%Y %H:%M",
    "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
    "%m.%d.%Y %H:%M:%S", "%m.%d.%Y %H:%M",
    "%m-%d-%YT%H:%M:%S", "%m-%d-%YT%H:%M",
    "%d %b %Y %H:%M:%S", "%d %b %Y %H:%M",
    "%d %B %Y %H:%M:%S", "%d %B %Y %H:%M",
    "%b %d %Y %H:%M:%S", "%b %d %Y %H:%M",
    "%B %d %Y %H:%M:%S", "%B %d %Y %H:%M",
    "%Y-%m-%d", "%Y/%m/%d", "%Y.%m.%d",
    "%d-%m-%Y", "%d/%m/%Y", "%d.%m.%Y",
    "%m-%d-%Y", "%m/%d/%Y", "%m.%d.%Y",
    "%d %b %Y", "%d %B %Y", "%b %d %Y", "%B %d %Y"
  )
  naive_idx <- which(is.na(parsed) & !has_timezone & !time_empty)
  if (length(naive_idx) > 0L) {
    naive_parsed <- .parse_time_formats(time[naive_idx], formats, timezone)
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
    stop("Could not parse time values. Sample: ",
         paste(utils::head(time_original[invalid], 3), collapse = ", "),
         ". Supported: ISO-8601 (with Z/offset), YYYY-MM-DD HH:MM[:SS], ",
         "DD-MM-YYYY, MM-DD-YYYY, month names, compact YYYYMMDDHHMMSS, ",
         "Unix timestamps. Use 'custom_format' otherwise.", call. = FALSE)
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
