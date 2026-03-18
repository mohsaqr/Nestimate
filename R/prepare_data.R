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
#' @param time_threshold Numeric. Maximum gap in seconds between consecutive
#'   events before a new session starts. Only used when \code{time} is
#'   provided. Default: 900 (15 minutes).
#' @param custom_format Character or NULL. Custom \code{strptime} format
#'   string for parsing timestamps. Default: NULL (auto-detect).
#' @param is_unix_time Logical. If TRUE, treat numeric time values as Unix
#'   timestamps. Default: FALSE (auto-detected for numeric columns).
#' @param unix_time_unit Character. Unit for Unix timestamps:
#'   \code{"seconds"}, \code{"milliseconds"}, or \code{"microseconds"}.
#'   Default: \code{"seconds"}.
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
#' \dontrun{
#' # Basic: actor + action + time
#' prepared <- prepare_data(df, actor = "student", action = "code",
#'                          time = "timestamp")
#' net <- build_network(prepared$sequence_data, method = "tna")
#'
#' # Multiple actors + multiple sessions
#' prepared <- prepare_data(df, actor = c("student", "group"),
#'                          action = "code", time = "timestamp",
#'                          session = c("course", "semester"))
#'
#' # Custom session threshold (5 minutes)
#' prepared <- prepare_data(df, actor = "student", action = "code",
#'                          time = "timestamp", time_threshold = 300)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{prepare_onehot}}
#'
#' @export
prepare_data <- function(data,
                         actor,
                         action,
                         time = NULL,
                         order = NULL,
                         session = NULL,
                         time_threshold = 900,
                         custom_format = NULL,
                         is_unix_time = FALSE,
                         unix_time_unit = c("seconds", "milliseconds",
                                            "microseconds")) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(action), length(action) == 1, action %in% names(data))
  stopifnot(is.numeric(time_threshold), length(time_threshold) == 1,
            time_threshold > 0)
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
      df$.actor <- interaction(df[, actor, drop = FALSE], sep = "-",
                               drop = TRUE)
      actor_col <- ".actor"
    } else {
      actor_col <- actor
    }
  }

  # ---- Session (explicit grouping) ----
  if (!is.null(session)) {
    stopifnot(is.character(session), all(session %in% names(df)))
    if (length(session) > 1L) {
      df$.session_explicit <- interaction(
        df[, session, drop = FALSE], sep = "-", drop = TRUE
      )
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
  if (!is.null(session)) {
    df$.base_group <- interaction(
      df[[actor_col]], df$.session_explicit, sep = " | ", drop = TRUE
    )
  } else {
    df$.base_group <- df[[actor_col]]
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
                          unix_time_unit = unix_time_unit)
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

    df$.session_id <- paste0(df$.base_group, " s", df$.inferred_nr)

  } else {
    # No time: sort by base_group + order
    df <- df[base::order(df$.base_group, df[[order_col]]), ]

    # Each base group = one session
    df$.session_id <- as.character(df$.base_group)
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
    time_data <- as.data.frame(time_mat)
    names(time_data) <- paste0("time_T", seq_len(max_len))
    for (j in seq_len(ncol(time_data))) {
      time_data[[j]] <- as.POSIXct(time_data[[j]], origin = "1970-01-01")
    }
  } else {
    time_data <- NULL
  }

  # Meta data
  meta_data <- data.frame(.session_id = sessions, stringsAsFactors = FALSE)
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

#' Parse time values from various formats
#'
#' Tries 40+ date/time formats, Unix timestamps, and custom formats.
#' Matches tna::parse_time logic.
#'
#' @param time Character or numeric vector of time values.
#' @param custom_format Character or NULL. Custom strptime format.
#' @param is_unix_time Logical. Treat numeric as Unix timestamp.
#' @param unix_time_unit Character. "seconds", "milliseconds", "microseconds".
#' @return POSIXct vector.
#' @noRd
.parse_time <- function(time, custom_format = NULL, is_unix_time = FALSE,
                        unix_time_unit = "seconds") {
  # Already POSIXct
  if (inherits(time, c("POSIXct", "POSIXlt"))) return(time)

  # Numeric Unix timestamps
  if (is.numeric(time) && is_unix_time) {
    divisor <- switch(unix_time_unit,
      seconds = 1,
      milliseconds = 1000,
      microseconds = 1e6
    )
    return(as.POSIXct(time / divisor, origin = "1970-01-01"))
  }

  # Character parsing
  time <- trimws(as.character(time))
  # Strip fractional seconds and timezone suffixes for format matching
  time_clean <- gsub("(\\.\\d{1,3})?[A-Za-z ]*$", "", time)

  # Custom format first
  if (!is.null(custom_format)) {
    parsed <- as.POSIXct(strptime(time_clean, format = custom_format))
    if (!all(is.na(parsed))) return(parsed)
  }

  # Try standard formats
  formats <- c(
    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
    "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M",
    "%Y.%m.%d %H:%M:%S", "%Y.%m.%d %H:%M",
    "%Y-%m-%dT%H:%M:%S", "%Y-%m-%dT%H:%M", "%Y-%m-%dT%H:%M:%OS",
    "%Y-%m-%d %H:%M:%S%z", "%Y-%m-%d %H:%M%z",
    "%Y-%m-%d %H:%M:%S %z", "%Y-%m-%d %H:%M %z",
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

  for (fmt in formats) {
    parsed <- as.POSIXct(strptime(time_clean, format = fmt))
    if (!all(is.na(parsed))) return(parsed)
  }

  # Last resort: try as numeric Unix timestamp
  numeric_time <- suppressWarnings(as.numeric(time))
  if (!any(is.na(numeric_time))) {
    divisor <- switch(unix_time_unit,
      seconds = 1, milliseconds = 1000, microseconds = 1e6
    )
    return(as.POSIXct(numeric_time / divisor, origin = "1970-01-01"))
  }

  stop("Could not parse time values. Sample: ",
       paste(utils::head(time, 3), collapse = ", "),
       ". Use 'custom_format' argument.", call. = FALSE)
}


#' Aggregate extra columns per session
#'
#' Numeric columns: mean. Character/factor columns: mode (most frequent value);
#' if tied, uses the first occurring value and emits a message.
#' @noRd
.aggregate_metadata <- function(df, session_col, extra_cols) {
  sessions <- unique(df[[session_col]])

  agg_list <- lapply(extra_cols, function(col) {
    vals <- df[[col]]
    if (is.numeric(vals)) {
      tapply(vals, df[[session_col]], mean, na.rm = TRUE)[sessions]
    } else {
      tapply(vals, df[[session_col]], function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0L) return(NA_character_)
        tab <- table(v)
        max_count <- max(tab)
        modes <- names(tab)[tab == max_count]
        if (length(modes) > 1L) {
          # Tie: use first occurring value
          chosen <- v[v %in% modes][1L]
          message(sprintf(
            "Column '%s': tied mode (%s), using first value '%s'",
            col, paste(modes, collapse = ", "), chosen
          ))
          chosen
        } else {
          modes
        }
      })[sessions]
    }
  })
  names(agg_list) <- extra_cols

  result <- data.frame(
    .session_id = sessions,
    stringsAsFactors = FALSE
  )
  for (col in extra_cols) {
    result[[col]] <- unname(agg_list[[col]])
  }
  result
}
