#' Convert Wide Sequences to Long Format
#'
#' @description
#' `wide_to_long()` turns sequences in wide format (one row per sequence,
#' one column per time point) into long format, with one row per sequence
#' and time point. The time point columns are those named `time_prefix`
#' followed by a number (`V1`, `V2`, ...); every other column is kept and
#' repeated on each row of its sequence.
#'
#' @param data A data frame in wide format, one row per sequence.
#' @param id_col The column identifying each sequence. Default `NULL`, which
#'   adds a column `id` holding the row number of the wide input (replacing
#'   any existing column of that name). A name that is not a column raises
#'   an error.
#' @param time_prefix The prefix of the time point columns: with `"V"`
#'   (default), the columns `V1`, `V2`, ... are read. No column matching it
#'   raises an error.
#' @param action_col The name of the state column in the output. Default
#'   `"Action"`.
#' @param time_col The name of the time column in the output. Default
#'   `"Time"`.
#' @param drop_na Whether to drop the rows whose state is missing, such as
#'   the padding after a shorter sequence ends. Default `TRUE`.
#'
#' @return A data frame in long format, one row per sequence and time point,
#'   sorted by identifier and then time. Its columns are the identifier, the
#'   other kept columns, the time and the state:
#' \describe{
#'   \item{id}{The sequence identifier, named by `id_col`; with
#'     `id_col = NULL` it is called `id` and holds the row number of the wide
#'     input (integer).}
#'   \item{Time}{The time point within the sequence (integer), taken from
#'     the number in the wide column name; named by `time_col`.}
#'   \item{Action}{The state at that time point; named by `action_col`.}
#' }
#'
#' @details
#' The long format is the event-log layout read by [build_network()] with
#' `action` (and [prepare()]): give `action = "Action"`, `actor = "id"` and
#' `time = "Time"` (or the names chosen here).
#'
#' @examples
#' wide_data <- data.frame(
#'   V1 = c("A", "B", "C"), V2 = c("B", "C", "A"), V3 = c("C", "A", "B")
#' )
#' long_data <- wide_to_long(wide_data)
#' head(long_data)
#'
#' @seealso \code{\link{long_to_wide}} for the reverse conversion,
#'   \code{\link{prepare_for_tna}} for preparing data for TNA analysis.
#'
#' @export
wide_to_long <- function(data,
                         id_col = NULL,
                         time_prefix = "V",
                         action_col = "Action",
                         time_col = "Time",
                         drop_na = TRUE) {
  # Input validation
  .assert("`data` must be a data frame" = is.data.frame(data))
  .assert("`time_prefix` must be a single character string" =
            is.character(time_prefix) && length(time_prefix) == 1)
  .assert("`action_col` must be a single character string" =
            is.character(action_col) && length(action_col) == 1)
  .assert("`time_col` must be a single character string" =
            is.character(time_col) && length(time_col) == 1)
  .assert("`drop_na` must be TRUE or FALSE" =
            is.logical(drop_na) && length(drop_na) == 1)

  # Identify time columns
  time_cols <- grep(paste0("^", time_prefix, "[0-9]+$"), names(data), value = TRUE)
  if (length(time_cols) == 0) {
    stop(errorCondition(
      paste0("No columns matching time prefix '", time_prefix, "' found."),
      class = "nestimate_bad_input", call = NULL))
  }

  # Add ID column if needed
  if (is.null(id_col)) {
    data$id <- seq_len(nrow(data))
    id_col <- "id"
  } else if (!id_col %in% names(data)) {
    stop(errorCondition(paste0("ID column '", id_col, "' not found in data."),
                        class = "nestimate_bad_input", call = NULL))
  }

  # Identify other columns to preserve
  other_cols <- setdiff(names(data), c(time_cols, id_col))

  # Extract time indices from column names
  time_nums <- as.integer(gsub(paste0("^", time_prefix), "", time_cols))

  # Manual stack: one data frame per time column, then rbind
  parts <- lapply(seq_along(time_cols), function(i) {
    row <- data[, c(id_col, other_cols), drop = FALSE]
    row[[time_col]] <- time_nums[i]
    row[[action_col]] <- data[[time_cols[i]]]
    row
  })
  result <- do.call(rbind, parts)

  # Drop NA values if requested
  if (drop_na) {
    result <- result[!is.na(result[[action_col]]), ]
  }

  # Sort by ID and time
  result <- result[order(result[[id_col]], result[[time_col]]), ]

  # Reset row names
  rownames(result) <- NULL

  return(result)
}


#' Convert Long Format to Wide Sequences
#'
#' @description
#' `long_to_wide()` turns a long event table (one row per action) into wide
#' sequences, with one row per sequence and one column per position
#' (`V1`, `V2`, ...). Each sequence's actions are ordered by `time_col`, and
#' the k-th action fills column k.
#'
#' @param data A data frame in long format, one row per action.
#' @param id_col The column identifying each sequence. Default `"Actor"`;
#'   a missing column raises an error.
#' @param time_col The column that orders the actions within a sequence.
#'   Default `"Time"`. Its values only order the actions; they do not set
#'   the column numbers, so gaps between time points are not kept. Actions
#'   with the same time keep their row order. When the column is absent, the
#'   row order is used.
#' @param action_col The column holding the actions (states). Default
#'   `"Action"`; a missing column raises an error.
#' @param time_prefix The prefix of the position columns in the output.
#'   Default `"V"`.
#' @param fill_na How to handle sequences of different lengths. With `TRUE`
#'   (default), shorter sequences are padded with `NA`, so every row has the
#'   columns of the longest sequence. With `FALSE`, only the positions present
#'   in every sequence are kept, so longer sequences are cut to the length of
#'   the shortest.
#'
#' @return A data frame in wide format, one row per sequence in the sorted
#'   order of `id_col`: the `id_col` column followed by the position columns
#'   (`V1`, `V2`, ..., named with `time_prefix`) holding the action at each
#'   position. Other columns of `data` are dropped. A missing action stays a
#'   missing cell.
#'
#' @details
#' The wide format is the one [build_network()], `tna::tna()` and related
#' functions read without column declarations, since the position columns
#' form a numbered run; [build_network()] sets the `id_col` column aside as
#' metadata, with a message.
#'
#' @examples
#' long_data <- data.frame(
#'   Actor = rep(1:3, each = 4),
#'   Time = rep(1:4, 3),
#'   Action = sample(c("A", "B", "C"), 12, replace = TRUE)
#' )
#' wide_data <- long_to_wide(long_data, id_col = "Actor")
#' head(wide_data)
#'
#' @seealso \code{\link{wide_to_long}} for the reverse conversion,
#'   \code{\link{prepare_for_tna}} for preparing data for TNA analysis.
#'
#' @export
long_to_wide <- function(data,
                         id_col = "Actor",
                         time_col = "Time",
                         action_col = "Action",
                         time_prefix = "V",
                         fill_na = TRUE) {
  # Input validation
  .assert("`data` must be a data frame" = is.data.frame(data))
  .assert("`id_col` must be a single character string" =
            is.character(id_col) && length(id_col) == 1)
  .assert("`time_col` must be a single character string" =
            is.character(time_col) && length(time_col) == 1)
  .assert("`action_col` must be a single character string" =
            is.character(action_col) && length(action_col) == 1)
  .assert("`time_prefix` must be a single character string" =
            is.character(time_prefix) && length(time_prefix) == 1)
  .assert("`fill_na` must be TRUE or FALSE" =
            is.logical(fill_na) && length(fill_na) == 1)

  # Check required columns exist
  required_cols <- c(id_col, action_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(errorCondition(
      paste0(
        "Missing required columns: ",
        paste(missing_cols, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  # Sort by id and time (stable sort preserves original order for ties)
  data$.orig_order <- seq_len(nrow(data))
  if (time_col %in% names(data)) {
    data <- data[order(data[[id_col]], data[[time_col]], data$.orig_order), ]
  } else {
    data <- data[order(data[[id_col]], data$.orig_order), ]
  }
  data$.orig_order <- NULL

  # Compute within-group sequential index
  data[[".time_idx"]] <- ave(
    seq_len(nrow(data)), data[[id_col]], FUN = seq_along
  )

  # Create column names
  data[[".time_name"]] <- paste0(time_prefix, data[[".time_idx"]])

  # Spread manually: split by id, build one row per id
  groups <- split(data, data[[id_col]])
  wide_list <- lapply(groups, function(g) {
    row <- data.frame(g[[id_col]][1], stringsAsFactors = FALSE)
    names(row) <- id_col
    # Loop kept: one column is added per event, in event order, and a
    # repeated time name overwrites, as `[[<-` does.
    for (i in seq_len(nrow(g))) {
      row[[g[[".time_name"]][i]]] <- g[[action_col]][i]
    }
    row
  })
  # Align columns (sequences may differ in length).
  # fill_na = TRUE  : pad shorter sequences' missing time points with NA so
  #                    every row shares the same V1..Vn columns.
  # fill_na = FALSE : do not fill missing time points with NA; keep only the
  #                    time-point columns present in every sequence (truncate
  #                    ragged tails) so no NA-padded cells are emitted.
  if (fill_na) {
    all_cols <- unique(unlist(lapply(wide_list, names)))
    wide_list <- lapply(wide_list, function(row) {
      missing <- setdiff(all_cols, names(row))
      row[missing] <- rep(list(NA_character_), length(missing))
      row[, all_cols, drop = FALSE]
    })
  } else {
    common_cols <- Reduce(intersect, lapply(wide_list, names))
    wide_list <- lapply(wide_list, function(row) {
      row[, common_cols, drop = FALSE]
    })
  }
  result <- do.call(rbind, wide_list)

  # Reorder columns to ensure V1, V2, V3... order
  time_cols <- grep(paste0("^", time_prefix, "[0-9]+$"), names(result), value = TRUE)
  time_nums <- as.integer(gsub(time_prefix, "", time_cols))
  time_cols <- time_cols[order(time_nums)]
  other_cols <- setdiff(names(result), time_cols)
  result <- result[, c(other_cols, time_cols), drop = FALSE]

  # Convert to data frame
  result <- as.data.frame(result)
  rownames(result) <- NULL

  return(result)
}


#' Prepare Data for TNA Analysis
#'
#' @description
#' `prepare_for_tna()` returns the sequence columns of a data set as a wide
#' data frame for `tna::tna()` and related functions: one row per sequence
#' and the position columns `V1`, `V2`, ... in numeric order, with factors
#' turned into character. Long input is converted with [long_to_wide()]
#' first, padding shorter sequences with `NA`. Every other column (an id,
#' metadata) is dropped.
#'
#' @param data A data frame of sequences: wide, with position columns named
#'   `V1`, `V2`, ..., or long, with one row per action.
#' @param type How to read `data`. `"sequences"` (default) reads wide
#'   position columns; `"long"` converts a long table with [long_to_wide()],
#'   which needs `id_col` and `action_col`; `"auto"` reads long data when
#'   `action_col` is a column and no `V<number>` column exists, wide data in
#'   the opposite case, and, when both exist, long data only if `id_col` is a
#'   column and there are more than twice as many rows as ids. Data with
#'   neither raises an error. Whatever the type, the result needs `V<number>`
#'   columns; others, such as `T1`, `T2`, are not renamed and raise an error.
#' @param state_names The expected states, used with `validate`. Default
#'   `NULL`, which accepts the states found in the data.
#' @param id_col The column identifying sequences in long data. Default
#'   `"Actor"`.
#' @param time_col The column ordering the actions in long data. Default
#'   `"Time"`.
#' @param action_col The column holding the actions in long data. Default
#'   `"Action"`.
#' @param validate Whether to check the states against `state_names`; states
#'   not listed are named in a warning. Ignored when `state_names` is `NULL`.
#'   Default `TRUE`.
#'
#' @return A data frame with one row per sequence and only the position
#'   columns `V1`, `V2`, ..., in numeric order, holding character states (or
#'   the original type for non-factor columns), with no row names.
#'
#' @examples
#' # From wide format sequences
#' sequences <- data.frame(
#'   V1 = c("A","B","C","A"), V2 = c("B","C","A","B"),
#'   V3 = c("C","A","B","C"), V4 = c("A","B","A","B")
#' )
#' tna_data <- prepare_for_tna(sequences, type = "sequences")
#'
#' @seealso \code{\link{wide_to_long}}, \code{\link{long_to_wide}} for
#'   format conversions.
#'
#' @export
prepare_for_tna <- function(data,
                            type = c("sequences", "long", "auto"),
                            state_names = NULL,
                            id_col = "Actor",
                            time_col = "Time",
                            action_col = "Action",
                            validate = TRUE) {
  # Input validation
  .assert("`data` must be a data frame" = is.data.frame(data))
  type <- match.arg(type)

  # Auto-detect format
  if (type == "auto") {
    has_time_cols <- any(grepl("^V[0-9]+$", names(data)))
    has_action_col <- action_col %in% names(data)

    if (has_time_cols && !has_action_col) {
      type <- "sequences"
    } else if (has_action_col && !has_time_cols) {
      type <- "long"
    } else if (has_time_cols && has_action_col) {
      # If both, use long if nrow >> number of unique IDs
      if (id_col %in% names(data)) {
        n_ids <- length(unique(data[[id_col]]))
        if (nrow(data) > n_ids * 2) {
          type <- "long"
        } else {
          type <- "sequences"
        }
      } else {
        type <- "sequences"
      }
    } else {
      stop(errorCondition(
        "Cannot auto-detect data format. Please specify 'type' explicitly.",
        class = "nestimate_bad_input", call = NULL))
    }
  }

  # Convert long to wide if needed
  if (type == "long") {
    if (!action_col %in% names(data)) {
      stop(errorCondition(
        paste0("Action column '", action_col, "' not found in data."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (!id_col %in% names(data)) {
      stop(errorCondition(paste0("ID column '", id_col, "' not found in data."),
                          class = "nestimate_bad_input", call = NULL))
    }

    data <- long_to_wide(
      data,
      id_col = id_col,
      time_col = time_col,
      action_col = action_col
    )
  }

  # Identify sequence columns (V1, V2, ...)
  seq_cols <- grep("^V[0-9]+$", names(data), value = TRUE)
  if (length(seq_cols) == 0) {
    stop(errorCondition("No sequence columns (V1, V2, ...) found in data.",
                        class = "nestimate_bad_input", call = NULL))
  }

  # Extract just the sequence columns
  result <- data[, seq_cols, drop = FALSE]

  # Convert factors to characters
  is_fct <- vapply(result, is.factor, logical(1L))
  result[is_fct] <- lapply(result[is_fct], as.character)

  # Validate state names if requested
  if (validate && !is.null(state_names)) {
    all_states <- unique(unlist(result, use.names = FALSE))
    all_states <- all_states[!is.na(all_states)]
    unknown_states <- setdiff(all_states, state_names)
    if (length(unknown_states) > 0) {
      warning(warningCondition(
        paste0("Unknown states found: ",
               paste(unknown_states, collapse = ", ")),
        class = "nestimate_unknown_states", call = NULL))
    }
  }

  # Ensure proper column ordering
  col_nums <- as.integer(gsub("V", "", seq_cols))
  result <- result[, order(col_nums), drop = FALSE]

  # Convert to data frame
  result <- as.data.frame(result)
  rownames(result) <- NULL

  return(result)
}


#' Convert Action Column to One-Hot Encoding
#'
#' @description
#' `action_to_onehot()` adds one 0/1 indicator column per state to a data
#' frame with a categorical action column: the column of a row's state is 1
#' and the others are 0. The result is the one-hot layout that
#' [build_network()] reads with `codes` or `format = "onehot"`, and [wtna()]
#' reads directly.
#'
#' @param data A data frame with an action column.
#' @param action_col The action column. Default `"Action"`; a name that is
#'   not a column raises an error.
#' @param states The states to make columns for, in that order. Default
#'   `NULL`, every non-missing value of the action column in order of first
#'   appearance. A listed state that never occurs gets a column of zeros, and
#'   a row whose state is not listed is 0 in every column.
#' @param drop_action Whether to remove the action column. Default `TRUE`.
#' @param sort_states Whether to sort the state columns alphabetically (in
#'   level order for a factor action column). Default `FALSE`.
#' @param prefix A prefix for the state column names. Without one, a state
#'   named like an existing column replaces that column. Default `""`.
#'
#' @return The input data frame with one integer column per state appended,
#'   named `paste0(prefix, state)`. All other columns are kept; the action
#'   column is removed unless `drop_action = FALSE`. A row with a missing
#'   action is `NA` in every state column.
#'
#' @examples
#' long_data <- data.frame(
#'   Actor = rep(1:3, each = 4),
#'   Time = rep(1:4, 3),
#'   Action = sample(c("A", "B", "C"), 12, replace = TRUE)
#' )
#' onehot_data <- action_to_onehot(long_data)
#' head(onehot_data)
#'
#' @export
action_to_onehot <- function(data, action_col = "Action", states = NULL,
                             drop_action = TRUE, sort_states = FALSE,
                             prefix = "") {
  .assert("`data` must be a data frame" = is.data.frame(data))
  .assert("`action_col` must name one column of `data`" =
            is.character(action_col) && length(action_col) == 1 &&
            action_col %in% names(data))

  if (is.null(states)) {
    states <- unique(data[[action_col]])
    states <- states[!is.na(states)]
  }

  if (sort_states) {
    states <- sort(states)
  }

  # Loop kept: columns are added one state at a time, and a state named like
  # an existing column replaces it before the next state is read.
  for (state in states) {
    col_name <- paste0(prefix, state)
    data[[col_name]] <- as.integer(data[[action_col]] == state)
  }

  if (drop_action) {
    data[[action_col]] <- NULL
  }

  return(data)
}


#' Import One-Hot Encoded Data into Sequence Format
#'
#' @description
#' Converts binary indicator (one-hot) data into the wide sequence format
#' expected by \code{\link{build_network}} and \code{tna::tna()}. Each binary
#' column represents a state; rows where the value is 1 are marked with the
#' column name. Supports optional windowed aggregation.
#'
#' Simultaneous active states are preserved using the same window-span
#' representation as \code{tna::import_onehot()}: each input row/window is
#' expanded to one sequence slot per code and transition counting occurs between
#' windows, not between simultaneous states inside the same row.
#'
#' @param data Data frame with binary (0/1) indicator columns.
#' @param cols Character vector. Names of the one-hot columns to use.
#' @param actor Character or NULL. Name of the actor/ID column. If NULL,
#'   all rows are treated as a single sequence. Default: NULL.
#' @param session Character or NULL. Name of the session column for
#'   sub-grouping within actors. Default: NULL.
#' @param interval Integer or NULL. Number of rows per time point in the
#'   output. If NULL, all rows become a single time point group.
#'   Default: NULL.
#' @param window_size Integer (>= 1). Number of consecutive rows to aggregate
#'   into each window. Default: 3. Set \code{window_size = 1} for no
#'   windowing (each row is its own time point).
#' @param window_type Character. \code{"non-overlapping"} (fixed, separate
#'   windows) or \code{"overlapping"} (rolling, step = 1: each window covers
#'   exactly \code{window_size} rows, so a group of n rows gives
#'   n - window_size + 1 windows, and a shorter group gives none).
#'   Default: \code{"non-overlapping"}.
#' @param aggregate Logical. If TRUE, aggregate within each window by
#'   taking the first non-NA indicator per column. Default: FALSE.
#'
#' @return A data frame in wide format, one row per actor/session sequence,
#'   with columns named \code{W<window>_T<slot>} where each cell contains a
#'   state name or \code{NA}. Attributes \code{windowed} (always
#'   \code{TRUE}), \code{window_size}, \code{window_span} (the number of
#'   \code{cols}) and \code{codes} (the \code{cols} themselves) are set on
#'   the result.
#'
#' @noRd
.prepare_onehot <- function(data,
                            cols,
                            actor = NULL,
                            session = NULL,
                            interval = NULL,
                            window_size = 3L,
                            window_type = c("non-overlapping", "overlapping"),
                            aggregate = FALSE) {
  .assert("`data` must be a data frame" = is.data.frame(data))
  .assert("`cols` must be a character vector of column names" =
            is.character(cols) && length(cols) >= 1)
  .assert("`cols` must name columns of `data`" = all(cols %in% names(data)))
  window_size <- as.integer(window_size)
  .assert("`window_size` must be a whole number of at least 1" =
            window_size >= 1L)
  window_type <- match.arg(window_type)

  interval <- if (is.null(interval)) ncol(data) * nrow(data) else interval
  tna_window_type <- switch(window_type,
                            "non-overlapping" = "tumbling",
                            "overlapping" = "sliding")

  work <- data
  code_values <- lapply(cols, function(col) {
    ifelse(work[[col]] == 1, col, NA_character_)
  })
  names(code_values) <- cols
  work[cols] <- code_values

  actor_col <- if (is.null(actor)) ".actor" else actor
  session_col <- if (is.null(session)) ".session" else session
  if (!is.null(actor)) {
    .assert("`actor` must name a column of `data`" = actor %in% names(data))
  }
  if (!is.null(session)) {
    .assert("`session` must name a column of `data`" =
              session %in% names(data))
  }
  if (is.null(actor)) work[[actor_col]] <- 1L
  if (is.null(session)) work[[session_col]] <- 1L

  group_ids <- .observed_group_id(
    work[, c(actor_col, session_col), drop = FALSE],
    context = "actor/session"
  )
  groups <- lapply(.split_first_seen(group_ids),
                   function(i) work[i, , drop = FALSE])

  lag_chr <- function(x, n) {
    if (n == 0L) return(x)
    if (n >= length(x)) return(rep(NA_character_, length(x)))
    c(rep(NA_character_, n), head(x, -n))
  }

  first_non_missing <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) NA_character_ else x[1L]
  }

  work <- do.call(rbind, lapply(groups, function(g) {
    rownames(g) <- NULL
    if (tna_window_type == "sliding") {
      # Row t marks the states active in rows t - w + 1, ..., t of the data.
      # Lags are of the original rows: lagging the result would widen the
      # window with each step. The first w - 1 rows are incomplete windows.
      lags <- seq_len(window_size - 1L)
      g[cols] <- lapply(cols, function(col) {
        x <- g[[col]]
        active <- Reduce(`|`, lapply(lags, \(k) !is.na(lag_chr(x, k))), !is.na(x))
        ifelse(active, col, NA_character_)
      })
      g <- g[seq_len(nrow(g)) >= window_size, , drop = FALSE]
      g$.window <- floor(seq_len(nrow(g)) - 1L)
    } else {
      g$.window <- floor((seq_len(nrow(g)) - 1L) / window_size)
    }
    g
  }))

  if (aggregate && nrow(work) > 0L) {
    ord <- do.call(order, work[, c(actor_col, session_col, ".window"),
                                drop = FALSE])
    work <- work[ord, , drop = FALSE]
    agg_ids <- .observed_group_id(
      work[, c(actor_col, session_col, ".window"), drop = FALSE],
      context = "actor/session/window"
    )
    work <- do.call(rbind, lapply(.split_first_seen(agg_ids), function(i) {
      g <- work[i, , drop = FALSE]
      out <- g[1L, , drop = FALSE]
      out[cols] <- lapply(cols, function(col) first_non_missing(g[[col]]))
      out
    }))
    rownames(work) <- NULL
  }

  if (is.null(work) || nrow(work) == 0L) {
    result <- as.data.frame(setNames(replicate(0, character(0), simplify = FALSE),
                                     character(0)))
  } else {
    work$.window_grp <- work$.window %/% interval
    work$.window_idx <- work$.window %% interval

    # One row of `work` contributes one entry per code column, in `cols` order.
    # Built in a single stacked allocation rather than a data.frame per row.
    n_code <- length(cols)
    long <- data.frame(
      .actor      = rep(work[[actor_col]], each = n_code),
      .session    = rep(work[[session_col]], each = n_code),
      .window_grp = rep(work$.window_grp, each = n_code),
      .window_idx = rep(work$.window_idx, each = n_code),
      value       = as.character(t(as.matrix(work[, cols, drop = FALSE]))),
      stringsAsFactors = FALSE
    )

    obs_ids <- .observed_group_id(
      long[, c(".actor", ".session", ".window_grp", ".window_idx"),
           drop = FALSE],
      context = "actor/session/window"
    )
    long$.obs <- ave(seq_len(nrow(long)), obs_ids, FUN = seq_along)
    id_ids <- .observed_group_id(
      long[, c(".actor", ".session", ".window_grp"), drop = FALSE],
      context = "actor/session"
    )
    # `.obs` counts within (id group, window_idx), so within one id group each
    # (window_idx, obs) pair -- and therefore each cell name -- occurs once.
    cell_names <- sprintf("W%d_T%d", long$.window_idx, long$.obs)
    wide_list <- lapply(.split_first_seen(id_ids), function(i) {
      row <- as.list(long$value[i])
      names(row) <- cell_names[i]
      as.data.frame(row, stringsAsFactors = FALSE)
    })

    all_cols <- unique(unlist(lapply(wide_list, names), use.names = FALSE))
    wide_list <- lapply(wide_list, function(row) {
      missing <- setdiff(all_cols, names(row))
      row[missing] <- rep(list(NA_character_), length(missing))
      row[, all_cols, drop = FALSE]
    })
    result <- do.call(rbind, wide_list)
    rownames(result) <- NULL
  }

  if (ncol(result) > 0L) {
    result[] <- lapply(result, as.character)
  }

  attr(result, "windowed") <- TRUE
  attr(result, "window_size") <- window_size^(!aggregate)
  attr(result, "window_span") <- length(cols)
  attr(result, "codes") <- cols

  result
}
