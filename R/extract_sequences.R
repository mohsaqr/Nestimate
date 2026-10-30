#' Extract Sequences in Wide Format
#'
#' Returns the sequences behind any Nestimate input as a wide data frame, one
#' row per sequence and one column per position. Accepts the output of
#' \code{\link{prepare}()}, a long event log, a wide data frame or character
#' matrix, a \code{netobject}, a grouped \code{netobject_group}, or a
#' \code{tna} model. Model objects return their stored, decoded states.
#'
#' @param x Sequence data: a \code{nestimate_data} object from
#'   \code{\link{prepare}()}; a long event log (supply \code{action});
#'   a wide \code{data.frame} or character matrix; a \code{netobject},
#'   \code{netobject_group}, \code{cograph_network} or \code{tna} built from
#'   sequence data.
#' @param cols Character vector of state columns for wide data frames. When
#'   \code{NULL}, numeric columns of a mixed data frame are treated as
#'   covariates and the remaining columns as states. Ignored for model
#'   objects.
#' @param actor,action,time Column names of a long event log, passed to
#'   \code{\link{prepare}()}. Supplying \code{action} marks \code{x} as long.
#' @param meta Logical. When \code{TRUE}, the covariate columns that belong
#'   to each sequence (actor, session, group, metadata) are placed before the
#'   state columns. Default \code{FALSE}.
#' @param ... Further arguments passed to \code{\link{prepare}()} for long
#'   input (e.g. \code{time_threshold}, \code{order}).
#'
#' @return A \code{data.frame} with one row per sequence. State columns are
#'   character, with \code{NA} after a sequence ends; their names are the
#'   position names of the input (\code{T1}, \code{T2}, ... for prepared or
#'   long data). With \code{meta = TRUE} the covariate columns come first.
#'   A grouped network adds its grouping column. Input that has no sequence
#'   data raises \code{nestimate_bad_input}.
#'
#' @examples
#' extract_sequences(trajectories)
#'
#' extract_sequences(group_regulation_long, actor = "Actor",
#'                   action = "Action", time = "Time", meta = TRUE)
#'
#' net <- build_tna(trajectories)
#' extract_sequences(net)
#'
#' @seealso \code{\link{prepare}}, \code{\link{sequence_patterns}}
#' @export
extract_sequences <- function(x, cols = NULL, actor = NULL, action = NULL,
                              time = NULL, meta = FALSE, ...) {
  stopifnot("`meta` must be TRUE or FALSE" = isTRUE(meta) || isFALSE(meta))
  seqs <- .as_sequences(x, cols = cols, actor = actor, action = action,
                        time = time, ...)
  states <- as.data.frame(seqs$states, stringsAsFactors = FALSE)
  if (!meta || ncol(seqs$units) == 0L) return(states)
  out <- cbind(seqs$units, states)
  rownames(out) <- NULL
  out
}


#' Resolve any sequence input to states plus a unit table
#'
#' @return A list with \code{states} (character matrix, one row per
#'   sequence), \code{units} (data.frame of per-sequence covariates, same
#'   rows; possibly zero columns) and \code{alphabet} (sorted states).
#' @noRd
# Methods whose $data are sequences of states.
.sequence_methods <- c("relative", "frequency", "co_occurrence", "attention",
                       "ngram", "gap", "reverse", "entropy")

# Name of the column that holds a grouped network's group labels: the
# grouping column, or for several the joined names ("Achiever-Course"),
# matching the joined group labels ("High-A").
.group_col_name <- function(x) {
  paste(attr(x, "group_col") %||% "group", collapse = "-")
}

.as_sequences <- function(x, cols = NULL, actor = NULL, action = NULL,
                          time = NULL, exclude = NULL, ...) {
  if (is.data.frame(x) && !inherits(x, "nestimate_data") && !is.null(action)) {
    x <- prepare(x, actor = actor, action = action, time = time, ...)
  }

  if (inherits(x, "nestimate_data")) {
    states <- x$sequence_data
    units <- x$meta_data
  } else if (inherits(x, "netobject_group")) {
    group_col <- .group_col_name(x)
    parts <- lapply(x, .as_sequences)
    width <- max(vapply(parts, \(p) ncol(p$states), integer(1)))
    states <- do.call(rbind, lapply(parts, \(p) .pad_states(p$states, width)))
    units <- do.call(rbind, Map(\(p, g) {
      u <- p$units
      u[[group_col]] <- rep(g, nrow(p$states))
      u
    }, parts, names(x)))
  } else if (inherits(x, c("netobject", "cograph_network", "tna"))) {
    # an association or window network carries numeric or one-hot data, not
    # states; reading it as sequences would give meaningless "patterns"
    if (inherits(x, "netobject") && !is.null(x$method) &&
        !x$method %in% .sequence_methods) {
      stop(errorCondition(sprintf(paste(
        "`x` is a \"%s\" network, whose data are not sequences of states;",
        "build it with a transition method (\"relative\", \"frequency\",",
        "...) or pass the sequences."), x$method),
        class = "nestimate_bad_input", call = NULL))
    }
    if (is.null(x$data)) {
      stop(errorCondition(
        "`x` holds no sequence data; build it from sequences, not a matrix.",
        class = "nestimate_bad_input", call = NULL))
    }
    states <- .coerce_sequence_input(x)
    if (!is.data.frame(states)) {
      stop(errorCondition("`x` holds no sequence data.",
                          class = "nestimate_bad_input", call = NULL))
    }
    units <- if (inherits(x, "netobject")) x$metadata else NULL
  } else if (is.data.frame(x) || is.matrix(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    missing_cols <- setdiff(c(cols, exclude), names(df))
    if (length(missing_cols) > 0L) {
      stop(errorCondition(
        sprintf("column(s) not found in `x`: %s.",
                paste(missing_cols, collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
    state_cols <- .build_sequence_cols(df, cols = cols, exclude = exclude) %||%
      names(df)
    states <- df[state_cols]
    units <- df[setdiff(names(df), state_cols)]
  } else {
    stop(errorCondition(
      paste("`x` must be prepared data, a long or wide data frame, a",
            "character matrix, a netobject, a netobject_group or a tna",
            "model."),
      class = "nestimate_bad_input", call = NULL))
  }

  states <- .states_matrix(states)
  units <- if (is.null(units)) {
    data.frame(row.names = seq_len(nrow(states)))
  } else {
    as.data.frame(units, stringsAsFactors = FALSE)
  }
  rownames(units) <- NULL
  stopifnot("sequence and covariate rows must align" =
              nrow(units) == nrow(states))
  if (nrow(states) == 0L || all(is.na(states))) {
    stop(errorCondition("`x` contains no observed states.",
                        class = "nestimate_bad_input", call = NULL))
  }
  list(states = states, units = units,
       alphabet = sort(unique(states[!is.na(states)])))
}


#' Character state matrix with void and missing markers as NA
#'
#' Same markers as \code{build_network()}: empty strings and TraMineR's void
#' (\code{\%}) and missing (\code{*}) symbols, so an \code{stslist} reads
#' with its padding removed.
#' @noRd
.states_matrix <- function(states) {
  df <- as.data.frame(states, stringsAsFactors = FALSE)
  # as.matrix() on mixed columns pads numbers with format(); convert first
  df[] <- lapply(df, as.character)
  m <- .clean_states(as.matrix(df))
  rownames(m) <- NULL
  m
}


#' Pad a state matrix with NA columns up to `width`
#' @noRd
.pad_states <- function(m, width) {
  if (ncol(m) >= width) return(m)
  pad <- matrix(NA_character_, nrow(m), width - ncol(m))
  extra <- setdiff(paste0("T", seq_len(width)), colnames(m))
  colnames(pad) <- utils::head(extra, width - ncol(m))
  cbind(m, pad)
}
