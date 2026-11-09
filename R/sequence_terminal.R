# Sequence-with-dropout helpers: surface the terminal-NA pattern in
# wide-format sequence data and turn it into either (a) a tidy
# per-actor endpoint summary, or (b) a chain with an extra absorbing
# state representing dropout / course-end / censoring.

#' Tidy per-actor endpoint summary of a wide-format sequence dataset
#'
#' @description
#' `actor_endpoints()` returns one row per sequence (actor) of a
#' wide-format sequence table, with its first and last observed states, the
#' positions where they occur, the number of observed cells and whether the
#' sequence ends early. A sequence ends early (`dropped_out`) when its last
#' observed state comes before the last column, so every cell after it is
#' missing: the *terminal-NA* pattern of drop-out, course end or censoring.
#' Internal gaps do not count.
#'
#' The flag records only where the observations stop; it cannot say why.
#' A sequence recorded over fewer positions than the others looks the same
#' as one that left. [mark_terminal_state()] turns the terminal-NA cells
#' into an explicit end state, and [chain_structure()] then describes the
#' chain with that state.
#'
#' @param data The sequences: a matrix or data frame with one row per actor
#'   and one column per position, holding state labels, with `NA` for a
#'   missing observation. Without `cols`, every column is read as a
#'   position, so the frame must be one column or a numbered run of
#'   positions (`T1`, `T2`, ...); a frame with other columns beside them,
#'   such as an id, raises an error of class `nestimate_bad_input`. A table
#'   from [mark_terminal_state()] or [mark_first_state()] is read as it is,
#'   because those functions mark their columns as positions.
#' @param cols The columns holding the positions, in order; the other
#'   columns are ignored. A name not in `data` raises an error of class
#'   `nestimate_bad_input`. Default `NULL` (every column).
#' @return A data frame with one row per actor and the columns:
#'   \describe{
#'     \item{`actor`}{The row name when `data` has row names other than
#'       `1`, `2`, ...; otherwise the row number.}
#'     \item{`first_state`, `last_state`}{The first and the last observed
#'       state.}
#'     \item{`first_step`, `last_step`}{The positions (column numbers among
#'       the position columns) of the first and the last observed state.}
#'     \item{`n_observed`}{The number of observed cells.}
#'     \item{`dropped_out`}{`TRUE` when `last_step` is before the last
#'       column, so the row ends in missing cells.}
#'   }
#'   A row with no observed state has `NA` states and steps, `n_observed`
#'   0 and `dropped_out` `FALSE`.
#' @examples
#' actor_endpoints(trajectories) |> head()
#'
#' @seealso [mark_terminal_state()], [chain_structure()]
#' @export
actor_endpoints <- function(data, cols = NULL) {
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  .assert("`data` must be a data frame or a matrix" = is.data.frame(data))
  .check_declared_columns(data, "actor_endpoints()", cols = cols)
  if (is.null(cols)) {
    .refuse_undeclared_roles(data, "actor_endpoints()", whole_frame = TRUE,
                             takes_action = FALSE)
  } else {
    data <- data[, cols, drop = FALSE]
  }

  M <- as.matrix(data)
  storage.mode(M) <- "character"
  n <- nrow(M); k <- ncol(M)

  endpoints <- lapply(seq_len(n), function(i) {
    row <- M[i, ]
    obs <- which(!is.na(row))
    if (length(obs) == 0L) {
      return(list(first_state = NA_character_, last_state = NA_character_,
                  first_step = NA_integer_, last_step = NA_integer_,
                  n_observed = 0L, dropped_out = FALSE))
    }
    first_step <- min(obs); last_step <- max(obs)
    after <- if (last_step < k) row[(last_step + 1L):k] else character(0)
    dropped <- last_step < k && all(is.na(after))
    list(first_state = row[first_step],
         last_state  = row[last_step],
         first_step  = first_step,
         last_step   = last_step,
         n_observed  = length(obs),
         dropped_out = dropped)
  })

  actor_id <- if (!is.null(rownames(M)) && !identical(rownames(M),
                                                       as.character(seq_len(n))))
    rownames(M) else seq_len(n)

  data.frame(
    actor       = actor_id,
    first_state = vapply(endpoints, `[[`, character(1), "first_state"),
    last_state  = vapply(endpoints, `[[`, character(1), "last_state"),
    first_step  = vapply(endpoints, `[[`, integer(1),   "first_step"),
    last_step   = vapply(endpoints, `[[`, integer(1),   "last_step"),
    n_observed  = vapply(endpoints, `[[`, integer(1),   "n_observed"),
    dropped_out = vapply(endpoints, `[[`, logical(1),   "dropped_out"),
    stringsAsFactors = FALSE
  )
}


#' Mark leading-NA cells with an explicit state label
#'
#' @description
#' `mark_first_state()` fills every missing cell before each sequence's
#' first observed state with an explicit state, `"Start"` by default, and
#' returns the wide table. It is the counterpart of [mark_terminal_state()]
#' for *left-censored* data, such as students who join a course at
#' different times: a network built from the result shows into
#' which states the late entrants arrive. Internal and trailing missing
#' cells are left as they are, and so are rows observed from the first
#' position or not at all.
#'
#' The marked state is a transient entry point, not an absorbing state:
#' once a sequence leaves it, it never returns. In a transition network
#' built from the result, its row holds the moves from it, which go to the
#' first observed states and, for sequences with two or more leading
#' missing cells, back to itself.
#'
#' @param data The sequences: a matrix or data frame with one row per actor
#'   and one column per position, holding state labels, with `NA` for a
#'   missing observation. Without `cols`, the frame must be one column or a
#'   numbered run of positions (`T1`, `T2`, ...); other columns beside them
#'   raise an error of class `nestimate_bad_input`.
#' @param state The label written into the leading missing cells, a single
#'   non-empty string. Default `"Start"`. If the label already occurs in
#'   `data`, a warning of class `nestimate_state_collision` is given and `_1`, `_2`, ... is appended until it is
#'   unique, so the marker is never confused with an observed state.
#' @param cols The columns holding the positions, in order; only they are
#'   returned. A name not in `data` raises an error of class
#'   `nestimate_bad_input`. Default `NULL` (every column).
#' @return A data frame of character columns with the shape of `data` (or of
#'   its `cols`), the leading missing cells filled with the marker. The label
#'   used is kept in the `"leading_state"` attribute, which differs from
#'   `state` when it had to be made unique. The columns are marked as
#'   positions, so [build_network()] reads them without `state_cols`.
#' @examples
#' M <- mark_first_state(trajectories, state = "Start")
#' # In a chain built from M, "Start" is a transient entry point.
#'
#' @seealso [mark_terminal_state()], [actor_endpoints()]
#' @export
mark_first_state <- function(data, state = "Start", cols = NULL) {
  .assert("`state` must be a single non-empty string" =
            is.character(state) && length(state) == 1L && nzchar(state))
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  .assert("`data` must be a data frame or a matrix" = is.data.frame(data))
  .check_declared_columns(data, "mark_first_state()", cols = cols)
  if (is.null(cols)) {
    .refuse_undeclared_roles(data, "mark_first_state()", whole_frame = TRUE,
                             takes_action = FALSE)
  } else {
    data <- data[, cols, drop = FALSE]
  }
  if (state %in% unlist(lapply(data, as.character))) {
    warning(warningCondition(
      sprintf("State '%s' already appears in `data`; using a unique replacement label.",
              state),
      class = "nestimate_state_collision", call = NULL))
    base <- state
    i <- 1L
    while (state %in% unlist(lapply(data, as.character))) {
      state <- paste0(base, "_", i); i <- i + 1L
    }
  }

  M <- as.matrix(data)
  storage.mode(M) <- "character"
  n <- nrow(M); k <- ncol(M)
  if (n > 0L && k > 0L) {
    # Cells before each row's first observed state; an all-missing row is left.
    observed   <- !is.na(M)
    first_step <- max.col(observed, ties.method = "first")
    M[rowSums(observed) > 0L & col(M) < first_step] <- state
  }
  out <- as.data.frame(M, stringsAsFactors = FALSE)
  attr(out, "leading_state") <- state
  .mark_positions(out)
}

#' Mark terminal-NA cells with an explicit state label
#'
#' @description
#' `mark_terminal_state()` fills every missing cell after each sequence's
#' last observed state with an explicit state, `"End"` by default, and
#' returns the wide table. This turns right-censored sequences (drop-out,
#' course end) into data for an absorbing-chain model: in a transition
#' network built from the result with [build_network()], the marked state
#' is entered once a sequence stops and never left. Internal and leading
#' missing cells are left as they are, and so are rows observed to the last
#' position or not at all.
#'
#' [chain_structure()] on that network gives the mean time to absorption,
#' the absorption probabilities and each state's classification.
#' [markov_stability()] is not the right summary for an absorbing chain:
#' its stationary distribution collapses onto the absorbing state.
#'
#' @details
#' The marked state is absorbing because each sequence that reaches it
#' stays there to the last column, so its only observed move is to itself
#' (`P[state, state] = 1` in a `"relative"` network). That needs at least
#' one sequence with two or more trailing missing cells: when every
#' sequence stops at most one cell early, the state has no outgoing
#' transitions and its row in the network is zero.
#'
#' Unlike `build_network(end = )`, which adds an end state only in the cell
#' after each sequence's last observed state, with no outgoing edges and no
#' self-loop, this function fills all the trailing cells.
#'
#' @param data The sequences: a matrix or data frame with one row per actor
#'   and one column per position, holding state labels, with `NA` for a
#'   missing observation. Without `cols`, the frame must be one column or a
#'   numbered run of positions (`T1`, `T2`, ...); other columns beside them
#'   raise an error of class `nestimate_bad_input`.
#' @param state The label written into the trailing missing cells, a single
#'   non-empty string. Default `"End"`. If the label already occurs in
#'   `data`, a warning of class `nestimate_state_collision` is given and `_1`, `_2`, ... is appended until it is
#'   unique, so the absorbing marker is never confused with an observed
#'   state.
#' @param cols The columns holding the positions, in order; only they are
#'   returned. A name not in `data` raises an error of class
#'   `nestimate_bad_input`. Default `NULL` (every column).
#' @return A data frame of character columns with the shape of `data` (or of
#'   its `cols`), the trailing missing cells filled with the marker. The
#'   label used is kept in the `"terminal_state"` attribute, which differs
#'   from `state` when it had to be made unique. The columns are marked as
#'   positions, so [build_network()] reads them without `state_cols`.
#' @examples
#' M <- mark_terminal_state(trajectories, state = "Dropout")
#' net <- build_network(M, method = "relative")
#' chain_structure(net)
#'
#' @seealso [actor_endpoints()], [chain_structure()], [build_network()]
#' @export
mark_terminal_state <- function(data, state = "End", cols = NULL) {
  .assert("`state` must be a single non-empty string" =
            is.character(state) && length(state) == 1L && nzchar(state))
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  .assert("`data` must be a data frame or a matrix" = is.data.frame(data))
  .check_declared_columns(data, "mark_terminal_state()", cols = cols)
  if (is.null(cols)) {
    .refuse_undeclared_roles(data, "mark_terminal_state()", whole_frame = TRUE,
                             takes_action = FALSE)
  } else {
    data <- data[, cols, drop = FALSE]
  }
  if (state %in% unlist(lapply(data, as.character))) {
    warning(warningCondition(
      sprintf("State '%s' already appears in `data`; using a unique replacement label.",
              state),
      class = "nestimate_state_collision", call = NULL))
    base <- state
    i <- 1L
    while (state %in% unlist(lapply(data, as.character))) {
      state <- paste0(base, "_", i); i <- i + 1L
    }
  }

  M <- as.matrix(data)
  storage.mode(M) <- "character"
  n <- nrow(M); k <- ncol(M)
  if (n > 0L && k > 0L) {
    # Cells after each row's last observed state; an all-missing row is left.
    observed  <- !is.na(M)
    last_step <- max.col(observed, ties.method = "last")
    M[rowSums(observed) > 0L & col(M) > last_step] <- state
  }
  out <- as.data.frame(M, stringsAsFactors = FALSE)
  attr(out, "terminal_state") <- state
  .mark_positions(out)
}
