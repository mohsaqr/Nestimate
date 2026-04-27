# Sequence-with-dropout helpers: surface the terminal-NA pattern in
# wide-format sequence data and turn it into either (a) a tidy
# per-actor endpoint summary, or (b) a chain with an extra absorbing
# state representing dropout / course-end / censoring.

#' Tidy per-actor endpoint summary of a wide-format sequence dataset.
#'
#' For each actor (row), reports the first and last observed states,
#' the time indices at which they appear, the number of observed
#' steps, and a `dropped_out` flag that is `TRUE` when the actor
#' has a *terminal-NA* pattern (after the final observed step,
#' every remaining cell is `NA`).
#'
#' @param data A wide-format matrix or data.frame where rows are
#'   actors and columns are time steps. Cells are state labels;
#'   `NA` represents missing observations. If `data` is a
#'   `data.frame`, non-character/-factor columns (e.g. an `id`
#'   column) are dropped via the `cols` argument.
#' @param cols Optional character vector of state-column names. If
#'   `NULL` (default) every column is treated as a state column.
#' @return A tidy `data.frame` with one row per actor and columns:
#'   \describe{
#'     \item{`actor`}{Row number (or row name if present).}
#'     \item{`first_state`}{First non-NA state.}
#'     \item{`last_state`}{Last non-NA state.}
#'     \item{`first_step`}{Column index of the first observed state.}
#'     \item{`last_step`}{Column index of the last observed state.}
#'     \item{`n_observed`}{Number of non-NA cells.}
#'     \item{`dropped_out`}{`TRUE` iff every cell after `last_step`
#'       is `NA` and `last_step < ncol(data)`.}
#'   }
#' @examples
#' actor_endpoints(trajectories) |> head()
#'
#' @seealso [mark_terminal_state()], [chain_structure()]
#' @export
actor_endpoints <- function(data, cols = NULL) {
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  stopifnot(is.data.frame(data))
  if (!is.null(cols)) data <- data[, cols, drop = FALSE]

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


#' Mark leading-NA cells with an explicit state label.
#'
#' Mirror of [mark_terminal_state()] for *left-censored* sequence data.
#' Replaces every cell *before* each row's first observed state with
#' the label given by `state`. The resulting chain has a structurally
#' *recurrent* "Start" state that everyone enters from — useful for
#' cohort-entry analyses where students join at different time points
#' and you want a uniform pre-observation marker.
#'
#' @param data A wide-format matrix or data.frame (rows = actors,
#'   cols = time steps) of state labels with `NA` for missing
#'   observations.
#' @param state Character. Label to insert in leading-NA cells.
#'   Default `"Start"`.
#' @param cols Optional state-column names; otherwise all columns.
#' @return A `data.frame` of the same shape as `data` with leading
#'   NAs filled by `state`.
#' @details
#'   Unlike [mark_terminal_state()], the marked state is **not**
#'   absorbing in the resulting transition matrix — every transition
#'   from "Start" goes to one of the original states (the actor's
#'   first observed state), and the "Start" row is row-stochastic
#'   exactly as the data dictates.
#' @examples
#' M <- mark_first_state(trajectories, state = "Start")
#' # In a chain built from M, "Start" is a transient entry point.
#'
#' @seealso [mark_terminal_state()], [actor_endpoints()]
#' @export
mark_first_state <- function(data, state = "Start", cols = NULL) {
  stopifnot(is.character(state), length(state) == 1L, nzchar(state))
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  stopifnot(is.data.frame(data))
  if (!is.null(cols)) data <- data[, cols, drop = FALSE]
  if (state %in% unlist(lapply(data, as.character))) {
    warning(sprintf("State '%s' already appears in `data`; using a unique replacement label.",
                    state), call. = FALSE)
    base <- state
    i <- 1L
    while (state %in% unlist(lapply(data, as.character))) {
      state <- paste0(base, "_", i); i <- i + 1L
    }
  }

  M <- as.matrix(data)
  storage.mode(M) <- "character"
  n <- nrow(M); k <- ncol(M)
  for (i in seq_len(n)) {
    row <- M[i, ]
    obs <- which(!is.na(row))
    if (length(obs) == 0L) next
    first_step <- min(obs)
    if (first_step > 1L) M[i, seq_len(first_step - 1L)] <- state
  }
  out <- as.data.frame(M, stringsAsFactors = FALSE)
  attr(out, "leading_state") <- state
  out
}

#' Mark terminal-NA cells with an explicit state label.
#'
#' Replaces every cell after each row's last observed state with the
#' label given by `state`, leaving non-terminal NAs untouched. The
#' result, passed to [build_network()], yields a Markov chain in
#' which the marked state is **absorbing** by construction
#' (`P[state, state] = 1`).
#'
#' @param data A wide-format matrix or data.frame (rows = actors,
#'   cols = time steps) of state labels with `NA` for missing
#'   observations.
#' @param state Character. Label to insert in terminal-NA cells.
#'   Default `"End"`.
#' @param cols Optional state-column names; otherwise all columns.
#' @return A `data.frame` of the same shape as `data` with terminal
#'   NAs filled by `state`.
#' @details
#'   This is the small piece of pre-processing required to turn
#'   right-censored sequence data into an absorbing-chain model. The
#'   chain on the resulting matrix has one extra state (`state`)
#'   which is structurally absorbing because every cell after the
#'   actor's last observed step has been set to `state` — the chain
#'   stays there forever once entered.
#'
#'   Use [chain_structure()] on the result to compute mean absorption
#'   time, absorption probabilities, and per-state classification.
#'   Note that `markov_stability()` is *not* the right summary for
#'   absorbing chains; its stationary distribution will collapse to
#'   the absorbing state.
#' @examples
#' M <- mark_terminal_state(trajectories, state = "Dropout")
#' net <- build_network(M, method = "relative")
#' chain_structure(net)
#'
#' @seealso [actor_endpoints()], [chain_structure()], [build_network()]
#' @export
mark_terminal_state <- function(data, state = "End", cols = NULL) {
  stopifnot(is.character(state), length(state) == 1L, nzchar(state))
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  stopifnot(is.data.frame(data))
  if (!is.null(cols)) data <- data[, cols, drop = FALSE]
  if (state %in% unlist(lapply(data, as.character))) {
    warning(sprintf("State '%s' already appears in `data`; using a unique replacement label.",
                    state), call. = FALSE)
    base <- state
    i <- 1L
    while (state %in% unlist(lapply(data, as.character))) {
      state <- paste0(base, "_", i); i <- i + 1L
    }
  }

  M <- as.matrix(data)
  storage.mode(M) <- "character"
  n <- nrow(M); k <- ncol(M)
  for (i in seq_len(n)) {
    row <- M[i, ]
    obs <- which(!is.na(row))
    if (length(obs) == 0L) next
    last_step <- max(obs)
    if (last_step < k) M[i, (last_step + 1L):k] <- state
  }
  out <- as.data.frame(M, stringsAsFactors = FALSE)
  attr(out, "terminal_state") <- state
  out
}
