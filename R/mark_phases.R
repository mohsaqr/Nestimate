#' Split each actor's sequence into phases
#'
#' @description
#' `mark_phases()` divides the events of each actor into `n` consecutive
#' phases of nearly equal length and records the phase of every event in a
#' new column. The result can be passed to [build_network()] with
#' `group = "phase"` to estimate one network per phase, for example to
#' compare the first and second halves of each session.
#'
#' The events of an actor are ordered by `time` (or kept in row order) and
#' numbered \eqn{k = 1, \dots, L}. Event \eqn{k} belongs to phase
#' \eqn{\lceil k n / L \rceil}, so the phases of one actor differ in length
#' by at most one event, and with an odd number of events and `n = 2` the
#' middle event belongs to the second phase. An actor with fewer than `n`
#' events has no events in some phases.
#'
#' @param data An event log: a data frame with one row per event.
#' @param actor The column identifying whose events form one sequence (a
#'   session, a student). Default `NULL`, which treats the whole log as one
#'   sequence.
#' @param time The column ordering the events within each actor: numeric, a
#'   `Date` or a `POSIXct`. Default `NULL`, which keeps the row order.
#' @param n The number of phases, a whole number of at least 2. Default `2`.
#' @param labels The labels of the phases, in order, `n` distinct non-empty
#'   strings. Default `NULL`, which uses `"1"`, `"2"`, ..., `"n"`.
#' @param name The name of the new column. Default `"phase"`. A column of
#'   that name already in `data` raises an error of class
#'   `nestimate_bad_input`.
#' @return `data` with one added character column, `name`, holding each
#'   event's phase label; the rows and the other columns are unchanged and
#'   keep their order. An absent `actor` or `time` column, missing values
#'   in either, a `time` column of another type, or an invalid `n` or `labels`
#'   raises an error of class `nestimate_bad_input`.
#' @examples
#' halves <- mark_phases(human_long, actor = "session_id",
#'                       time = "order_in_session", n = 2,
#'                       labels = c("early", "late"))
#' head(halves)
#' net <- build_network(halves, method = "relative", actor = "session_id",
#'                      action = "cluster", time = "order_in_session",
#'                      group = "phase")
#' net
#'
#' @seealso [build_network()] for the per-phase networks,
#'   [sequence_compare()] and [permutation()] for comparing them.
#' @export
mark_phases <- function(data, actor = NULL, time = NULL, n = 2L,
                        labels = NULL, name = "phase") {
  .assert(
    "`data` must be a data frame" = is.data.frame(data),
    "`n` must be a single whole number of at least 2" =
      is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 2 &&
      n == round(n),
    "`name` must be a single non-empty string" =
      is.character(name) && length(name) == 1L && nzchar(name),
    "`actor` must be NULL or one column name" =
      is.null(actor) || (is.character(actor) && length(actor) == 1L),
    "`time` must be NULL or one column name" =
      is.null(time) || (is.character(time) && length(time) == 1L)
  )
  n <- as.integer(n)
  labels <- labels %||% as.character(seq_len(n))
  .assert("`labels` must be `n` distinct non-empty strings" =
            is.character(labels) && length(labels) == n &&
            !anyNA(labels) && all(nzchar(labels)) && !anyDuplicated(labels))
  .check_declared_columns(data, "mark_phases()", actor = actor, time = time)
  if (name %in% names(data)) {
    stop(errorCondition(sprintf(paste0(
      "mark_phases(): `data` already has a column '%s'; choose another ",
      "`name`."), name), class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(actor)) {
    .assert("`actor` must have no missing values" = !anyNA(data[[actor]]))
  }
  if (!is.null(time)) {
    tv <- data[[time]]
    .assert("`time` must be numeric, a Date or a POSIXct" =
              is.numeric(tv) || inherits(tv, c("Date", "POSIXt")),
            "`time` must have no missing values" = !anyNA(tv))
  }

  n_rows <- nrow(data)
  actor_id <- if (is.null(actor)) rep(1L, n_rows) else
    match(data[[actor]], unique(data[[actor]]))
  time_vals <- if (is.null(time)) seq_len(n_rows) else data[[time]]
  # events of each actor in time order; ties keep the row order
  ord <- order(actor_id, time_vals, seq_len(n_rows))
  position <- integer(n_rows)
  position[ord] <- stats::ave(seq_len(n_rows), actor_id[ord], FUN = seq_along)
  length_of_actor <- tabulate(actor_id)[actor_id]
  phase <- as.integer(ceiling(position * n / length_of_actor))

  data[[name]] <- labels[phase]
  stopifnot("one phase per event" = !anyNA(data[[name]]))
  data
}
