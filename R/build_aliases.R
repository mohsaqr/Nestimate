# ---- The aliases of build_network(), one per method ----
#
# Each alias is build_network() with its method fixed. It names the arguments
# of its method's kind of network, in one order for every alias of that kind
# (the leading positions kept from earlier releases, so positional calls still
# bind as they did), and passes on only the arguments its caller set
# (.forward_to_build_network()).

# ---- Sequences and event logs ----

#' Estimate a transition network analysis (TNA) model
#'
#' @description
#' `build_network(method = "relative")`, or for short `build_tna()`,
#' estimates a transition network analysis (TNA) model from sequences of
#' states. The edge from state A to state B is the probability that B comes
#' next after A: the A-to-B transitions divided by all transitions out of A
#' (Saqr et al., 2025).
#'
#' A state's edges compare within its row, not across rows. The weights are
#' point estimates: check how far they move under resampling with
#' [bootstrap_network()], and compare groups with [permutation()].
#'
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, an event log (one row per event, with `action`
#'   set), one-hot code columns (declared with `codes`), or a square matrix
#'   of transition counts or weights, which is row-normalised (every row
#'   needs a positive sum).
#' @param start A start state put before every sequence, so the network
#'   shows where sequences begin (an edge from it to each sequence's first
#'   state). `FALSE` (default) adds none, `TRUE` adds one labelled
#'   `"Start"`, and a single string adds one with that label. One-hot data,
#'   whose rows are time points rather than sequences, raises an error.
#' @param end An end state put in the cell after each sequence's last
#'   observed state, so the network shows where sequences stop. It has no
#'   outgoing edges and no self-loop, unlike [mark_terminal_state()], which
#'   fills all trailing missing cells with an absorbing state. `FALSE`
#'   (default), `TRUE` (`"End"`) or a single string, as `start`.
#' @param actor The column naming who each row belongs to (a person, a
#'   team). In an event log, each actor's events form a sequence, split
#'   further by `session` or `time_threshold`; in one-hot data, each actor's
#'   rows. Default `NULL`.
#' @param action The column holding each event's state; setting it declares
#'   the data an event log. Column names give no role: `action`, `actor`,
#'   `time` and `session` are read only when set. A column called `Action`
#'   with nothing declared raises an error of class `nestimate_bad_input`,
#'   and an undeclared time or session column beside `action` is named in a
#'   message of class `nestimate_column_role`. Default `NULL`.
#' @param time The timestamp column of an event log. Events are ordered by
#'   it, and a gap longer than `time_threshold` starts a new sequence.
#'   Default `NULL`.
#' @param session The column identifying sessions in an event log; each
#'   session is its own sequence, also without `actor`. Default `NULL`.
#' @param order The column that orders events sharing a timestamp; without
#'   it, those events keep their row order. Default `NULL`.
#' @param time_threshold The longest gap, in seconds, between two events of
#'   the same sequence in an event log; a longer gap starts a new sequence.
#'   `FALSE` turns the split off, so each actor (or actor-session) is one
#'   sequence. Default `900`.
#' @param timezone The Olson time zone for timestamps without an offset
#'   (passed to [prepare()]); timestamps with an offset, such as `...Z` or
#'   `+02:00`, are converted from it. Default `"UTC"`.
#' @param state_cols The columns of a wide table that are sequence
#'   positions; every other column goes to the network's metadata. Declared
#'   columns are positions even when they hold only 0 and 1. With nothing
#'   declared, the frame must be one column or hold position columns named
#'   as a numbered run (a shared prefix and a number, `T1`, `T2`, ..., or
#'   `1`, `2`, ...); numeric columns beside the run and character or factor
#'   columns outside it are then covariates. Any other frame, and columns of
#'   only 0 and 1, raise an error of class `nestimate_bad_input`: nothing
#'   says whether they are positions or one-hot codes. Default `NULL`.
#' @param metadata_cols Columns to keep as the network's metadata
#'   (covariates), never as positions; they cannot overlap `state_cols`. A
#'   column left out of a numbered run of positions is named in a message of
#'   class `nestimate_column_role`, and listing it here or in `state_cols`
#'   settles it. Default `NULL`.
#' @param codes The one-hot code columns (one per state, 1 when the state is
#'   active); setting them declares one-hot data. Default `NULL`.
#' @param window_size For one-hot data: how many consecutive rows form one
#'   window. Default `3`. Windows on sequence data are set through `params`
#'   instead.
#' @param mode For one-hot data: `"non-overlapping"` (default) or
#'   `"overlapping"` windows; [build_wtna()] describes how overlapping windows
#'   are counted.
#' @param params Settings for the estimator, as a named list. Transitions
#'   are counted with tna's `weighted` and `concat` options, and
#'   `window_size` (with `window_span`, or `windowed = TRUE` for a size of 2)
#'   counts them between consecutive non-overlapping blocks of that many
#'   positions, as `tna::tna()` does; `mode = "overlapping"` then raises a
#'   warning of class `nestimate_ignored_mode`. The network keeps these
#'   settings, and the resampling verbs such as [bootstrap_network()]
#'   re-estimate with them. Default `list()`.
#' @param ... Further arguments for [build_network()]: `format` (`"auto"`,
#'   `"wide"`, `"long"` or `"onehot"`, to force a reading of `data`), and for
#'   event logs [prepare()]'s time options `unix_time_unit`, `is_unix_time`
#'   and `custom_format`. `method` cannot be given: each alias has its own
#'   (error of class `nestimate_bad_input`).
#' @inheritParams build_network
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()] and the weights with
#'   [extract_transition_matrix()], and draw it with `cograph::splot()`.
#' @references
#' Saqr, M., L\enc{ó}{o}pez-Pernas, S., T\enc{ö}{o}rm\enc{ä}{a}nen, T.,
#' Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition network
#' analysis: A novel framework for modeling, visualizing, and identifying the
#' temporal patterns of learners and learning processes. In
#' \emph{Proceedings of the 15th International Learning Analytics and
#' Knowledge Conference (LAK25)} (pp. 351--361). ACM.
#' \doi{10.1145/3706468.3706513}
#'
#' Tikka, S., L\enc{ó}{o}pez-Pernas, S., & Saqr, M. (2025). tna: An R
#' package for transition network analysis. \emph{Applied Psychological
#' Measurement}. \doi{10.1177/01466216251348840}
#' @seealso [build_network()] for every method; [build_ftna()],
#'   [build_atna()] and [build_cna()] for the other transition networks;
#'   [bootstrap_network()] and [permutation()] for inference.
#' @examples
#' build_network(trajectories, method = "relative")
#' build_tna(trajectories)
#'
#' build_tna(group_regulation_long, actor = "Actor", action = "Action",
#'           time = "Time", group = "Course", threshold = 0.05)
#' @export
build_tna <- function(data, start = FALSE, end = FALSE,
                      actor = NULL, action = NULL, time = NULL,
                      session = NULL, ..., order = NULL, time_threshold = 900,
                      timezone = "UTC", state_cols = NULL,
                      metadata_cols = NULL, codes = NULL, window_size = 3L,
                      mode = c("non-overlapping", "overlapping"),
                      group = NULL, scaling = NULL, threshold = 0,
                      labels = NULL, params = list()) {
  .forward_to_build_network("build_tna", "relative", ...names())
}

#' Estimate a frequency transition network analysis (FTNA) model
#'
#' @description
#' `build_network(method = "frequency")`, or for short `build_ftna()`,
#' estimates a frequency transition network analysis (FTNA) model from
#' sequences of states. The edge from state A to state B is the number of
#' times B comes right after A (Saqr et al., 2025).
#'
#' Counts grow with the number and length of sequences, so they compare
#' within one network; for transition probabilities use [build_tna()]. A
#' matrix given as `data` must hold non-negative whole counts. Check how far
#' the counts move under resampling with [bootstrap_network()].
#'
#' @inheritParams build_tna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()] and the counts with [extract_transition_matrix()],
#'   and draw it with `cograph::splot()`.
#' @references
#' Saqr, M., L\enc{ó}{o}pez-Pernas, S., T\enc{ö}{o}rm\enc{ä}{a}nen, T.,
#' Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition network
#' analysis: A novel framework for modeling, visualizing, and identifying the
#' temporal patterns of learners and learning processes. In
#' \emph{Proceedings of the 15th International Learning Analytics and
#' Knowledge Conference (LAK25)} (pp. 351--361). ACM.
#' \doi{10.1145/3706468.3706513}
#' @seealso [build_network()] for every method; [build_tna()],
#'   [build_atna()] and [build_cna()] for the other transition networks.
#' @examples
#' build_network(trajectories, method = "frequency")
#' build_ftna(trajectories)
#'
#' build_ftna(group_regulation_long, actor = "Actor", action = "Action",
#'            time = "Time", group = "Achiever")
#' @export
build_ftna <- function(data, start = FALSE, end = FALSE,
                       actor = NULL, action = NULL, time = NULL,
                       session = NULL, ..., order = NULL, time_threshold = 900,
                       timezone = "UTC", state_cols = NULL,
                       metadata_cols = NULL, codes = NULL, window_size = 3L,
                       mode = c("non-overlapping", "overlapping"),
                       group = NULL, scaling = NULL, threshold = 0,
                       labels = NULL, params = list()) {
  .forward_to_build_network("build_ftna", "frequency", ...names())
}

#' Estimate an attention transition network analysis (ATNA) model
#'
#' @description
#' `build_network(method = "attention")`, or for short `build_atna()`,
#' estimates an attention transition network analysis (ATNA) model, in which
#' every later state in a sequence counts, not only the next one. A pair of
#' positions d steps apart adds `exp(-d / lambda)` to the edge from the
#' earlier state to the later one, so near successors weigh most (Saqr et
#' al., 2025).
#'
#' The weights are summed decay weights, not counts or probabilities, so
#' they compare within one network; add `scaling = "normalize"` for rows
#' that sum to 1. Attention networks read sequences and event logs, not
#' one-hot data.
#'
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, or an event log (one row per event, with `action`
#'   set).
#' @param lambda The decay length: larger values let distant states count
#'   more. Default `1`.
#' @param direction `"forward"` (default) links each state to the later
#'   ones, `"backward"` to the earlier ones, and `"both"` to all.
#' @param decay A function of `(ti, tj, lambda)` returning one weight per
#'   pair of positions, in place of the exponential decay. Default `NULL`.
#' @param params Settings for the estimator, as a named list, as for
#'   [build_tna()]; also `duration` (one value per position, accumulated
#'   into times) or `time_matrix` (one row per sequence, one column per
#'   position), which measure the distance in time instead of positions.
#'   Default `list()`.
#' @inheritParams build_tna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' Saqr, M., L\enc{ó}{o}pez-Pernas, S., T\enc{ö}{o}rm\enc{ä}{a}nen, T.,
#' Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition network
#' analysis: A novel framework for modeling, visualizing, and identifying the
#' temporal patterns of learners and learning processes. In
#' \emph{Proceedings of the 15th International Learning Analytics and
#' Knowledge Conference (LAK25)} (pp. 351--361). ACM.
#' \doi{10.1145/3706468.3706513}
#' @seealso [build_network()] for every method; [build_tna()] for the
#'   transitions to the next state only.
#' @examples
#' build_network(trajectories, method = "attention")
#' build_atna(trajectories)
#'
#' build_atna(trajectories, lambda = 2, direction = "both")
#' @export
build_atna <- function(data, start = FALSE, end = FALSE,
                       actor = NULL, action = NULL, time = NULL,
                       session = NULL, ..., order = NULL, time_threshold = 900,
                       timezone = "UTC", state_cols = NULL,
                       metadata_cols = NULL, lambda = 1,
                       direction = c("forward", "backward", "both"),
                       decay = NULL, group = NULL, scaling = NULL,
                       threshold = 0, labels = NULL, params = list()) {
  .forward_to_build_network("build_atna", "attention", ...names())
}

#' Estimate a co-occurrence network analysis (CNA) model
#'
#' @description
#' `build_network(method = "co_occurrence")`, or for short `build_cna()`,
#' estimates a co-occurrence network analysis (CNA) model: an undirected
#' network of states that occur in the same sequence. The edge between state
#' A and state B counts the pairs of positions in one sequence that hold A
#' and B, in either order (Saqr et al., 2025).
#'
#' An edge says how often two states occur together, not which came first;
#' for order, use [build_tna()]. One-hot codes are counted in windows, as
#' [build_wtna()] does.
#'
#' @param type `"frequency"` (default) keeps the counts; `"relative"`
#'   divides each row by its sum, which makes the network directed.
#' @inheritParams build_tna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' Saqr, M., L\enc{ó}{o}pez-Pernas, S., T\enc{ö}{o}rm\enc{ä}{a}nen, T.,
#' Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition network
#' analysis: A novel framework for modeling, visualizing, and identifying the
#' temporal patterns of learners and learning processes. In
#' \emph{Proceedings of the 15th International Learning Analytics and
#' Knowledge Conference (LAK25)} (pp. 351--361). ACM.
#' \doi{10.1145/3706468.3706513}
#' @seealso [build_network()] for every method; [build_tna()] for the order
#'   of states.
#' @examples
#' build_network(trajectories, method = "co_occurrence")
#' build_cna(trajectories)
#'
#' build_cna(trajectories, scaling = "max", threshold = 0.1)
#' @export
build_cna <- function(data, start = FALSE, end = FALSE,
                      actor = NULL, action = NULL, time = NULL,
                      session = NULL, ..., order = NULL, time_threshold = 900,
                      timezone = "UTC", state_cols = NULL,
                      metadata_cols = NULL, codes = NULL, window_size = 3L,
                      mode = c("non-overlapping", "overlapping"),
                      type = c("frequency", "relative"), group = NULL,
                      scaling = NULL, threshold = 0, labels = NULL,
                      params = list()) {
  .forward_to_build_network("build_cna", "co_occurrence", ...names())
}

#' Estimate an n-gram transition network analysis model
#'
#' @description
#' `build_network(method = "ngram")`, or for short `build_ngram()`,
#' estimates an n-gram transition network analysis model from windows of
#' `n_gram` consecutive states. Each window adds each of its adjacent pairs
#' once, so a pair inside a sequence is counted up to `n_gram - 1` times, as
#' `tna::build_model(type = "n-gram")` counts it.
#'
#' The weights are raw counts; add `scaling = "normalize"` for row
#' probabilities. With `n_gram = 2` the network equals [build_ftna()].
#'
#' @param n_gram The window length, a whole number of at least 2. Default
#'   `2`.
#' @inheritParams build_atna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' Tikka, S., L\enc{ó}{o}pez-Pernas, S., & Saqr, M. (2025). tna: An R
#' package for transition network analysis. \emph{Applied Psychological
#' Measurement}. \doi{10.1177/01466216251348840}
#' @seealso [build_network()] for every method; [build_gap()] and
#'   [build_reverse()] for the other tna model types.
#' @examples
#' build_network(trajectories, method = "ngram", n_gram = 3)
#' build_ngram(trajectories, n_gram = 3)
#' @export
build_ngram <- function(data, start = FALSE, end = FALSE,
                        actor = NULL, action = NULL, time = NULL,
                        session = NULL, order = NULL, time_threshold = 900,
                        timezone = "UTC", state_cols = NULL,
                        metadata_cols = NULL, n_gram = 2L, group = NULL,
                        scaling = NULL, threshold = 0, labels = NULL,
                        params = list(), ...) {
  .forward_to_build_network("build_ngram", "ngram", ...names())
}

#' Estimate a gap transition network analysis model
#'
#' @description
#' `build_network(method = "gap")`, or for short `build_gap()`, estimates a
#' gap transition network analysis model, which links each state to the
#' states up to `max_gap + 1` positions later, each pair weighted by one
#' over its distance, as `tna::build_model(type = "gap")` does.
#'
#' The weights are raw weighted counts; add `scaling = "normalize"` for row
#' probabilities.
#'
#' @param max_gap How many states may lie between the two states of a pair,
#'   a whole number of at least 0. Default `1`.
#' @inheritParams build_atna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' Tikka, S., L\enc{ó}{o}pez-Pernas, S., & Saqr, M. (2025). tna: An R
#' package for transition network analysis. \emph{Applied Psychological
#' Measurement}. \doi{10.1177/01466216251348840}
#' @seealso [build_network()] for every method; [build_ngram()] and
#'   [build_reverse()] for the other tna model types.
#' @examples
#' build_network(trajectories, method = "gap", max_gap = 2)
#' build_gap(trajectories, max_gap = 2)
#' @export
build_gap <- function(data, start = FALSE, end = FALSE,
                      actor = NULL, action = NULL, time = NULL,
                      session = NULL, order = NULL, time_threshold = 900,
                      timezone = "UTC", state_cols = NULL,
                      metadata_cols = NULL, max_gap = 1L, group = NULL,
                      scaling = NULL, threshold = 0, labels = NULL,
                      params = list(), ...) {
  .forward_to_build_network("build_gap", "gap", ...names())
}

#' Estimate a reverse transition network analysis model
#'
#' @description
#' `build_network(method = "reverse")`, or for short `build_reverse()`,
#' estimates a reverse (reply) transition network analysis model: the edge
#' from state A to state B counts the times A comes right after B, so each
#' state points to the one it answers. It is the transpose of the frequency
#' network, as `tna::build_model(type = "reverse")` computes it.
#'
#' The weights are raw counts; add `scaling = "normalize"` for row
#' probabilities. [build_latency()] weights the same replies by their speed.
#'
#' @param weighted Whether each pair is weighted by one over the length of
#'   its sequence. Default `FALSE`.
#' @inheritParams build_atna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' Tikka, S., L\enc{ó}{o}pez-Pernas, S., & Saqr, M. (2025). tna: An R
#' package for transition network analysis. \emph{Applied Psychological
#' Measurement}. \doi{10.1177/01466216251348840}
#' @seealso [build_network()] for every method; [build_latency()] for
#'   replies weighted by speed.
#' @examples
#' build_network(trajectories, method = "reverse")
#' build_reverse(trajectories)
#' @export
build_reverse <- function(data, start = FALSE, end = FALSE,
                          actor = NULL, action = NULL, time = NULL,
                          session = NULL, order = NULL, time_threshold = 900,
                          timezone = "UTC", state_cols = NULL,
                          metadata_cols = NULL, weighted = FALSE,
                          group = NULL, scaling = NULL, threshold = 0,
                          labels = NULL, params = list(), ...) {
  .forward_to_build_network("build_reverse", "reverse", ...names())
}

# ---- Message logs ----

#' Estimate a latency network analysis model
#'
#' @description
#' `build_network(method = "latency")`, or for short `build_latency()`,
#' estimates a latency network analysis model from a message log. The edge
#' from sender A to sender B is the average weight of A's replies to B,
#' each reply weighted by its speed: 1 for an immediate reply, falling
#' linearly to `floor` at `cap` (Templeton et al., 2022).
#'
#' Its edges are the replies of [build_reverse()]: within each sequence (an
#' `actor`'s events, split by `session` or `time_threshold`), each event
#' points from its sender, the responder, to the sender of the event before
#' it. A reply that takes \eqn{t} seconds weighs
#' \eqn{1 - (1 - f) \min(t / c, 1)}{1 - (1 - f) * min(t / c, 1)}, with
#' \eqn{f} the `floor` and \eqn{c} the `cap` in seconds, and the weights of
#' a pair are averaged, or summed with `aggregate = "sum"`. A sender who posts twice in
#' a row is a self-loop on the diagonal, weighted by the same rule. Shifting
#' every time by a constant leaves the weights unchanged, and so does a
#' change of time unit when the cap is a percentile. [extract_edges()] adds
#' each edge's number of replies (`n_replies`) and their median latency in
#' seconds (`median_latency`). `start`, `end` and
#' `params = list(weighted = TRUE)` raise an error of class
#' `nestimate_bad_input`. [bootstrap_network()] resamples latency networks
#' by their sequences and judges stability on the speed score; the other
#' resampling verbs, such as [permutation()], raise an error of class
#' `nestimate_latency_resampling` on one.
#'
#' In [build_network()] the four settings are named `latency_floor`,
#' `latency_cap`, `latency_aggregate` and `latency_scope`.
#'
#' @param data A message log: one row per message, `action` the sender.
#' @param action The column holding each message's sender.
#' @param time The timestamp column. Events are ordered by it, and a gap
#'   longer than `time_threshold` starts a new sequence.
#' @param floor The weight of a reply that takes `cap` or longer, a single
#'   number from 0 to 1. `1` ignores time (with `aggregate = "sum"`, the
#'   counts of [build_reverse()]); `0` gives replies at or past the cap no
#'   weight.
#'   Default `0.5`.
#' @param cap The time at which a reply's weight reaches `floor`. A number in
#'   (0, 1] is a percentile of the latencies of the replies to another
#'   sender, pooled over all of `data` (with `group`, over every group, so
#'   the groups share one cap); a string such as `"120s"` is a number of
#'   seconds. The resolved cap, in seconds, is stored on the network.
#'   Default `0.95`.
#' @param aggregate How the replies of a pair of senders make one weight:
#'   `"mean"` (default) averages their weights, the pair's typical reply
#'   speed; `"sum"` adds them, so an edge grows with both the speed and the
#'   number of replies.
#' @param scope With a percentile `cap`: the units that get their own cap.
#'   `"global"` (default) takes one cap from all of `data`; `"actor"` one per
#'   value of `actor`; `"session"` one per sequence; `"responder"` one per
#'   responding sender. A unit with fewer than 5 replies to another sender
#'   takes the global cap, and a message of class
#'   `nestimate_latency_fallback` says how many units did. A cap in seconds
#'   applies to every reply, so a scope other than `"global"` is then
#'   ignored with a message of class `nestimate_ignored_argument`.
#' @param ... Further arguments for [build_network()]: [prepare()]'s time
#'   options `unix_time_unit`, `is_unix_time` and `custom_format`. `method`
#'   cannot be given (error of class `nestimate_bad_input`).
#' @inheritParams build_tna
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it. Beside the weights it keeps `initial`,
#'   `n_replies` and `median_latency` (matrices laid out as the weights),
#'   `latency_cap_seconds` (the resolved global cap), `latency_unit_caps`
#'   (the per-unit caps of a non-global `scope`) and `latency_time`, the
#'   event times in seconds.
#' @references
#' Templeton, E. M., Chang, L. J., Reynolds, E. A., Cone LeBeaumont, M. D.,
#' & Wheatley, T. (2022). Fast response times signal social connection in
#' conversation. \emph{Proceedings of the National Academy of Sciences},
#' 119(4), e2116915119. \doi{10.1073/pnas.2116915119}
#' @seealso [build_network()] for every method; [build_reverse()] for the
#'   unweighted replies.
#' @examples
#' build_network(group_regulation_long, method = "latency", actor = "Actor",
#'               action = "Action", time = "Time")
#' build_latency(group_regulation_long, action = "Action", time = "Time",
#'               actor = "Actor")
#' @export
build_latency <- function(data, action, time, actor = NULL, session = NULL,
                          order = NULL, time_threshold = 900,
                          timezone = "UTC", floor = 0.5, cap = 0.95,
                          aggregate = c("mean", "sum"),
                          scope = c("global", "actor", "session",
                                    "responder"),
                          group = NULL, scaling = NULL, threshold = 0,
                          labels = NULL, ...) {
  .forward_to_build_network("build_latency", "latency", ...names(),
                            rename = c(floor = "latency_floor",
                                       cap = "latency_cap",
                                       aggregate = "latency_aggregate",
                                       scope = "latency_scope"))
}

# ---- One-hot codes ----

#' Estimate a window transition network analysis (WTNA) model
#'
#' @description
#' `build_network(method = "wtna")`, or for short `build_wtna()`, estimates
#' a window transition network analysis (WTNA) model from one-hot codes. The
#' rows are cut into windows of `window_size` rows, and the edge from code A
#' to code B counts the times B is active in the window after one where A is
#' active.
#'
#' Windows do not cross actors (`actor`) or sessions (`session`).
#' `build_network(method = "wtna_cooccurrence")` counts the codes active in
#' the same window instead; [wtna()] returns both networks at once.
#'
#' @param data One-hot data: one row per time point and one column of 0 and
#'   1 per code.
#' @param codes The code columns (one per state, 1 when the state is
#'   active). Default `NULL`: every column outside the other roles.
#' @param window_size How many consecutive rows form one window. Default
#'   `3`.
#' @param mode `"non-overlapping"` (default) or `"overlapping"` windows:
#'   overlapping ones slide by one row, as
#'   `tna::import_onehot(window_type = "sliding")`, and are counted in
#'   blocks of `window_size`.
#' @param type `"frequency"` (default) keeps the counts; `"relative"`
#'   divides each row by its sum.
#' @param actor The column naming who each row belongs to; windows do not
#'   cross actors. Default `NULL`.
#' @param session The column identifying sessions; windows do not cross
#'   sessions. Default `NULL`.
#' @param metadata_cols Columns to keep as the network's metadata, never
#'   read as codes. Default `NULL`.
#' @param params Settings for the estimator, as a named list. The network
#'   keeps them. Default `list()`.
#' @param ... Further arguments for [build_network()]. `method` cannot be
#'   given (error of class `nestimate_bad_input`).
#' @inheritParams build_network
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @seealso [build_network()] for every method; [wtna()] for the transition
#'   and co-occurrence networks together; [action_to_onehot()] to turn an
#'   event log into one-hot codes.
#' @examples
#' codes <- c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")
#' build_network(learning_activities, method = "wtna", codes = codes,
#'               actor = "student")
#' build_wtna(learning_activities, codes = codes, actor = "student")
#' @export
build_wtna <- function(data, codes = NULL, window_size = 3L,
                       mode = c("non-overlapping", "overlapping"),
                       type = c("frequency", "relative"), actor = NULL,
                       session = NULL, metadata_cols = NULL, group = NULL,
                       scaling = NULL, threshold = 0, labels = NULL,
                       params = list(), ...) {
  .forward_to_build_network("build_wtna", "wtna", ...names())
}

# ---- Numeric variables ----

#' Estimate a correlation network
#'
#' @description
#' `build_network(method = "cor")`, or for short `build_cor()`, estimates a
#' correlation network of numeric variables. The edge between two variables
#' is their correlation (Epskamp & Fried, 2018).
#'
#' A correlation edge is a marginal association: it can come from a shared
#' link with a third variable. To separate direct from indirect links, use
#' [build_pcor()] or [build_glasso()]. Check the stability of the edges with
#' [bootstrap_network()].
#'
#' @param data The variables: a data frame with one row per observation and
#'   one column per variable, or a square symmetric correlation or
#'   covariance matrix, which needs the sample size in `params` (as in
#'   `params = list(n = 250)`). Non-numeric columns are set aside; a frame of
#'   two or more character or factor columns and fewer than two numeric ones
#'   is read as sequences and turned into per-row state counts first. A
#'   column named `rid` is set aside; columns with non-syntactic names,
#'   all-missing and constant columns, and rows with a missing value are
#'   dropped, each with a message; fewer than 2 variables or 3 complete rows
#'   raise an error.
#' @param actor The id column naming whose observation each row is (a
#'   person). With `level`, the decomposition groups by it; without, it is
#'   set aside, not modelled. Default `NULL`.
#' @param level Multilevel decomposition, with `actor`: `"between"`
#'   estimates the network of person means, `"within"` the network of
#'   deviations from them, and `"both"` returns the two as a
#'   `netobject_ml`. Default `NULL` (no decomposition).
#' @param cor_method `"pearson"` (default), `"spearman"` or `"kendall"`,
#'   passed to [stats::cor()].
#' @param predictability Whether to compute each node's predictability (the
#'   R-squared of the node on its neighbours). It is stored on the network,
#'   and `cograph::splot()` draws it as a ring around each node. Default
#'   `TRUE`.
#' @param metadata_cols Columns to keep as the network's metadata
#'   (covariates), never read as variables. Default `NULL`.
#' @param params Settings for the estimator, as a named list, such as `n`,
#'   the sample size when `data` is a correlation or covariance matrix. The
#'   network stores them, so [bootstrap_network()] re-estimates with the
#'   same settings. Default `list()`.
#' @param ... Further arguments for [build_network()]. `method` cannot be
#'   given (error of class `nestimate_bad_input`).
#' @inheritParams build_network
#' @return A `netobject` (a `netobject_group` with `group`, a `netobject_ml`
#'   with `level = "both"`), as [build_network()] returns it: print it for
#'   an overview, read the edges with [extract_edges()], and draw it with
#'   `cograph::splot()`.
#' @references
#' Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
#' correlation networks. \emph{Psychological Methods}, 23(4), 617--634.
#' \doi{10.1037/met0000167}
#'
#' Saqr, M., Beck, E., & L\enc{ó}{o}pez-Pernas, S. (2024). Psychological
#' networks: A modern approach to analysis of learning and complex learning
#' processes. In M. Saqr & S. L\enc{ó}{o}pez-Pernas (Eds.), \emph{Learning
#' Analytics Methods and Tutorials: A Practical Guide Using R}. Springer.
#' \doi{10.1007/978-3-031-54464-4_19}
#' @seealso [build_network()] for every method; [build_pcor()],
#'   [build_glasso()] and [build_ising()] for the other association
#'   networks; [bootstrap_network()] and [nct()] for inference.
#' @examples
#' build_network(srl_strategies, method = "cor")
#' build_cor(srl_strategies)
#'
#' build_cor(srl_strategies, cor_method = "spearman", threshold = 0.2)
#' @export
build_cor <- function(data, actor = NULL, level = NULL, ...,
                      cor_method = c("pearson", "spearman", "kendall"),
                      predictability = TRUE, metadata_cols = NULL,
                      group = NULL, scaling = NULL, threshold = 0,
                      labels = NULL, params = list()) {
  .forward_to_build_network("build_cor", "cor", ...names())
}

#' Estimate a partial correlation network
#'
#' @description
#' `build_network(method = "pcor")`, or for short `build_pcor()`, estimates
#' a partial correlation network of numeric variables. The edge between two
#' variables is their partial correlation given all the others, read from
#' the inverse of the correlation matrix, without regularisation (Epskamp &
#' Fried, 2018).
#'
#' Unregularised partial correlations are noisy when there are many
#' variables for the number of observations, and they cannot be estimated
#' when the correlation matrix is singular (more variables than
#' observations, or collinear ones): that raises an error, and a
#' near-singular matrix a warning. [build_glasso()] estimates a sparse
#' network in those cases.
#'
#' @inheritParams build_cor
#' @return A `netobject` (a `netobject_group` with `group`, a `netobject_ml`
#'   with `level = "both"`), as [build_network()] returns it.
#' @references
#' Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
#' correlation networks. \emph{Psychological Methods}, 23(4), 617--634.
#' \doi{10.1037/met0000167}
#' @seealso [build_network()] for every method; [build_cor()],
#'   [build_glasso()] and [build_ising()] for the other association
#'   networks.
#' @examples
#' build_network(srl_strategies, method = "pcor")
#' build_pcor(srl_strategies)
#' @export
build_pcor <- function(data, actor = NULL, level = NULL, ...,
                       cor_method = c("pearson", "spearman", "kendall"),
                       predictability = TRUE, metadata_cols = NULL,
                       group = NULL, scaling = NULL, threshold = 0,
                       labels = NULL, params = list()) {
  .forward_to_build_network("build_pcor", "pcor", ...names())
}

#' Estimate an EBIC graphical lasso network
#'
#' @description
#' `build_network(method = "glasso")`, or for short `build_glasso()`,
#' estimates an EBIC graphical lasso network, a Gaussian graphical model,
#' from numeric variables. The edge between two variables is their partial
#' correlation given all the others, shrunk by the graphical lasso, with the
#' penalty chosen by EBIC (Friedman et al., 2008; Foygel & Drton, 2010).
#'
#' The penalty is chosen along a path of `nlambda` values; an L1 penalty
#' sets small partial correlations to exactly zero. An absent edge means the
#' penalty removed it at the selected lambda, not that the two variables are
#' shown to be independent. Check the stability of the edges and
#' centralities with [boot_glasso()].
#'
#' @param gamma The EBIC hyperparameter: larger values give sparser
#'   networks. Default `0.5`.
#' @param nlambda The number of penalties tried. Default `100`.
#' @param params Settings for the estimator, as a named list:
#'   `lambda.min.ratio` (default `0.01`), `penalize.diagonal` (default
#'   `FALSE`), `refit` (default `FALSE`; `TRUE` re-estimates the selected
#'   edges without the penalty), and `n`, the sample size when `data` is a
#'   correlation or covariance matrix. The network stores them, so the
#'   resampling verbs re-estimate with the same settings. Default `list()`.
#' @inheritParams build_cor
#' @return A `netobject` (a `netobject_group` with `group`, a `netobject_ml`
#'   with `level = "both"`), as [build_network()] returns it.
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
#' covariance estimation with the graphical lasso. \emph{Biostatistics},
#' 9(3), 432--441. \doi{10.1093/biostatistics/kxm045}
#'
#' Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
#' for Gaussian graphical models. \emph{Advances in Neural Information
#' Processing Systems}, 23, 604--612.
#' @seealso [build_network()] for every method; [build_pcor()] for the
#'   unregularised network; [boot_glasso()] and [nct()] for inference.
#' @examples
#' build_network(srl_strategies, method = "glasso")
#' build_glasso(srl_strategies)
#'
#' build_glasso(srl_strategies, gamma = 0.75)
#' @export
build_glasso <- function(data, actor = NULL, level = NULL, ..., gamma = 0.5,
                         nlambda = 100L,
                         cor_method = c("pearson", "spearman", "kendall"),
                         predictability = TRUE, metadata_cols = NULL,
                         group = NULL, scaling = NULL, threshold = 0,
                         labels = NULL, params = list()) {
  .forward_to_build_network("build_glasso", "glasso", ...names())
}

#' Estimate an Ising model
#'
#' @description
#' `build_network(method = "ising")`, or for short `build_ising()`,
#' estimates an Ising model of binary variables. The edge between two
#' variables is their interaction, fitted by L1-penalised logistic
#' regression of each variable on the others, with the penalty chosen by
#' EBIC, and the two estimates of each edge combined into one undirected
#' weight (van Borkulo et al., 2014). It needs the glmnet package.
#'
#' An edge is a conditional association between two variables given all the
#' others; an absent edge means the penalty removed it.
#'
#' @param data The binary variables: a data frame with one row per
#'   observation and one column per variable, holding only 0 and 1.
#'   Non-numeric columns and a column named `rid` are set aside; columns
#'   with non-syntactic names, all-missing and constant columns, and rows
#'   with a missing value are dropped, each with a message; fewer than 2
#'   variables or 3 complete rows, or any value other than 0 and 1, raise an
#'   error.
#' @param actor An id column naming whose observation each row is; it is
#'   set aside, not modelled. Default `NULL`.
#' @param gamma The EBIC hyperparameter. Default `0.25`.
#' @param nlambda The number of penalties tried. Default `100`.
#' @param rule `"AND"` (default) keeps an edge when both regressions select
#'   it, `"OR"` when either does.
#' @param params Settings for the estimator, as a named list. The network
#'   stores them, so the resampling verbs re-estimate with the same settings.
#'   Default `list()`.
#' @inheritParams build_cor
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
#' L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
#' constructing networks from binary data. \emph{Scientific Reports}, 4,
#' 5918. \doi{10.1038/srep05918}
#' @seealso [build_network()] for every method; [build_mgm()] for mixed
#'   continuous and categorical variables.
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1)
#'   binary <- data.frame(matrix(rbinom(200, 1, 0.5), ncol = 5))
#'   build_network(binary, method = "ising")
#'   build_ising(binary, rule = "OR")
#' }
#' @export
build_ising <- function(data, actor = NULL, ..., gamma = 0.25, nlambda = 100L,
                        rule = c("AND", "OR"), metadata_cols = NULL,
                        group = NULL, scaling = NULL, threshold = 0,
                        labels = NULL, params = list()) {
  .forward_to_build_network("build_ising", "ising", ...names())
}

#' Estimate a mixed graphical model (MGM)
#'
#' @description
#' `build_network(method = "mgm")`, or for short `build_mgm()`, estimates a
#' mixed graphical model (MGM) of continuous and categorical variables. The
#' edge between two variables is the strength of their interaction, fitted
#' by penalised node-wise regression with the penalty chosen by EBIC
#' (Haslbeck & Waldorp, 2020).
#'
#' Factor and character columns are categorical, and an integer column with
#' at most 10 values is categorical with a warning. Missing continuous
#' values are filled with the observed column mean and missing categorical
#' values with the modal level; this keeps the rows and is \pkg{psychnets}'
#' documented `na_method = "pairwise"` policy, not pairwise deletion.
#'
#' @param data The variables: a data frame with one row per observation and
#'   one column per variable, continuous or categorical.
#' @param level With `actor`, a multilevel decomposition (`"between"`,
#'   `"within"` or `"both"`), as for [build_cor()]. A numeric vector instead
#'   gives one category count per modelled variable (`1` for continuous
#'   ones); it is checked and kept for compatibility, and the fit reads the
#'   categories from the data. Default `NULL`.
#' @param lambdaGam The EBIC hyperparameter. Default `0.25`.
#' @param ruleReg `"AND"` (default) keeps an edge when both regressions
#'   select it, `"OR"` when either does.
#' @param threshold Weights whose absolute value is below this are set to
#'   zero (default `0`); or `"LW"` or `"none"`, mgm's rule for thresholding
#'   the coefficients (mgm's default is `"LW"`).
#' @param params Settings for the estimator, as a named list. The network
#'   stores them. Default `list()`.
#' @inheritParams build_cor
#' @return A `netobject` (a `netobject_group` with `group`), as
#'   [build_network()] returns it.
#' @references
#' Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating time-varying
#' mixed graphical models in high-dimensional data. \emph{Journal of
#' Statistical Software}, 93(8), 1--46. \doi{10.18637/jss.v093.i08}
#' @seealso [build_network()] for every method; [build_ising()] for binary
#'   variables, [build_glasso()] for continuous ones.
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1)
#'   mixed <- data.frame(x = rnorm(100), y = rnorm(100),
#'                       g = factor(sample(c("a", "b"), 100, TRUE)))
#'   build_network(mixed, method = "mgm")
#'   build_mgm(mixed)
#' }
#' @export
build_mgm <- function(data, actor = NULL, level = NULL, lambdaGam = 0.25,
                      ruleReg = c("AND", "OR"), metadata_cols = NULL,
                      group = NULL, scaling = NULL, threshold = 0,
                      labels = NULL, params = list(), ...) {
  .forward_to_build_network("build_mgm", "mgm", ...names())
}
