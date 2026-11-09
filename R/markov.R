# ---- Markov Chain Analysis: Passage Times and Stability ----

# Columns of the passage-time tile table used in ggplot aesthetics
utils::globalVariables(c("fill", "label", "txt_col"))

#' @noRd
.mpt_not_a_chain <- function(what) {
  stop(errorCondition(sprintf(paste0(
    "%s is not a transition matrix. Markov-chain measures need ",
    "non-negative transition weights (method \"relative\" or ",
    "\"frequency\")."), what), class = "nestimate_not_markov", call = NULL))
}

#' @noRd
.mpt_extract_P <- function(x) {
  if (is.matrix(x) && is.numeric(x)) {
    if (any(x < 0, na.rm = TRUE)) .mpt_not_a_chain("A matrix with negative entries")
    return(x)
  }
  if (inherits(x, "netobject") || inherits(x, "cograph_network") ||
      inherits(x, "tna")) {
    P <- x$weights
    if (is.null(P) || !is.matrix(P) || !is.numeric(P))
      stop(errorCondition("Object has no numeric weight matrix.",
                          class = "nestimate_bad_input", call = NULL))
    # Correlation-type networks (signed, or positive by chance) are not
    # transition weights; row-normalising them makes a meaningless "chain".
    # A co-occurrence network row-normalises to its random walk, a valid one.
    assoc <- c("cor", "pcor", "glasso", "ising", "mgm")
    if (isTRUE(x$method %in% assoc) || any(P < 0, na.rm = TRUE)) {
      .mpt_not_a_chain(sprintf("A \"%s\" network",
                               x$method %||% "undirected"))
    }
    return(P)
  }
  if (is.matrix(x) && (is.character(x) || is.logical(x))) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (is.data.frame(x)) {
    .refuse_undeclared_roles(x, "Reading a chain", via_network = TRUE,
                           reads_onehot = TRUE)
    net <- build_network(x, method = "relative")
    return(net$weights)
  }
  stop(errorCondition(
    paste0("'x' must be a numeric matrix (transition matrix), a netobject / ",
    "cograph_network / tna object, or a wide sequence data.frame / ",
    "character matrix."),
    class = "nestimate_bad_input", call = NULL)
  )
}

#' @noRd
.mpt_normalize_rows <- function(P, state_names, normalize = TRUE) {
  row_sums <- rowSums(P)

  zero_rows <- which(row_sums <= 0)
  if (length(zero_rows)) {
    bad <- state_names[zero_rows]
    stop(errorCondition(
      paste0("Transition matrix has zero-sum row(s) for state(s): ",
      paste(bad, collapse = ", "),
      ". These states have no outgoing transitions, so their rows hold no ",
      "transition probabilities. Remove the state(s), or make an end state ",
      "absorbing with mark_terminal_state() instead of `end = TRUE`."),
      class = "nestimate_bad_input", call = NULL)
    )
  }

  if (any(abs(row_sums - 1) > 1e-6)) {
    if (!normalize)
      stop(errorCondition("Transition matrix rows must sum to 1.",
                          class = "nestimate_bad_input", call = NULL))
    warning(warningCondition("Rows do not sum to 1; normalizing.",
                             class = "nestimate_rows_normalized",
                             call = NULL))
    P <- P / row_sums
  }

  P
}

#' @noRd
.mpt_stationary <- function(P) {
  ev  <- eigen(t(P))
  idx <- which.min(abs(Re(ev$values) - 1))
  pi  <- abs(Re(ev$vectors[, idx]))
  pi / sum(pi)
}

#' @noRd
.mpt_full <- function(P, pi) {
  n     <- nrow(P)
  PI    <- matrix(pi, nrow = n, ncol = n, byrow = TRUE)
  Z     <- solve(diag(n) - P + PI)
  # Kemeny-Snell: M[i,j] = (Z[j,j] - Z[i,j]) / pi[j]
  Zd    <- matrix(diag(Z), nrow = n, ncol = n, byrow = TRUE)
  pimat <- matrix(pi,      nrow = n, ncol = n, byrow = TRUE)
  M     <- (Zd - Z) / pimat
  diag(M) <- 1 / pi
  M
}

# Mean first passage times of a chain that is not irreducible. m[i, j] =
# E[T_j | X_0 = i] is finite only when j is reached from i with probability
# 1; it is then the minimal non-negative solution of
# m_i = 1 + sum_{k != j} P[i, k] m_k over those states (Norris 1997,
# Thm 1.3.5), and Inf otherwise. The diagonal is the mean return time
# 1 + sum_{k != j} P[j, k] m[k, j] (Inf when a step leaves j for a state
# that may never come back). On an irreducible chain this equals the
# Kemeny-Snell matrix of .mpt_full().
#' @noRd
.mpt_hitting_times <- function(P) {
  n <- nrow(P)
  A <- (P > 0) * 1L
  cols <- lapply(seq_len(n), function(j) {
    # Support graph in which the walk stops at the target j: only paths
    # that have not yet hit j decide whether j is certain.
    A_stop <- A
    A_stop[j, ] <- 0L
    R <- .cs_reachability(A_stop)
    diag(R) <- 1L
    # j is hit surely from i iff every state i can reach before j can reach j.
    sure <- vapply(seq_len(n), function(i) all(R[R[i, ] > 0, j] > 0),
                   logical(1))
    m <- rep(Inf, n)
    S <- setdiff(which(sure), j)
    if (length(S) > 0L) {
      m[S] <- drop(solve(diag(length(S)) - P[S, S, drop = FALSE],
                         rep(1, length(S))))
    }
    leave <- setdiff(which(P[j, ] > 0), j)
    m[j] <- if (any(is.infinite(m[leave]))) Inf else
      1 + sum(P[j, leave] * m[leave])
    m
  })
  M <- do.call(cbind, cols)
  dimnames(M) <- dimnames(P)
  M
}

# ---- passage_time() ----

#' Mean First Passage Times
#'
#' @description
#' `passage_time()` computes the mean first passage times of a Markov chain:
#' the entry in row \eqn{i} and column \eqn{j} is the expected number of
#' steps a chain started in state \eqn{i} takes to reach state \eqn{j} for
#' the first time, and the diagonal entry is the mean recurrence time,
#' the expected number of steps to return to state \eqn{i}. On an
#' irreducible chain it equals \eqn{1/\pi_i}.
#' The stationary distribution \eqn{\pi} is returned with them. It reads a
#' transition network, a `tna` model, a transition matrix or wide sequences,
#' and analyses a grouped network one network at a time.
#'
#' The times describe the estimated chain, in which a transition never seen
#' in the data counts as impossible. In an irreducible chain every state can
#' be reached from every other and every time is finite; in any other chain
#' a time into a state that is not reached with certainty is `Inf` (see
#' *Details*), so check the structure with [chain_structure()] first. [markov_stability()] condenses the matrix into
#' per-state summaries beside each state's persistence and sojourn time.
#'
#' @param x The chain: a `netobject`, bare `cograph_network` or `tna` model
#'   with non-negative transition weights, typically method `"relative"` or
#'   `"frequency"` (a co-occurrence network is read as its random walk), a
#'   non-negative numeric transition matrix (states are named `S1`, `S2`,
#'   ... when it has no names), or a wide sequence data frame or character
#'   matrix, from which a `"relative"` network is estimated. A data frame
#'   needs its position columns named as a numbered run, as in
#'   [build_network()]; build the network first for any other layout. A
#'   `netobject_group` gives one result per group. An association network
#'   (`cor`, `pcor`, `glasso`, `ising`, `mgm`) or negative weights raise an
#'   error of class `nestimate_not_markov`. For the `print()` methods: a
#'   `net_mpt`, a `net_mpt_group`, or the result of `summary()`; for the
#'   `plot()` methods: a `net_mpt` or `net_mpt_group`.
#' @param object For the `summary()` method: a `net_mpt`.
#' @param states The states to keep in the result. The times are computed
#'   on the whole chain and then subset, so they still count paths through
#'   the states left out. A name that is not a state raises an error.
#'   Default `NULL` (all states).
#' @param normalize Whether to rescale rows that do not sum to 1 so that
#'   they do, with a warning, as for a `"frequency"` network. With `FALSE`,
#'   such rows raise an error. A state with no outgoing transitions (a zero
#'   row) raises an error either way. Default `TRUE`.
#' @param ... For `print()` of a `net_mpt_group`: passed to `print()` of
#'   each group's result (for example `digits`). Otherwise not used.
#' @return A `net_mpt` object; for a `netobject_group`, a `net_mpt_group`,
#'   a named list holding one `net_mpt` per group. `print()` shows the
#'   passage-time matrix and the stationary distribution. `summary()`
#'   returns a `summary.net_mpt` object, a list holding the per-state table
#'   (`table`) and the summarised `net_mpt` (`object`); its `print()` shows
#'   the table, one row per state: `state`, `return_time` (mean recurrence time),
#'   `stationary` (\eqn{\pi_i}), `mean_out` (the mean number of steps from
#'   the state to the others) and
#'   `mean_in` (the mean number of steps from the others to the state),
#'   rounded to 2 decimals (`stationary` to 4). `plot()` returns a `ggplot`
#'   heatmap of the matrix, rows the starting state and columns the target;
#'   for a `net_mpt_group`, one panel per group on a shared colour scale.
#'   The `print()` methods return their input invisibly. A `net_mpt` holds:
#' \describe{
#'   \item{matrix}{The \eqn{n \times n} matrix of mean first passage times:
#'     row \eqn{i}, column \eqn{j} is the expected number of steps from
#'     state \eqn{i} to state \eqn{j}; the diagonal is the mean recurrence
#'     time, equal to \eqn{1/\pi_i} on an irreducible chain.}
#'   \item{stationary}{The stationary distribution \eqn{\pi}, a named
#'     numeric vector; zero outside the closed class of a chain with one
#'     closed class, and `NA` for a chain with several, whose stationary
#'     distribution is not unique.}
#'   \item{return_times}{The mean recurrence time of each state (the
#'     diagonal of `matrix`), a named numeric vector: \eqn{1/\pi_i} in an
#'     irreducible chain, `Inf` for a transient state.}
#'   \item{states}{The state names.}
#' }
#'
#' @details
#' The times follow the Kemeny-Snell fundamental matrix formula:
#' \deqn{M_{ij} = \frac{Z_{jj} - Z_{ij}}{\pi_j}, \quad
#'       Z = (I - P + \Pi)^{-1}}
#' where \eqn{\Pi_{ij} = \pi_j}, and \eqn{\pi} is the left eigenvector
#' of \eqn{P} for eigenvalue 1, scaled to sum to 1. The fundamental matrix
#' exists for every irreducible chain (an ergodic chain in Kemeny and
#' Snell's terms), periodic ones included.
#'
#' A chain that is not irreducible (for example one with an absorbing state,
#' or with several closed classes) has no fundamental matrix. Its times are
#' the expected hitting times: the time from \eqn{i} to \eqn{j} is `Inf`
#' when the chain started in \eqn{i} may never reach \eqn{j}, and
#' otherwise the minimal non-negative solution of
#' \eqn{m_i = 1 + \sum_{k \ne j} P_{ik} m_k} (Norris, 1997). Times into
#' the states outside a closed class are therefore `Inf` from inside it, and
#' times into the closed class are finite from every state that cannot
#' avoid it. The diagonal is the mean return time, `Inf` for a transient
#' state. On an irreducible chain both routes give the same matrix.
#' Nestimate 1.0.9 and earlier raised the zero stationary probabilities to
#' machine precision instead, giving times near \eqn{10^{15}} with a
#' warning, and stopped with an error for a chain with several closed
#' classes.
#'
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' pt  <- passage_time(net)
#' print(pt)
#' \donttest{
#' plot(pt)
#' }
#'
#' @seealso \code{\link{markov_stability}}, \code{\link{build_network}}
#'
#' @references
#' Kemeny, J.G. and Snell, J.L. (1976). \emph{Finite Markov Chains}.
#' Springer-Verlag.
#'
#' Norris, J. R. (1997). \emph{Markov Chains}. Cambridge University Press.
#'
#' @export
passage_time <- function(x, states = NULL, normalize = TRUE) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      passage_time(net, states = states, normalize = normalize)
    })
    class(out) <- c("net_mpt_group", "list")
    return(out)
  }
  P           <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names     <- paste0("S", seq_len(nrow(P)))
    colnames(P)     <- rownames(P) <- state_names
  }

  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)

  classes <- .cs_scc((P > 0) * 1L)
  if (length(classes) == 1L) {
    # Irreducible: Kemeny-Snell fundamental matrix, pi > 0 everywhere.
    pi <- .mpt_stationary(P)
    M  <- .mpt_full(P, pi)
  } else {
    # Reducible: a time into a state that is not reached with certainty is
    # Inf (formerly pi was floored at machine epsilon, giving ~1e15). pi is
    # unique only with one closed class, where it is zero off that class.
    M <- .mpt_hitting_times(P)
    closed <- Filter(function(cl) {
      all(P[cl, setdiff(seq_len(nrow(P)), cl), drop = FALSE] == 0)
    }, classes)
    pi <- rep(if (length(closed) == 1L) 0 else NA_real_, nrow(P))
    if (length(closed) == 1L) {
      cl <- closed[[1L]]
      pi[cl] <- if (length(cl) == 1L) 1 else
        .mpt_stationary(P[cl, cl, drop = FALSE])
    }
  }
  names(pi)   <- state_names
  dimnames(M) <- list(state_names, state_names)
  return_times <- diag(M)

  if (!is.null(states)) {
    bad <- setdiff(states, state_names)
    if (length(bad))
      stop(errorCondition(
        paste0("Unknown states: ", paste(bad, collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    M  <- M[states, states, drop = FALSE]
    pi <- pi[states]
    return_times <- return_times[states]
  }

  structure(
    list(
      matrix       = M,
      stationary   = pi,
      return_times = return_times,
      states       = rownames(M)
    ),
    class = "net_mpt"
  )
}

#' @rdname passage_time
#' @param log_scale For the `plot()` methods: whether to colour the tiles by
#'   the logarithm of the passage time, so that short and long times can
#'   both be told apart; the legend and the tile labels still show steps.
#'   Default `TRUE`.
#' @param digits Decimal places of the printed matrix (`print()`) and of
#'   the tile labels (`plot()`). The stationary distribution is always
#'   printed with 4. Default `1`.
#' @param title For the `plot()` methods: the plot title. Default
#'   `"Mean First Passage Times"`.
#' @param low,high For the `plot()` methods: the fill colours for the
#'   shortest and the longest passage times. Default `"#004d00"` (dark
#'   green) and `"#ccffcc"` (pale green).
#' @export
print.net_mpt <- function(x, digits = 1, ...) {
  n <- length(x$states)
  cat(sprintf("Mean First Passage Times (%d states)\n\n", n))
  print(round(x$matrix, digits))
  cat("\nStationary distribution:\n")
  print(round(x$stationary, 4))
  invisible(x)
}

#' @rdname passage_time
#' @export
print.net_mpt_group <- function(x, ...) {
  cat(sprintf("Mean First Passage Times -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  invisible(lapply(names(x), function(nm) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }))
  invisible(x)
}

#' @rdname passage_time
#' @export
summary.net_mpt <- function(object, ...) {
  M  <- object$matrix
  n  <- nrow(M)
  st <- object$states

  df <- data.frame(
    state       = st,
    return_time = round(object$return_times, 2),
    stationary  = round(object$stationary,   4),
    mean_out    = vapply(seq_len(n),
                         function(i) round(mean(M[i, -i]), 2), numeric(1)),
    mean_in     = vapply(seq_len(n),
                         function(i) round(mean(M[-i, i]), 2), numeric(1)),
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  structure(list(table = df, object = object), class = "summary.net_mpt")
}

#' @rdname passage_time
#' @export
print.summary.net_mpt <- function(x, ...) {
  cat("Mean First Passage Times - Summary\n\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' @rdname passage_time
#' @export
plot.net_mpt <- function(x,
                          log_scale = TRUE,
                          digits    = 1,
                          title     = "Mean First Passage Times",
                          low       = "#004d00",
                          high      = "#ccffcc",
                          ...) {
  .mpt_check_plot_args(log_scale, digits)
  .mpt_heatmap(.mpt_tiles(x, log_scale, digits), log_scale, title, low, high)
}

#' @rdname passage_time
#' @export
plot.net_mpt_group <- function(x,
                                log_scale = TRUE,
                                digits    = 1,
                                title     = "Mean First Passage Times",
                                low       = "#004d00",
                                high      = "#ccffcc",
                                ...) {
  .mpt_check_plot_args(log_scale, digits)
  tiles <- do.call(rbind, Map(\(g, m) {
    cbind(group = g, .mpt_tiles(m, log_scale, digits))
  }, names(x), x))
  tiles$group <- factor(tiles$group, levels = names(x))
  .mpt_heatmap(tiles, log_scale, title, low, high) +
    ggplot2::facet_wrap(~ group)
}

#' Arguments of the passage-time plots, checked before any drawing
#' @noRd
.mpt_check_plot_args <- function(log_scale, digits) {
  .assert(
    "`log_scale` must be TRUE or FALSE" =
      (is.logical(log_scale) || is.numeric(log_scale)) &&
      length(log_scale) == 1L && !is.na(log_scale),
    "`digits` must be numeric" = is.numeric(digits)
  )
}

#' Long tile table of one passage-time matrix
#' @noRd
.mpt_tiles <- function(x, log_scale, digits) {
  states <- x$states
  n <- length(states)
  # row-major, to match from = rep(each), to = rep(times): row = from
  vals <- as.vector(t(x$matrix))
  data.frame(
    from  = factor(rep(states, each = n), levels = rev(states)),
    to    = factor(rep(states, times = n), levels = states),
    fill  = if (log_scale) log(pmax(vals, .Machine$double.eps)) else vals,
    label = formatC(round(vals, digits), format = "f", digits = digits,
                    flag = " "),
    stringsAsFactors = FALSE
  )
}

#' Passage-time heatmap from a tile table (text colour from the shared range)
#' @noRd
.mpt_heatmap <- function(df, log_scale, title, low, high) {
  mid <- mean(range(df$fill, na.rm = TRUE))
  df$txt_col <- ifelse(df$fill <= mid, "white", "#003300")
  ggplot2::ggplot(df, ggplot2::aes(x = to, y = from, fill = fill)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label, color = txt_col),
                       size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(
      low  = low, high = high,
      name = if (log_scale) "log(Steps)" else "Steps",
      labels = if (log_scale) function(b) round(exp(b), 1) else ggplot2::waiver()
    ) +
    ggplot2::labs(x = "To", y = "From", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


# ---- markov_stability() ----

#' Markov Stability Analysis
#'
#' @description
#' `markov_stability()` describes how stable each state of a Markov chain
#' is and returns a `net_markov_stability` object holding one row per state:
#' its persistence (the probability of staying, \eqn{P_{ii}}), stationary
#' probability (\eqn{\pi_i}, the long-run share of steps spent in it), mean
#' recurrence time, sojourn time (the expected number of
#' consecutive steps in it, \eqn{1/(1 - P_{ii})}), and the mean of its first
#' passage times to the other states and of theirs to it, taken from
#' [passage_time()]. It reads the same inputs as [passage_time()] and
#' analyses a grouped network one network at a time.
#'
#' Persistence and sojourn time describe staying in a state; the passage
#' times describe how long it takes to get from it to the other states, or
#' from them to it. All of them are computed from the estimated transition
#' matrix. On an irreducible chain, the mean recurrence time is
#' \eqn{1/\pi_i}. On a reducible chain, recurrence times are computed within
#' each closed class and are `Inf` for transient states; passage times to
#' states not reached with certainty are `Inf`. The stationary probabilities
#' are `NA` when there is more than one closed class, because the stationary
#' distribution is not unique. [chain_structure()] identifies the classes;
#' how far the underlying transition probabilities move under resampling is
#' shown by [bootstrap_network()].
#'
#' @param x The chain, as in [passage_time()]: a `netobject`, bare
#'   `cograph_network` or `tna` model with non-negative transition weights,
#'   a non-negative numeric transition matrix, or a wide sequence data frame
#'   or character matrix (a `"relative"` network is estimated from it). A
#'   `netobject_group` gives one result per group. An association network or
#'   negative weights raise an error of class `nestimate_not_markov`. For
#'   the `print()` and `plot()` methods: a `net_markov_stability` or
#'   `net_markov_stability_group`.
#' @param normalize Whether to rescale rows that do not sum to 1 so that
#'   they do, with a warning, as for a `"frequency"` network. With `FALSE`,
#'   such rows raise an error. A state with no outgoing transitions (a zero
#'   row) raises an error either way. Default `TRUE`.
#' @param ... Not used (`print()` of a group passes it on to each group's
#'   `print()`, which ignores it).
#' @return A `net_markov_stability` object; for a `netobject_group`, a
#'   `net_markov_stability_group`, a named list holding one such object per
#'   group. `print()` shows the per-state table and returns its input
#'   invisibly. `summary()` prints the most accessible state (the lowest
#'   mean passage time from the other states) and the most persistent one
#'   (the longest sojourn time), and returns the per-state table as a data
#'   frame. `plot()` returns a `ggplot` with one panel of bars per metric and
#'   the states ordered by stationary probability; with `combined = FALSE`,
#'   invisibly, a list of one-panel plots named by metric. For a group,
#'   `plot()` draws one panel per metric with one row per state and one
#'   point per group (told apart by colour and shape) joined by a line, the
#'   states ordered by their mean stationary probability; more than 81
#'   groups (9 colours times 9 shapes) raise an error of class
#'   `nestimate_bad_input`. A `net_markov_stability` holds:
#' \describe{
#'   \item{stability}{The per-state table, one row per state, with the
#'     columns `state`, `persistence` (\eqn{P_{ii}}), `stationary_prob`
#'     (\eqn{\pi_i}), `return_time` (the mean recurrence time), `sojourn_time`
#'     (\eqn{1/(1-P_{ii})}), `avg_time_to_others` (the mean first passage
#'     time from the state to the others) and `avg_time_from_others` (the
#'     mean first passage time from the others to the state).
#'     `persistence` and `stationary_prob` are rounded to 4 decimals, the
#'     times to 2.}
#'   \item{mpt}{The `net_mpt` from [passage_time()] the times come from.}
#' }
#'
#' @details
#' \strong{Sojourn time} is the expected number of consecutive steps spent
#' in a state before leaving it, \eqn{1/(1-P_{ii})}, the mean of a geometric
#' distribution. A state with `persistence = 1` never leaves and has
#' `sojourn_time = Inf`.
#'
#' \strong{avg_time_to_others} is the mean passage time from the state to
#' each other state: a high value means the chain takes long to get
#' anywhere else from it.
#'
#' \strong{avg_time_from_others} is the mean passage time from each other
#' state to this one: a low value means the state is reached quickly from
#' elsewhere.
#'
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' ms  <- markov_stability(net)
#' print(ms)
#' \donttest{
#' plot(ms)
#' }
#'
#' @seealso \code{\link{passage_time}}
#'
#' @references
#' Kemeny, J.G. and Snell, J.L. (1976). \emph{Finite Markov Chains}.
#' Springer-Verlag.
#'
#' @export
markov_stability <- function(x, normalize = TRUE) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      markov_stability(net, normalize = normalize)
    })
    class(out) <- c("net_markov_stability_group", "list")
    return(out)
  }
  P           <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names     <- paste0("S", seq_len(nrow(P)))
    colnames(P)     <- rownames(P) <- state_names
  }

  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)

  mpt  <- passage_time(P, normalize = FALSE)
  M    <- mpt$matrix
  pi   <- mpt$stationary
  n    <- nrow(P)

  persistence <- diag(P)
  sojourn     <- ifelse(persistence >= 1 - .Machine$double.eps,
                        Inf, 1 / (1 - persistence))
  avg_to   <- vapply(seq_len(n), function(i) mean(M[i, -i]), numeric(1))
  avg_from <- vapply(seq_len(n), function(i) mean(M[-i, i]), numeric(1))

  stability_df <- data.frame(
    state                = state_names,
    persistence          = round(persistence,  4),
    stationary_prob      = round(pi,            4),
    return_time          = round(mpt$return_times, 2),
    sojourn_time         = round(sojourn,       2),
    avg_time_to_others   = round(avg_to,        2),
    avg_time_from_others = round(avg_from,      2),
    stringsAsFactors     = FALSE
  )
  rownames(stability_df) <- NULL

  structure(list(stability = stability_df, mpt = mpt),
            class = "net_markov_stability")
}

#' @rdname markov_stability
#' @param metrics For the `plot()` methods: the metrics to draw, one panel
#'   each, any of `"persistence"`, `"stationary_prob"`, `"return_time"`,
#'   `"sojourn_time"`, `"avg_time_to_others"` and `"avg_time_from_others"`.
#'   Default all six.
#' @param combined For `plot()` of a single network: `TRUE` (default) draws
#'   the metrics as panels of one plot; `FALSE` returns a named list of
#'   one-panel plots, one per metric, to print or save separately. Ignored
#'   for a group, whose `plot()` always draws one plot, with a warning of
#'   class `nestimate_ignored_argument`.
#' @export
print.net_markov_stability <- function(x, ...) {
  cat("Markov Stability Analysis\n\n")
  print(x$stability, row.names = FALSE)
  invisible(x)
}

#' @rdname markov_stability
#' @export
print.net_markov_stability_group <- function(x, ...) {
  cat(sprintf("Markov Stability -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  invisible(lapply(names(x), function(nm) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }))
  invisible(x)
}

#' @rdname markov_stability
#' @param object For the `summary()` method: a `net_markov_stability`.
#' @export
summary.net_markov_stability <- function(object, ...) {
  df      <- object$stability
  attract <- df$state[which.min(df$avg_time_from_others)]
  sticky  <- df$state[which.max(df$sojourn_time)]
  cat(sprintf("Most accessible state (attractor): %s\n", attract))
  cat(sprintf("Most persistent state (stickiest): %s\n\n", sticky))
  df
}

#' @rdname markov_stability
#' @export
plot.net_markov_stability <- function(x,
                                       metrics = c("persistence",
                                                   "stationary_prob",
                                                   "return_time",
                                                   "sojourn_time",
                                                   "avg_time_to_others",
                                                   "avg_time_from_others"),
                                       combined = TRUE,
                                       ...) {
  df <- x$stability
  metrics <- match.arg(metrics, several.ok = TRUE)
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))

  labels <- c(
    persistence          = "Persistence",
    stationary_prob      = "Stationary Prob.",
    return_time          = "Return Time",
    sojourn_time         = "Sojourn Time",
    avg_time_to_others   = "Avg. Steps to Others",
    avg_time_from_others = "Avg. Steps from Others"
  )

  st_ord <- factor(df$state, levels = df$state[order(df$stationary_prob)])


  plot_rows <- lapply(metrics, function(m) {
    data.frame(state  = st_ord,
               metric = labels[m],
               value  = df[[m]],
               stringsAsFactors = FALSE,
               row.names = NULL)
  })
  plot_df <- do.call(rbind, plot_rows)
  plot_df$metric <- factor(plot_df$metric, levels = labels[metrics])

  base_plot <- function(d, ttl) {
    # One fill: the state is labelled on the axis, and a per-state palette
    # runs out of colours beyond eight states
    ggplot2::ggplot(d, ggplot2::aes(x = state, y = value)) +
      ggplot2::geom_col(fill = "#0072B2") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL, title = ttl) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(face = "bold")
      )
  }
  if (!combined) {
    levs <- levels(plot_df$metric)
    plots <- lapply(levs, function(lv) {
      sub <- plot_df[plot_df$metric == lv, , drop = FALSE]
      base_plot(sub, lv)
    })
    names(plots) <- as.character(metrics)
    return(invisible(plots))
  }
  base_plot(plot_df, "Markov Stability Metrics") +
    ggplot2::facet_wrap(~ metric, scales = "free_x", ncol = 2)
}

#' @rdname markov_stability
#' @export
plot.net_markov_stability_group <- function(x,
                                             metrics = c("persistence",
                                                         "stationary_prob",
                                                         "return_time",
                                                         "sojourn_time",
                                                         "avg_time_to_others",
                                                         "avg_time_from_others"),
                                             ...) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  if ("combined" %in% names(list(...))) {
    warning(warningCondition(paste0(
      "`combined` is ignored for a group: its plot() always draws one plot, ",
      "one panel per metric."),
      class = "nestimate_ignored_argument", call = NULL))
  }
  labels <- c(
    persistence          = "Persistence",
    stationary_prob      = "Stationary Prob.",
    return_time          = "Return Time",
    sojourn_time         = "Sojourn Time",
    avg_time_to_others   = "Avg. Steps to Others",
    avg_time_from_others = "Avg. Steps from Others"
  )
  d <- do.call(rbind, lapply(names(x), \(g) {
    st <- x[[g]]$stability
    do.call(rbind, lapply(metrics, \(m) {
      data.frame(group = g, state = st$state, metric = labels[[m]],
                 value = st[[m]], stringsAsFactors = FALSE)
    }))
  }))
  # States ordered by their mean stationary probability across groups
  stationary <- do.call(rbind, lapply(x, \(s) s$stability))
  mean_pi <- tapply(stationary$stationary_prob, stationary$state, mean)
  d$state <- factor(d$state, levels = names(sort(mean_pi)))
  d$metric <- factor(d$metric, levels = labels[metrics])
  d$group <- factor(d$group, levels = names(x))
  marks <- .group_marks(nlevels(d$group))
  ggplot2::ggplot(d, ggplot2::aes(x = value, y = state)) +
    ggplot2::geom_line(ggplot2::aes(group = state), colour = "grey75",
                       linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = group, shape = group),
                        size = 2.6) +
    ggplot2::scale_colour_manual(values = marks$colour, name = NULL) +
    ggplot2::scale_shape_manual(values = marks$shape, name = NULL) +
    ggplot2::facet_wrap(~ metric, scales = "free_x", ncol = 2) +
    ggplot2::labs(x = NULL, y = NULL, title = "Markov Stability Metrics") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom"
    )
}


#' Colour and shape per group: a distinct pair for up to 81 groups
#'
#' Okabe-Ito colours cycle every 9 groups and the shape shifts by one at each
#' cycle, so no two groups share both colour and shape and groups never
#' differ by colour alone.
#' @param k Number of groups.
#' @return List with \code{colour} and \code{shape}, each of length \code{k}.
#' @noRd
.group_marks <- function(k) {
  shapes <- c(16, 17, 15, 18, 3, 4, 8, 1, 2)
  if (k > length(.okabe_ito) * length(shapes)) {
    stop(errorCondition(
      sprintf("%d groups cannot be told apart in one plot (at most %d).",
              k, length(.okabe_ito) * length(shapes)),
      class = "nestimate_bad_input", call = NULL))
  }
  j <- seq_len(k) - 1L
  list(colour = .okabe_ito[j %% 9L + 1L],
       shape = shapes[(j + j %/% 9L) %% 9L + 1L])
}
