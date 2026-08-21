# ---- Transition Entropy ----

#' @noRd
.te_row_entropy <- function(P, base) {
  ## 0 * log 0 := 0; vectorised, no loops
  L      <- ifelse(P > 0, log(P, base = base), 0)
  -rowSums(P * L)
}

# ---- transition_entropy() ----

#' Transition Entropy of a Markov Chain
#'
#' @description
#' Computes per-state branching entropy, stationary entropy, and the
#' chain-level entropy rate of a Markov transition process. The entropy rate
#' \eqn{H = -\sum_i \pi_i \sum_j P_{ij} \log_b P_{ij}} (with \eqn{\pi} the
#' stationary distribution from the eigendecomposition of \eqn{P^\top} at
#' \eqn{\lambda = 1}) is the Shannon-McMillan-Breiman per-step uncertainty
#' of trajectories - the canonical information-theoretic summary of a
#' transition matrix, introduced to behavioral research as gaze transition
#' entropy by Krejtz et al. (2015) and tracked in real time as mobile
#' transition matrix entropy by Krejtz et al. (2025). The normalized fields
#' (\code{*_norm}, division by \eqn{\log_b n}) are the scale-free variants
#' those papers report.
#'
#' @param x A \code{netobject}, \code{cograph_network}, \code{tna} object,
#'   row-stochastic numeric transition matrix, or a wide sequence data.frame
#'   (rows = actors, columns = time-steps; a relative transition network is
#'   built automatically). Group dispatch on \code{netobject_group}.
#' @param base Numeric. Logarithm base. \code{2} (default) for bits,
#'   \code{exp(1)} for nats, \code{10} for hartleys.
#' @param normalize Logical. If \code{TRUE} (default), rows that do not sum
#'   to 1 are normalised automatically (with a warning).
#'
#' @return An object of class \code{"net_transition_entropy"} with:
#' \describe{
#'   \item{row_entropy}{Named numeric vector, length \eqn{n}. Per-state
#'     branching entropy \eqn{H(P_{i\cdot}) = -\sum_j P_{ij} \log P_{ij}}.}
#'   \item{row_entropy_norm}{Named numeric vector. \code{row_entropy}
#'     divided by the ceiling \eqn{\log_b n} (in \eqn{[0, 1]}; all zeros
#'     when \eqn{n = 1}).}
#'   \item{stationary}{Named numeric vector. Stationary distribution
#'     \eqn{\pi}.}
#'   \item{stationary_entropy}{Scalar. \eqn{H(\pi) = -\sum_i \pi_i \log \pi_i}
#'     - the entropy of \eqn{\pi} treated as an i.i.d. distribution. Upper
#'     bound on the entropy rate.}
#'   \item{stationary_entropy_norm}{Scalar. \code{stationary_entropy}
#'     divided by the ceiling \eqn{\log_b n}.}
#'   \item{entropy_rate}{Scalar. \eqn{h(P) = \sum_i \pi_i H(P_{i\cdot})} -
#'     the Shannon-McMillan-Breiman entropy rate.}
#'   \item{entropy_rate_norm}{Scalar. \code{entropy_rate} divided by the
#'     ceiling \eqn{\log_b n}.}
#'   \item{redundancy}{Scalar. \eqn{H(\pi) - h(P)}, the entropy deficit
#'     attributable to serial dependence; zero for an i.i.d. chain (rows of
#'     \eqn{P} all equal \eqn{\pi}).}
#'   \item{redundancy_norm}{Scalar. The \emph{relative} redundancy
#'     \eqn{(H(\pi) - h(P)) / H(\pi)} (the fraction of the stationary
#'     entropy removed by order-1 memory), \strong{not}
#'     \code{redundancy} divided by \eqn{\log_b n}; \code{0} when
#'     \eqn{H(\pi) = 0}.}
#'   \item{max_entropy}{Scalar. The normalising ceiling \eqn{\log_b n}.}
#'   \item{base}{Logarithm base used.}
#'   \item{states}{Character vector of state names.}
#' }
#'
#' @details
#' Convention \eqn{0 \log 0 := 0} is applied, so absorbing or
#' deterministic rows contribute zero per-row entropy. The chain need not be
#' irreducible; \eqn{\pi} is computed from the eigendecomposition of
#' \eqn{P^\top} as elsewhere in the package. For non-ergodic chains the
#' returned \eqn{\pi} is one stationary distribution among many - interpret
#' with the help of \code{\link{chain_structure}}.
#'
#' The relation \eqn{h(P) \leq H(\pi)} holds with equality iff successive
#' states are independent. The deficit \eqn{H(\pi) - h(P)} is reported as
#' \code{redundancy} - a measure of how much memory the chain has at order 1.
#'
#' @examples
#' \donttest{
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' te  <- transition_entropy(net)
#' print(te)
#' summary(te)
#' plot(te)
#' }
#'
#' @seealso \code{\link{entropy_network}} for the edge-level decomposition,
#'   \code{\link{entropy_trajectory}} for the sliding-window version,
#'   \code{\link{entropy_bayes}} for credible intervals;
#'   \code{\link{markov_stability}}, \code{\link{passage_time}},
#'   \code{\link{markov_order_test}}, \code{\link{chain_structure}}
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information Theory},
#' 2nd ed., chapter 4. Wiley.
#'
#' Krejtz, K., Duchowski, A., Szmidt, T., Krejtz, I., Gonzalez Perilli, F.,
#' Pires, A., Vilaro, A., & Villalobos, N. (2015). Gaze transition entropy.
#' \emph{ACM Transactions on Applied Perception}, 13(1), 4:1-4:20.
#' \doi{10.1145/2834121}
#'
#' Krejtz, K., Hughes, C.J., Stasiak, I., Duchowski, A., & Krejtz, I.
#' (2025). Real-time mobile transition matrix entropy based on eye and head
#' movements. \emph{Proceedings of ETRA '25}. \doi{10.1145/3715669.3723128}
#'
#' Shannon, C.E. (1948). A mathematical theory of communication.
#' \emph{Bell System Technical Journal}, 27, 379-423.
#'
#' @export
transition_entropy <- function(x, base = 2, normalize = TRUE) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      transition_entropy(net, base = base, normalize = normalize)
    })
    class(out) <- c("net_transition_entropy_group", "list")
    return(out)
  }

  stopifnot(is.numeric(base), length(base) == 1L, base > 0, base != 1)

  P           <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }

  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)

  pi <- .mpt_stationary(P)
  names(pi) <- state_names
  if (any(pi <= 0)) {
    warning("Non-positive stationary probabilities; chain may not be ergodic.",
            call. = FALSE)
    pi <- pmax(pi, .Machine$double.eps)
    pi <- pi / sum(pi)
  }

  H_row <- .te_row_entropy(P, base = base)
  names(H_row) <- state_names

  H_pi  <- {
    p     <- pi[pi > 0]
    -sum(p * log(p, base = base))
  }
  h_P   <- sum(pi * H_row)
  H_max <- log(length(state_names), base = base)

  ## Normalised forms (unit-free, in [0, 1]); H_max == 0 only when n == 1
  norm <- function(x) if (H_max > 0) x / H_max else x * 0
  H_row_norm <- norm(H_row)
  H_pi_norm  <- norm(H_pi)
  h_P_norm   <- norm(h_P)
  redund_norm <- if (H_pi > 0) (H_pi - h_P) / H_pi else 0

  structure(
    list(
      row_entropy             = H_row,
      row_entropy_norm        = H_row_norm,
      stationary              = pi,
      stationary_entropy      = H_pi,
      stationary_entropy_norm = H_pi_norm,
      entropy_rate            = h_P,
      entropy_rate_norm       = h_P_norm,
      redundancy              = H_pi - h_P,
      redundancy_norm         = redund_norm,
      max_entropy             = H_max,
      base                    = base,
      states                  = state_names
    ),
    class = "net_transition_entropy"
  )
}

# ---- entropy_network() ----

#' Transition Entropy Network
#'
#' @description
#' Decomposes the entropy rate of a Markov transition process edge by edge
#' and returns the decomposition as a network. The entropy rate
#' \eqn{H = -\sum_{ij} \pi_i P_{ij} \log P_{ij}}
#' (\code{\link{transition_entropy}}; Krejtz et al. 2015, 2025) is an
#' additive sum over transitions, so every edge \eqn{i \to j} owns the exact
#' term \eqn{\pi_i P_{ij} \log(1/P_{ij})} of the chain-level uncertainty.
#' \strong{No new quantity is estimated}: the network displays the summands
#' of the entropy-rate equation on the transition graph, locating the
#' process's uncertainty spatially. The returned object is a regular
#' \code{netobject} / \code{cograph_network}, so it prints, summarises, and
#' plots (\code{cograph::splot()}) like any other Nestimate network.
#'
#' @param x A \code{netobject}, \code{cograph_network}, \code{tna} object,
#'   row-stochastic numeric transition matrix, or a wide sequence data.frame
#'   (rows = actors, columns = time-steps; a relative transition network is
#'   built automatically).
#' @param base Numeric. Logarithm base. \code{2} (default) for bits,
#'   \code{exp(1)} for nats, \code{10} for hartleys.
#' @param weight Character. Edge weight definition:
#'   \describe{
#'     \item{\code{"contribution"}}{(default) \eqn{\pi_i P_{ij} \log_b(1/P_{ij})}
#'       - the edge's additive share of the entropy rate; all edge weights
#'       sum exactly to \eqn{h(P)}. Thick edges are transitions that are both
#'       frequent \emph{and} unpredictable.}
#'     \item{\code{"surprisal"}}{\eqn{\log_b(1/P_{ij})} - the information
#'       content of observing the transition, ignoring how often it occurs.
#'       Thick edges are rare, surprising transitions. Low surprisal =
#'       expected transition; this is the edge-level predictability reading.}
#'     \item{\code{"production"}}{\eqn{F_{ij} \log_b(F_{ij}/F_{ji})} with
#'       \eqn{F_{ij} = \pi_i P_{ij}} - the edge's contribution to the entropy
#'       production rate (irreversibility). Positive on the dominant
#'       direction of a pair, negative on the reverse; the two always sum to
#'       \eqn{(F_{ij}-F_{ji}) \log_b(F_{ij}/F_{ji}) \ge 0}. Pairs with flow
#'       in only one direction have infinite production and are excluded
#'       (weight 0); their count is reported in
#'       \code{$params$n_oneway_pairs}. Total over included pairs is
#'       \code{$params$production_rate} - \code{0} iff the chain is
#'       reversible (detailed balance).}
#'   }
#' @param scaling Character. \code{"none"} (default) keeps weights in
#'   \code{base}-units. \code{"share"} (contribution only) expresses each
#'   edge as its \emph{percentage} of the entropy rate - weights sum to
#'   100 and each label reads "this transition holds x% of the process's
#'   unpredictability"; comparable across networks regardless of state
#'   count or base. \code{"chance"} (surprisal only) divides by
#'   \eqn{\log_b n} - the surprisal of a chance-level transition when all
#'   \eqn{n} next states are equally likely: \code{1} = chance level,
#'   below 1 = expected transition, above 1 = rarer than chance.
#' @param normalize Logical. If \code{TRUE} (default), rows that do not sum
#'   to 1 are normalised automatically (with a warning).
#'
#' @return A \code{netobject} (also class \code{cograph_network}) whose
#'   \code{$weights} hold the per-edge entropy quantities. When \code{x} is a
#'   fitted network (or sequence data, from which a relative network is
#'   built), the result is that network with entropy weights swapped in -
#'   \code{$inits}, \code{$meta}, \code{$node_groups}, and node coordinates
#'   are inherited, so it plots with the same TNA styling and layout as its
#'   source. \code{$method} is \code{"entropy"}. \code{$params} carries
#'   \code{base}, \code{weight}, \code{entropy_rate}, and the stationary
#'   distribution \code{stationary}.
#'
#'   The object also declares the entropy house style through the
#'   \code{$meta$splot} producer contract (honoured by cograph >= 2.4.4):
#'   no minimum-weight pruning (bit values are smaller than probabilities),
#'   2-digit edge labels, vermilion edges, node rings showing the stationary
#'   distribution, and TNA rather than psychometric geometry. Because the
#'   contract states the styling outright, cograph needs no knowledge of the
#'   \code{"entropy"} method. \code{cograph::splot(ent)} therefore renders
#'   correctly with no arguments; any user argument overrides the contract.
#'
#' @details
#' Impossible transitions (\eqn{P_{ij} = 0}) get weight \code{0} under both
#' definitions (the \eqn{0 \log 0 := 0} convention), so the entropy network
#' has the same support as the transition network. Self-loops are retained
#' like any other edge. Deterministic transitions (\eqn{P_{ij} = 1}) also get
#' weight \code{0}: observing the inevitable carries no information.
#'
#' @examples
#' net <- build_network(group_regulation_long,
#'                      method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' ent <- entropy_network(net)
#' ent
#' \donttest{
#' if (requireNamespace("cograph", quietly = TRUE)) {
#'   cograph::splot(ent)
#' }
#' }
#'
#' @seealso \code{\link{transition_entropy}} for the chain- and state-level
#'   summary, \code{\link{entropy_bayes}} for credible intervals on the
#'   decomposition, \code{\link{build_network}}.
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information Theory},
#' 2nd ed., chapter 4. Wiley.
#'
#' Krejtz, K., Duchowski, A., Szmidt, T., Krejtz, I., Gonzalez Perilli, F.,
#' Pires, A., Vilaro, A., & Villalobos, N. (2015). Gaze transition entropy.
#' \emph{ACM Transactions on Applied Perception}, 13(1), 4:1-4:20.
#' \doi{10.1145/2834121}
#'
#' Krejtz, K., Hughes, C.J., Stasiak, I., Duchowski, A., & Krejtz, I.
#' (2025). Real-time mobile transition matrix entropy based on eye and head
#' movements. \emph{Proceedings of ETRA '25}. \doi{10.1145/3715669.3723128}
#'
#' @export
entropy_network <- function(x, base = 2,
                            weight = c("contribution", "surprisal",
                                       "production"),
                            scaling = c("none", "share", "chance"),
                            normalize = TRUE) {
  weight  <- match.arg(weight)
  scaling <- match.arg(scaling)
  stopifnot(is.numeric(base), length(base) == 1L, base > 0, base != 1)
  if (scaling == "chance" && weight != "surprisal") {
    stop("scaling = 'chance' applies to weight = 'surprisal' only.",
         call. = FALSE)
  }
  if (scaling == "share" && weight != "contribution") {
    stop("scaling = 'share' applies to weight = 'contribution' only.",
         call. = FALSE)
  }

  ## Resolve a source netobject so the result inherits its styling metadata
  ## (inits, meta, node_groups, coordinates) and renders identically to the
  ## transition network it decomposes.
  src <- NULL
  if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    src <- x
  } else if (is.data.frame(x) ||
             (is.matrix(x) && (is.character(x) || is.logical(x)))) {
    src <- build_network(as.data.frame(x, stringsAsFactors = FALSE),
                         method = "relative")
  }

  P           <- .mpt_extract_P(if (is.null(src)) x else src)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }

  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)

  pi <- .mpt_stationary(P)
  names(pi) <- state_names

  ## surprisal -log P, with 0 log 0 := 0; pi * (...) scales row i by pi[i]
  S <- ifelse(P > 0, -log(P, base = base), 0)
  n_oneway <- NULL
  W <- switch(weight,
    surprisal    = S,
    contribution = pi * P * S,
    production   = {
      Fm <- pi * P
      Rv <- t(Fm)
      two_way  <- Fm > 0 & Rv > 0
      n_oneway <- sum(Fm > 0 & Rv == 0)
      ifelse(two_way, Fm * log(ifelse(two_way, Fm / Rv, 1), base = base), 0)
    }
  )
  if (scaling == "chance") {
    W <- W / log(length(state_names), base = base)
  }
  if (scaling == "share") {
    W <- 100 * W / sum(W)
  }
  dimnames(W) <- list(state_names, state_names)

  if (is.null(src)) {
    inits <- if (inherits(x, "tna")) x$inits else NULL
    net   <- .wrap_netobject(W, data = NULL, method = "entropy",
                             directed = TRUE, inits = inits)
  } else {
    net           <- src
    net$weights   <- W
    net$edges     <- .extract_edges_from_matrix(W, directed = TRUE)
    net$n_edges   <- nrow(net$edges)
    net$directed  <- TRUE
    net$method    <- "entropy"
    net$scaling   <- NULL
    net$threshold <- 0
  }
  net$params <- list(
    base         = base,
    weight       = weight,
    scaling      = scaling,
    entropy_rate = sum(pi * P * S),
    stationary   = pi
  )
  if (!identical(scaling, "none")) net$scaling <- scaling
  if (identical(weight, "production")) {
    net$params$production_rate <- sum(W)
    net$params$n_oneway_pairs  <- n_oneway
  }

  ## Producer-side styling contract (honoured by cograph >= 2.4.4, silently
  ## ignored by older versions): the entropy house style. TNA geometry is
  ## requested explicitly -- `splot.netobject` picks tna vs psych styling from
  ## a method-name whitelist that predates "entropy", and its else branch sets
  ## `psych_styling <- TRUE` only when the argument is still NULL. Declaring
  ## both flags here fills them in before that branch runs, so cograph needs
  ## no knowledge of the "entropy" method at all. Otherwise, bit-valued weights
  ## are never pruned by the probability-calibrated minimum, edge labels get
  ## 2 digits, node rings show the stationary distribution (the row weight
  ## mu_i of the decomposition), and edges use Okabe-Ito vermilion so bits
  ## are not misread as transition probabilities. The production variant has
  ## signed weights, so it keeps cograph's positive/negative edge colouring
  ## instead of the flat vermilion.
  style <- list(
    minimum           = 0,
    edge_label_digits = if (identical(scaling, "share")) 1 else 2,
    donut_fill        = as.numeric(pi),
    donut_empty       = FALSE,
    tna_styling       = TRUE,
    psych_styling     = FALSE
  )
  if (!identical(weight, "production")) style$edge_color <- "#D55E00"
  net$meta$splot <- list(defaults = style)

  ## Chance-scaled surprisal ships a fully-styled rendering: the drawn edge
  ## set and widths come from the familiar transition probabilities (the
  ## contract's `weight` override), while surprisal is encoded visually -
  ## sequential blues for transitions likelier than chance (darker =
  ## stronger), dashed grey for those at or below chance (P <= 1/n).
  if (identical(weight, "surprisal") && identical(scaling, "chance")) {
    ramp     <- grDevices::colorRampPalette(c("#C6DBEF", "#08306B"))(100)
    strength <- pmin(pmax(1 - W, 0), 1)
    colm     <- matrix(ramp[1L + floor(99 * strength)], nrow(W), ncol(W),
                       dimnames = dimnames(W))
    colm[W >= 1] <- "grey72"
    ltym     <- ifelse(W >= 1, 2, 1)
    net$prob <- P
    net$meta$splot <- list(
      weight   = "prob",
      defaults = list(
        minimum           = 0,
        edge_label_digits = 2,
        edge_color        = colm,
        edge_style        = ltym,
        donut_fill        = as.numeric(pi),
        donut_empty       = FALSE,
        tna_styling       = TRUE,
        psych_styling     = FALSE
      )
    )
  }
  net
}

# ---- entropy_trajectory() ----

#' Sliding-Window Transition Entropy Trajectory
#'
#' @description
#' Tracks transition entropy over time: transitions are built within each
#' actor's sequence, pooled in temporal order, and a window of fixed size
#' slides across the stream, yielding one entropy estimate per window. This
#' turns transition entropy from a snapshot into a process measure -
#' declining entropy signals routinization, rising entropy exploration, and
#' level shifts mark phase changes (cf. Krejtz et al., 2025, who track gaze
#' transition entropy through task phases this way).
#'
#' @param data A long-format data.frame of timestamped events.
#' @param action Character. Name of the column holding the state/action.
#' @param actor Character or NULL. Column identifying sequences; transitions
#'   are only formed between consecutive events of the same actor. NULL
#'   treats the data as one sequence.
#' @param time Character or NULL. Timestamp column used to order events
#'   within actors and to place windows on a real time axis. NULL keeps the
#'   row order and uses the transition index as the axis.
#' @param group Character or NULL. Column splitting the data into parallel
#'   trajectories (e.g. condition, achievement level).
#' @param window Integer. Number of transitions per window (default
#'   \code{500}).
#' @param step Integer. Stride between window starts (default
#'   \code{window \%/\% 5}).
#' @param base Numeric. Logarithm base (default \code{2}, bits).
#'
#' @return An object of class \code{"net_entropy_trajectory"} with:
#' \describe{
#'   \item{trajectory}{Tidy data.frame, one row per window: \code{group},
#'     \code{window} (index), \code{time} (window midpoint; transition
#'     index when no \code{time} column), \code{time_start},
#'     \code{time_end}, \code{n_transitions}, \code{n_states} (distinct
#'     states in the window), \code{entropy} (bits per transition) and
#'     \code{entropy_norm} (divided by \eqn{\log_b} of the window's
#'     active-state count; in \eqn{[0, 1]}).}
#'   \item{window, step, base, states}{Call metadata; \code{states} is the
#'     global state set.}
#' }
#'
#' @details
#' Per-window entropy is the empirical conditional entropy
#' \eqn{-\sum_{ij} (n_{ij}/N) \log_b(n_{ij}/n_{i\cdot})} - rows weighted by
#' observed occupancy rather than the eigenvector stationary distribution.
#' Within a short window the chain is routinely non-ergodic (absorbing
#' fragments, unvisited states), where the eigenvector is undefined or
#' misleading; the empirical estimator is the standard windowed choice and
#' converges to the stationary entropy rate for long stationary stretches.
#'
#' Windows shorter than \code{window} at the tail are dropped; if the whole
#' stream is shorter than \code{window}, one window covering everything is
#' returned with a warning.
#'
#' @examples
#' \donttest{
#' tr <- entropy_trajectory(group_regulation_long,
#'                          action = "Action", actor = "Actor",
#'                          time = "Time", group = "Achiever")
#' tr
#' summary(tr)
#' plot(tr)
#' }
#'
#' @seealso \code{\link{transition_entropy}} for the whole-process snapshot,
#'   \code{\link{entropy_bayes}} for credible intervals on it.
#'
#' @references
#' Krejtz, K., Hughes, C.J., Stasiak, I., Duchowski, A., & Krejtz, I.
#' (2025). Real-Time Mobile Transition Matrix Entropy Based on Eye and Head
#' Movements. \emph{ETRA '25}. doi:10.1145/3715669.3723128
#'
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information Theory},
#' 2nd ed. Wiley.
#'
#' @export
entropy_trajectory <- function(data, action, actor = NULL, time = NULL,
                               group = NULL, window = 500L, step = NULL,
                               base = 2) {
  stopifnot(
    is.data.frame(data),
    is.character(action), length(action) == 1L, action %in% names(data),
    is.null(actor) || (is.character(actor) && length(actor) == 1L &&
                         actor %in% names(data)),
    is.null(time) || (is.character(time) && length(time) == 1L &&
                        time %in% names(data)),
    is.null(group) || (is.character(group) && length(group) == 1L &&
                         group %in% names(data)),
    is.numeric(window), length(window) == 1L, window >= 10,
    is.numeric(base), length(base) == 1L, base > 0, base != 1
  )
  window <- as.integer(window)
  if (is.null(step)) step <- max(1L, window %/% 5L)
  step <- as.integer(step)
  stopifnot(step >= 1L)

  ## Order events within actor (and group) by time, then build transitions
  ## between consecutive events of the same actor.
  grp_vals   <- if (is.null(group)) rep("all", nrow(data)) else
    as.character(data[[group]])
  actor_vals <- if (is.null(actor)) rep(1L, nrow(data)) else
    .observed_group_id(data[, c(if (!is.null(group)) group else NULL,
                                actor), drop = FALSE])
  time_vals  <- if (is.null(time)) seq_len(nrow(data)) else data[[time]]

  ord <- order(grp_vals, actor_vals, time_vals)
  a   <- as.character(data[[action]])[ord]
  g   <- grp_vals[ord]
  act <- actor_vals[ord]
  tv  <- time_vals[ord]

  states <- sort(unique(a[!is.na(a)]))
  S      <- length(states)
  ai     <- match(a, states)

  n    <- length(ai)
  keep <- act[-n] == act[-1] & g[-n] == g[-1] &
    !is.na(ai[-n]) & !is.na(ai[-1])
  from_i <- ai[-n][keep]
  to_i   <- ai[-1][keep]
  tr_t   <- tv[-1][keep]
  tr_g   <- g[-1][keep]

  if (length(from_i) < 2L) {
    stop("Fewer than two transitions; cannot build a trajectory.",
         call. = FALSE)
  }

  pair_enc <- (from_i - 1L) * S + to_i

  ## One window's empirical conditional entropy + bookkeeping
  win_stats <- function(idx) {
    cnt <- tabulate(pair_enc[idx], nbins = S * S)
    ni  <- tabulate(from_i[idx],   nbins = S)
    nz  <- which(cnt > 0)
    fi  <- ((nz - 1L) %/% S) + 1L
    N   <- length(idx)
    h   <- -sum(cnt[nz] / N * log(cnt[nz] / ni[fi], base = base))
    n_states <- length(unique(c(from_i[idx], to_i[idx])))
    ceiling_ <- log(max(n_states, 2L), base = base)
    c(h, h / ceiling_, N, n_states)
  }

  per_group <- lapply(split(seq_along(pair_enc), tr_g), function(gi) {
    n_tr <- length(gi)
    if (n_tr < window) {
      warning(sprintf(
        "Group stream (%d transitions) shorter than window (%d); using one window.",
        n_tr, window), call. = FALSE)
      starts <- 1L
      w_len  <- n_tr
    } else {
      starts <- seq.int(1L, n_tr - window + 1L, by = step)
      w_len  <- window
    }
    stats <- vapply(starts,
                    function(s) win_stats(gi[s:(s + w_len - 1L)]),
                    numeric(4))
    t_start <- tr_t[gi[starts]]
    t_end   <- tr_t[gi[starts + w_len - 1L]]
    data.frame(
      group         = tr_g[gi[1L]],
      window        = seq_along(starts),
      time          = t_start + (t_end - t_start) / 2,
      time_start    = t_start,
      time_end      = t_end,
      n_transitions = as.integer(stats[3, ]),
      n_states      = as.integer(stats[4, ]),
      entropy       = stats[1, ],
      entropy_norm  = stats[2, ],
      stringsAsFactors = FALSE
    )
  })
  trajectory <- do.call(rbind, per_group)
  rownames(trajectory) <- NULL

  structure(
    list(
      trajectory = trajectory,
      window     = window,
      step       = step,
      base       = base,
      states     = states
    ),
    class = "net_entropy_trajectory"
  )
}

#' Print method for `net_entropy_trajectory`
#'
#' @param x A `net_entropy_trajectory` object.
#' @param digits Integer. Digits to round numeric output. Default `3`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.net_entropy_trajectory <- function(x, digits = 3, ...) {
  tr <- x$trajectory
  cat(sprintf(
    "Entropy Trajectory (%d states; window %d transitions, step %d)\n\n",
    length(x$states), x$window, x$step))
  agg <- vapply(split(tr$entropy, tr$group), function(v) {
    c(length(v), mean(v), min(v), max(v))
  }, numeric(4))
  cat(sprintf("  %-12s %3.0f windows   mean %.*f   range [%.*f, %.*f]\n",
              paste0(colnames(agg), ":"), agg[1, ],
              digits, agg[2, ], digits, agg[3, ], digits, agg[4, ]))
  cat("\nUse plot() for the trajectory, $trajectory for the tidy table.\n")
  invisible(x)
}

#' Summary method for `net_entropy_trajectory`
#'
#' @param object A `net_entropy_trajectory` object.
#' @param ... Ignored.
#' @return Tidy per-group data.frame: windows, mean/sd/min/max entropy,
#'   entropy at the first and last window, and their difference (negative =
#'   routinization).
#' @export
summary.net_entropy_trajectory <- function(object, ...) {
  tr <- object$trajectory
  agg <- lapply(split(tr, tr$group), function(d) {
    data.frame(
      group     = d$group[1L],
      windows   = nrow(d),
      mean      = mean(d$entropy),
      sd        = stats::sd(d$entropy),
      min       = min(d$entropy),
      max       = max(d$entropy),
      first     = d$entropy[1L],
      last      = d$entropy[nrow(d)],
      change    = d$entropy[nrow(d)] - d$entropy[1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, agg)
  rownames(out) <- NULL
  out
}

#' Plot method for `net_entropy_trajectory`
#'
#' @description
#' Raw per-window entropy as faint lines with a loess-smoothed trend per
#' group, Okabe-Ito coloured. Declining trend = routinization; level shifts
#' = phase changes.
#'
#' @param x A `net_entropy_trajectory` object.
#' @param normalized Logical. Plot `entropy_norm` instead of raw bits
#'   (default `FALSE`).
#' @param span Numeric. Loess span (default `0.4`).
#' @param title Character. Plot title.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.net_entropy_trajectory <- function(x, normalized = FALSE, span = 0.4,
                                        title = "Transition entropy over time",
                                        ...) {
  tr <- x$trajectory
  tr$y <- if (normalized) tr$entropy_norm else tr$entropy
  ylab <- if (normalized) "Normalized entropy" else
    sprintf("Entropy (%s per transition)",
            if (x$base == 2) "bits" else sprintf("log_%g", x$base))

  single <- length(unique(tr$group)) == 1L
  p <- ggplot2::ggplot(tr, ggplot2::aes(x = .data$time, y = .data$y,
                                        colour = .data$group)) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.35) +
    ggplot2::geom_smooth(method = "loess", span = span, se = TRUE,
                         linewidth = 0.9, alpha = 0.15) +
    ggplot2::scale_colour_manual(values = .okabe_ito, name = NULL) +
    ggplot2::labs(x = if (inherits(tr$time, "POSIXct")) "Time" else
      "Transition index",
      y = ylab, title = title,
      subtitle = sprintf("Sliding window: %d transitions, step %d",
                         x$window, x$step)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = if (single) "none" else "bottom")
  p
}

# ---- entropy_bayes() ----

#' @noRd
.te_extract_counts <- function(x) {
  if (is.matrix(x) && is.numeric(x)) {
    if (any(x < 0) || any(x != floor(x))) {
      stop("Bayesian entropy estimation needs transition COUNTS ",
           "(non-negative whole numbers), not probabilities. Supply a ",
           "count matrix, a frequency network ",
           "(build_network(method = 'frequency')), or sequence data.",
           call. = FALSE)
    }
    if (is.null(rownames(x))) {
      dimnames(x) <- list(paste0("S", seq_len(nrow(x))),
                          paste0("S", seq_len(nrow(x))))
    }
    return(x)
  }
  if (inherits(x, "netobject") || inherits(x, "cograph_network") ||
      inherits(x, "tna")) {
    w <- x$weights
    if (is.matrix(w) && is.numeric(w) && all(w >= 0) &&
        all(w == floor(w)) && sum(w) > 0) {
      return(w)  # already a frequency (count) network
    }
    if (!is.null(x$data) && is.data.frame(x$data)) {
      cw <- build_network(x$data, method = "frequency")$weights
      if (!identical(rownames(cw), rownames(w)) &&
          setequal(rownames(cw), rownames(w))) {
        cw <- cw[rownames(w), rownames(w)]
      }
      return(cw)
    }
    stop("Cannot recover transition counts from this object. Supply a ",
         "frequency network (build_network(method = 'frequency')), a count ",
         "matrix, or the sequence data.", call. = FALSE)
  }
  if (is.data.frame(x) ||
      (is.matrix(x) && (is.character(x) || is.logical(x)))) {
    return(build_network(as.data.frame(x, stringsAsFactors = FALSE),
                         method = "frequency")$weights)
  }
  stop("'x' must be a count matrix, a netobject / cograph_network / tna ",
       "object, or a sequence data.frame.", call. = FALSE)
}

#' Bayesian Transition Entropy
#'
#' @description
#' Bayesian estimation of the transition entropy quantities of
#' \code{\link{transition_entropy}} and the edge-level decomposition of
#' \code{\link{entropy_network}}. Each row of the transition matrix gets an
#' independent Dirichlet posterior (counts + \code{prior}); Monte Carlo
#' draws propagate count uncertainty into the entropy rate, the per-state
#' branching entropies, and every edge's entropy contribution, yielding
#' posterior means and credible intervals.
#'
#' The practical purpose is to \strong{exclude unstable estimates}: an edge
#' whose contribution rests on a handful of observations has a wide
#' posterior, and is flagged non-credible unless it credibly accounts for at
#' least \code{min_share} of the process entropy. \code{$model} is the
#' entropy network with non-credible edges zeroed - the stable entropy
#' skeleton.
#'
#' @param x A frequency \code{netobject}
#'   (\code{build_network(method = "frequency")}), any \code{netobject} that
#'   carries its \code{$data} (counts are rebuilt automatically), a
#'   count matrix, or a wide sequence data.frame. Group dispatch on
#'   \code{netobject_group}.
#' @param prior Numeric. Dirichlet prior concentration added to every cell
#'   of the count matrix (default \code{0.5}, Jeffreys).
#' @param draws Integer. Number of Monte Carlo posterior draws
#'   (default \code{4000}).
#' @param ci Numeric in (0, 1). Credible interval mass (default \code{0.95}).
#' @param min_share Numeric in [0, 1). An edge is credible when the lower
#'   bound of the credible interval of its \emph{share} of the entropy rate
#'   exceeds this value (default \code{0.01}: the edge credibly accounts for
#'   at least 1% of process entropy).
#' @param base Numeric. Logarithm base (default \code{2}, bits).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"net_entropy_bayes"} with:
#' \describe{
#'   \item{summary}{Tidy data.frame - one row per chain-level quantity
#'     (\code{entropy_rate}, \code{stationary_entropy}, \code{redundancy}),
#'     with posterior \code{mean}, \code{sd}, \code{ci_lower},
#'     \code{ci_upper}.}
#'   \item{states}{Tidy data.frame - one row per state: posterior mean/CI of
#'     the row entropy and of the stationary probability.}
#'   \item{edges}{Tidy data.frame - one row per observed transition:
#'     posterior mean/sd/CI of the contribution (bits), posterior mean/CI of
#'     its share of the entropy rate, and the \code{credible} flag.}
#'   \item{network}{\code{netobject} - posterior-mean entropy network
#'     (all edges), carrying the entropy house style.}
#'   \item{model}{\code{netobject} - the pruned entropy network: posterior
#'     means where \code{credible}, \code{0} elsewhere.}
#'   \item{draws_entropy_rate}{Numeric vector of posterior entropy-rate
#'     draws (for further analysis or plotting).}
#'   \item{prior, draws, ci, min_share, base, states_names}{Call metadata.}
#' }
#'
#' @details
#' With \code{prior > 0} the posterior puts mass on every transition, so
#' contribution draws are strictly positive and a naive "CI excludes zero"
#' rule would flag every edge as credible. The share criterion is used
#' instead: an edge is stable when it credibly carries at least
#' \code{min_share} of \eqn{h(P)}. Unobserved transitions (count 0) get
#' only prior mass and are never credible under any sensible
#' \code{min_share}.
#'
#' The posterior mean entropy rate is typically slightly below the plug-in
#' estimate on sparse data (Dirichlet smoothing pulls rows toward uniform
#' but averages over uncertainty); the difference vanishes as counts grow.
#'
#' @examples
#' \donttest{
#' net <- build_network(group_regulation_long, method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' eb <- entropy_bayes(net, seed = 1)
#' eb
#' summary(eb)
#' plot(eb)
#' }
#'
#' @seealso \code{\link{transition_entropy}}, \code{\link{entropy_network}},
#'   \code{\link{bayes_compare}}
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information Theory},
#' 2nd ed. Wiley.
#'
#' @export
entropy_bayes <- function(x, prior = 0.5, draws = 4000, ci = 0.95,
                          min_share = 0.01, base = 2, seed = NULL) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      entropy_bayes(net, prior = prior, draws = draws, ci = ci,
                    min_share = min_share, base = base, seed = seed)
    })
    class(out) <- c("net_entropy_bayes_group", "list")
    return(out)
  }

  stopifnot(
    is.numeric(prior), length(prior) == 1L, prior > 0,
    is.numeric(draws), length(draws) == 1L, draws >= 100,
    is.numeric(ci), length(ci) == 1L, ci > 0, ci < 1,
    is.numeric(min_share), length(min_share) == 1L,
    min_share >= 0, min_share < 1,
    is.numeric(base), length(base) == 1L, base > 0, base != 1
  )
  if (!is.null(seed)) set.seed(seed)

  counts <- .te_extract_counts(x)
  states <- rownames(counts)
  n      <- nrow(counts)
  draws  <- as.integer(draws)

  ## Posterior: independent Dirichlet(counts + prior) per row, sampled via
  ## normalized Gamma draws; array is n x n x draws, filled column-major so
  ## rep(alpha, draws) lays one full matrix per slice.
  alpha <- counts + prior
  A  <- array(stats::rgamma(n * n * draws, shape = rep(as.vector(alpha), draws)),
              dim = c(n, n, draws))
  rs <- apply(A, c(1, 3), sum)                                # n x draws
  P  <- A / aperm(array(rs, c(n, draws, n)), c(1, 3, 2))
  S  <- -log(P, base = base)                                  # P > 0 a.s.
  H  <- apply(P * S, c(1, 3), sum)                            # row entropy
  pi_mat <- vapply(seq_len(draws),
                   function(d) .mpt_stationary(P[, , d]), numeric(n))
  h    <- colSums(pi_mat * H)                                 # entropy rate
  H_pi <- colSums(-pi_mat * log(pi_mat, base = base))
  red  <- H_pi - h
  Cn   <- P * S * aperm(array(pi_mat, c(n, draws, n)), c(1, 3, 2))
  Shr  <- Cn / aperm(array(matrix(h, n, draws, byrow = TRUE),
                           c(n, draws, n)), c(1, 3, 2))       # share of h

  lo <- (1 - ci) / 2
  hi <- 1 - lo
  qs <- function(v) c(mean(v), stats::sd(v),
                      stats::quantile(v, c(lo, hi), names = FALSE))

  chain <- t(vapply(list(entropy_rate = h, stationary_entropy = H_pi,
                         redundancy = red), qs, numeric(4)))
  chain_df <- data.frame(
    quantity = rownames(chain),
    mean = chain[, 1], sd = chain[, 2],
    ci_lower = chain[, 3], ci_upper = chain[, 4],
    row.names = NULL, stringsAsFactors = FALSE
  )

  Hq  <- apply(H, 1, qs)                                      # 4 x n
  piq <- apply(pi_mat, 1, qs)
  states_df <- data.frame(
    state = states,
    row_entropy = Hq[1, ], row_entropy_sd = Hq[2, ],
    row_entropy_lower = Hq[3, ], row_entropy_upper = Hq[4, ],
    stationary = piq[1, ],
    stationary_lower = piq[3, ], stationary_upper = piq[4, ],
    row.names = NULL, stringsAsFactors = FALSE
  )

  ## Edge table restricted to observed transitions - unobserved cells carry
  ## only prior mass and are reported nowhere.
  obs <- which(counts > 0, arr.ind = TRUE)
  edge_draws  <- matrix(Cn[cbind(rep(obs[, 1], draws), rep(obs[, 2], draws),
                                 rep(seq_len(draws), each = nrow(obs)))],
                        nrow = nrow(obs))
  share_draws <- matrix(Shr[cbind(rep(obs[, 1], draws), rep(obs[, 2], draws),
                                  rep(seq_len(draws), each = nrow(obs)))],
                        nrow = nrow(obs))
  eq <- apply(edge_draws, 1, qs)
  sq <- apply(share_draws, 1, qs)
  edges_df <- data.frame(
    from = states[obs[, 1]],
    to   = states[obs[, 2]],
    count = counts[obs],
    contribution = eq[1, ], sd = eq[2, ],
    ci_lower = eq[3, ], ci_upper = eq[4, ],
    share = sq[1, ], share_lower = sq[3, ], share_upper = sq[4, ],
    row.names = NULL, stringsAsFactors = FALSE
  )
  edges_df$credible <- edges_df$share_lower > min_share
  edges_df <- edges_df[order(-edges_df$contribution), ]
  rownames(edges_df) <- NULL

  ## Posterior-mean entropy network (all observed edges) and the pruned
  ## stable model, both carrying the entropy house style.
  M <- apply(Cn, c(1, 2), mean)
  M[counts == 0] <- 0
  dimnames(M) <- list(states, states)
  network <- .wrap_netobject(M, data = NULL, method = "entropy",
                             directed = TRUE)
  cred_mat <- matrix(FALSE, n, n, dimnames = dimnames(M))
  cred_mat[cbind(match(edges_df$from, states),
                 match(edges_df$to, states))] <- edges_df$credible
  Mp <- M * cred_mat
  model <- .wrap_netobject(Mp, data = NULL, method = "entropy",
                           directed = TRUE)
  pi_mean <- rowMeans(pi_mat)
  style <- list(minimum = 0, edge_label_digits = 2, edge_color = "#D55E00",
                donut_fill = as.numeric(pi_mean), donut_empty = FALSE)
  ep <- list(base = base, weight = "contribution",
             entropy_rate = mean(h), stationary = stats::setNames(pi_mean, states))
  network$params <- ep
  model$params   <- ep
  network$meta$splot <- list(defaults = style)
  model$meta$splot   <- list(defaults = style)

  structure(
    list(
      summary            = chain_df,
      states             = states_df,
      edges              = edges_df,
      network            = network,
      model              = model,
      draws_entropy_rate = h,
      prior              = prior,
      draws              = draws,
      ci                 = ci,
      min_share          = min_share,
      base               = base,
      states_names       = states
    ),
    class = "net_entropy_bayes"
  )
}

#' Print method for `net_entropy_bayes`
#'
#' @param x A `net_entropy_bayes` object.
#' @param digits Integer. Digits to round numeric output. Default `3`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.net_entropy_bayes <- function(x, digits = 3, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf(
    "Bayesian Transition Entropy (%d states, %s; Dirichlet prior %.2g, %d draws)\n\n",
    length(x$states_names), unit, x$prior, x$draws))
  s <- x$summary
  cat(sprintf("  %-19s %.*f [%.*f, %.*f]\n",
              paste0(s$quantity, ":"),
              digits, s$mean, digits, s$ci_lower, digits, s$ci_upper))
  n_cred <- sum(x$edges$credible)
  cat(sprintf(
    "\nEdges: %d observed; %d credible (share of h(P) credibly > %g%%).\n",
    nrow(x$edges), n_cred, 100 * x$min_share))
  cat("Use summary() for the edge table, plot() for the posterior, and\n")
  cat("$model for the pruned stable entropy network.\n")
  invisible(x)
}

#' Print method for `net_entropy_bayes_group`
#'
#' @param x A `net_entropy_bayes_group`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.net_entropy_bayes_group <- function(x, ...) {
  cat(sprintf("Bayesian Transition Entropy - %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  hs <- vapply(x, function(o) {
    s <- o$summary
    sprintf("%.3f [%.3f, %.3f]",
            s$mean[1], s$ci_lower[1], s$ci_upper[1])
  }, character(1))
  cat("Entropy rate per group (posterior mean [CI]):\n")
  cat(sprintf("  %-12s %s\n", paste0(names(hs), ":"), hs))
  invisible(x)
}

#' Summary method for `net_entropy_bayes`
#'
#' @param object A `net_entropy_bayes` object.
#' @param ... Ignored.
#' @return The tidy edge table (data.frame), sorted by posterior mean
#'   contribution, invisibly printed with the chain-level summary.
#' @export
summary.net_entropy_bayes <- function(object, ...) {
  cat("Chain-level posterior:\n")
  print(cbind(object$summary[1],
              round(object$summary[-1], 4)), row.names = FALSE)
  cat("\nPer-edge posterior (sorted by contribution):\n")
  tbl <- object$edges
  num <- vapply(tbl, is.numeric, logical(1))
  tbl[num] <- lapply(tbl[num], round, 4)
  print(tbl, row.names = FALSE)
  invisible(object$edges)
}

#' Plot method for `net_entropy_bayes`
#'
#' @description
#' Forest plot of the per-edge entropy contributions: posterior mean and
#' credible interval, credible edges in Okabe-Ito blue, unstable
#' (non-credible) edges in grey. The dashed line marks \code{min_share} of
#' the posterior-mean entropy rate - the stability criterion.
#'
#' @param x A `net_entropy_bayes` object.
#' @param top Integer. Show at most this many edges, by posterior mean
#'   contribution (default `25`).
#' @param title Character. Plot title.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.net_entropy_bayes <- function(x, top = 25,
                                   title = "Bayesian edge entropy contributions",
                                   ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  df <- utils::head(x$edges, top)
  ## ASCII arrow: a graphics device in a non-UTF8 locale cannot convert
  ## U+2192 and errors outright in grid's text measurement (mbcsToSbcs),
  ## rather than degrading. Matches the "A -> B" arrow notation the HON
  ## family already uses. Console-only arrows elsewhere are unaffected.
  df$edge <- factor(paste(df$from, "->", df$to),
                    levels = rev(paste(df$from, "->", df$to)))
  h_mean <- x$summary$mean[x$summary$quantity == "entropy_rate"]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$contribution, y = .data$edge,
                                   colour = .data$credible)) +
    ggplot2::geom_vline(xintercept = x$min_share * h_mean,
                        linetype = "dashed", colour = "grey40",
                        linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_lower,
                                        xmax = .data$ci_upper),
                           orientation = "y", width = 0.25, linewidth = 0.5) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#0072B2", "FALSE" = "grey60"),
      labels = c("TRUE" = "credible", "FALSE" = "unstable"),
      name = NULL
    ) +
    ggplot2::labs(
      x = sprintf("Contribution to h(P) (%s)", unit),
      y = NULL, title = title,
      subtitle = sprintf(
        "Posterior mean and %g%% CI; dashed line = %g%% of h(P)",
        100 * x$ci, 100 * x$min_share)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

#' Print method for `net_transition_entropy`
#'
#' @param x A `net_transition_entropy` object.
#' @param digits Integer. Digits to round numeric output. Default `3`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.net_transition_entropy <- function(x, digits = 3, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf("Transition Entropy (%d states, %s; ceiling = %.3f)\n\n",
              length(x$states), unit, x$max_entropy))
  cat(sprintf("                          raw            normalised\n"))
  cat(sprintf("  Entropy rate    h(P)  = %.*f %-6s  %.3f\n",
              digits, x$entropy_rate, unit, x$entropy_rate_norm))
  cat(sprintf("  Stationary    H(pi)  = %.*f %-6s  %.3f\n",
              digits, x$stationary_entropy, unit, x$stationary_entropy_norm))
  cat(sprintf("  Redundancy   H(pi)-h = %.*f %-6s  %.3f\n",
              digits, x$redundancy, unit, x$redundancy_norm))
  cat(sprintf(paste0(
    "\nNormalised: h(P) and H(pi) are raw / log_%g(n_states) ",
    "(0 = deterministic, 1 = uniform);\n",
    "  redundancy is the relative redundancy (H(pi) - h(P)) / H(pi), ",
    "not raw / log_%g(n_states).\n"),
    x$base, x$base))
  invisible(x)
}

#' Print method for `net_transition_entropy_group`
#'
#' @param x A `net_transition_entropy_group`.
#' @param ... Forwarded to `print.net_transition_entropy`.
#' @return `x` invisibly.
#' @export
print.net_transition_entropy_group <- function(x, ...) {
  cat(sprintf("Transition Entropy - %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  rates <- vapply(x, function(o) o$entropy_rate, numeric(1))
  cat("Entropy rate per group:\n")
  print(round(rates, 4))
  invisible(x)
}

#' Summary method for `net_transition_entropy`
#'
#' @description
#' Returns a tidy per-state contribution table sorted by share of the
#' chain-level entropy rate (largest first), so the dominant contributors
#' to \eqn{h(P)} are visible at a glance. Each row contains the
#' stationary mass, the raw and normalised row entropy, the additive
#' contribution \eqn{\pi_i H(P_{i\cdot})}, and that contribution as a
#' percentage of \eqn{h(P)}.
#'
#' @param object A `net_transition_entropy` object.
#' @param ... Ignored.
#' @return A `summary.net_transition_entropy` containing
#'   \describe{
#'     \item{table}{tidy per-state data.frame, sorted by `contribution_pct`
#'       descending}
#'     \item{chain}{tidy chain-level data.frame with raw and normalised
#'       \eqn{h(P)}, \eqn{H(\pi)}, redundancy, and ceiling}
#'     \item{base}{logarithm base used}
#'   }
#' @export
summary.net_transition_entropy <- function(object, ...) {
  contrib  <- object$stationary * object$row_entropy
  total    <- object$entropy_rate
  pct      <- if (total > 0) 100 * contrib / total else contrib * 0

  per_state <- data.frame(
    state            = object$states,
    stationary       = unname(object$stationary),
    row_entropy      = unname(object$row_entropy),
    row_entropy_norm = unname(object$row_entropy_norm),
    contribution     = unname(contrib),
    contribution_pct = unname(pct),
    stringsAsFactors = FALSE
  )
  per_state <- per_state[order(-per_state$contribution_pct), ]
  rownames(per_state) <- NULL

  unit <- switch(as.character(object$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", object$base))
  if (isTRUE(all.equal(object$base, exp(1)))) unit <- "nats"

  chain <- data.frame(
    quantity = c("entropy_rate h(P)",
                 "stationary H(pi)",
                 "redundancy H(pi)-h(P)",
                 sprintf("ceiling log_%g(n)", object$base)),
    raw        = c(object$entropy_rate,
                   object$stationary_entropy,
                   object$redundancy,
                   object$max_entropy),
    normalised = c(object$entropy_rate_norm,
                   object$stationary_entropy_norm,
                   object$redundancy_norm,
                   1),
    unit       = unit,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      table = per_state,
      chain = chain,
      base  = object$base
    ),
    class = "summary.net_transition_entropy"
  )
}

#' Print method for `summary.net_transition_entropy`
#'
#' @param x A `summary.net_transition_entropy`.
#' @param digits Integer. Digits to round numeric output. Default `3`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.summary.net_transition_entropy <- function(x, digits = 3, ...) {
  unit <- x$chain$unit[1L]

  cat(sprintf("Transition Entropy Summary (%s)\n\n", unit))

  cat("Per-state contribution to h(P):\n")
  tbl <- x$table
  tbl$stationary       <- round(tbl$stationary,       digits)
  tbl$row_entropy      <- round(tbl$row_entropy,      digits)
  tbl$row_entropy_norm <- round(tbl$row_entropy_norm, digits)
  tbl$contribution     <- round(tbl$contribution,     digits)
  tbl$contribution_pct <- round(tbl$contribution_pct, 1)
  print(tbl, row.names = FALSE)

  cat("\nChain-level summary:\n")
  ch <- x$chain
  ch$raw        <- round(ch$raw,        digits)
  ch$normalised <- round(ch$normalised, digits)
  ch$unit       <- NULL
  print(ch, row.names = FALSE)

  cat(sprintf(paste0(
    "\nNormalised: row entropy, h(P) and H(pi) are raw / log_%g(n_states), ",
    "in [0, 1];\n",
    "  redundancy is the relative redundancy (H(pi) - h(P)) / H(pi).\n"),
    x$base))
  invisible(x)
}

#' Plot method for `net_transition_entropy`
#'
#' @description
#' Bar chart of per-state row entropy with overlaid horizontal lines at
#' the entropy rate \eqn{h(P)} (chain-level summary) and the maximum row
#' entropy \eqn{\log_b n} (uniform branching). Bar widths are proportional
#' to the stationary probability so the visual area sums to the entropy
#' rate.
#'
#' @param x A `net_transition_entropy` object.
#' @param title Character. Plot title.
#' @param fill Character. Bar fill colour. Default Okabe-Ito blue.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.net_transition_entropy <- function(x,
                                        title = "Transition Entropy",
                                        fill  = "#0072B2",
                                        ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  n      <- length(x$states)
  H_max  <- log(n, base = x$base)

  df <- data.frame(
    state       = factor(x$states, levels = x$states),
    row_entropy = x$row_entropy,
    stationary  = x$stationary,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$state,
                                   y = .data$row_entropy,
                                   width = .data$stationary /
                                     max(.data$stationary))) +
    ggplot2::geom_col(fill = fill, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = x$entropy_rate,
                        linetype = "dashed", colour = "#D55E00",
                        linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = H_max,
                        linetype = "dotted", colour = "grey30",
                        linewidth = 0.4) +
    ggplot2::annotate("text",
                      x = 0.6, y = x$entropy_rate,
                      label = sprintf("h(P) = %.3f", x$entropy_rate),
                      hjust = 0, vjust = -0.4,
                      colour = "#D55E00", size = 3.4) +
    ggplot2::annotate("text",
                      x = 0.6, y = H_max,
                      label = sprintf("log_%g(n) = %.3f", x$base, H_max),
                      hjust = 0, vjust = -0.4,
                      colour = "grey30", size = 3.2) +
    ggplot2::labs(x = "State",
                  y = sprintf("Row entropy (%s)", unit),
                  title    = title,
                  subtitle = "Bar width proportional to stationary probability") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title         = ggplot2::element_text(face = "bold")
    )
}
