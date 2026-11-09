# ---- Transition Entropy ----

#' @noRd
# A chain with a zero stationary probability is not ergodic; the entropy
# verbs report it the same way.
# Name of the unit for a logarithm base: bits, nats, hartleys, or log_b.
.te_unit <- function(base) {
  if (isTRUE(all.equal(base, 2))) return("bits")
  if (isTRUE(all.equal(base, exp(1)))) return("nats")
  if (isTRUE(all.equal(base, 10))) return("hartleys")
  sprintf("log_%g", base)
}

.te_warn_not_ergodic <- function() {
  warning(warningCondition(
    "Non-positive stationary probabilities; chain may not be ergodic.",
    class = "nestimate_not_ergodic", call = NULL))
}

.te_row_entropy <- function(P, base) {
  ## 0 * log 0 := 0; vectorised, no loops
  L      <- ifelse(P > 0, log(P, base = base), 0)
  -rowSums(P * L)
}

# ---- transition_entropy() ----

#' Transition Entropy of a Markov Chain
#'
#' @description
#' `transition_entropy()` measures how uncertain a Markov chain's next step
#' is, at three levels: the branching entropy of each state (how evenly its
#' outgoing transitions spread), the entropy of the stationary distribution,
#' and the chain's entropy rate, the average uncertainty per step. It takes a
#' transition network, a `tna` model, a transition matrix or wide sequence
#' data; a `netobject_group` gives one result per group. The entropy rate
#' \deqn{h(P) = -\sum_i \pi_i \sum_j P_{ij} \log_b P_{ij},}
#' with \eqn{\pi} the stationary distribution, is the per-step uncertainty of
#' the chain's trajectories (Cover & Thomas, 2006, ch. 4). It was introduced
#' to behavioural research as gaze transition entropy by Krejtz et al. (2015)
#' and tracked in real time as mobile transition matrix entropy by Krejtz et
#' al. (2025). The normalised values divide by \eqn{\log_b n}, the entropy of
#' uniform branching among \eqn{n} states; they are the scale-free variants
#' those papers report.
#'
#' An entropy rate near 0 means the current state almost fixes the next one;
#' near the ceiling \eqn{\log_b n} (normalised 1), every next state is about
#' equally likely. The redundancy \eqn{H(\pi) - h(P)} is the uncertainty that
#' knowing the current state removes, zero when successive states are
#' independent. All of these are point estimates from one transition matrix:
#' [entropy_bayes()] gives their credible intervals, [entropy_network()]
#' shows which transitions carry the entropy, and [entropy_trajectory()] how
#' it changes over time. A chain with more than one closed class has more
#' than one stationary distribution; check its classes with
#' [chain_structure()].
#'
#' @param x The chain: a `netobject` or `cograph_network` with non-negative
#'   transition weights, a `tna` model, a numeric transition matrix, or wide
#'   sequence data (a data frame or a character matrix, one row per sequence)
#'   that holds a numbered run of positions (T1, T2, ...), from which a
#'   relative transition network is built; other columns are set aside with
#'   a message. For a frame without such a run, build the network with
#'   [build_network()] (`state_cols =` or `action =`) and pass that. A `netobject_group` gives one result per group. An association
#'   network (`cor`, `pcor`, `glasso`, `ising`, `mgm`) or a matrix with
#'   negative entries raises an error of class `nestimate_not_markov`. For the
#'   `print()` and `plot()` methods: a `transition_entropy()` result;
#'   `print()` also takes a group result or the result of `summary()`.
#' @param base The logarithm base, which sets the unit: `2` (default) for
#'   bits, `exp(1)` for nats, `10` for hartleys. Any positive number other
#'   than 1.
#' @param normalize Whether rows that do not sum to 1 (a frequency network, a
#'   count matrix) are divided by their sums, with a warning. With `FALSE`,
#'   such rows raise an error. Default `TRUE`. A state with no outgoing
#'   transitions (a row summing to zero) is an error either way.
#'
#' @return A `net_transition_entropy` object, or for a `netobject_group` a
#'   `net_transition_entropy_group`: a named list of them, one per group.
#'   `print()` shows the entropy rate, stationary entropy and redundancy, raw
#'   and normalised, and for a group the entropy rate of each group. For a
#'   group, `summary()` returns a data frame with one row per group (`group`,
#'   `n_states`, the raw and normalised entropy rate and stationary entropy,
#'   `redundancy` and `unit`), and `plot()` draws the bar chart below in one
#'   panel per group.
#'   `summary()` returns a `summary.net_transition_entropy` object, printed by
#'   its own `print()` method, that holds `table`, one row per state
#'   (stationary probability, raw and normalised branching entropy, the
#'   state's additive contribution \eqn{\pi_i H(P_{i\cdot})} to the entropy
#'   rate and that contribution as a percentage of it), sorted by contribution,
#'   largest first; `chain`, one row per chain-level quantity (entropy rate,
#'   stationary entropy, redundancy and the ceiling, raw and normalised),
#'   with the unit; and the logarithm `base`.
#'   `plot()` returns a `ggplot` bar chart of each state's branching entropy,
#'   with bar widths proportional to the stationary probability (so a bar's
#'   area is proportional to the state's contribution), a dashed line at the
#'   entropy rate \eqn{h(P)} and a dotted line at the ceiling \eqn{\log_b n}
#'   (uniform branching). The `print()` methods return their input
#'   invisibly. The object's fields:
#' \describe{
#'   \item{row_entropy}{Named numeric vector of length \eqn{n}: each state's
#'     branching entropy \eqn{H(P_{i\cdot}) = -\sum_j P_{ij} \log_b P_{ij}}.}
#'   \item{row_entropy_norm}{`row_entropy` divided by the ceiling
#'     \eqn{\log_b n}, in \eqn{[0, 1]}; all zeros when \eqn{n = 1}.}
#'   \item{stationary}{Named numeric vector: the stationary distribution
#'     \eqn{\pi}.}
#'   \item{stationary_entropy}{\eqn{H(\pi) = -\sum_i \pi_i \log_b \pi_i}, the
#'     entropy of \eqn{\pi} read as an independent draw each step; an upper
#'     bound on the entropy rate.}
#'   \item{stationary_entropy_norm}{`stationary_entropy` divided by the
#'     ceiling \eqn{\log_b n}.}
#'   \item{entropy_rate}{\eqn{h(P) = \sum_i \pi_i H(P_{i\cdot})}, the
#'     Shannon-McMillan-Breiman entropy rate.}
#'   \item{entropy_rate_norm}{`entropy_rate` divided by the ceiling
#'     \eqn{\log_b n}.}
#'   \item{redundancy}{\eqn{H(\pi) - h(P)}, the entropy deficit due to serial
#'     dependence; zero for an independent chain (every row of \eqn{P} equal
#'     to \eqn{\pi}).}
#'   \item{redundancy_norm}{The relative redundancy
#'     \eqn{(H(\pi) - h(P)) / H(\pi)}, the share of the stationary entropy
#'     that order-1 memory removes; \strong{not} `redundancy` divided by
#'     \eqn{\log_b n}. `0` when \eqn{H(\pi) = 0}.}
#'   \item{max_entropy}{The normalising ceiling \eqn{\log_b n}.}
#'   \item{base}{The logarithm base used.}
#'   \item{states}{Character vector of state names.}
#' }
#'
#' @details
#' The convention \eqn{0 \log 0 = 0} applies, so an absorbing or
#' deterministic row has zero branching entropy. The chain need not be
#' irreducible: \eqn{\pi} is the eigenvector of \eqn{P^\top} for the
#' eigenvalue closest to 1, as elsewhere in the package. For a chain with
#' more than one closed class this is one stationary distribution among
#' many; when some of
#' its probabilities are zero, a warning says so, and they are set to machine
#' epsilon and the distribution renormalised before the entropies are
#' computed.
#'
#' \eqn{h(P) \le H(\pi)} holds, with equality if and only if successive
#' states are independent. The deficit \eqn{H(\pi) - h(P)} is reported as
#' `redundancy`, a measure of how much memory the chain has at order 1.
#'
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' te  <- transition_entropy(net)
#' print(te)
#' summary(te)
#' plot(te)
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

  .assert(is.numeric(base), length(base) == 1L, base > 0, base != 1)

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
    .te_warn_not_ergodic()
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
#' `entropy_network()` splits a chain's entropy rate into one term per
#' transition and returns the terms as a network. The entropy rate
#' \eqn{h(P) = -\sum_{ij} \pi_i P_{ij} \log_b P_{ij}}
#' ([transition_entropy()]; Krejtz et al., 2015, 2025) is a sum over
#' transitions, so every edge \eqn{i \to j} owns the exact term
#' \eqn{\pi_i P_{ij} \log_b(1/P_{ij})} of the chain-level uncertainty. No new
#' quantity is estimated: the network places the summands of the
#' entropy-rate equation on the transition graph, which shows where in the
#' process the uncertainty lies. Two other edge weights are available: each
#' transition's surprisal, and its contribution to the entropy production
#' rate (irreversibility). The result is a regular `netobject`, so it prints,
#' summarises and plots (`cograph::splot()`) like any other Nestimate
#' network.
#'
#' Under the default weight, a heavy edge is a transition that is both
#' frequent and unpredictable; under `"surprisal"`, a heavy edge is a rare
#' one, whatever its frequency. The weights are computed from one estimated
#' transition matrix, so an edge that rests on a few observations can carry a
#' large weight by chance: [entropy_bayes()] gives each contribution a
#' credible interval and keeps only the edges that credibly hold a share of
#' the entropy rate.
#'
#' @param x The chain: a `netobject` or `cograph_network` with non-negative
#'   transition weights, a `tna` model, a numeric transition matrix, or wide
#'   sequence data (a data frame or a character matrix, one row per sequence)
#'   that holds a numbered run of positions (T1, T2, ...), from which a
#'   relative transition network is built; other columns are set aside with
#'   a message. For a frame without such a run, build the network with
#'   [build_network()] (`state_cols =` or `action =`) and pass that. An association network (`cor`, `pcor`, `glasso`, `ising`, `mgm`)
#'   or a matrix with negative entries raises an error of class
#'   `nestimate_not_markov`. A `netobject_group` is not accepted; pass one
#'   network at a time.
#' @param base The logarithm base, which sets the unit: `2` (default) for
#'   bits, `exp(1)` for nats, `10` for hartleys. Any positive number other
#'   than 1.
#' @param weight What each edge's weight measures:
#'   \describe{
#'     \item{`"contribution"`}{(default) \eqn{\pi_i P_{ij} \log_b(1/P_{ij})},
#'       the edge's additive share of the entropy rate; all edge weights sum
#'       exactly to \eqn{h(P)}. Heavy edges are transitions that are both
#'       frequent \emph{and} unpredictable.}
#'     \item{`"surprisal"`}{\eqn{\log_b(1/P_{ij})}, the information content of
#'       observing the transition, ignoring how often it occurs. Heavy edges
#'       are rare, surprising transitions; a low surprisal marks an expected
#'       transition, which is the edge-level reading of predictability.}
#'     \item{`"production"`}{\eqn{F_{ij} \log_b(F_{ij}/F_{ji})} with
#'       \eqn{F_{ij} = \pi_i P_{ij}}, the edge's contribution to the entropy
#'       production rate (irreversibility). It is positive on the dominant
#'       direction of a pair and negative on the reverse, and the two always
#'       sum to \eqn{(F_{ij}-F_{ji}) \log_b(F_{ij}/F_{ji}) \ge 0}. A pair with
#'       flow in one direction only has infinite production and is left out
#'       (weight 0); the number of such pairs is stored with the network's
#'       parameters. The total over the included pairs is the production rate.
#'       It is `0` for a reversible chain (detailed balance); when one-way
#'       pairs were left out, a total of 0 does not show reversibility, since
#'       such pairs mean the chain is not reversible.}
#'   }
#' @param scaling How the weights are expressed. `"none"` (default) keeps
#'   them in `base` units. `"share"` (contribution only) gives each edge as a
#'   \emph{percentage} of the entropy rate: the weights sum to 100, each one
#'   reads "this transition holds x% of the process's unpredictability", and
#'   they compare across networks whatever the number of states or the base.
#'   `"chance"` (surprisal only) divides by \eqn{\log_b n}, the surprisal of a
#'   transition when all \eqn{n} next states are equally likely: `1` is chance
#'   level, below 1 an expected transition, above 1 one rarer than chance.
#'   `"share"` with another weight, or `"chance"` with another weight, raises
#'   an error.
#' @param normalize Whether rows that do not sum to 1 (a frequency network, a
#'   count matrix) are divided by their sums, with a warning. With `FALSE`,
#'   such rows raise an error. Default `TRUE`. A state with no outgoing
#'   transitions (a row summing to zero) is an error either way.
#'
#' @return A `netobject` (class `c("netobject", "cograph_network")`; for a
#'   bare `cograph_network` input, that network with its own class) with
#'   method `"entropy"`, whose weights are the per-edge entropy quantities.
#'   Print it for an overview; `summary()` (see [build_network()]) gives its
#'   network metrics,
#'   [extract_edges()] one row per edge, and `cograph::splot()` draws it.
#'   When `x` is a network (or sequence data, from which a relative network
#'   is built), the result is that network with the entropy weights swapped
#'   in: its initial probabilities, metadata, node groups and node
#'   coordinates are kept, so it plots with the same TNA styling and layout
#'   as its source. The stored parameters are `base`, `weight`, `scaling`,
#'   the entropy rate `entropy_rate` (whatever the weight) and the stationary
#'   distribution `stationary`, plus `production_rate` and `n_oneway_pairs`
#'   when `weight = "production"`.
#'
#'   The object also declares the entropy house style through its
#'   `meta$splot` producer contract (honoured by cograph >= 2.4.4): no
#'   minimum-weight pruning (bit values are smaller than probabilities), edge
#'   labels with 2 digits (1 for `"share"`), vermilion edges, node rings
#'   showing the stationary distribution, and TNA rather than psychometric
#'   geometry. The production weights are signed, so they keep cograph's
#'   positive and negative edge colours instead of the vermilion. With
#'   `weight = "surprisal"` and `scaling = "chance"`, the edges are drawn and
#'   sized by the transition probabilities, transitions likelier than chance
#'   in blues (darker for likelier) and those at or below chance dashed grey.
#'   Because the contract states the styling outright, cograph needs no
#'   knowledge of the `"entropy"` method: `cograph::splot(ent)` renders it
#'   with no arguments, and any argument given overrides the contract.
#'
#' @details
#' An impossible transition (\eqn{P_{ij} = 0}) gets weight `0` under every
#' definition (the \eqn{0 \log 0 = 0} convention), so the entropy network has
#' no edge the transition network lacks. Self-loops are kept like any other
#' edge; their production weight is always 0. A deterministic transition
#' (\eqn{P_{ij} = 1}) also gets contribution and surprisal `0`: observing the
#' inevitable carries no information.
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
  .assert(is.numeric(base), length(base) == 1L, base > 0, base != 1)
  if (scaling == "chance" && weight != "surprisal") {
    stop(errorCondition(
      "scaling = 'chance' applies to weight = 'surprisal' only.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (scaling == "share" && weight != "contribution") {
    stop(errorCondition(
      "scaling = 'share' applies to weight = 'contribution' only.",
      class = "nestimate_bad_input", call = NULL))
  }

  ## Resolve a source netobject so the result inherits its styling metadata
  ## (inits, meta, node_groups, coordinates) and renders identically to the
  ## transition network it decomposes.
  src <- NULL
  if (inherits(x, "netobject")) {
    src <- x
  } else if (inherits(x, "cograph_network")) {
    src <- .as_netobject(x)
  } else if (is.data.frame(x) ||
             (is.matrix(x) && (is.character(x) || is.logical(x)))) {
    if (is.data.frame(x)) {
      .refuse_undeclared_roles(x, "Reading sequences", via_network = TRUE,
                               reads_onehot = TRUE)
    }
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
  if (any(pi <= 0)) .te_warn_not_ergodic()

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
    if (length(state_names) < 2L) {
      stop(errorCondition(paste0(
        "scaling = 'chance' needs at least 2 states: with one, chance is ",
        "log(1) = 0."),
        class = "nestimate_bad_input", call = NULL))
    }
    W <- W / log(length(state_names), base = base)
  }
  if (scaling == "share") {
    if (!(sum(W) > 0)) {
      stop(errorCondition(paste0(
        "scaling = 'share' needs a positive entropy rate; this chain is ",
        "deterministic (entropy rate 0), so there is no total to share."),
        class = "nestimate_bad_input", call = NULL))
    }
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
#' `entropy_trajectory()` tracks transition entropy over time: it forms the
#' transitions within each actor's sequence of an event log, pools them in
#' time order, and slides a window of a fixed number of transitions across
#' the stream, giving one entropy estimate per window. A `group` column gives
#' one trajectory per group. This turns transition entropy from a snapshot
#' into a process measure, as Krejtz et al. (2025) do when they track gaze
#' transition entropy through the phases of a task.
#'
#' A falling trajectory means the transitions become more predictable
#' (routinisation), a rising one less predictable (exploration), and a level
#' shift marks a change of phase. Each window's estimate rests on `window`
#' transitions, and with the default `step` neighbouring windows share four
#' fifths of them, so their estimates are not independent: read the smoothed trend of `plot()` rather
#' than single windows. For the whole-process value with credible intervals,
#' use [transition_entropy()] and [entropy_bayes()].
#'
#' @param data An event log: a data frame with one row per event.
#' @param action The column holding each event's state.
#' @param actor The column identifying sequences: transitions are formed only
#'   between consecutive events of the same actor (within the same group).
#'   Default `NULL`, which reads the data as one sequence.
#' @param time The timestamp column, which orders the events within each
#'   actor and places the windows on a time axis. Default `NULL`, which keeps
#'   the row order and uses the row number of the event that ends each
#'   transition as the axis. The column must be numeric, a `Date` or a
#'   `POSIXct`; a character timestamp fails, so parse it first.
#' @param group A column splitting the data into parallel trajectories (for
#'   example a condition or an achievement level), each with its own windows.
#'   Default `NULL` (one trajectory, labelled `"all"`).
#' @param window The number of transitions per window, at least 10. Default
#'   `500`.
#' @param step The number of transitions between the starts of consecutive
#'   windows. Default `NULL`, which uses `window %/% 5` (at least 1).
#' @param base The logarithm base, which sets the unit: `2` (default) for
#'   bits. Any positive number other than 1.
#'
#' @return A `net_entropy_trajectory` object. `print()` shows each group's
#'   number of windows and its mean and range of entropy, and returns its
#'   input invisibly. `summary()` returns a data frame with one row per
#'   group: `group`, `windows`, the `mean`, `sd`, `min` and `max` of the
#'   window entropies, the entropy of the `first` and `last` windows, and
#'   their difference `change` (negative when the transitions became more
#'   predictable). `plot()` returns a `ggplot` with the per-window entropy as
#'   faint lines and a loess-smoothed trend (with its standard-error band) per
#'   group, in Okabe-Ito colours. The object's fields:
#' \describe{
#'   \item{trajectory}{Data frame, one row per window: `group`, `window`
#'     (index within the group), `time` (the window's midpoint; a row number
#'     when there is no `time` column), `time_start` and `time_end` (the
#'     time, or row number, of the window's first and last transition),
#'     `n_transitions`, `n_states` (distinct states in the window), `entropy`
#'     (in `base` units per transition) and `entropy_norm` (divided by
#'     \eqn{\log_b} of the window's number of states, floored at 2 states so
#'     the ceiling is never zero; in \eqn{[0, 1]}).}
#'   \item{window, step, base}{The settings used.}
#'   \item{states}{The states of the whole log.}
#' }
#'
#' @details
#' Each window's entropy is the empirical conditional entropy
#' \eqn{-\sum_{ij} (n_{ij}/N) \log_b(n_{ij}/n_{i\cdot})}, which weights rows
#' by their observed occupancy rather than by the stationary distribution
#' from the eigenvector. Within a short window the chain is often not ergodic
#' (absorbing fragments, unvisited states), where that eigenvector is
#' undefined or misleading; the empirical estimator is the standard windowed
#' choice and converges to the stationary entropy rate over long stationary
#' stretches.
#'
#' The transitions of all actors are pooled in time order within each group
#' (ties broken by actor, then by input order). A window shorter than
#' `window` at the end of a stream is dropped; a group whose stream is
#' shorter than `window` gets one window covering all of it, with a warning.
#' Fewer than two transitions in all raise an error.
#'
#' @examples
#' tr <- entropy_trajectory(group_regulation_long,
#'                          action = "Action", actor = "Actor",
#'                          time = "Time", group = "Achiever")
#' tr
#' summary(tr)
#' plot(tr)
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
  .assert(
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
  if (!is.null(time) && !(is.numeric(data[[time]]) ||
                           inherits(data[[time]], c("Date", "POSIXt")))) {
    stop(errorCondition(sprintf(paste0(
      "`time` column '%s' must be numeric, a Date or a POSIXct, not %s; ",
      "parse it first (for example with prepare())."),
      time, class(data[[time]])[1L]),
      class = "nestimate_bad_input", call = NULL))
  }
  window <- as.integer(window)
  if (is.null(step)) step <- max(1L, window %/% 5L)
  step <- as.integer(step)
  .assert(step >= 1L)

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
  tr_a   <- act[-1][keep]

  ## Transitions were built actor by actor; the trajectory runs through
  ## time, so pool them in time order within each group (ties: actor, then
  ## input order).
  by_time <- order(tr_g, tr_t, tr_a, seq_along(tr_t))
  from_i <- from_i[by_time]
  to_i   <- to_i[by_time]
  tr_t   <- tr_t[by_time]
  tr_g   <- tr_g[by_time]

  if (length(from_i) < 2L) {
    stop(errorCondition(
      "Fewer than two transitions; cannot build a trajectory.",
      class = "nestimate_bad_input", call = NULL))
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
      warning(warningCondition(sprintf(
        "Group stream (%d transitions) shorter than window (%d); using one window.",
        n_tr, window), class = "nestimate_short_stream", call = NULL))
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
      states     = states,
      time_col   = time
    ),
    class = "net_entropy_trajectory"
  )
}

#' @rdname entropy_trajectory
#' @param x For the `print()` and `plot()` methods: an
#'   `entropy_trajectory()` result.
#' @param digits For `print()`: decimal places of the entropies. Default
#'   `3`.
#' @param ... For the `print()`, `summary()` and `plot()` methods: not used.
#' @param object For `summary()`: an `entropy_trajectory()` result.
#' @param normalized For `plot()`: whether to draw `entropy_norm` instead of
#'   the raw entropy. Default `FALSE`.
#' @param span For `plot()`: the loess span of the trend line. Default
#'   `0.4`.
#' @param title For `plot()`: the plot title. Default
#'   `"Transition entropy over time"`.
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
              digits, agg[2, ], digits, agg[3, ], digits, agg[4, ]),
      sep = "")
  cat("\nUse plot() for the trajectory, summary() for per-group statistics.\n")
  invisible(x)
}

#' @rdname entropy_trajectory
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

#' @rdname entropy_trajectory
#' @export
plot.net_entropy_trajectory <- function(x, normalized = FALSE, span = 0.4,
                                        title = "Transition entropy over time",
                                        ...) {
  tr <- x$trajectory
  tr$y <- if (normalized) tr$entropy_norm else tr$entropy
  ylab <- if (normalized) "Normalised entropy" else
    sprintf("Entropy (%s per transition)", .te_unit(x$base))
  xlab <- if (inherits(tr$time, "POSIXct")) "Time" else
    if (inherits(tr$time, "Date")) "Date" else
      if (is.null(x$time_col)) "Row of the event ending the transition" else
        x$time_col

  n_groups <- length(unique(tr$group))
  single <- n_groups == 1L
  # More groups than palette colours: one panel per group, all in the
  # single-group colour, instead of a colour scale that cannot cover them.
  faceted <- n_groups > length(.okabe_ito)
  tr$colour_key <- if (faceted) tr$group[1L] else tr$group
  p <- ggplot2::ggplot(tr, ggplot2::aes(x = .data$time, y = .data$y,
                                        colour = .data$colour_key)) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.35) +
    ggplot2::geom_smooth(method = "loess", span = span, se = TRUE,
                         linewidth = 0.9, alpha = 0.15) +
    ggplot2::scale_colour_manual(values = .okabe_ito, name = NULL) +
    ggplot2::labs(x = xlab,
      y = ylab, title = title,
      subtitle = sprintf("Sliding window: %d transitions, step %d",
                         x$window, x$step)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = if (single || faceted) "none" else
                     "bottom")
  if (faceted) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$group))
  p
}

# ---- entropy_bayes() ----

#' @noRd
.te_extract_counts <- function(x) {
  if (is.matrix(x) && is.numeric(x)) {
    if (any(x < 0) || any(x != floor(x))) {
      stop(errorCondition(
        paste0(
          "Bayesian entropy estimation needs transition COUNTS ",
          "(non-negative whole numbers), not probabilities. Supply a ",
          "count matrix, a frequency network ",
          "(build_network(method = 'frequency')), or sequence data."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (is.null(rownames(x))) {
      dimnames(x) <- list(paste0("S", seq_len(nrow(x))),
                          paste0("S", seq_len(nrow(x))))
    }
    return(x)
  }
  if (inherits(x, "tna")) return(.te_tna_counts(x))
  if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    return(.te_network_counts(x))
  }
  if (is.data.frame(x) ||
      (is.matrix(x) && (is.character(x) || is.logical(x)))) {
    if (is.data.frame(x)) {
      .refuse_undeclared_roles(x, "Reading sequences", via_network = TRUE,
                               reads_onehot = TRUE)
    }
    return(build_network(as.data.frame(x, stringsAsFactors = FALSE),
                         method = "frequency")$weights)
  }
  stop(errorCondition(
    paste0(
      "'x' must be a count matrix, a netobject / cograph_network / tna ",
      "object, or a sequence data.frame."),
    class = "nestimate_bad_input", call = NULL))
}

# Counts of a network, decided by its declared method, never by its values:
# whole-number weights can be probabilities (a deterministic chain) or scaled
# counts. A sequence network with its data is recounted by the estimator's own
# per-sequence counting, which honours start/end states, concatenation, windows
# and one-hot units. Weights are read as counts only from an unscaled count
# network without data, or a network with no declared method.
#' @noRd
.te_network_counts <- function(x) {
  w <- x$weights
  method <- x$method %||% x$meta$tna$method
  if (isTRUE(method %in% .association_methods)) {
    stop(errorCondition(sprintf(paste(
      "This \"%s\" network holds associations between variables, not",
      "transitions. Pass sequence data or a frequency network."), method),
      class = "nestimate_bad_input", call = NULL))
  }
  declared <- is.character(method) && length(method) == 1L
  if (!declared && identical(x$meta$source, "tna")) {
    # a cograph conversion of a tna model: its sequences count only when they
    # reproduce its weights (windows and settings are not kept)
    counts <- .tna_conversion_counts(x)
    if (is.null(counts)) {
      stop(errorCondition(paste(
        "This network was converted from a tna model, and its sequences do not",
        "reproduce its weights (windows, scaling or other settings were lost).",
        "Pass the tna model itself or the sequences."),
        class = "nestimate_bad_input", call = NULL))
    }
    return(counts)
  }
  params <- x$params
  count_method <- if (declared) {
    .resampling_method(.resolve_method_alias(method), params)
  }
  states <- rownames(w)
  if (isTRUE(count_method %in% c("relative", "frequency", "co_occurrence")) &&
      is.data.frame(x$data) && nrow(x$data) > 0L && !is.null(states)) {
    # raw counts: neither row normalisation nor 1 / length weighting
    params$weighted <- NULL
    counting <- if (count_method == "co_occurrence") "co_occurrence" else
      "frequency"
    data <- .resampling_transition_data(x$data, x$metadata, params)
    per <- .precompute_per_sequence(data, counting, params, states)
    return(matrix(colSums(per), nrow = length(states), byrow = TRUE,
                  dimnames = list(states, states)))
  }
  unscaled <- is.null(x$scaling) || identical(x$scaling, "none")
  counts_held <- !declared ||
    (isTRUE(count_method %in% c("frequency", "co_occurrence", "reverse")) &&
       !.counts_are_relative(count_method, params) &&
       !isTRUE(.param_get(params, "weighted", FALSE)))
  if (counts_held && unscaled && is.matrix(w) && is.numeric(w) &&
      all(w >= 0) && all(w == floor(w)) && sum(w) > 0) {
    return(w)
  }
  stop(errorCondition(paste0(
    "Cannot recover transition counts from this ",
    if (declared) sprintf("\"%s\" ", method) else "",
    "network: its weights are not counts, and only a relative or frequency ",
    "network is recounted from its sequences. Supply a frequency network (build_network(method = 'frequency')), ",
    "a count matrix, or the sequence data."),
    class = "nestimate_bad_input", call = NULL))
}

# A tna model declares its type. An unscaled model of a counting type holds
# its counts. A relative or scaled model with no other settings is recounted
# from its sequences, decoded through its labels; one with windows or other
# settings cannot be replayed here and is refused, as is attention or gap.
#' @noRd
.te_tna_counts <- function(x) {
  w <- x$weights
  type <- attr(x, "type") %||% "relative"
  scaling <- attr(x, "scaling")
  unscaled <- length(scaling) == 0L || identical(scaling, "none")
  params <- attr(x, "params")
  # weighted models divide each sequence's counts by its length
  weighted <- isTRUE(params$weighted)
  if (unscaled && !weighted &&
      type %in% c("frequency", "co-occurrence", "reverse", "n-gram") &&
      is.matrix(w) && all(w >= 0) && all(w == floor(w)) && sum(w) > 0) {
    return(w)
  }
  # weighting changes the estimate, not the sequences, so a weighted model is
  # recounted from them like a relative one
  plain <- !isTRUE(params$windowed) &&
    length(setdiff(names(params), c("windowed", "weighted"))) == 0L
  labels <- x$labels %||% rownames(w)
  if (!type %in% c("relative", "frequency") || !plain || is.null(x$data) ||
      is.null(labels)) {
    stop(errorCondition(sprintf(paste(
      "Cannot recover transition counts from this tna model (type \"%s\"%s):",
      "its weights are not counts and its settings cannot be replayed. Supply",
      "an unscaled frequency model (tna::ftna()) or the sequence data."),
      type, if (plain) "" else ", windowed or with settings"),
      class = "nestimate_bad_input", call = NULL))
  }
  data <- unclass(x$data)
  decoded <- matrix(labels[data], nrow = nrow(data),
                    dimnames = list(NULL, colnames(data) %||%
                                      paste0("T", seq_len(ncol(data)))))
  decoded <- as.data.frame(decoded, stringsAsFactors = FALSE)
  cw <- build_network(decoded, method = "frequency",
                      state_cols = names(decoded))$weights
  missing <- setdiff(labels, rownames(cw))
  if (length(missing) > 0L) {
    out <- matrix(0, length(labels), length(labels),
                  dimnames = list(labels, labels))
    out[rownames(cw), colnames(cw)] <- cw
    cw <- out
  }
  cw[labels, labels, drop = FALSE]
}

#' Bayesian Transition Entropy
#'
#' @description
#' `entropy_bayes()` puts credible intervals on the entropy quantities of
#' [transition_entropy()] and on the edge-level decomposition of
#' [entropy_network()]. Each row of the transition count matrix gets an
#' independent Dirichlet posterior (counts + `prior`), and Monte Carlo draws
#' from it carry the uncertainty of the counts into the entropy rate, the
#' stationary entropy, the redundancy, each state's branching entropy and
#' stationary probability, and every edge's entropy contribution. Each is
#' reported as a posterior mean with an equal-tailed credible interval. The
#' counts are recovered from a network, a `tna` model, a count matrix or
#' sequence data (see `x`); a `netobject_group` gives one result per group.
#'
#' The practical purpose is to \strong{exclude unstable estimates}: an edge
#' whose contribution rests on a handful of observations has a wide
#' posterior, and it is flagged as not credible unless it credibly accounts
#' for at least `min_share` of the entropy rate. The pruned network keeps the
#' credible edges only, the stable entropy skeleton. The posterior treats
#' every transition as an independent draw from its row of a first-order
#' chain; whether a first-order chain describes the sequences is a separate
#' question, which [markov_order_test()] tests.
#'
#' @param x The transition counts: a frequency `netobject`
#'   (`build_network(method = "frequency")`); a relative, frequency or
#'   co-occurrence `netobject` that carries its data, recounted with the
#'   network's own start and end states, concatenation, windows and one-hot
#'   units; a one-hot network built with `build_network(method = "wtna")` or
#'   `"wtna_cooccurrence"`, with its data, recounted the same way (a [wtna()]
#'   result, whose method is `"wtna_transition"`, currently raises
#'   `nestimate_bad_input`); an unscaled frequency, co-occurrence or reverse
#'   network without data, whose weights are the counts; a network without a
#'   declared method whose weights are non-negative whole numbers, read as
#'   counts; a `tna` model; a count matrix (non-negative whole numbers); or
#'   wide sequence data (a data frame or a character matrix) that holds a
#'   numbered run of positions (T1, T2, ...), from which a frequency network
#'   is built, other columns set aside with a message. For a frame without
#'   such a run, pass a network built with [build_network()]. When a
#'   network declares its method, counts are read from that method, never
#'   from the values: a scaled or relative network is recounted from its
#'   data. An association network, or any other network
#'   whose counts cannot be recovered, raises an error of class
#'   `nestimate_bad_input`; a matrix of probabilities raises an error. A
#'   `netobject_group` gives one result per group. For the `print()` and
#'   `plot()` methods: an `entropy_bayes()` result; `print()` also takes a
#'   group result.
#' @param prior The Dirichlet concentration added to every cell of the count
#'   matrix, a positive number. Default `0.5` (Jeffreys).
#' @param draws The number of Monte Carlo posterior draws, at least 100.
#'   Default `4000`.
#' @param ci The mass of the credible intervals, between 0 and 1. Default
#'   `0.95`.
#' @param min_share The credibility threshold, a proportion in \eqn{[0, 1)}:
#'   an edge is credible when the lower limit of the credible interval of its
#'   \emph{share} of the entropy rate exceeds it. Default `0.01`: the edge
#'   credibly accounts for at least 1% of the entropy rate.
#' @param base The logarithm base, which sets the unit: `2` (default) for
#'   bits. Any positive number other than 1.
#' @param seed An integer that makes the draws reproducible. Default `NULL`.
#'   The session's own random stream is left as it was; a group uses the
#'   same seed for every
#'   group.
#'
#' @return A `net_entropy_bayes` object, or for a `netobject_group` a
#'   `net_entropy_bayes_group`: a named list of them, one per group.
#'   `print()` shows the posterior mean and credible interval of the entropy
#'   rate, stationary entropy and redundancy, and how many edges are
#'   credible; for a group, each group's entropy rate. For a group,
#'   `summary()` returns the chain-level table of every group stacked, with a
#'   `group` column, and `plot()` draws the forest plot below in one panel
#'   per group. For one network, `summary()` prints the
#'   chain-level table and the per-edge table and returns the per-edge data
#'   frame invisibly. `plot()` returns a `ggplot` forest plot of the edges
#'   with the largest posterior mean contributions: the mean and credible
#'   interval of each, credible edges in Okabe-Ito blue and unstable ones in
#'   grey, with a dashed line at `min_share` of the posterior mean entropy
#'   rate, the stability criterion. The `print()` methods return their input
#'   invisibly. The object's fields:
#' \describe{
#'   \item{summary}{Data frame, one row per chain-level quantity
#'     (`entropy_rate`, `stationary_entropy`, `redundancy`), with the
#'     posterior `mean`, `sd`, `ci_lower` and `ci_upper`.}
#'   \item{states}{Data frame, one row per state: the posterior mean, sd and
#'     credible interval of its branching entropy, and the posterior mean and
#'     credible interval of its stationary probability.}
#'   \item{edges}{Data frame, one row per observed transition (`from`, `to`,
#'     `count`), sorted by contribution, largest first: the posterior mean,
#'     sd and credible interval of its contribution (in `base` units), the
#'     posterior mean and credible interval of its share of the entropy rate
#'     (a proportion), and the `credible` flag.}
#'   \item{network}{A `netobject`: the posterior-mean entropy network, every
#'     observed edge included, with the entropy edge style (vermilion edges,
#'     2-digit labels, no minimum-weight pruning, node rings showing the
#'     posterior mean stationary distribution).}
#'   \item{model}{A `netobject`: the pruned entropy network, posterior means
#'     where `credible` and `0` elsewhere, styled the same way.}
#'   \item{draws_entropy_rate}{The posterior draws of the entropy rate, for
#'     further analysis or plotting.}
#'   \item{prior, draws, ci, min_share, base}{The settings used.}
#'   \item{states_names}{Character vector of state names.}
#' }
#'
#' @details
#' With `prior > 0` the posterior puts mass on every transition, so the
#' contribution draws are strictly positive and a rule such as "the interval
#' excludes zero" would flag every edge as credible. The share criterion is
#' used instead: an edge is stable when it credibly carries at least
#' `min_share` of \eqn{h(P)}. An unobserved transition (count 0) gets prior
#' mass only; it is left out of the edge table and of both networks.
#'
#' The posterior mean entropy rate differs from the plug-in estimate of
#' [transition_entropy()] in two opposite ways: the prior pulls every row
#' toward uniform, which raises the entropy, and averaging over the posterior
#' lowers it. Which effect wins depends on the data (on the first 200 events of
#' `group_regulation_long` the posterior mean is above the plug-in value); the
#' difference is largest on sparse counts and vanishes as the counts grow.
#'
#' @examples
#' net <- build_network(group_regulation_long, method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' eb <- entropy_bayes(net, draws = 1000, seed = 1)
#' eb
#' summary(eb)
#' plot(eb)
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

  .assert(
    is.numeric(prior), length(prior) == 1L, prior > 0,
    is.numeric(draws), length(draws) == 1L, draws >= 100,
    is.numeric(ci), length(ci) == 1L, ci > 0, ci < 1,
    is.numeric(min_share), length(min_share) == 1L,
    min_share >= 0, min_share < 1,
    is.numeric(base), length(base) == 1L, base > 0, base != 1
  )
  if (!is.null(seed)) .local_seed(seed)

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
                donut_fill = as.numeric(pi_mean), donut_empty = FALSE,
                tna_styling = TRUE, psych_styling = FALSE)
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

#' @rdname entropy_bayes
#' @param digits For `print()` of a result: decimal places of the
#'   chain-level estimates. Default `3`. A group result prints with three.
#' @param ... For the `print()`, `summary()` and `plot()` methods: not used.
#' @param object For `summary()`: an `entropy_bayes()` result.
#' @param top For `plot()`: the most edges to show, those with the largest
#'   posterior mean contribution. Default `25`.
#' @param title For `plot()`: the plot title. Default
#'   `"Bayesian edge entropy contributions"`.
#' @export
print.net_entropy_bayes <- function(x, digits = 3, ...) {
  unit <- .te_unit(x$base)

  cat(sprintf(
    "Bayesian Transition Entropy (%d states, %s; Dirichlet prior %.2g, %d draws)\n\n",
    length(x$states_names), unit, x$prior, x$draws))
  s <- x$summary
  cat(sprintf("  %-19s %.*f [%.*f, %.*f]\n",
              paste0(s$quantity, ":"),
              digits, s$mean, digits, s$ci_lower, digits, s$ci_upper),
      sep = "")
  n_cred <- sum(x$edges$credible)
  cat(sprintf(
    "\nEdges: %d observed; %d credible (share of h(P) credibly > %g%%).\n",
    nrow(x$edges), n_cred, 100 * x$min_share))
  cat("Use summary() for the edge table and plot() for the posterior.\n")
  invisible(x)
}

#' @rdname entropy_bayes
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
  cat(sprintf("  %-12s %s\n", paste0(names(hs), ":"), hs), sep = "")
  invisible(x)
}

#' @rdname entropy_bayes
#' @export
summary.net_entropy_bayes_group <- function(object, ...) {
  rows <- Map(function(o, g) {
    cbind(group = g, o$summary, stringsAsFactors = FALSE)
  }, object, names(object))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname entropy_bayes
#' @export
plot.net_entropy_bayes_group <- function(x, top = 25,
                                         title = paste(
                                           "Bayesian edge entropy",
                                           "contributions"),
                                         ...) {
  first <- x[[1L]]
  df <- do.call(rbind, Map(function(o, g) {
    e <- utils::head(o$edges, top)
    data.frame(group = g, edge = paste(e$from, "->", e$to),
               contribution = e$contribution, ci_lower = e$ci_lower,
               ci_upper = e$ci_upper, credible = e$credible,
               stringsAsFactors = FALSE)
  }, x, names(x)))
  # one y order per panel: prefix each edge with its group, strip it on the
  # axis (assigning duplicate levels instead would merge the panels' orders)
  key <- paste0(df$group, "\r", df$edge)
  df$row <- factor(key, levels = rev(unique(key)))
  thresholds <- do.call(rbind, Map(function(o, g) {
    h <- o$summary$mean[o$summary$quantity == "entropy_rate"]
    data.frame(group = g, cut = o$min_share * h, stringsAsFactors = FALSE)
  }, x, names(x)))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$contribution, y = .data$row,
                                   colour = .data$credible,
                                   shape = .data$credible)) +
    ggplot2::geom_vline(data = thresholds,
                        ggplot2::aes(xintercept = .data$cut),
                        linetype = "dashed", colour = "grey40",
                        linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_lower,
                                        xmax = .data$ci_upper),
                           orientation = "y", width = 0.25, linewidth = 0.5) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_y_discrete(labels = function(l) sub("^.*\r", "", l)) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#0072B2", "FALSE" = "grey60"),
      labels = c("TRUE" = "credible", "FALSE" = "unstable"), name = NULL) +
    ggplot2::scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "credible", "FALSE" = "unstable"), name = NULL) +
    ggplot2::facet_wrap(ggplot2::vars(.data$group), scales = "free_y") +
    ggplot2::labs(
      x = sprintf("Contribution to h(P) (%s)", .te_unit(first$base)),
      y = NULL, title = title,
      subtitle = sprintf(
        "Posterior mean and %g%% CI; dashed line = %g%% of h(P)",
        100 * first$ci, 100 * first$min_share)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

#' @rdname entropy_bayes
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

#' @rdname entropy_bayes
#' @export
plot.net_entropy_bayes <- function(x, top = 25,
                                   title = "Bayesian edge entropy contributions",
                                   ...) {
  unit <- .te_unit(x$base)

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

#' @rdname transition_entropy
#' @param digits For the `print()` methods of a result and of its summary:
#'   decimal places of the values. Default `3`. A result's normalised values
#'   always print with three decimals, the summary's percentages with one.
#' @param ... For the `print()`, `summary()` and `plot()` methods: not used.
#' @param object For `summary()`: a `transition_entropy()` result.
#' @param title For `plot()`: the plot title. Default `"Transition Entropy"`.
#' @param fill For `plot()`: the bar fill colour. Default `"#0072B2"`, the
#'   Okabe-Ito blue.
#' @export
print.net_transition_entropy <- function(x, digits = 3, ...) {
  unit <- .te_unit(x$base)

  cat(sprintf("Transition Entropy (%d states, %s; ceiling = %.3f)\n\n",
              length(x$states), unit, x$max_entropy))
  cat(sprintf("                          raw            normalised\n"))
  cat(sprintf("  Entropy rate    h(P)  = %.*f %-6s  %.3f\n",
              digits, round(x$entropy_rate, digits) + 0, unit,
              round(x$entropy_rate_norm, 3) + 0))
  cat(sprintf("  Stationary    H(pi)  = %.*f %-6s  %.3f\n",
              digits, round(x$stationary_entropy, digits) + 0, unit,
              round(x$stationary_entropy_norm, 3) + 0))
  cat(sprintf("  Redundancy   H(pi)-h = %.*f %-6s  %.3f\n",
              digits, round(x$redundancy, digits) + 0, unit,
              round(x$redundancy_norm, 3) + 0))
  cat(sprintf(paste0(
    "\nNormalised: h(P) and H(pi) are raw / log_%g(n_states) ",
    "(0 = deterministic, 1 = uniform);\n",
    "  redundancy is the relative redundancy (H(pi) - h(P)) / H(pi), ",
    "not raw / log_%g(n_states).\n"),
    x$base, x$base))
  invisible(x)
}

#' @rdname transition_entropy
#' @export
print.net_transition_entropy_group <- function(x, ...) {
  cat(sprintf("Transition Entropy - %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  rates <- vapply(x, function(o) o$entropy_rate, numeric(1))
  cat("Entropy rate per group:\n")
  print(round(rates, 4))
  invisible(x)
}

#' @rdname transition_entropy
#' @export
summary.net_transition_entropy_group <- function(object, ...) {
  rows <- Map(function(o, g) {
    data.frame(
      group                   = g,
      n_states                = length(o$states),
      entropy_rate            = o$entropy_rate,
      entropy_rate_norm       = o$entropy_rate_norm,
      stationary_entropy      = o$stationary_entropy,
      stationary_entropy_norm = o$stationary_entropy_norm,
      redundancy              = o$redundancy,
      unit                    = .te_unit(o$base),
      stringsAsFactors = FALSE
    )
  }, object, names(object))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname transition_entropy
#' @export
plot.net_transition_entropy_group <- function(x,
                                              title = "Transition Entropy",
                                              fill  = "#0072B2",
                                              ...) {
  base <- x[[1L]]$base
  bars <- do.call(rbind, Map(function(o, g) {
    data.frame(group = g, state = o$states, row_entropy = o$row_entropy,
               width = o$stationary / max(o$stationary),
               stringsAsFactors = FALSE)
  }, x, names(x)))
  bars$state <- factor(bars$state, levels = unique(bars$state))
  lines <- do.call(rbind, Map(function(o, g) {
    data.frame(group = g, entropy_rate = o$entropy_rate,
               ceiling = log(length(o$states), base = o$base),
               stringsAsFactors = FALSE)
  }, x, names(x)))
  lines$rate_label <- sprintf("h(P) = %.3f", lines$entropy_rate)

  ggplot2::ggplot(bars, ggplot2::aes(x = .data$state, y = .data$row_entropy,
                                     width = .data$width)) +
    ggplot2::geom_col(fill = fill, alpha = 0.85) +
    ggplot2::geom_hline(data = lines,
                        ggplot2::aes(yintercept = .data$entropy_rate),
                        linetype = "dashed", colour = "#D55E00",
                        linewidth = 0.6) +
    ggplot2::geom_hline(data = lines,
                        ggplot2::aes(yintercept = .data$ceiling),
                        linetype = "dotted", colour = "grey30",
                        linewidth = 0.4) +
    ggplot2::geom_text(data = lines,
                       ggplot2::aes(x = -Inf, y = .data$entropy_rate,
                                    label = .data$rate_label),
                       inherit.aes = FALSE, hjust = -0.05, vjust = -0.4,
                       colour = "#D55E00", size = 3.4) +
    ggplot2::facet_wrap(ggplot2::vars(.data$group)) +
    ggplot2::labs(x = "State",
                  y = sprintf("Row entropy (%s)", .te_unit(base)),
                  title = title,
                  subtitle = paste0(
                    "Bar width proportional to stationary probability; ",
                    "dotted line = ceiling log(n)")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title         = ggplot2::element_text(face = "bold")
    )
}

#' @rdname transition_entropy
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

  unit <- .te_unit(object$base)

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

#' @rdname transition_entropy
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

#' @rdname transition_entropy
#' @export
plot.net_transition_entropy <- function(x,
                                        title = "Transition Entropy",
                                        fill  = "#0072B2",
                                        ...) {
  unit <- .te_unit(x$base)

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
