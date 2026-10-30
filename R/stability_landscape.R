# ---- Stability landscapes (ball-in-cup) of a transition network ----

utils::globalVariables(c("u", "y", "ymin", "ymax", "label", "state", "group", "value", "depth", "h",
                         "label_y"))

#' @noRd
.ls_bad <- function(msg, class = "nestimate_bad_input") {
  stop(errorCondition(msg, class = class, call = NULL))
}

# Per-sequence transition counts: one row per sequence, k * k columns, where
# column f * k + t + 1 counts from-state f to to-state t (0-based codes). A
# pair counts when both positions are observed, as in build_network().
#' @noRd
.ls_counts <- function(states, alphabet) {
  k <- length(alphabet)
  codes <- matrix(match(states, alphabet) - 1L, nrow = nrow(states))
  if (ncol(codes) < 2L) return(matrix(0L, nrow(codes), k * k))
  cell <- codes[, -ncol(codes), drop = FALSE] * k +
    codes[, -1L, drop = FALSE] + 1L
  t(apply(cell, 1L, \(v) tabulate(v[!is.na(v)], nbins = k * k)))
}

# Cup depth of every state from a pooled count matrix: log of the expected
# stay, log(1 / (1 - P_ii)). NA for a state with no outgoing transition,
# Inf for one that is never left.
#' @noRd
.ls_depth <- function(n) {
  out_n <- rowSums(n)
  stay <- ifelse(out_n > 0, diag(n) / out_n, NA_real_)
  -log1p(-stay)
}

#' @noRd
.ls_irreducible <- function(P) {
  # repeated squaring of the one-step reachability matrix
  reach <- Reduce(\(r, i) (r %*% r) > 0,
                  seq_len(ceiling(log2(nrow(P))) + 1L),
                  (P > 0) | diag(nrow(P)) > 0)
  all(reach)
}

# Net circulation of a three-state chain: the net flux pi_1 P_12 - pi_2 P_21,
# which stationarity makes equal on all three edges of the cycle. NA when the
# chain has no unique stationary distribution.
#' @noRd
.ls_current <- function(n) {
  out_n <- rowSums(n)
  if (any(out_n == 0)) return(NA_real_)
  P <- n / out_n
  if (!.ls_irreducible(P)) return(NA_real_)
  pi <- .mpt_stationary(P)
  pi[1L] * P[1L, 2L] - pi[2L] * P[2L, 1L]
}

#' @noRd
.ls_sequences <- function(x) {
  if (inherits(x, c("netobject", "cograph_network")) &&
      !isTRUE(x$method %in% c("relative", "frequency"))) {
    .mpt_not_a_chain(sprintf("A \"%s\" network", x$method %||% "unknown"))
  }
  .as_sequences(x)
}

#' @noRd
.ls_states <- function(states, alphabet) {
  if (is.null(states)) return(alphabet)
  if (!is.character(states) || anyDuplicated(states) ||
      !setequal(states, alphabet)) {
    .ls_bad(sprintf(
      "`states` must name every state exactly once: %s.",
      paste(alphabet, collapse = ", ")))
  }
  states
}

# Percentile interval that tolerates Inf (a resample in which the state is
# never left): type 1 returns order statistics, so no Inf - Inf arithmetic.
#' @noRd
.ls_interval <- function(v, ci_level) {
  v <- v[!is.na(v)]
  if (length(v) < 2L) return(c(NA_real_, NA_real_))
  alpha <- (1 - ci_level) / 2
  stats::quantile(v, c(alpha, 1 - alpha), names = FALSE,
                  type = if (any(is.infinite(v))) 1L else 7L)
}

#' @noRd
.ls_fit <- function(states_mat, states, iter, ci_level, name) {
  k <- length(states)
  counts <- .ls_counts(states_mat, states)
  as_matrix <- \(v) matrix(v, k, k, byrow = TRUE,
                           dimnames = list(states, states))
  pooled <- as_matrix(colSums(counts))
  if (sum(pooled) == 0) {
    .ls_bad(sprintf(paste(
      "'%s': no transition is observed (no sequence has two adjacent",
      "observed states), so there is no landscape."), name),
      class = "nestimate_landscape_unestimable")
  }
  depth <- .ls_depth(pooled)
  absorbing <- states[is.infinite(depth)]
  if (length(absorbing)) {
    .ls_bad(sprintf(paste(
      "'%s': %s is never left, so its expected stay (and cup depth) is",
      "infinite."), name, paste(absorbing, collapse = ", ")),
      class = "nestimate_landscape_unestimable")
  }
  boot <- lapply(seq_len(iter), \(b) as_matrix(colSums(
    counts[sample.int(nrow(counts), replace = TRUE), , drop = FALSE])))
  boot_depth <- if (iter > 0L) vapply(boot, .ls_depth, numeric(k)) else
    matrix(NA_real_, k, 0L)
  boot_depth <- matrix(boot_depth, nrow = k)
  ci <- t(apply(boot_depth, 1L, .ls_interval, ci_level = ci_level))
  occupancy <- tabulate(match(states_mat, states), nbins = k)
  out_n <- rowSums(pooled)
  cups <- data.frame(
    state = states, depth = depth, ci_lower = ci[, 1], ci_upper = ci[, 2],
    dwell = exp(depth), stay = diag(pooled) / ifelse(out_n > 0, out_n, NA),
    share = occupancy / sum(occupancy), n_from = out_n,
    resamples = rowSums(!is.na(boot_depth)), row.names = NULL)
  off <- which(row(pooled) != col(pooled), arr.ind = TRUE)
  off <- off[order(off[, 1], off[, 2]), , drop = FALSE]
  exits <- data.frame(
    from = states[off[, 1]], to = states[off[, 2]],
    probability = pooled[off] / ifelse(out_n[off[, 1]] > 0,
                                       out_n[off[, 1]], NA),
    count = pooled[off], row.names = NULL)
  circulation <- NULL
  if (k == 3L) {
    cur <- .ls_current(pooled)
    bc <- if (iter > 0L) vapply(boot, .ls_current, numeric(1)) else numeric(0)
    cv <- .ls_interval(bc, ci_level)
    flip <- isTRUE(cur < 0)
    cyc <- if (flip) rev(states) else states
    circulation <- data.frame(
      cycle = paste(c(cyc, cyc[1]), collapse = " -> "),
      current = abs(cur),
      ci_lower = if (flip) -cv[2] else cv[1],
      ci_upper = if (flip) -cv[1] else cv[2],
      resamples = sum(!is.na(bc)))
  }
  list(cups = cups, exits = exits, circulation = circulation,
       n_sequences = nrow(states_mat), n_transitions = sum(pooled))
}

#' Stability Landscape of a Transition Network
#'
#' @description
#' Draws a transition network as a ball-in-cup landscape (Scheffer et al.
#' 2009) in which the depth of each cup is the stability of the state: how
#' long the process stays once it is there. A state with stay probability
#' \eqn{P_{ii}} is left after \eqn{1 / (1 - P_{ii})} steps on average (its
#' dwell time), and its cup depth is the log of that time,
#' \eqn{-\log(1 - P_{ii})}, so that dwell time = \eqn{e^{depth}}: the
#' escape-time law of a particle in a potential well (Kramers 1940). A state
#' that is left at once has depth 0, no cup at all. Intervals come from
#' resampling whole sequences.
#'
#' @param x A transition network from \code{build_network()} (method
#'   \code{"relative"} or \code{"frequency"}), a grouped one
#'   (\code{netobject_group}, for example \code{build_network()} of a
#'   \code{build_clusters()} result) or a wide sequence data frame or matrix.
#'   The network must carry its sequences. For the \code{print()} and
#'   \code{plot()} methods: an object of class \code{net_landscape}.
#' @param object A \code{net_landscape} object (for \code{summary()}).
#' @param states Character vector naming every state once, in the order to
#'   draw them (for example from least to most engaged). Default: the
#'   alphabetical order.
#' @param iter Integer. Bootstrap resamples of sequences. Default 1000.
#' @param ci_level Numeric in (0, 1). Interval level. Default 0.95.
#' @param seed Integer or \code{NULL}. Seed for the resampling, scoped to the
#'   call (the global random state is restored).
#' @param table One of \code{"cups"}, \code{"exits"}, \code{"circulation"}.
#'   In \code{summary()}: which table to return.
#' @param digits Integer. Decimal places in \code{print()}. Default 2.
#' @param ... Ignored.
#'
#' @details
#' Depth is a property of the state alone; how often the process enters the
#' state is its \code{share} of the observed time steps, reported beside the
#' depth. Where the process goes when it leaves is in the exit probabilities
#' (\code{summary(x, table = "exits")}), which the plot writes under each
#' cup. With three states the plot also reports the net circulation, the net
#' probability flux \eqn{\pi_i P_{ij} - \pi_j P_{ji}} around the cycle, which
#' stationarity makes equal on all three edges (Schnakenberg 1976): zero for
#' a chain whose every transition is as common as its reverse.
#'
#' In a resample where a state is never left, its depth is infinite, so an
#' upper interval limit can be \code{Inf}; a resample in which a state never
#' occurs is left out of that state's interval (\code{resamples} counts the
#' rest). Only the cup floors carry data in the plot; the curve between them
#' is a smooth join, and every rim is drawn at 0.
#'
#' @return An object of class \code{"net_landscape"}: one fit per network
#'   (\code{"network"} for a single one, the group names for a grouped one),
#'   with \code{states}, \code{iter} and \code{ci_level}.
#' @return In \code{summary()}: a data frame. \code{table = "cups"} (default):
#'   one row per network and state, with \code{group}, \code{state},
#'   \code{depth}, \code{ci_lower}, \code{ci_upper}, \code{dwell} (expected
#'   consecutive steps in the state), \code{stay} (\eqn{P_{ii}}),
#'   \code{share} (share of observed time steps in the state), \code{n_from}
#'   (transitions out of the state, self-transitions included) and
#'   \code{resamples}. \code{table = "exits"}: one row per network and pair
#'   of distinct states, with \code{probability} and \code{count}.
#'   \code{table = "circulation"}: one row per network (three states only),
#'   with \code{cycle} (the direction of the net flow), \code{current},
#'   \code{ci_lower}, \code{ci_upper} and \code{resamples} (resamples with a
#'   stationary distribution).
#' @return In \code{print()}: \code{x}, invisibly.
#' @return In \code{plot()}: a \code{ggplot} object, one panel per network.
#'
#' @section Errors:
#' \code{nestimate_not_markov} for a network that is not a transition
#' network; \code{nestimate_bad_input} for \code{states} that do not name
#' every state once, or input without sequences;
#' \code{nestimate_landscape_unestimable} for a state that is never left or
#' input with no observed transition.
#'
#' @examples
#' net <- build_network(trajectories, method = "relative")
#' land <- stability_landscape(net,
#'   states = c("Disengaged", "Average", "Active"), iter = 50, seed = 1)
#' land
#' summary(land, table = "exits")
#' \donttest{
#' plot(land)
#' }
#'
#' @seealso \code{\link{landscape_surface}}, \code{\link{passage_time}},
#'   \code{\link{markov_stability}}
#'
#' @references
#' Kramers, H. A. (1940). Brownian motion in a field of force and the
#' diffusion model of chemical reactions. \emph{Physica}, 7, 284--304.
#'
#' Scheffer, M., Bascompte, J., Brock, W. A., et al. (2009). Early-warning
#' signals for critical transitions. \emph{Nature}, 461, 53--59.
#'
#' Schnakenberg, J. (1976). Network theory of microscopic and macroscopic
#' behavior of master equation systems. \emph{Reviews of Modern Physics},
#' 48, 571--585.
#'
#' @export
stability_landscape <- function(x, states = NULL, iter = 1000L,
                                ci_level = 0.95, seed = NULL) {
  stopifnot(
    "`iter` must be a single non-negative integer" =
      length(iter) == 1L && is.numeric(iter) && is.finite(iter) && iter >= 0,
    "`ci_level` must be a single number in (0, 1)" =
      length(ci_level) == 1L && is.numeric(ci_level) && ci_level > 0 &&
      ci_level < 1
  )
  parts <- if (inherits(x, "netobject_group")) as.list(x) else
    list(network = x)
  seqs <- lapply(parts, .ls_sequences)
  alphabet <- sort(unique(unlist(lapply(seqs, \(s) s$alphabet))))
  if (length(alphabet) < 2L) .ls_bad("A landscape needs at least two states.")
  states <- .ls_states(states, alphabet)
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  fits <- Map(\(s, nm) .ls_fit(s$states, states, as.integer(iter), ci_level,
                               nm), seqs, names(parts))
  structure(list(fits = fits, states = states, iter = as.integer(iter),
                 ci_level = ci_level),
            class = "net_landscape")
}

#' @rdname stability_landscape
#' @export
summary.net_landscape <- function(object,
                                  table = c("cups", "exits", "circulation"),
                                  ...) {
  table <- match.arg(table)
  parts <- Map(\(f, g) {
    tab <- f[[table]]
    if (is.null(tab)) return(NULL)
    cbind(group = g, tab)
  }, object$fits, names(object$fits))
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    .ls_bad("The net circulation is defined for three states only.")
  }
  do.call(rbind, c(parts, make.row.names = FALSE))
}

#' @rdname stability_landscape
#' @export
print.net_landscape <- function(x, digits = 2, ...) {
  n_net <- length(x$fits)
  cat(sprintf("Stability landscape: %d states, %d network%s, %d resamples\n\n",
              length(x$states), n_net, if (n_net == 1L) "" else "s", x$iter))
  cups <- summary(x)
  cups <- cups[, c("group", "state", "depth", "ci_lower", "ci_upper",
                   "dwell", "stay", "share")]
  num <- vapply(cups, is.numeric, logical(1))
  cups[num] <- lapply(cups[num], round, digits = digits)
  if (n_net == 1L) cups$group <- NULL
  cat("Cups (depth = log of the expected stay; deeper = more stable):\n")
  print(cups, row.names = FALSE)
  circ <- Filter(Negate(is.null), lapply(x$fits, \(f) f$circulation))
  if (length(circ)) {
    cat("\nNet circulation (per step):\n")
    invisible(Map(\(cc, g) cat(sprintf(
      "  %s%s: %s [%s, %s]%s\n", if (n_net > 1L) paste0(g, ", ") else "",
      cc$cycle, format(round(cc$current, digits + 2L)),
      format(round(cc$ci_lower, digits + 2L)),
      format(round(cc$ci_upper, digits + 2L)),
      if (cc$resamples < x$iter)
        sprintf(" (%d of %d resamples have a stationary distribution)",
                cc$resamples, x$iter) else "")),
      circ, names(circ)))
  }
  cat("\nExit probabilities: summary(x, table = \"exits\")\n")
  invisible(x)
}

# One cup per state on [i - 0.5, i + 0.5], rim at 0, floor at -depth; a
# cosine ease, flat at the floor and at the rims. Only the floors carry data.
#' @noRd
.ls_cup_profile <- function(depth, n = 60L) {
  s <- seq(0, 1, length.out = n)
  d <- ifelse(is.na(depth), 0, depth)
  out <- do.call(rbind, lapply(seq_along(d), \(i) data.frame(
    x = i - 1.5 + s, u = -d[i] * (1 - cos(2 * pi * s)) / 2)))
  out[!duplicated(out$x), ]
}

#' @rdname stability_landscape
#' @export
plot.net_landscape <- function(x, ...) {
  states <- x$states
  k <- length(states)
  short <- \(s) substr(s, 1L, 3L)
  cups <- summary(x)
  cups$group <- factor(cups$group, levels = names(x$fits))
  exits <- summary(x, table = "exits")
  cups$exit_text <- vapply(seq_len(nrow(cups)), \(i) {
    e <- exits[exits$group == cups$group[i] & exits$from == cups$state[i], ]
    paste(sprintf("to %s %.2f", short(e$to), e$probability), collapse = ", ")
  }, character(1))
  cups$x <- match(cups$state, states) - 1
  cups$y <- -ifelse(is.na(cups$depth), 0, cups$depth)
  cups$ymin <- -cups$ci_upper
  cups$ymax <- -cups$ci_lower
  curve <- do.call(rbind, lapply(names(x$fits), \(g) transform(
    .ls_cup_profile(cups$depth[cups$group == g]), group = g)))
  curve$group <- factor(curve$group, levels = names(x$fits))
  deepest <- max(c(cups$depth, cups$ci_upper[is.finite(cups$ci_upper)]),
                 na.rm = TRUE)
  span <- max(deepest, 0.5)
  cups$label <- sprintf("dwell %.1f, share %.2f\n%s", cups$dwell, cups$share,
                        cups$exit_text)
  cups$label_y <- pmin(cups$y, -ifelse(is.finite(cups$ci_upper),
                                       cups$ci_upper, deepest),
                       na.rm = TRUE) - 0.08 * span
  cups$state <- factor(cups$state, levels = states)
  circ <- do.call(rbind, Map(\(f, g) if (is.null(f$circulation)) NULL else
    transform(f$circulation, group = g), x$fits, names(x$fits)))
  p <- ggplot2::ggplot(curve, ggplot2::aes(x, u)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 1, colour = "#000000") +
    ggplot2::geom_errorbar(
      data = cups, ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, width = 0.07, colour = "#999999", na.rm = TRUE) +
    ggplot2::geom_point(
      data = cups, ggplot2::aes(x = x, y = y + 0.05 * span, fill = state,
                                shape = state),
      size = 7, stroke = 0.6, colour = "#000000") +
    ggplot2::geom_text(
      data = cups, ggplot2::aes(x = x, y = label_y, label = label),
      size = 3.1, vjust = 1, lineheight = 0.95) +
    ggplot2::scale_fill_manual(values = rep_len(.okabe_ito, k)) +
    ggplot2::scale_shape_manual(values = rep_len(c(21, 22, 24, 23, 25), k)) +
    ggplot2::scale_x_continuous(breaks = seq(0, k - 1), labels = states,
                                expand = ggplot2::expansion(add = 0.05)) +
    ggplot2::scale_y_continuous(
      name = "Cup depth, log of the expected stay",
      labels = \(b) format(-b),
      expand = ggplot2::expansion(mult = c(0.45, 0.2)),
      sec.axis = ggplot2::sec_axis(
        ~ exp(-.), name = "Expected stay (steps)",
        breaks = c(1, 1.5, 2, 3, 5, 10, 20, 50, 100))) +
    ggplot2::labs(
      x = NULL, fill = NULL, shape = NULL,
      caption = sprintf(paste(
        "Bars: %d%% interval, %d resamples of sequences. Under each cup:",
        "expected stay, share of time steps, exit probabilities."),
        round(100 * x$ci_level), x$iter)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())
  if (!is.null(circ)) {
    circ$group <- factor(circ$group, levels = names(x$fits))
    circ$label <- sprintf("net circulation %s: %.3f [%.3f, %.3f] per step",
                          circ$cycle, circ$current, circ$ci_lower,
                          circ$ci_upper)
    p <- p + ggplot2::geom_text(
      data = circ, ggplot2::aes(x = (k - 1) / 2, y = Inf, label = label),
      inherit.aes = FALSE, vjust = 1.5, size = 3.4)
  }
  if (length(x$fits) > 1L) p <- p + ggplot2::facet_wrap(~ group, ncol = 1L)
  p
}

# ---- landscape_surface() ----

.ls_measure_labels <- c(
  self_loop_tendency        = "Self-loop tendency (higher = more stable)",
  transition_rate           = "Transition rate (higher = less stable)",
  longitudinal_entropy      = "Longitudinal entropy (higher = less stable)",
  simpson_diversity         = "Simpson diversity (higher = less stable)",
  complexity_index          = "Complexity index (higher = less stable)",
  mean_spell_duration       = "Mean spell, steps (higher = more stable)",
  max_spell_duration        = "Longest spell, steps (higher = more stable)",
  cyclic_feedback_strength  = "Recurrence (higher = more returns)",
  initial_state_persistence = "Initial-state persistence (higher = more inertia)",
  dominant_proportion       = "Dominant-state share (higher = more stable)"
)
.ls_higher_is_stable <- c("self_loop_tendency", "mean_spell_duration",
                          "max_spell_duration", "initial_state_persistence",
                          "dominant_proportion")
.ls_higher_is_unstable <- c("transition_rate", "longitudinal_entropy",
                            "simpson_diversity", "complexity_index")

#' @noRd
.ls_surface_one <- function(states_mat, states, measure, k, n_grid, name) {
  met <- sequence_metrics(states_mat)
  keep <- which(met$valid_n >= 2L & !is.na(met[[measure]]))
  if (length(keep) < 5L) {
    .ls_bad(sprintf(
      "'%s': fewer than 5 sequences with two or more observed states.", name),
      class = "nestimate_landscape_unestimable")
  }
  m <- met[[measure]][keep]
  counts <- .ls_counts(states_mat[keep, , drop = FALSE], states)
  kk <- length(states)
  k_eff <- min(k, length(keep))
  grid <- seq(stats::quantile(m, 0.05, names = FALSE),
              stats::quantile(m, 0.95, names = FALSE), length.out = n_grid)
  if (diff(range(m)) == 0) {
    .ls_bad(sprintf("'%s': %s is the same for every sequence, so there is no axis to trace.",
                    name, sQuote(measure, FALSE)))
  }
  floor_h <- diff(range(m)) / 50
  rows <- lapply(grid, \(y) {
    d <- abs(m - y)
    h <- max(sort(d)[k_eff], floor_h) * 1.01
    w <- ifelse(d < h, (1 - (d / h)^3)^3, 0)
    n <- matrix(colSums(counts * w), kk, kk, byrow = TRUE)
    depth <- .ls_depth(n)
    data.frame(group = name, measure = measure, value = y, state = states,
               depth = depth, dwell = exp(depth),
               effective_n = sum(w)^2 / sum(w^2))
  })
  do.call(rbind, rows)
}

#' Stability Landscape Across a Stability Measure
#'
#' @description
#' The cups of \code{\link{stability_landscape}()} drawn as a surface over a
#' second axis: a per-sequence stability measure from
#' \code{\link{sequence_metrics}()}, such as self-loop tendency or
#' longitudinal entropy. At each value of the measure, sequences are weighted
#' by a tricube kernel on their distance to that value over the \code{k}
#' nearest (Cleveland 1979), and each state's cup depth,
#' \eqn{-\log(1 - P_{ii})} (the log of its expected stay), is computed from
#' their weighted pooled transitions. Moving along the axis goes from the
#' most stable to the most volatile sequences.
#'
#' @param x A transition network from \code{build_network()} (method
#'   \code{"relative"} or \code{"frequency"}), a grouped one, or a wide
#'   sequence data frame or matrix. For the \code{print()} and \code{plot()}
#'   methods: an object of class \code{net_landscape_surface}.
#' @param object A \code{net_landscape_surface} object (for
#'   \code{summary()}).
#' @param measure Character vector of numeric \code{sequence_metrics()}
#'   columns, one surface each. Default \code{"longitudinal_entropy"}.
#' @param states Character vector naming every state once, in drawing order.
#'   Default: alphabetical.
#' @param k Integer. Nearest sequences per slice. Default 30 (capped at the
#'   number of sequences).
#' @param n_grid Integer. Slices between the 5th and 95th percentile of the
#'   measure. Default 25.
#' @param ends Logical. In \code{summary()}: \code{TRUE} (default) compares
#'   the low and the high end of each measure; \code{FALSE} returns every
#'   slice.
#' @param max_stay Numeric. In \code{plot()}: cups of states whose expected
#'   stay exceeds this many steps (including states never left in a slice)
#'   are drawn at the depth \code{log(max_stay)}, a display floor only.
#'   Default 50.
#' @param theta,phi Numeric. In \code{plot()}: viewing angles passed to
#'   \code{persp()}.
#' @param ncol Integer or \code{NULL}. In \code{plot()}: panels per row.
#'   Default: one row per network, one column per measure.
#' @param digits Integer. Decimal places in \code{print()}. Default 2.
#' @param ... Ignored.
#'
#' @details
#' The measures are computed from the same sequences as the depths, so the
#' cups of the most stable sequences are deep by construction (self-loop
#' tendency in particular is closely tied to the stay probabilities). What
#' the surface adds is which cup deepens as stability rises. A slice's
#' \code{effective_n} (\eqn{(\sum w)^2 / \sum w^2}) shows how many sequences
#' it rests on. A state that does not occur in a slice has an \code{NA}
#' depth and leaves a gap in the surface.
#'
#' @return An object of class \code{"net_landscape_surface"} holding the
#'   slices as one data frame (one row per network, measure, slice and
#'   state), with \code{states}, \code{measure} and \code{k}.
#' @return In \code{summary()}: with \code{ends = TRUE}, a data frame with
#'   one row per network, measure and state: \code{group}, \code{measure},
#'   \code{state}, \code{low_value}, \code{depth_at_low}, \code{high_value},
#'   \code{depth_at_high}; with \code{ends = FALSE}, every slice
#'   (\code{group}, \code{measure}, \code{value}, \code{state}, \code{depth},
#'   \code{dwell}, \code{effective_n}).
#' @return In \code{print()}: \code{x}, invisibly.
#' @return In \code{plot()}: \code{x}, invisibly; draws one \code{persp()}
#'   surface per network and measure on a shared depth scale.
#'
#' @section Errors:
#' \code{nestimate_not_markov}, \code{nestimate_bad_input} (unknown,
#' non-numeric or constant \code{measure}, bad \code{states}),
#' \code{nestimate_landscape_unestimable} (fewer than five usable sequences).
#' Warning \code{nestimate_thin_slice} when slices rest on fewer than five
#' effective sequences.
#'
#' @examples
#' net <- build_network(trajectories, method = "relative")
#' surf <- landscape_surface(net, measure = "longitudinal_entropy",
#'   states = c("Disengaged", "Average", "Active"))
#' surf
#' \donttest{
#' plot(surf)
#' }
#'
#' @seealso \code{\link{stability_landscape}}, \code{\link{sequence_metrics}}
#'
#' @references
#' Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. \emph{Journal of the American Statistical Association}, 74,
#' 829--836.
#'
#' Kramers, H. A. (1940). Brownian motion in a field of force and the
#' diffusion model of chemical reactions. \emph{Physica}, 7, 284--304.
#'
#' @export
landscape_surface <- function(x, measure = "longitudinal_entropy",
                              states = NULL, k = 30L, n_grid = 25L) {
  stopifnot(
    "`k` must be a single integer >= 5" =
      length(k) == 1L && is.numeric(k) && k >= 5,
    "`n_grid` must be a single integer >= 2" =
      length(n_grid) == 1L && is.numeric(n_grid) && n_grid >= 2
  )
  if (!is.character(measure) || length(measure) == 0L) {
    .ls_bad("`measure` must name one or more sequence_metrics() columns.")
  }
  parts <- if (inherits(x, "netobject_group")) as.list(x) else
    list(network = x)
  seqs <- lapply(parts, .ls_sequences)
  alphabet <- sort(unique(unlist(lapply(seqs, \(s) s$alphabet))))
  states <- .ls_states(states, alphabet)
  probe <- sequence_metrics(seqs[[1]]$states)
  numeric_cols <- names(probe)[vapply(probe, is.numeric, logical(1))]
  bad <- setdiff(measure, numeric_cols)
  if (length(bad)) {
    .ls_bad(sprintf(
      "`measure` must be numeric sequence_metrics() columns; not: %s.",
      paste(bad, collapse = ", ")))
  }
  slices <- do.call(rbind, unlist(lapply(names(parts), \(g)
    lapply(measure, \(m) .ls_surface_one(seqs[[g]]$states, states, m,
                                         as.integer(k), as.integer(n_grid),
                                         g))), recursive = FALSE))
  rownames(slices) <- NULL
  # tied measure values can leave a slice resting on one or two sequences
  thin <- slices$effective_n < 5
  if (any(thin)) {
    warning(warningCondition(sprintf(paste(
      "%d of %d slices rest on fewer than 5 effective sequences (minimum",
      "%.1f), often because many sequences share one value of the measure;",
      "their depths are unstable."), sum(thin) %/% length(states),
      nrow(slices) %/% length(states), min(slices$effective_n)),
      class = "nestimate_thin_slice"))
  }
  structure(list(slices = slices, states = states, measure = measure,
                 groups = names(parts), k = as.integer(k)),
            class = "net_landscape_surface")
}

#' @rdname landscape_surface
#' @export
summary.net_landscape_surface <- function(object, ends = TRUE, ...) {
  s <- object$slices
  if (!isTRUE(ends)) return(s)
  keys <- unique(s[, c("group", "measure", "state")])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), \(i) {
    r <- s[s$group == keys$group[i] & s$measure == keys$measure[i] &
             s$state == keys$state[i], ]
    lo <- which.min(r$value)
    hi <- which.max(r$value)
    data.frame(keys[i, ], low_value = r$value[lo], depth_at_low = r$depth[lo],
               high_value = r$value[hi], depth_at_high = r$depth[hi],
               row.names = NULL)
  }))
  rownames(out) <- NULL
  out
}

#' @rdname landscape_surface
#' @export
print.net_landscape_surface <- function(x, digits = 2, ...) {
  cat(sprintf(paste0("Landscape surface: %d states, %d network%s x %d measure%s,",
                     " %d nearest sequences per slice\n\n"),
              length(x$states), length(x$groups),
              if (length(x$groups) == 1L) "" else "s", length(x$measure),
              if (length(x$measure) == 1L) "" else "s", x$k))
  tab <- summary(x)
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], round, digits = digits)
  names(tab) <- c("group", "measure", "state", "low", "depth_at_low", "high",
                  "depth_at_high")
  if (length(x$groups) == 1L) tab$group <- NULL
  cat("Cup depth (log of the expected stay) at the low and the high end",
      "of each measure:\n")
  print(tab, row.names = FALSE)
  eff <- vapply(x$groups, \(g) {
    e <- x$slices$effective_n[x$slices$group == g]
    sprintf("%s %.1f-%.1f", g, min(e), max(e))
  }, character(1))
  cat(sprintf("\nEffective sequences per slice: %s\n",
              paste(eff, collapse = ", ")))
  invisible(x)
}

#' @rdname landscape_surface
#' @export
plot.net_landscape_surface <- function(x, max_stay = 50, theta = 10, phi = 30,
                                       ncol = NULL, ...) {
  stopifnot("`max_stay` must be a single number > 1" =
              length(max_stay) == 1L && is.numeric(max_stay) && max_stay > 1)
  states <- x$states
  kk <- length(states)
  floor_depth <- log(max_stay)
  panels <- unique(x$slices[, c("group", "measure")])
  surf_of <- function(g, m) {
    s <- x$slices[x$slices$group == g & x$slices$measure == m, ]
    values <- sort(unique(s$value))
    # the stable end of every measure faces the viewer, so its (deep) cups
    # are not hidden behind the rims of the shallow ones
    if (m %in% .ls_higher_is_stable) values <- rev(values)
    profs <- lapply(values, \(v) {
      d <- s$depth[s$value == v]
      d <- ifelse(is.finite(d) | is.na(d), pmin(d, floor_depth), floor_depth)
      p <- .ls_cup_profile(d, n = 16L)
      # a state absent from the slice leaves a gap, not a flat floor
      p$u[is.na(d[pmin(floor(p$x + 0.5) + 1, kk)])] <- NA_real_
      p
    })
    list(x = profs[[1]]$x, y = seq_along(values), values = values,
         z = vapply(profs, \(p) p$u, numeric(nrow(profs[[1]]))))
  }
  surfs <- Map(surf_of, panels$group, panels$measure)
  zlim <- range(c(0, unlist(lapply(surfs, \(s) s$z))), finite = TRUE)
  if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)   # a flat surface
  depth_ticks <- seq(0, floor(-zlim[1]))
  ramp <- grDevices::colorRampPalette(c("#009E73", "#F0E442", "#E69F00"))(100L)
  ncol <- ncol %||% length(x$measure)
  old_par <- graphics::par(mfrow = c(ceiling(nrow(panels) / ncol), ncol),
                           mar = c(1, 2, 3, 6), oma = c(1.5, 0, 0, 0),
                           xpd = NA)
  on.exit(graphics::par(old_par), add = TRUE)
  invisible(Map(\(s, g, m) {
    z <- s$z
    zf <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] +
             z[-nrow(z), -ncol(z)]) / 4
    col <- ramp[cut(zf, breaks = seq(zlim[1], zlim[2], length.out = 101L),
                    include.lowest = TRUE)]
    label <- unname(.ls_measure_labels[m])
    if (is.na(label)) label <- m
    main_title <- if (length(x$groups) > 1L) paste0(g, ": ", label) else label
    pm <- graphics::persp(s$x, s$y, z, zlim = zlim, theta = theta, phi = phi,
                          col = col, border = "grey25", lwd = 0.3,
                          shade = 0.35, box = FALSE, expand = 0.7,
                          main = main_title, cex.main = 0.9)
    front <- min(s$y)
    back <- max(s$y)
    left <- min(s$x)
    right <- max(s$x)
    graphics::text(grDevices::trans3d(seq(0, kk - 1), front, zlim[1], pm),
                   labels = states, cex = 0.75, font = 2, adj = c(0.5, 1.6))
    # depth scale on the front-left edge, in log-stay units
    graphics::segments(
      grDevices::trans3d(left, front, 0, pm)$x,
      grDevices::trans3d(left, front, 0, pm)$y,
      grDevices::trans3d(left, front, zlim[1], pm)$x,
      grDevices::trans3d(left, front, zlim[1], pm)$y, col = "grey40")
    graphics::text(grDevices::trans3d(left, front, -depth_ticks, pm),
                   labels = depth_ticks, cex = 0.65, pos = 2)
    graphics::text(grDevices::trans3d(left, front, zlim[1] / 2, pm),
                   labels = "depth", cex = 0.7, pos = 2, offset = 1.6)
    # the two ends of the measure, with their values
    graphics::text(grDevices::trans3d(right, c(front, back), 0, pm),
                   labels = sprintf("%s: %s",
                                    if (m %in% c(.ls_higher_is_stable,
                                                 .ls_higher_is_unstable))
                                      c("most stable", "least stable") else
                                      c("low", "high"),
                                    format(round(c(s$values[1],
                                                   s$values[length(s$values)]),
                                                 2))),
                   cex = 0.65, pos = 4)
  }, surfs, panels$group, panels$measure))
  graphics::mtext(sprintf(paste(
    "Front: the most stable end of each measure. Cup depth = log expected",
    "stay; stays above %g steps drawn at log(%g)."),
    max_stay, max_stay), side = 1, outer = TRUE, line = 0.3, cex = 0.7)
  invisible(x)
}


# ---- landscape_trend() ----

# Pooled k x k count matrix from per-sequence counts and sequence weights.
#' @noRd
.ls_pool <- function(counts, w, k) {
  matrix(colSums(counts * w), k, k, byrow = TRUE)
}

# Occupancy landscape of one pooled count matrix: cups at -log(pi) (pi the
# stationary distribution) and ridges at -log of the symmetrised flux
# (pi_i P_ij + pi_j P_ji) / 2 between neighbours in state order and between
# the last and the first state. All NA when the chain has no unique
# stationary distribution; a ridge never crossed is NA.
#' @noRd
.ls_n_ridges <- function(k) if (k == 2L) 1L else k

.ls_occupancy <- function(n) {
  k <- nrow(n)
  out_n <- rowSums(n)
  n_out <- k + .ls_n_ridges(k)
  if (any(out_n == 0)) return(rep(NA_real_, n_out))
  P <- n / out_n
  if (!.ls_irreducible(P)) return(rep(NA_real_, n_out))
  pi <- .mpt_stationary(P)
  nxt <- c(seq_len(k)[-1L], 1L)
  sym <- (pi * P[cbind(seq_len(k), nxt)] + pi[nxt] * P[cbind(nxt, seq_len(k))]) / 2
  # with two states the last-first ridge is the one pair already counted
  sym <- sym[seq_len(.ls_n_ridges(k))]
  c(-log(pi), ifelse(sym > 0, -log(sym), NA_real_))
}

# Smooth join through a sequence of extrema (cups and ridges): a cosine ease,
# monotone and flat at each extremum. Only the extrema carry data.
#' @noRd
.ls_profile <- function(ext_x, ext_u, n = 8L) {
  s <- seq(0, 1, length.out = n)
  out <- do.call(rbind, lapply(seq_len(length(ext_x) - 1L), \(i) data.frame(
    x = ext_x[i] + s * (ext_x[i + 1L] - ext_x[i]),
    u = ext_u[i] + (ext_u[i + 1L] - ext_u[i]) * (1 - cos(pi * s)) / 2)))
  out[!duplicated(out$x), ]
}

#' @noRd
.ls_trend_slices <- function(states_mat, covariate, states, along, window, k,
                             n_grid) {
  kk <- length(states)
  if (identical(along, "phase")) {
    n_pos <- ncol(states_mat)
    if (window > n_pos) {
      .ls_bad(sprintf("`window` (%d) is longer than the sequences (%d).",
                      window, n_pos))
    }
    starts <- seq_len(n_pos - window + 1L)
    counts <- lapply(starts, \(t) .ls_counts(
      states_mat[, t:(t + window - 1L), drop = FALSE], states))
    return(list(
      value = starts + (window - 1) / 2,
      label = sprintf("%d-%d", starts, starts + window - 1L),
      pooled_of = function(idx) lapply(counts, \(cn)
        .ls_pool(cn[idx, , drop = FALSE], 1, kk)),
      n_used = nrow(states_mat)))
  }
  keep <- which(!is.na(covariate) & rowSums(!is.na(states_mat)) >= 2L)
  if (length(keep) < 5L) {
    .ls_bad("Fewer than 5 sequences with the covariate and two or more observed states.",
            class = "nestimate_landscape_unestimable")
  }
  m <- covariate[keep]
  counts <- .ls_counts(states_mat[keep, , drop = FALSE], states)
  grid <- seq(stats::quantile(m, 0.05, names = FALSE),
              stats::quantile(m, 0.95, names = FALSE), length.out = n_grid)
  if (diff(range(m)) == 0) {
    .ls_bad(sprintf("`%s` is the same for every sequence, so there is no axis to trace.",
                    along))
  }
  floor_h <- diff(range(m)) / 50
  weights_at <- function(mm, y) {
    d <- abs(mm - y)
    h <- max(sort(d)[min(k, length(d))], floor_h) * 1.01
    ifelse(d < h, (1 - (d / h)^3)^3, 0)
  }
  list(
    value = grid, label = format(round(grid, 2)),
    pooled_of = function(idx) lapply(grid, \(y)
      .ls_pool(counts[idx, , drop = FALSE], weights_at(m[idx], y), kk)),
    n_used = length(keep))
}

#' Stability Landscape Along the Program or a Covariate
#'
#' @description
#' The landscape of \code{\link{stability_landscape}()} traced along a second
#' axis. The axis is either the program phase, \code{along = "phase"}: the
#' landscape is estimated from the transitions inside a sliding window of
#' \code{window} time steps; or a numeric per-sequence covariate such as a
#' grade, \code{along = "grade"}: at each value of the covariate it comes
#' from the \code{k} nearest sequences, weighted by a tricube kernel
#' (Cleveland 1979). Intervals come from resampling whole sequences.
#'
#' Two heights are available. \code{depth = "stay"} (default) is the depth
#' of \code{stability_landscape()}, the log of each state's expected stay,
#' \eqn{-\log(1 - P_{ii})}. \code{depth = "occupancy"} is the quasi-potential
#' \eqn{-\log \pi_i} of the stationary distribution, with ridges at
#' \eqn{-\log} of the symmetrised probability flux
#' \eqn{(\pi_i P_{ij} + \pi_j P_{ji}) / 2} between neighbouring states
#' (Wang et al. 2008): a low cup is a state the process spends much of its
#' time in, a high ridge a switch it rarely makes.
#'
#' @param x A transition network from \code{build_network()} (method
#'   \code{"relative"} or \code{"frequency"}), a grouped one, or a wide
#'   sequence data frame. A covariate is read from the network's metadata,
#'   which \code{build_network()} fills with the numeric columns of a wide
#'   data frame. For the \code{print()} and \code{plot()} methods: an object
#'   of class \code{net_landscape_trend}.
#' @param object A \code{net_landscape_trend} object (for \code{summary()}).
#' @param along \code{"phase"} (default), or the name of a numeric
#'   per-sequence covariate.
#' @param depth \code{"stay"} (default) or \code{"occupancy"}, see above.
#' @param states Character vector naming every state once, in drawing order.
#'   Default: alphabetical.
#' @param window Integer. With \code{along = "phase"}: time steps per window.
#'   Default 5.
#' @param k Integer. With a covariate: nearest sequences per slice. Default
#'   30 (capped at the number of sequences).
#' @param n_grid Integer. With a covariate: slices between its 5th and 95th
#'   percentile. Default 20.
#' @param iter Integer. Bootstrap resamples of sequences. Default 500.
#' @param ci_level Numeric in (0, 1). Interval level. Default 0.95.
#' @param seed Integer or \code{NULL}. Seed for the resampling, scoped to the
#'   call.
#' @param type In \code{plot()}: \code{"lines"} (default), one line per state
#'   with its interval band; or \code{"surface"}, a 3D \code{persp()}
#'   landscape with the states across and the second axis in depth, one
#'   panel per network on a shared height scale.
#' @param theta,phi Numeric. In \code{plot(type = "surface")}: viewing
#'   angles.
#' @param digits Integer. Decimal places in \code{print()}. Default 2.
#' @param ... Ignored.
#'
#' @details
#' With \code{depth = "stay"}, a window or slice in which a state is never
#' left gives an infinite depth, and one in which it does not occur gives
#' \code{NA}. With \code{depth = "occupancy"}, a window whose chain has no
#' unique stationary distribution is \code{NA} throughout. Such points are
#' kept in \code{summary()}, left out of the lines and the surface, and
#' counted in the plot caption. Interval bands are clipped at the panel
#' edge, which is set by the estimates. Consecutive windows share all but
#' one step and neighbouring slices most of their sequences, so neighbouring
#' points are not independent.
#'
#' @return An object of class \code{"net_landscape_trend"}.
#' @return In \code{summary()}: a data frame with one row per network,
#'   position and feature: \code{group}, \code{along}, \code{value} (window
#'   centre or covariate value), \code{label} (window steps or value),
#'   \code{feature} (\code{"well"}, or \code{"ridge"} with
#'   \code{depth = "occupancy"}), \code{state} (a state, or two for a
#'   ridge), the height (\code{depth} for \code{"stay"}, \code{potential}
#'   for \code{"occupancy"}), \code{ci_lower}, \code{ci_upper}, \code{dwell}
#'   (\code{"stay"} only) and \code{n_from} (transitions out of the state in
#'   the window, or their weighted sum in the slice; wells only).
#' @return In \code{print()}: \code{x}, invisibly.
#' @return In \code{plot()}: with \code{type = "lines"} a \code{ggplot}
#'   object; with \code{type = "surface"} \code{x}, invisibly, after drawing.
#'
#' @section Errors:
#' \code{nestimate_not_markov}; \code{nestimate_bad_input} for an unknown,
#' non-numeric or constant covariate, bad \code{states}, or a \code{window} longer than
#' the sequences; \code{nestimate_landscape_unestimable} for fewer than five
#' usable sequences.
#'
#' @examples
#' net <- build_network(trajectories, method = "relative")
#' trend <- landscape_trend(net, along = "phase",
#'   states = c("Disengaged", "Average", "Active"), iter = 50, seed = 1)
#' trend
#' \donttest{
#' plot(trend)
#' occupancy <- landscape_trend(net, depth = "occupancy", window = 6,
#'   states = c("Disengaged", "Average", "Active"), iter = 50, seed = 1)
#' plot(occupancy, type = "surface")
#' }
#'
#' @seealso \code{\link{stability_landscape}}, \code{\link{landscape_surface}}
#'
#' @references
#' Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. \emph{Journal of the American Statistical Association}, 74,
#' 829--836.
#'
#' Kramers, H. A. (1940). Brownian motion in a field of force and the
#' diffusion model of chemical reactions. \emph{Physica}, 7, 284--304.
#'
#' Wang, J., Xu, L., & Wang, E. (2008). Potential landscape and flux
#' framework of nonequilibrium networks. \emph{PNAS}, 105, 12271--12276.
#'
#' @export
landscape_trend <- function(x, along = "phase",
                            depth = c("stay", "occupancy"), states = NULL,
                            window = 5L, k = 30L, n_grid = 20L, iter = 500L,
                            ci_level = 0.95, seed = NULL) {
  depth <- match.arg(depth)
  stopifnot(
    "`along` must be a single string" =
      is.character(along) && length(along) == 1L,
    "`window` must be a single integer >= 2" =
      length(window) == 1L && is.numeric(window) && window >= 2,
    "`k` must be a single integer >= 5" =
      length(k) == 1L && is.numeric(k) && k >= 5,
    "`n_grid` must be a single integer >= 2" =
      length(n_grid) == 1L && is.numeric(n_grid) && n_grid >= 2,
    "`iter` must be a single non-negative integer" =
      length(iter) == 1L && is.numeric(iter) && iter >= 0,
    "`ci_level` must be a single number in (0, 1)" =
      length(ci_level) == 1L && ci_level > 0 && ci_level < 1
  )
  parts <- if (inherits(x, "netobject_group")) as.list(x) else
    list(network = x)
  seqs <- lapply(parts, .ls_sequences)
  alphabet <- sort(unique(unlist(lapply(seqs, \(s) s$alphabet))))
  states <- .ls_states(states, alphabet)
  covariates <- lapply(seqs, \(s) {
    if (identical(along, "phase")) return(NULL)
    v <- s$units[[along]]
    if (is.null(v) || !is.numeric(v)) {
      .ls_bad(sprintf(paste(
        "`along` must be \"phase\" or a numeric covariate of the sequences;",
        "found: %s."), paste(c(names(s$units), "phase"), collapse = ", ")))
    }
    v
  })
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  kk <- length(states)
  nxt <- c(seq_len(kk)[-1L], 1L)
  metric <- if (depth == "stay") .ls_depth else .ls_occupancy
  n_metric <- if (depth == "stay") kk else kk + .ls_n_ridges(kk)
  tabs <- Map(\(s, cv, g) {
    sl <- .ls_trend_slices(s$states, cv, states, along, as.integer(window),
                           as.integer(k), as.integer(n_grid))
    pooled <- sl$pooled_of(seq_len(sl$n_used))
    est <- vapply(pooled, metric, numeric(n_metric))
    boot <- lapply(seq_len(iter), \(b) vapply(
      sl$pooled_of(sample.int(sl$n_used, replace = TRUE)), metric,
      numeric(n_metric)))
    n_pos <- length(sl$value)
    ci <- vapply(seq_len(n_metric * n_pos), \(cell) .ls_interval(
      vapply(boot, \(bm) bm[cell], numeric(1)), ci_level), numeric(2))
    n_from <- vapply(pooled, rowSums, numeric(kk))
    out <- data.frame(
      group = g, along = along,
      value = rep(sl$value, each = n_metric),
      label = rep(sl$label, each = n_metric),
      feature = rep(c(rep("well", kk), rep("ridge", n_metric - kk)),
                    times = n_pos),
      state = rep(c(states, if (n_metric > kk)
        paste(states, states[nxt], sep = " | ")[seq_len(n_metric - kk)]),
        times = n_pos),
      height = as.vector(est), ci_lower = ci[1, ], ci_upper = ci[2, ],
      n_from = as.vector(rbind(n_from, matrix(NA_real_, n_metric - kk,
                                              n_pos))))
    if (depth == "stay") {
      out$dwell <- exp(out$height)
      names(out)[names(out) == "height"] <- "depth"
    } else {
      names(out)[names(out) == "height"] <- "potential"
    }
    out
  }, seqs, covariates, names(parts))
  structure(list(table = do.call(rbind, c(tabs, make.row.names = FALSE)),
                 states = states, along = along, depth = depth,
                 window = as.integer(window), k = as.integer(k),
                 iter = as.integer(iter), ci_level = ci_level,
                 groups = names(parts)),
            class = "net_landscape_trend")
}

#' @noRd
.ls_height_col <- function(x) if (x$depth == "stay") "depth" else "potential"

#' @rdname landscape_trend
#' @export
summary.net_landscape_trend <- function(object, ...) object$table

#' @rdname landscape_trend
#' @export
print.net_landscape_trend <- function(x, digits = 2, ...) {
  axis <- if (identical(x$along, "phase"))
    sprintf("program phase (%d-step windows)", x$window) else
    sprintf("%s (%d nearest sequences per slice)", x$along, x$k)
  cat(sprintf("Stability landscape along %s: %d states, %d network%s, %d resamples\n\n",
              axis, length(x$states), length(x$groups),
              if (length(x$groups) == 1L) "" else "s", x$iter))
  hc <- .ls_height_col(x)
  tab <- x$table[x$table$feature == "well", ]
  ends <- do.call(rbind, lapply(split(tab, list(tab$group, tab$state),
                                      drop = TRUE), \(r) {
    r <- r[order(r$value), ]
    h <- r[[hc]]
    fin <- h[is.finite(h)]
    data.frame(group = r$group[1], state = r$state[1],
               first = r$label[1], at_first = h[1],
               last = r$label[nrow(r)], at_last = h[nrow(r)],
               min = if (length(fin)) min(fin) else NA_real_,
               max = if (length(fin)) max(fin) else NA_real_)
  }))
  ends <- ends[order(match(ends$group, x$groups),
                     match(ends$state, x$states)), ]
  num <- vapply(ends, is.numeric, logical(1))
  ends[num] <- lapply(ends[num], round, digits = digits)
  if (length(x$groups) == 1L) ends$group <- NULL
  cat(if (x$depth == "stay")
    "Cup depth (log of the expected stay; higher = more stable) at the two ends and its range:\n" else
    "Cup potential (-log pi; lower = more time in the state) at the two ends and its range:\n")
  print(ends, row.names = FALSE)
  bad <- sum(!is.finite(x$table[[hc]]))
  if (bad > 0L) {
    cat(sprintf("\n%d value%s not estimable (state never left, absent, or no stationary distribution).\n",
                bad, if (bad == 1L) "" else "s"))
  }
  invisible(x)
}

#' @rdname landscape_trend
#' @export
plot.net_landscape_trend <- function(x, type = c("lines", "surface"),
                                     theta = 10, phi = 30, ...) {
  type <- match.arg(type)
  if (type == "surface") return(.ls_trend_surface(x, theta, phi))
  hc <- .ls_height_col(x)
  tab <- x$table[x$table$feature == "well", ]
  tab$h <- tab[[hc]]
  bad <- sum(!is.finite(tab$h))
  tab <- tab[is.finite(tab$h), ]
  tab$ci_upper[!is.finite(tab$ci_upper)] <- NA_real_
  tab$state <- factor(tab$state, levels = x$states)
  tab$group <- factor(tab$group, levels = x$groups)
  kk <- length(x$states)
  # the panel follows the estimates; bands wider than that are clipped
  lims <- range(tab$h)
  lims <- lims + c(-0.15, 0.15) * max(diff(lims), 0.5)
  if (x$depth == "stay") lims[1] <- max(lims[1], 0)
  last <- do.call(rbind, lapply(split(tab, list(tab$group, tab$state),
                                      drop = TRUE), \(r) r[which.max(r$value), ]))
  xlab <- if (identical(x$along, "phase"))
    sprintf("Program phase (centre of a %d-step window)", x$window) else
    x$along
  y_scale <- if (x$depth == "stay") {
    ggplot2::scale_y_continuous(
      name = "Cup depth, log of the expected stay",
      sec.axis = ggplot2::sec_axis(~ exp(.), name = "Expected stay (steps)",
                                   breaks = c(1, 2, 5, 10, 20, 50)))
  } else {
    ggplot2::scale_y_continuous(
      name = "Cup potential, -log pi (lower = more time in the state)")
  }
  p <- ggplot2::ggplot(tab, ggplot2::aes(value, h, colour = state,
                                         fill = state, shape = state,
                                         linetype = state)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                         alpha = 0.15, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.2, fill = "white") +
    ggplot2::geom_text(data = last, ggplot2::aes(label = state), hjust = -0.15,
                       size = 3.4, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = rep_len(.okabe_ito, kk)) +
    ggplot2::scale_fill_manual(values = rep_len(.okabe_ito, kk)) +
    ggplot2::scale_shape_manual(values = rep_len(c(21, 22, 24, 23, 25), kk)) +
    ggplot2::scale_linetype_manual(
      values = rep_len(c("solid", "dashed", "dotdash", "longdash", "twodash"),
                       kk)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.15))) +
    y_scale +
    ggplot2::coord_cartesian(ylim = lims) +
    ggplot2::labs(x = xlab, colour = NULL, fill = NULL, shape = NULL,
                  linetype = NULL,
                  caption = paste0(sprintf(
                    "Bands: %d%% interval, %d resamples of sequences, clipped at the panel edge.",
                    round(100 * x$ci_level), x$iter),
                    if (bad > 0L) sprintf(" %d point%s not estimable, not drawn.",
                                          bad, if (bad == 1L) "" else "s")
                    else "")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
  if (length(x$groups) > 1L) p <- p + ggplot2::facet_wrap(~ group, ncol = 1L)
  p
}

# 3D landscape: states across (each once; with occupancy the last-first
# ridge is the outer wall at both ends), the second axis in depth.
#' @noRd
.ls_trend_surface <- function(x, theta, phi) {
  hc <- .ls_height_col(x)
  kk <- length(x$states)
  surf_of <- function(g) {
    s <- x$table[x$table$group == g, ]
    values <- sort(unique(s$value))
    profs <- lapply(values, \(v) {
      r <- s[s$value == v, ]
      well <- r[[hc]][r$feature == "well"]
      if (x$depth == "stay") {
        d <- ifelse(is.finite(well), well, NA_real_)
        p <- .ls_cup_profile(d, n = 16L)
        p$u[is.na(d[pmin(floor(p$x + 0.5) + 1, kk)])] <- NA_real_
        return(p)
      }
      ridge <- rep_len(r[[hc]][r$feature == "ridge"], kk)
      .ls_profile(seq(-0.5, kk - 0.5, by = 0.5),
                  c(ridge[kk], c(rbind(well, ridge))))
    })
    list(x = profs[[1]]$x, y = values,
         z = vapply(profs, \(p) p$u, numeric(nrow(profs[[1]]))))
  }
  surfs <- lapply(x$groups, surf_of)
  zlim <- range(unlist(lapply(surfs, \(s) s$z)), finite = TRUE)
  if (x$depth == "stay") zlim[2] <- max(zlim[2], 0)
  if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)   # a flat surface
  ramp <- grDevices::colorRampPalette(c("#009E73", "#F0E442", "#E69F00"))(100L)
  old_par <- graphics::par(mfrow = c(1L, length(surfs)), mar = c(2, 3, 3, 3))
  on.exit(graphics::par(old_par), add = TRUE)
  ylab <- if (identical(x$along, "phase")) "Course (window centre)" else x$along
  zlab <- if (x$depth == "stay") "-depth (log stay)" else "Potential"
  invisible(Map(\(s, g) {
    z <- s$z
    zf <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] +
             z[-nrow(z), -ncol(z)]) / 4
    if (all(is.na(zf))) {
      # fewer than two neighbouring estimable positions: nothing to draw
      graphics::plot.new()
      graphics::title(main = if (length(x$groups) > 1L) g else "")
      graphics::text(0.5, 0.5, cex = 0.9, labels = sprintf(
        "Not estimable:\nfewer than two neighbouring %s\nwith a value for every state.",
        if (identical(x$along, "phase")) "windows" else "slices"))
      return(invisible(NULL))
    }
    col <- ramp[cut(zf, breaks = seq(zlim[1], zlim[2], length.out = 101L),
                    include.lowest = TRUE)]
    pm <- graphics::persp(s$x, s$y, z, zlim = zlim, theta = theta, phi = phi,
                          col = col, border = "grey25", lwd = 0.3,
                          shade = 0.35, ticktype = "detailed", nticks = 4,
                          expand = 0.7, xlab = "", ylab = ylab, zlab = zlab,
                          main = if (length(x$groups) > 1L) g else "")
    graphics::text(grDevices::trans3d(seq(0, kk - 1), min(s$y), zlim[1], pm),
                   labels = x$states, cex = 0.8, font = 2, adj = c(0.5, 1.6))
  }, surfs, x$groups))
  invisible(x)
}
