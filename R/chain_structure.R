# chain_structure() -- qualitative analysis of a discrete-time Markov chain.
#
# Computes properties that are independent of any starting distribution:
# state classification (absorbing / recurrent / transient), communicating
# classes (SCCs of the support graph), per-state period, hitting probability
# matrix H[i,j] = P(ever reach j | start at i), and absorption analysis when
# absorbing states exist. Used as a diagnostic before trusting passage_time()
# / markov_stability(): if a chain is not regular, the stationary distribution
# those report can mix multiple behavioural phases.
#
# All linear algebra; no new dependencies.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Reachability matrix of a directed adjacency matrix.
#'
#' `R[i, j] = 1` iff there is a path of any positive length from i to j (note:
#' a self-loop or a cycle through i marks `R[i, i] = 1`; a state with no
#' outgoing edges has `R[i, ]` all zero except itself which is also zero).
#' Computed via repeated boolean multiplication; O(n^4) but n is small for
#' Markov chains in this package.
#' @noRd
.cs_reachability <- function(A) {
  n <- nrow(A)
  R <- A > 0
  for (k in seq_len(n - 1L)) {
    R_new <- (R | (R %*% A) > 0)
    if (identical(R_new, R)) break
    R <- R_new
  }
  R * 1L
}

#' Strongly connected components (communicating classes) of A.
#'
#' Two states are in the same SCC iff each is reachable from the other.
#' Returns a list of integer index vectors. Self-loop is counted as
#' "reachable from self".
#' @noRd
.cs_scc <- function(A) {
  n <- nrow(A)
  if (n == 0L) return(list())
  R <- .cs_reachability(A)
  diag(R) <- 1L
  M <- R & t(R)
  keys <- apply(M, 1L, function(row) paste(which(row), collapse = ","))
  unique_keys <- unique(keys)
  lapply(unique_keys, function(k) which(keys == k))
}

#' GCD of a non-empty integer vector via Euclid (no external deps).
#' @noRd
.cs_gcd <- function(x) {
  x <- x[x > 0L]
  if (length(x) == 0L) return(NA_integer_)
  out <- x[1L]
  for (v in x[-1L]) {
    while (v != 0L) { tmp <- out %% v; out <- v; v <- tmp }
  }
  as.integer(out)
}

#' Period of a recurrent class given its support submatrix.
#'
#' Period of state i = gcd{n >= 1 : P^n\eqn{[i, i]} > 0}. For an irreducible
#' aperiodic class on a self-looped state, returns 1. Returns NA if no
#' return is possible (which shouldn't happen for a recurrent class).
#' @noRd
.cs_class_period <- function(A_sub) {
  n <- nrow(A_sub)
  if (n == 0L) return(NA_integer_)
  if (n == 1L) return(if (A_sub[1L, 1L] > 0) 1L else NA_integer_)
  An <- A_sub > 0
  An_running <- An
  cycle_lengths <- integer()
  if (any(diag(An_running))) cycle_lengths <- c(cycle_lengths, 1L)
  for (k in seq_len(2L * n)) {
    An_running <- ((An_running %*% An) > 0)
    if (any(diag(An_running))) {
      cycle_lengths <- c(cycle_lengths, k + 1L)
      if (length(cycle_lengths) >= 3L && .cs_gcd(cycle_lengths) == 1L) break
    }
  }
  .cs_gcd(cycle_lengths)
}

#' Hitting probability matrix `H[i, j] = P(T_j < infty | X_0 = i)`.
#'
#' For i != j: T_j = inf{n >= 0 : X_n = j}, so the off-diagonal hitting
#' probability is the standard "eventually reach j" probability. For i == j:
#' uses T_j = inf{n >= 1 : X_n = j} (return-time convention, matching
#' `markovchain::hittingProbabilities`), so `H[j, j] = sum_k P[j, k] H[k, j]`.
#'
#' For each column j, the off-diagonal system is solved over the states from
#' which j is reachable; states that cannot reach j get `H[i, j] = 0` (the
#' standard *minimal non-negative solution*; cf. Norris, *Markov Chains*,
#' Theorem 1.3.2). Without the reachability restriction the system is
#' rank-deficient whenever a closed class disjoint from j exists.
#' @noRd
.cs_hitting <- function(P, reach = NULL) {
  n <- nrow(P)
  state_names <- rownames(P)
  H <- matrix(0, n, n, dimnames = list(state_names, state_names))
  if (is.null(reach)) {
    A <- (P > 0) * 1L
    reach <- .cs_reachability(A)
    diag(reach) <- 1L
  }
  for (j in seq_len(n)) {
    if (n == 1L) {
      H[j, j] <- if (P[1L, 1L] > 0) 1 else 0
      next
    }
    can_reach <- setdiff(which(reach[, j] > 0), j)
    if (length(can_reach) > 0L) {
      P_sub <- P[can_reach, can_reach, drop = FALSE]
      rhs   <- P[can_reach, j]
      I_sub <- diag(length(can_reach))
      h <- tryCatch(
        drop(solve(I_sub - P_sub, rhs)),
        error = function(e) rep(NA_real_, length(can_reach))
      )
      H[can_reach, j] <- h
    }
  }
  # Diagonal via return-time recursion. H[j, j] is the probability of
  # returning to j in >= 1 steps: the self-loop contributes P[j, j],
  # and any path leaving j must eventually hit j again with probability
  # H[k, j] (already computed off-diagonal).
  for (j in seq_len(n)) {
    off_idx <- setdiff(seq_len(n), j)
    H[j, j] <- P[j, j] + sum(P[j, off_idx] * H[off_idx, j])
  }
  H[!is.na(H)] <- pmin(pmax(H[!is.na(H)], 0), 1)
  H
}

#' Absorption analysis (transient -> closed classes) via fundamental matrix.
#'
#' Kemeny & Snell (1976, ch. III): lump every closed (recurrent) class into
#' one absorbing state. With Q the transient-to-transient block and R the
#' transient-to-class block (row sums of P over each class), N = (I - Q)^-1,
#' B = N R gives the probability of ending in each closed class (rows sum to
#' 1) and t = N 1 the expected steps to enter one. A singleton absorbing
#' state is its own class, labelled by its name; a class of several states is
#' labelled by its states joined with "+". Using only the absorbing states as
#' columns (as before) dropped the mass that flows into non-absorbing closed
#' classes, so rows did not sum to 1. I - Q is non-singular because Q holds
#' only transient states.
#'
#' Returns NULL fields when there is no transient state.
#' @noRd
.cs_absorption <- function(P, transient_idx, closed_classes, state_names) {
  if (length(transient_idx) == 0L || length(closed_classes) == 0L) {
    return(list(probabilities = NULL, mean_time = NULL))
  }
  closed_classes <- closed_classes[order(vapply(closed_classes, min,
                                                integer(1)))]
  Q <- P[transient_idx, transient_idx, drop = FALSE]
  R <- vapply(closed_classes, function(cl) {
    rowSums(P[transient_idx, cl, drop = FALSE])
  }, numeric(length(transient_idx)))
  R <- matrix(R, nrow = length(transient_idx))
  N <- solve(diag(nrow(Q)) - Q)
  abs_probs <- N %*% R
  rownames(abs_probs) <- state_names[transient_idx]
  colnames(abs_probs) <- vapply(closed_classes, function(cl) {
    paste(state_names[sort(cl)], collapse = "+")
  }, character(1))
  abs_time <- as.numeric(N %*% rep(1, nrow(N)))
  names(abs_time) <- state_names[transient_idx]
  list(probabilities = abs_probs, mean_time = abs_time)
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Qualitative structure of a discrete-time Markov chain
#'
#' @description
#' `chain_structure()` describes the qualitative structure of a
#' discrete-time Markov chain: the properties that depend only on which
#' transitions are possible (the support of the transition matrix), not on
#' any starting distribution. It classifies every state as absorbing,
#' recurrent or transient, finds the communicating classes and the period of
#' each recurrent state, says whether the chain is irreducible, aperiodic,
#' regular and reversible, and computes the hitting probabilities and, when
#' transient states exist, the probability that each ends in each closed
#' class (an absorbing state or a closed class of several states) and the
#' expected number of steps it takes to enter one. It reads a transition network, a `tna` model, a transition
#' matrix or sequences, and analyses a grouped network one network at a
#' time.
#'
#' It is a diagnostic to run before trusting [passage_time()] or
#' [markov_stability()]. Both assume a regular chain (irreducible and
#' aperiodic), whose stationary distribution is unique; `is_regular` says
#' whether that holds. The structure follows the estimated matrix: a
#' transition never observed in the data counts as impossible, so a chain
#' estimated from few sequences can look reducible because a rare move was
#' not seen.
#'
#' @details
#' The fundamental-matrix absorption analysis follows Kemeny & Snell (1976),
#' with every closed class lumped into one absorbing state, so each
#' transient state's absorption probabilities sum to 1;
#' the hitting probabilities are the minimal non-negative solution of the
#' linear system in Norris (1997), solved over the states that can reach the
#' target. The diagonal of the hitting matrix uses the return-time
#' convention, the probability of returning to the state in one or more
#' steps, as `markovchain::hittingProbabilities()` does.
#'
#' @param x The chain: a `netobject`, bare `cograph_network` or `tna` model
#'   with transition weights, a numeric transition matrix, or a wide
#'   sequence data frame or character matrix (estimated with
#'   `build_network(method = "relative")`). A `netobject_group` is analysed
#'   network by network. An association network (`cor`, `pcor`, `glasso`,
#'   `ising`, `mgm`) or a matrix with negative entries raises an error of
#'   class `nestimate_not_markov`. For the `print()` and `plot()` methods: a
#'   `chain_structure` object (`print()` also takes a
#'   `chain_structure_group` or a `summary_chain_structure`).
#' @param normalize Whether to rescale rows that do not sum to 1 so that
#'   they do, with a warning, as [passage_time()] does. Default `TRUE`. With
#'   `FALSE`, such a matrix raises an error. A state with no outgoing
#'   transitions (a row of zeros) raises an error either way.
#' @param tol The tolerance below which a transition probability counts as
#'   zero when the support graph is built (it drives `classification`,
#'   `communicating_classes`, `period` and `hitting_probabilities`), and the
#'   tolerance of the reversibility check (detailed balance). Default
#'   `1e-10`. It does **not** govern the absorbing-state test: a state is
#'   absorbing only when `P[i, i]` equals 1 to a fixed internal tolerance of
#'   `.Machine$double.eps^0.5`, so raising `tol` to ignore tiny transition
#'   probabilities never reclassifies a near-deterministic state as
#'   absorbing.
#' @return A `chain_structure` object, or for a `netobject_group` a
#'   `chain_structure_group` (class `c("chain_structure_group", "list")`), a
#'   named list holding one `chain_structure` per network. `print()` shows
#'   the chain-level properties (for a group, one block per network);
#'   `summary()` returns a data frame (see *Summary columns*): for a
#'   `chain_structure`, one row per state, of class
#'   `c("summary_chain_structure", "data.frame")`, carrying the chain-level
#'   flags `is_regular`, `is_irreducible`, `is_aperiodic`, `is_reversible`,
#'   `n_classes` and `absorbing_states` as attributes that its `print()`
#'   method shows as a header above the table; for a
#'   `chain_structure_group`, a plain data frame with one row per group and
#'   state. `plot()` returns a `ggplot` heatmap of the hitting probabilities
#'   (see *Plot colours*). The `print()` methods return their input
#'   invisibly.
#'
#'   The fields of a `chain_structure`:
#'   \describe{
#'     \item{`states`}{Character vector of state names.}
#'     \item{`classification`}{Named character vector, one of
#'       `"absorbing"`, `"recurrent"` or `"transient"` per state.}
#'     \item{`communicating_classes`}{List of state-name vectors, one per
#'       strongly connected component of the support graph.}
#'     \item{`recurrent_classes`}{The communicating classes that are closed
#'       (no transition leaves the class).}
#'     \item{`transient_classes`}{The communicating classes that are not
#'       closed.}
#'     \item{`absorbing_states`}{Character vector of states with
#'       `P[i, i] = 1`, tested to within `.Machine$double.eps^0.5`; `tol`
#'       does not relax this.}
#'     \item{`period`}{Named integer vector: the period of each recurrent
#'       state, `NA` for transient states.}
#'     \item{`is_irreducible`}{`TRUE` when there is exactly one
#'       communicating class.}
#'     \item{`is_aperiodic`}{`TRUE` when every recurrent state has period
#'       1.}
#'     \item{`is_regular`}{`is_irreducible && is_aperiodic`.}
#'     \item{`is_reversible`}{`TRUE` when the chain satisfies detailed
#'       balance against its stationary distribution; `NA` for a chain that
#'       is not irreducible (no unique stationary distribution).}
#'     \item{`hitting_probabilities`}{`n x n` matrix whose `[i, j]` entry is
#'       the probability of ever reaching `j` from `i`, computed over the
#'       same `tol`-thresholded support graph as `classification`, so the
#'       two agree: a state in a closed class has hitting probability 0 to
#'       every state outside it.}
#'     \item{`absorption_probabilities`}{`n_transient x n_closed` matrix
#'       whose `[i, j]` entry is the probability that the chain started in
#'       transient state `i` ends in closed class `j`; each row sums to 1.
#'       An absorbing state is a column named after it; a closed class of
#'       several states is one column named by its states joined with `+`
#'       (for example `B+C`). `NULL` when the chain has no transient state.}
#'     \item{`mean_absorption_time`}{Named numeric vector of the expected
#'       number of steps each transient state takes to enter a closed class,
#'       whichever it is. `NULL` when `absorption_probabilities` is.}
#'     \item{`P`}{The transition matrix analysed (after any row
#'       normalisation).}
#'   }
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' cs  <- chain_structure(net)
#' print(cs)
#' \donttest{
#' summary(cs)
#' }
#'
#' @seealso [passage_time()], [markov_stability()], [build_network()]
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). \emph{Finite Markov Chains}.
#' Springer-Verlag.
#'
#' Norris, J. R. (1997). \emph{Markov Chains}. Cambridge University Press.
#'
#' @export
chain_structure <- function(x, normalize = TRUE, tol = 1e-10) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      chain_structure(net, normalize = normalize, tol = tol)
    })
    class(out) <- c("chain_structure_group", "list")
    return(out)
  }
  P <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }
  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)
  n <- nrow(P)

  # ---- Support graph and SCCs ----
  A <- (P > tol) * 1L
  classes <- .cs_scc(A)

  # A class is closed (recurrent) iff no outgoing edges leave it.
  is_closed <- vapply(classes, function(cl) {
    if (n == length(cl)) return(TRUE)
    out_rows <- A[cl, , drop = FALSE]
    !any(out_rows[, setdiff(seq_len(n), cl), drop = FALSE] > 0)
  }, logical(1))
  recurrent_classes <- classes[is_closed]
  transient_classes <- classes[!is_closed]

  recurrent_idx <- if (length(recurrent_classes) > 0L)
    sort(unlist(recurrent_classes)) else integer(0)
  transient_idx <- if (length(transient_classes) > 0L)
    sort(unlist(transient_classes)) else integer(0)
  # A state is absorbing iff P[i, i] == 1 exactly. The @return contract
  # promises P[i, i] = 1, so this test must NOT use the user-facing `tol`
  # (which is documented only for the support graph + reversibility);
  # otherwise raising `tol` would silently relabel a near-deterministic
  # recurrent state (e.g. P[a,a] = 0.99) as "absorbing". Use a fixed,
  # tiny machine-precision tolerance independent of `tol`.
  absorbing_tol <- .Machine$double.eps^0.5
  absorbing_mask <- vapply(seq_len(n), function(i)
    abs(P[i, i] - 1) < absorbing_tol, logical(1))
  absorbing_idx <- which(absorbing_mask)

  classification <- character(n)
  classification[recurrent_idx] <- "recurrent"
  classification[transient_idx] <- "transient"
  classification[absorbing_idx] <- "absorbing"
  names(classification) <- state_names

  # ---- Period per state (only meaningful for recurrent states) ----
  period <- rep(NA_integer_, n)
  for (cl in recurrent_classes) {
    cl_period <- .cs_class_period(A[cl, cl, drop = FALSE])
    period[cl] <- cl_period
  }
  names(period) <- state_names

  # ---- Chain-level properties ----
  is_irreducible <- length(classes) == 1L
  is_aperiodic <- length(recurrent_classes) > 0L &&
                  all(period[recurrent_idx] == 1L, na.rm = TRUE)
  is_regular <- is_irreducible && is_aperiodic

  is_reversible <- if (is_irreducible) {
    pi <- tryCatch(.mpt_stationary(P), error = function(e) NULL)
    if (is.null(pi)) NA else {
      flux <- pi * P
      max_dev <- max(abs(flux - t(flux)))
      max_dev < tol * max(abs(flux), 1)
    }
  } else NA

  # ---- Hitting probabilities and absorption analysis ----
  # Drive hitting probabilities off the SAME tol-thresholded support graph
  # `A` that classification uses, so the returned object is internally
  # consistent. Without this, .cs_hitting() would rebuild its own support
  # with a hard-coded `P > 0`, and a state classified "absorbing" via a
  # sub-`tol` escape edge could still show hitting probability 1 to other
  # states (a self-contradictory object).
  reach <- .cs_reachability(A)
  diag(reach) <- 1L
  H <- .cs_hitting(P, reach = reach)
  abs_info <- .cs_absorption(P, transient_idx, recurrent_classes, state_names)

  structure(list(
    states                   = state_names,
    classification           = classification,
    communicating_classes    = lapply(classes, function(cl) state_names[cl]),
    recurrent_classes        = lapply(recurrent_classes, function(cl) state_names[cl]),
    transient_classes        = lapply(transient_classes, function(cl) state_names[cl]),
    absorbing_states         = state_names[absorbing_idx],
    period                   = period,
    is_irreducible           = is_irreducible,
    is_aperiodic             = is_aperiodic,
    is_regular               = is_regular,
    is_reversible            = is_reversible,
    hitting_probabilities    = H,
    absorption_probabilities = abs_info$probabilities,
    mean_absorption_time     = abs_info$mean_time,
    P                        = P
  ), class = "chain_structure")
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname chain_structure
#' @param ... Ignored, except by the `print()` methods of a
#'   `chain_structure_group`, which pass it to each network's `print()`, and
#'   of a `summary_chain_structure`, which pass it to `print.data.frame()`.
#' @param show_values For the `plot()` method: whether to print each cell's
#'   probability inside it. Default `TRUE`. Values are never printed for more
#'   than 20 states; set `FALSE` earlier when the labels overlap (from about
#'   10 states).
#' @param digits For the `plot()` method: decimal places of the cell values.
#'   Default `2`.
#' @param object For the `summary()` methods: a `chain_structure` or
#'   `chain_structure_group` object.
#' @export
print.chain_structure <- function(x, ...) {
  rev <- x$is_reversible
  cat(sprintf("Chain structure  [%d states, %d communicating classes]\n",
              length(x$states), length(x$communicating_classes)))
  cat(sprintf(
    "  irreducible: %s   aperiodic: %s   regular: %s   reversible: %s\n",
    x$is_irreducible, x$is_aperiodic, x$is_regular,
    if (is.na(rev)) "NA" else rev))
  cat(sprintf("  recurrent classes: %d   transient classes: %d\n",
              length(x$recurrent_classes), length(x$transient_classes)))
  if (length(x$absorbing_states) > 0L) {
    cat(sprintf("  absorbing states: %s\n",
                paste(x$absorbing_states, collapse = ", ")))
  }
  cat("\nUse summary(x) for the per-state table, plot(x) for the heatmap.\n")
  invisible(x)
}

#' @rdname chain_structure
#' @section Plot colours:
#' The `plot()` method draws the hitting-probability matrix as a heatmap,
#' with rows and columns ordered by communicating class (closed classes
#' first, larger classes before smaller ones) so the block structure is
#' visible. Cell colour runs from white (0) to dark blue (1) and encodes
#' `P(ever reach j | start at i)`; the diagonal uses the return-time
#' convention, `P(return to j in >= 1 steps)`. The state labels on both axes
#' name each state's class, as in `Active (recurrent)`, and the subtitle
#' lists the chain-level properties. A
#' chain that is not irreducible shows zero entries between classes that
#' cannot reach each other; when every state ends in one absorbing state,
#' that state's column is all 1s.
#' @export
plot.chain_structure <- function(x, show_values = TRUE, digits = 2L, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting.", call. = FALSE)
  }

  # Order states: members of each communicating class kept contiguous,
  # classes ordered by [recurrent first, then transient], and within each
  # group classes ordered by size (descending) for stable layout.
  classes <- x$communicating_classes
  is_rec  <- vapply(classes, function(cl) {
    any(x$classification[cl] %in% c("recurrent", "absorbing"))
  }, logical(1))
  rec_classes <- classes[is_rec]
  tra_classes <- classes[!is_rec]
  rec_classes <- rec_classes[order(-vapply(rec_classes, length, integer(1)))]
  tra_classes <- tra_classes[order(-vapply(tra_classes, length, integer(1)))]
  ordered_states <- unlist(c(rec_classes, tra_classes), use.names = FALSE)
  if (length(ordered_states) != length(x$states)) {
    ordered_states <- x$states  # fallback
  }

  H <- x$hitting_probabilities[ordered_states, ordered_states, drop = FALSE]
  cls <- x$classification[ordered_states]

  # Long form for ggplot
  df <- expand.grid(from = ordered_states, to = ordered_states,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$prob <- as.numeric(H)
  df$from <- factor(df$from, levels = ordered_states)
  df$to   <- factor(df$to,   levels = ordered_states)

  # Each axis label names the state's class in text: ggplot2 does not
  # support one colour per tick label, and a colour alone would be the only
  # cue for the class.
  class_labels <- stats::setNames(sprintf("%s (%s)", ordered_states, cls),
                                  ordered_states)

  subtitle <- sprintf(
    "irreducible: %s   aperiodic: %s   regular: %s   reversible: %s   classes: %d",
    x$is_irreducible, x$is_aperiodic, x$is_regular,
    if (is.na(x$is_reversible)) "NA" else x$is_reversible,
    length(x$communicating_classes)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$to, y = .data$from,
                                         fill = .data$prob)) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.4) +
    ggplot2::scale_fill_gradient(low = "white", high = "#08306B",
                                  limits = c(0, 1),
                                  name = "P(reach)",
                                  na.value = "grey90") +
    ggplot2::scale_x_discrete(position = "top", labels = class_labels) +
    ggplot2::scale_y_discrete(limits = rev(ordered_states),
                              labels = class_labels) +
    ggplot2::labs(
      title = "Hitting probabilities",
      subtitle = subtitle,
      x = "to (target state)",
      y = "from (source state)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(angle = 45, hjust = 0,
                                                vjust = 0),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30")
    )

  if (show_values && nrow(df) <= 400L) {
    df$label <- ifelse(df$prob > 0,
                       formatC(df$prob, digits = digits, format = "f"),
                       "")
    p <- p + ggplot2::geom_text(
      data = df,
      ggplot2::aes(label = .data$label,
                    colour = .data$prob > 0.5),
      size = 3, show.legend = FALSE
    ) +
      ggplot2::scale_colour_manual(
        values = c(`TRUE` = "white", `FALSE` = "black")
      )
  }

  p
}

#' @rdname chain_structure
#' @section Summary columns:
#' `summary()` of a `chain_structure` returns one row per state, combining
#' every per-state result: `state`, `classification`, `period`,
#' `persistence` (the diagonal of the transition matrix),
#' `return_probability` (the diagonal of the hitting matrix) and
#' `sojourn_steps` (`1 / (1 - persistence)`, `Inf` for an absorbing state).
#' An irreducible chain adds `stationary_probability`. A chain with transient
#' states adds `absorption_probability` (one closed class) or one
#' `absorbed_in_<class>` column per closed class (named as the columns of
#' `absorption_probabilities`), and `mean_absorption_time` (steps to enter a
#' closed class); these are `NA` for states that are not transient. Probabilities are rounded to 4 decimals and step counts to 2.
#' Columns are ordered for readability: identifiers first, classification
#' second, per-state results last.
#'
#' `summary()` of a `chain_structure_group` stacks these tables into one
#' data frame with a `group` column first, one row per group and state; a
#' column that only some groups have is `NA` for the others.
#' @export
summary.chain_structure <- function(object, ...) {
  state_names <- object$states
  P <- object$P
  H <- object$hitting_probabilities

  persistence <- diag(P)
  sojourn <- ifelse(persistence < 1 - .Machine$double.eps,
                    1 / (1 - persistence), Inf)

  out <- data.frame(
    state              = state_names,
    classification     = unname(object$classification),
    period             = unname(object$period),
    persistence        = round(unname(persistence), 4),
    return_probability = round(unname(diag(H)), 4),
    sojourn_steps      = round(unname(sojourn), 2),
    stringsAsFactors   = FALSE
  )

  if (object$is_irreducible) {
    pi <- tryCatch(.mpt_stationary(P), error = function(e) NULL)
    if (!is.null(pi)) {
      out$stationary_probability <- round(unname(pi), 4)
    }
  }

  if (!is.null(object$absorption_probabilities)) {
    abs_states <- colnames(object$absorption_probabilities)
    abs_probs <- object$absorption_probabilities
    for (a in abs_states) {
      colname <- if (length(abs_states) == 1L) "absorption_probability"
                 else sprintf("absorbed_in_%s", a)
      vals <- rep(NA_real_, length(state_names))
      vals[match(rownames(abs_probs), state_names)] <- round(abs_probs[, a], 4)
      out[[colname]] <- vals
    }
    abs_time <- rep(NA_real_, length(state_names))
    abs_time[match(names(object$mean_absorption_time), state_names)] <-
      round(object$mean_absorption_time, 2)
    out$mean_absorption_time <- abs_time
  }

  rownames(out) <- NULL
  class(out) <- c("summary_chain_structure", "data.frame")
  attr(out, "is_regular") <- object$is_regular
  attr(out, "is_irreducible") <- object$is_irreducible
  attr(out, "is_aperiodic") <- object$is_aperiodic
  attr(out, "is_reversible") <- object$is_reversible
  attr(out, "n_classes") <- length(object$communicating_classes)
  attr(out, "absorbing_states") <- object$absorbing_states
  out
}

#' @rdname chain_structure
#' @export
print.chain_structure_group <- function(x, ...) {
  cat(sprintf("Chain structure -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }
  invisible(x)
}

#' @rdname chain_structure
#' @export
summary.chain_structure_group <- function(object, ...) {
  parts <- lapply(names(object), function(nm) {
    s <- summary(object[[nm]])
    s$group <- nm
    s
  })
  # Union of columns across groups, preserving order
  all_cols <- unique(unlist(lapply(parts, colnames)))
  parts <- lapply(parts, function(d) {
    miss <- setdiff(all_cols, colnames(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, parts)
  out <- out[, c("group", setdiff(all_cols, "group")), drop = FALSE]
  rownames(out) <- NULL
  # The merged frame is a plain comparison table -- drop the
  # per-chain `summary_chain_structure` class (and its attributes) so
  # `print.data.frame` is dispatched, not the per-chain pretty-printer.
  class(out) <- "data.frame"
  out
}

#' @rdname chain_structure
#' @export
print.summary_chain_structure <- function(x, ...) {
  rev <- attr(x, "is_reversible")
  cat(sprintf(
    "Chain structure summary  [%d states, %d classes]\n",
    nrow(x), attr(x, "n_classes")))
  cat(sprintf(
    "  irreducible: %s   aperiodic: %s   regular: %s   reversible: %s\n",
    attr(x, "is_irreducible"), attr(x, "is_aperiodic"),
    attr(x, "is_regular"), if (is.na(rev)) "NA" else rev))
  abs_states <- attr(x, "absorbing_states")
  if (length(abs_states) > 0L) {
    cat(sprintf("  absorbing states: %s\n",
                paste(abs_states, collapse = ", ")))
  }
  cat("\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}
