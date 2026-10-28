# ---- HYPA: Hypothesis Testing for Path Anomalies ----
#
# Implements HYPA (LaRock et al. 2020) for detecting anomalous paths in
# sequential data. Uses a multi-hypergeometric null model on k-th order
# De Bruijn graphs to identify over/under-represented paths.

# ---------------------------------------------------------------------------
# Internal: possible edges and Xi fitting (LaRock et al. 2020; tlarock/hypa)
# ---------------------------------------------------------------------------

#' Possible edges of the order-k De Bruijn graph
#'
#' An edge v -> w is possible when v has out-weight, w has in-weight, w
#' continues v (the last k - 1 states of v are the first k - 1 states of w)
#' and the step from the last state of v to the last state of w is an
#' observed first-order transition. Follows \code{computeXiHigherOrder()} in
#' tlarock/hypa. At order 1 the possible edges are the observed transitions.
#'
#' @param adj Square adjacency matrix of the order-k De Bruijn graph; node
#'   names are k-grams joined by \code{.HON_SEP}.
#' @param first_order Character vector of observed first-order transitions,
#'   each \code{paste(from, to, sep = .HON_SEP)}.
#' @param k Integer order.
#' @return Logical matrix with the dimensions of \code{adj}.
#' @noRd
.hypa_possible_edges <- function(adj, first_order, k) {
  if (k <= 1L) return(adj > 0)
  nodes <- rownames(adj)
  parts <- strsplit(nodes, .HON_SEP, fixed = TRUE)
  src <- nodes[rowSums(adj) > 0]
  tgt <- nodes[colSums(adj) > 0]
  key_of <- function(x, drop_first) {
    vapply(parts[match(x, nodes)], function(p) {
      paste(if (drop_first) p[-1L] else p[-length(p)], collapse = .HON_SEP)
    }, character(1L))
  }
  last_of <- function(x) {
    vapply(parts[match(x, nodes)], function(p) p[length(p)], character(1L))
  }
  pairs <- merge(
    data.frame(from = src, key = key_of(src, TRUE), from_last = last_of(src),
               stringsAsFactors = FALSE),
    data.frame(to = tgt, key = key_of(tgt, FALSE), to_last = last_of(tgt),
               stringsAsFactors = FALSE),
    by = "key"
  )
  pairs <- pairs[paste(pairs$from_last, pairs$to_last, sep = .HON_SEP) %in%
                   first_order, , drop = FALSE]
  mask <- matrix(FALSE, nrow(adj), ncol(adj), dimnames = dimnames(adj))
  mask[cbind(match(pairs$from, nodes), match(pairs$to, nodes))] <- TRUE
  stopifnot("every observed path must be a possible edge" = all(mask[adj > 0]))
  mask
}

#' Fit the propensity matrix Xi
#'
#' Starts from \eqn{\Xi_{vw} = s_{out}(v) s_{in}(w)} on the possible edges
#' and alternately rescales rows and columns (iterative proportional
#' fitting, the continuous form of \code{fitXi()} / \code{xifix_row()} in
#' tlarock/hypa) until the expected out- and in-strengths
#' \eqn{m \sum_w \Xi_{vw} / \sum \Xi} match the observed strengths, so that
#' the null model preserves every node's strength. The observed adjacency
#' is itself a solution supported on the possible edges, so the fit always
#' has a target. No rounding happens inside the loop; the fitted matrix is
#' scaled to \eqn{\sum \Xi = m^2} at the end. The original Python code rounds
#' Xi to integers inside every sweep, which shifts expected counts by up to
#' about 0.002 on corpora of a few thousand paths.
#'
#' @param adj Square adjacency matrix of the De Bruijn graph.
#' @param mask Logical matrix of possible edges (default: observed edges).
#' @param tol Relative tolerance: the fit has converged when the largest
#'   absolute strength error is below \code{tol * m}.
#' @param max_iter Maximum number of row/column sweeps.
#' @return Matrix Xi with the dimensions of \code{adj}, summing to
#'   \eqn{m^2}; attribute \code{"fit"} holds \code{max_error} (largest
#'   absolute strength error of the expected counts), \code{iterations} and
#'   \code{converged}. A classed warning \code{nestimate_hypa_no_converge}
#'   is raised when the tolerance is not reached.
#' @noRd
.hypa_fit_xi <- function(adj, mask = adj > 0, tol = 1e-10, max_iter = 10000L) {
  s_out <- rowSums(adj)
  s_in <- colSums(adj)
  m <- sum(adj)
  xi0 <- outer(s_out, s_in) * mask
  if (m == 0 || sum(xi0) == 0) {
    return(structure(xi0, fit = list(max_error = 0, iterations = 0L,
                                     converged = TRUE)))
  }
  safe_ratio <- function(target, current) {
    ifelse(current > 0, target / current, 0)
  }
  # Expected counts are the fitted matrix itself (it sums to m after each
  # column step), so the error is read directly from its margins.
  strength_error <- function(e) {
    max(abs(rowSums(e) - s_out), abs(colSums(e) - s_in))
  }
  sweep_once <- function(e) {
    e <- e * safe_ratio(s_out, rowSums(e))
    t(t(e) * safe_ratio(s_in, colSums(e)))
  }
  start <- xi0 * m / sum(xi0)
  fit <- Reduce(function(state, i) {
    if (state$done) return(state)
    e <- sweep_once(state$e)
    err <- strength_error(e)
    list(e = e, err = err, iter = state$iter + 1L, done = err < tol * m)
  }, seq_len(max_iter),
  list(e = start, err = strength_error(start), iter = 0L,
       done = strength_error(start) < tol * m))
  if (!fit$done) {
    warning(warningCondition(
      sprintf(paste0("HYPA Xi fit did not converge in %d sweeps ",
                     "(largest strength error %.3g of %d paths)"),
              max_iter, fit$err, as.integer(m)),
      class = "nestimate_hypa_no_converge", call = NULL))
  }
  structure(fit$e * m, fit = list(max_error = fit$err, iterations = fit$iter,
                                  converged = fit$done))
}

# ---------------------------------------------------------------------------
# Internal: Compute HYPA scores
# ---------------------------------------------------------------------------

#' Compute hypergeometric p-values for each edge
#'
#' For each possible edge (v,w), i.e. each cell with Xi > 0, observed or
#' not, with observed weight f (possibly 0), computes:
#'   \code{p_under = P(X <= f)} and \code{p_over = P(X > f) = 1 - p_under}
#'   where \code{X ~ Hypergeometric(N, K, n)},
#'   \code{N = round(sum(Xi))}, \code{K = round(Xi[v,w])}, \code{n = sum(adj)}
#'
#' @param adj Adjacency matrix (edge weights = path frequencies).
#' @param xi Fitted propensity matrix.
#' @return Data frame with path, from, to, observed, expected, ratio,
#'   p_value, p_under, p_over, and anomaly.
#' @noRd
.hypa_compute_scores <- function(adj, xi) {
  n <- nrow(adj)
  nodes <- rownames(adj)

  # Total pool
  xi_total <- sum(xi)
  N_total <- round(xi_total)
  n_draws <- sum(adj)

  # Score every possible edge (Xi > 0), including unobserved paths
  edge_idx <- which(xi > 0, arr.ind = TRUE)
  if (nrow(edge_idx) == 0L) {
    return(data.frame(path = character(0L), from = character(0L),
                      to = character(0L), observed = integer(0L),
                      expected = numeric(0L), ratio = numeric(0L),
                      p_value = numeric(0L), p_under = numeric(0L),
                      p_over = numeric(0L), anomaly = character(0L),
                      stringsAsFactors = FALSE))
  }

  results <- lapply(seq_len(nrow(edge_idx)), function(idx) {
    i <- edge_idx[idx, 1L]
    j <- edge_idx[idx, 2L]
    f_obs <- adj[i, j]
    K <- round(xi[i, j])

    # Clamp parameters to valid range
    K <- min(K, N_total)
    K <- max(K, 0L)
    n_clamp <- min(n_draws, N_total)

    # LaRock et al. (2020): HYPA = P(X <= f); a path is over-represented when
    # HYPA >= 1 - alpha, i.e. when the strict upper tail P(X > f) <= alpha.
    p_under <- stats::phyper(f_obs, K, N_total - K, n_clamp)
    p_over <- stats::phyper(f_obs, K, N_total - K, n_clamp,
                            lower.tail = FALSE)

    # Expected count from the continuous fit; the tails use the integer
    # pool sizes the hypergeometric distribution requires.
    expected <- if (xi_total > 0) n_draws * xi[i, j] / xi_total else 0

    # Reconstruct full path; from = context, to = next state
    from_parts <- strsplit(nodes[i], .HON_SEP, fixed = TRUE)[[1L]]
    to_parts <- strsplit(nodes[j], .HON_SEP, fixed = TRUE)[[1L]]
    next_state <- to_parts[length(to_parts)]
    path <- paste(c(from_parts, next_state), collapse = " -> ")

    ratio <- if (expected > 0) f_obs / expected else Inf

    data.frame(
      path = path,
      from = paste(from_parts, collapse = " -> "),
      to = next_state,
      observed = as.integer(f_obs),
      expected = expected,
      ratio = ratio,
      p_value = p_under,
      p_under = p_under,
      p_over = p_over,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, results)
  # Classify anomalies
  result$anomaly <- ifelse(result$p_under < 0.05, "under",
                           ifelse(result$p_over < 0.05, "over", "normal"))

  result
}

#' Split trajectories at internal missing states
#'
#' A path must not span a missing observation, so each run of non-missing
#' states becomes its own trajectory; runs shorter than two states carry no
#' transition and are dropped.
#' @param trajectories List of character vectors.
#' @return List of character vectors without \code{NA}, each of length >= 2.
#' @noRd
.hypa_split_gaps <- function(trajectories) {
  pieces <- lapply(trajectories, function(tr) {
    present <- !is.na(tr)
    unname(split(tr[present], cumsum(!present)[present]))
  })
  pieces <- unlist(pieces, recursive = FALSE)
  Filter(function(tr) length(tr) >= 2L, pieces)
}

# ---------------------------------------------------------------------------
# Main function: build_hypa
# ---------------------------------------------------------------------------

#' HYPA scores for a single De Bruijn order
#'
#' Internal helper. Builds the order-\code{ord} De Bruijn graph from already
#' parsed \code{trajectories}, fits the propensity matrix, computes
#' hypergeometric tail probabilities, applies multiple-testing correction and
#' classifies each edge as over-/under-represented or normal.
#'
#' @param trajectories List of character vectors (parsed trajectories), as
#'   produced by \code{.hon_parse_input()}.
#' @param ord Integer scalar. Order of the De Bruijn graph for this layer.
#' @param alpha Numeric. Significance threshold for anomaly classification.
#' @param min_count Integer. Minimum observed count for an edge to be eligible
#'   for an anomaly classification.
#' @param p_adjust Character. Multiple-testing correction method (a value from
#'   \code{stats::p.adjust.methods} or \code{"none"}).
#' @return A plain \code{list} (no S3 class) with elements \code{scores},
#'   \code{over}, \code{under}, \code{adjacency}, \code{weights}, \code{xi},
#'   \code{edges}, \code{nodes}, \code{meta}, \code{n_over}, \code{n_under},
#'   \code{n_edges}, \code{order}; or \code{NULL} when the order produces no
#'   edges.
#' @noRd
.hypa_one_order <- function(trajectories, ord, alpha, min_count, p_adjust) {
  kg <- .mogen_count_kgrams(trajectories, ord)
  if (nrow(kg$edges) == 0L) {
    return(NULL)
  }
  nodes <- kg$nodes
  n <- length(nodes)
  adj <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
  idx <- cbind(match(kg$edges$from, nodes), match(kg$edges$to, nodes))
  adj[idx] <- kg$edges$weight

  first_order <- unique(unlist(lapply(trajectories, function(tr) {
    if (length(tr) < 2L) return(character(0L))
    paste(tr[-length(tr)], tr[-1L], sep = .HON_SEP)
  })))
  mask <- .hypa_possible_edges(adj, first_order, ord)
  xi <- .hypa_fit_xi(adj, mask)
  scores <- .hypa_compute_scores(adj, xi)

  if (p_adjust != "none" && nrow(scores) > 0L) {
    scores$p_adjusted_under <- stats::p.adjust(scores$p_under, method = p_adjust)
    scores$p_adjusted_over  <- stats::p.adjust(scores$p_over,  method = p_adjust)
  } else {
    scores$p_adjusted_under <- scores$p_under
    scores$p_adjusted_over  <- scores$p_over
  }

  # Over-representation needs at least min_count observations; a deficit needs
  # at least min_count expected occurrences, so a path that never occurs can
  # still be under-represented.
  is_under <- scores$expected >= min_count & scores$p_adjusted_under < alpha
  # HYPA >= 1 - alpha (LaRock et al. 2020, Sec. 3) <=> P(X > f) <= alpha
  is_over  <- scores$observed >= min_count & scores$p_adjusted_over <= alpha
  scores$anomaly <- ifelse(is_under, "under",
                           ifelse(is_over, "over", "normal"))
  scores$order <- ord

  anom_over  <- scores[scores$anomaly == "over",  , drop = FALSE]
  anom_under <- scores[scores$anomaly == "under", , drop = FALSE]
  normal     <- scores[scores$anomaly == "normal", , drop = FALSE]
  anom_over  <- anom_over[order(-anom_over$ratio),  , drop = FALSE]
  anom_under <- anom_under[order(anom_under$ratio), , drop = FALSE]
  scores <- rbind(anom_over, anom_under, normal)
  rownames(scores) <- NULL

  cg <- .ho_cograph_fields(adj, nodes, method = "hypa")
  list(
    scores = scores, over = anom_over, under = anom_under,
    adjacency = adj, weights = cg$weights, xi = xi,
    edges = cg$edges, nodes = cg$nodes, meta = cg$meta,
    n_over = nrow(anom_over), n_under = nrow(anom_under),
    n_edges = nrow(scores), order = ord
  )
}

#' Detect Path Anomalies via HYPA
#'
#' Constructs a k-th order De Bruijn graph from sequential trajectory data and
#' uses a hypergeometric null model to detect paths with anomalous frequencies.
#' Paths occurring more or less often than expected under the null model are
#' flagged as over- or under-represented.
#'
#' Internal missing states split a trajectory, so no path spans a gap;
#' pieces shorter than two states are dropped. State labels must not
#' contain the control characters \code{"\\x01"} or \code{"\\x02"}, which
#' are used internally as separators (error class
#' \code{nestimate_hypa_bad_label}). The propensity matrix is fitted by
#' continuous iterative proportional fitting; when the fit does not reach
#' its tolerance a warning of class \code{nestimate_hypa_no_converge} is
#' raised. The original Python implementation rounds the propensities
#' inside every sweep, so its expected counts can differ from these by a few
#' thousandths.
#'
#' @param data A data.frame (rows = trajectories), list of character vectors,
#'   \code{tna} object, or \code{netobject} with sequence data. For
#'   \code{tna}/\code{netobject}, numeric state IDs are automatically
#'   converted to label names.
#' @param order Integer scalar or integer vector. Order(s) of the De Bruijn
#'   graph (default \code{2L}). An order of \code{k} tests paths of
#'   \code{k + 1} states (\code{k} transitions), each conditioned on the
#'   \code{k}-state context that precedes its last state. Orders must be
#'   unique. When a vector is supplied, one De Bruijn layer
#'   is built per order and the per-order results are stored in
#'   \code{$by_order} (named by order). The orders are sorted ascending
#'   internally, so the \code{cograph_network} slots
#'   (\code{$weights}, \code{$edges}, \code{$adjacency}, \code{$xi},
#'   \code{$nodes}, \code{$meta}) always describe the network of the
#'   \emph{lowest order that produced a layer} (a requested order with no
#'   edges is dropped from \code{$by_order}, \code{$order} and \code{$k}),
#'   regardless of the order in which the vector is given; \code{$scores}
#'   aggregates every built order.
#' @param alpha Numeric in `(0, 0.5)`. Significance threshold for anomaly
#'   classification (default 0.05). A path is labelled \code{"under"} when its
#'   adjusted lower-tail p-value \code{p_adjusted_under} is below
#'   \code{alpha}, \code{"over"} when its adjusted upper-tail p-value
#'   \code{p_adjusted_over} is below \code{alpha}, and \code{"normal"}
#'   otherwise (the under-representation test is applied first).
#' @param min_count Integer. Minimum count for a path to be classified as
#'   anomalous (default 5). An over-represented path must be observed at
#'   least \code{min_count} times; an under-represented path must be
#'   expected at least \code{min_count} times under the null, so a path
#'   that never occurs in the data can still be flagged as
#'   under-represented. Other paths are classified as \code{"normal"}.
#' @param p_adjust Character. Optional multiple testing correction of
#'   p-values. Default \code{"none"}: paths are classified against the
#'   threshold \code{alpha} directly, as in LaRock et al. (2020), which
#'   applies no correction. Accepts any method from
#'   \code{\link[stats]{p.adjust.methods}}, e.g. \code{"BH"} for
#'   Benjamini-Hochberg FDR control. Under- and over-representation
#'   p-values are adjusted separately, within each order, over every
#'   possible path. With small counts the paper's upper-tail score
#'   \code{P(X > f)} rejects more often than \code{alpha} under the null, so
#'   neither the uncorrected nor the corrected classification carries an
#'   exact type-I or FDR guarantee; \code{min_count} limits the problem.
#' @param k Deprecated. Former name of \code{order}; if supplied it overrides
#'   \code{order} and emits a deprecation message. Use \code{order} instead.
#' @return An object of class \code{c("net_hypa", "cograph_network")} with
#'   components:
#'   \describe{
#'     \item{scores}{Data frame with path, from, to, observed, expected,
#'       ratio, p_value, p_under, p_over, p_adjusted_under,
#'       p_adjusted_over, anomaly, order
#'       columns (one block of rows per requested order). One row per
#'       possible path of the De Bruijn graph: every observed path plus every
#'       unobserved path that the lower-order transitions permit
#'       (\code{observed = 0}). \code{expected} comes from a propensity
#'       matrix fitted so that the expected counts preserve the observed in-
#'       and out-strength of every De Bruijn node (LaRock et al. 2020); at
#'       order 2 this equals the first-order Markov expectation
#'       \eqn{n(a b \cdot) n(\cdot b c) / n(\cdot b \cdot)}. The \code{path}
#'       column shows the full state sequence
#'       (e.g., "A -> B -> C"); \code{from} is the context (conditioning
#'       states); \code{to} is the next state; \code{ratio} is
#'       observed / expected; \code{p_value} is retained as an alias for
#'       \code{p_under}, the raw lower-tail hypergeometric CDF value;
#'       \code{p_over} is the strict upper-tail probability
#'       \code{P(X > observed) = 1 - p_under}, so that a path is
#'       over-represented when its HYPA score is at least
#'       \code{1 - alpha}, i.e. \code{p_over <= alpha} (LaRock et al.,
#'       2020); under-represented when \code{p_under < alpha}; \code{p_adjusted_under} and
#'       \code{p_adjusted_over} are the corrected p-values for under- and
#'       over-representation tests respectively.}
#'     \item{ho_edges}{Alias for \code{scores} (all orders, arrow notation).}
#'     \item{over}{Subset of \code{scores} classified as over-represented.}
#'     \item{under}{Subset of \code{scores} classified as under-represented.}
#'     \item{adjacency}{Weighted adjacency matrix of the lowest-order De
#'       Bruijn graph.}
#'     \item{weights}{cograph weight matrix of the lowest-order graph.}
#'     \item{xi}{Fitted propensity matrix of the lowest-order graph.}
#'     \item{edges}{cograph edge data.frame of the lowest-order graph.}
#'     \item{by_order}{Named list of per-order result lists.}
#'     \item{order}{Integer vector of orders actually built (sorted ascending).}
#'     \item{k}{Back-compatibility alias for \code{order}.}
#'     \item{alpha}{Significance threshold used.}
#'     \item{p_adjust}{Multiple testing correction method used.}
#'     \item{n_anomalous}{Number of anomalous paths detected (all orders).}
#'     \item{n_over}{Number of over-represented paths (all orders).}
#'     \item{n_under}{Number of under-represented paths (all orders).}
#'     \item{n_edges}{Total number of tested (possible) paths, all orders.}
#'     \item{nodes}{data.frame (\code{id}, \code{label}, \code{name}) of the
#'       lowest-order De Bruijn graph nodes (arrow notation).}
#'     \item{directed}{Logical. Always \code{TRUE}.}
#'     \item{meta}{cograph meta list of the lowest-order graph.}
#'     \item{node_groups}{Always \code{NULL}.}
#'   }
#'
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' & Schweitzer, F. (2020). HYPA: Efficient Detection of Path Anomalies in
#' Time Series Data on Networks. \emph{SDM 2020}, 460-468.
#'
#' @examples
#' seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
#' hyp <- build_hypa(seqs, order = 2)
#'
#' \donttest{
#' trajs <- list(c("A","B","C"), c("A","B","C"), c("A","B","C"),
#'               c("A","B","D"), c("C","B","D"), c("C","B","A"))
#' h <- build_hypa(trajs, order = 2)
#' print(h)
#' }
#'
#' @export
build_hypa <- function(data, order = 2L, alpha = 0.05, min_count = 5L,
                       p_adjust = "none", k = NULL) {
  if (!is.null(k)) {
    .Deprecated(msg = "'k' is deprecated; use 'order' instead.")
    order <- k
  }
  data <- .coerce_sequence_input(data)

  valid_methods <- c(stats::p.adjust.methods, "none")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% valid_methods) {
    stop(sprintf("'p_adjust' must be one of: %s",
                 paste(valid_methods, collapse = ", ")), call. = FALSE)
  }
  stopifnot(
    "'data' must be a data.frame or list" =
      is.data.frame(data) || is.list(data),
    "'order' must contain only whole numbers" =
      is.numeric(order) && length(order) >= 1L && all(is.finite(order)) &&
      all(order == round(order)),
    "'order' must contain only integers >= 1" = all(order >= 1L),
    "'order' must not contain duplicates" = !anyDuplicated(order),
    "'alpha' must be a single number in (0, 0.5)" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
      alpha > 0 && alpha < 0.5,
    "'min_count' must be a single whole number >= 1" =
      is.numeric(min_count) && length(min_count) == 1L &&
      is.finite(min_count) && min_count == round(min_count) && min_count >= 1
  )
  # Sort ascending so the primary cograph slot (built from per_order[[1L]])
  # is deterministically the lowest order that produces a layer, independent
  # of how the `order` vector is arranged (see @param order).
  order <- sort(as.integer(order))
  min_count <- as.integer(min_count)

  trajectories <- .hypa_split_gaps(.hon_parse_input(data, collapse_repeats = FALSE))
  states <- unlist(trajectories, use.names = FALSE)
  if (any(grepl("[\x01\x02]", states, useBytes = TRUE))) {
    stop(errorCondition(
      "state labels must not contain the control characters \\x01 or \\x02",
      class = "nestimate_hypa_bad_label", call = NULL))
  }
  if (length(trajectories) == 0L) {
    stop("No valid trajectories (each must have at least 2 states)")
  }

  per_order <- lapply(order, function(o) {
    .hypa_one_order(trajectories, o, alpha, min_count, p_adjust)
  })
  names(per_order) <- as.character(order)
  per_order <- Filter(Negate(is.null), per_order)
  if (length(per_order) == 0L) {
    stop("No edges at any requested order (paths too short or too few)")
  }

  scores_all <- do.call(rbind, lapply(per_order, `[[`, "scores"))
  over_all   <- do.call(rbind, lapply(per_order, `[[`, "over"))
  under_all  <- do.call(rbind, lapply(per_order, `[[`, "under"))
  rownames(scores_all) <- NULL
  if (!is.null(over_all))  rownames(over_all)  <- NULL
  if (!is.null(under_all)) rownames(under_all) <- NULL

  # Orders that actually produced a layer. Filter() above dropped any
  # requested order with no edges (e.g. an order higher than the sequence
  # lengths allow), so $order/$k must report only the built layers to stay
  # consistent with $by_order and $scores$order (Codex review).
  built_order <- as.integer(names(per_order))
  # Primary cograph slot = network of the lowest built order.
  # `order` was sorted ascending above, so per_order[[1L]] is the lowest.
  primary <- per_order[[1L]]
  result <- list(
    scores = scores_all,
    ho_edges = scores_all,
    edges = primary$edges,
    over = over_all,
    under = under_all,
    adjacency = primary$adjacency,
    weights = primary$weights,
    xi = primary$xi,
    by_order = per_order,
    order = built_order,
    k = built_order,  # back-compat alias; prefer $order
    alpha = alpha,
    p_adjust = p_adjust,
    n_anomalous = sum(vapply(per_order, function(x) x$n_over + x$n_under, integer(1))),
    n_over  = sum(vapply(per_order, `[[`, integer(1), "n_over")),
    n_under = sum(vapply(per_order, `[[`, integer(1), "n_under")),
    n_edges = sum(vapply(per_order, `[[`, integer(1), "n_edges")),
    nodes = primary$nodes,
    directed = TRUE,
    meta = primary$meta,
    node_groups = NULL
  )
  class(result) <- c("net_hypa", "cograph_network")
  result
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname build_hypa
#' @param x For the \code{print()} method: an object of class \code{net_hypa}.
#' @param ... In \code{print.net_hypa()} and \code{summary.net_hypa()}: Additional arguments (ignored).
#' @param object For the \code{summary()} method: an object of class \code{net_hypa}.
#' @param n Integer. Maximum number of paths to display per category (default: 10).
#' @param type Character. Which anomalies to show: \code{"all"} (default), \code{"over"}, or \code{"under"}.
#' @param order_by Character. Ranking used within each anomaly direction: \code{"sig"} (default) ranks by the active tail probability, \code{"freq"} (or its alias \code{"frequency"}) by observed count, \code{"ratio"} by observed/expected ratio, and \code{"path"} alphabetically.
#' @return In \code{print.net_hypa()}: The input object, invisibly.
#' @return In \code{summary.net_hypa()}: A data frame of the reported anomalies, at most \code{n} rows per direction, with columns \code{order} (the De Bruijn order the path was found at), \code{path}, \code{observed}, \code{expected}, \code{ratio}, \code{p_tail} (the raw tail probability in the reported direction: \code{p_over} for over-represented paths, \code{p_under} for under-represented ones) and \code{direction} (\code{"over"}/\code{"under"}). Returned invisibly; the summary text and the top-\code{n} tables are printed as a side effect. When no anomalies were detected, a zero-row data frame with the same columns except \code{order} is returned.
#' @export
print.net_hypa <- function(x, ...) {
  ord_str <- paste(x$order, collapse = ", ")
  cat("HYPA: Path Anomaly Detection\n")
  cat(sprintf("  Order(s):     %s\n", ord_str))
  cat(sprintf("  Edges:        %d\n", x$n_edges))
  cat(sprintf("  Anomalous:    %d (alpha=%.2f, p_adjust=%s)\n",
              x$n_anomalous, x$alpha, x$p_adjust %||% "none"))
  cat(sprintf("    Over-repr:  %d\n", x$n_over))
  cat(sprintf("    Under-repr: %d\n", x$n_under))
  if (length(x$order) > 1L) {
    per <- vapply(x$by_order, function(po)
      sprintf("order %d: %d edges (%d over, %d under)",
              po$order, po$n_edges, po$n_over, po$n_under),
      character(1))
    cat("  Per-order:\n")
    for (line in per) cat(sprintf("    %s\n", line))
  }
  invisible(x)
}

#' @rdname build_hypa
#' @export
summary.net_hypa <- function(object, n = 10L,
                             type = c("all", "over", "under"),
                             order_by = c("sig", "freq", "frequency",
                                          "ratio", "path"),
                             ...) {
  type <- match.arg(type)
  order_by <- match.arg(order_by)
  if (order_by == "frequency") order_by <- "freq"

  reorder_df <- function(df, dir) {
    if (is.null(df) || nrow(df) == 0L) return(df)
    # p_value is the hypergeometric CDF P(X <= obs): close to 1 for
    # over-represented paths, close to 0 for under-represented. So "by
    # sig" picks the tail extremeness in each direction.
    idx <- switch(order_by,
      ratio = if (dir == "over") order(-df$ratio) else order(df$ratio),
      sig   = if (dir == "over") order(df$p_over) else order(df$p_under),
      freq  = order(-df$observed),
      path  = order(df$path)
    )
    df[idx, , drop = FALSE]
  }

  over_sorted  <- reorder_df(object$over,  "over")
  under_sorted <- reorder_df(object$under, "under")

  cat("HYPA Summary\n\n")
  cat(sprintf("  Order(s): %s | Edges: %d\n",
              paste(object$order, collapse = ", "), object$n_edges))
  cat(sprintf("  Alpha: %.2f | p_adjust: %s\n",
              object$alpha, object$p_adjust %||% "none"))
  cat(sprintf("  Anomalous: %d (over: %d, under: %d) | order_by: %s\n\n",
              object$n_anomalous, object$n_over, object$n_under, order_by))

  display_df <- function(df, dir) {
    if (is.null(df) || nrow(df) == 0L) return(df)
    out <- df
    out$p_tail <- if (dir == "over") out$p_over else out$p_under
    out
  }

  show_cols <- c("order", "path", "observed", "expected", "ratio", "p_tail")
  if (length(object$order) == 1L) show_cols <- setdiff(show_cols, "order")

  if (type %in% c("all", "over") && object$n_over > 0L) {
    cat("  Over-represented (top", min(n, object$n_over), "):\n")
    over_display <- display_df(over_sorted, "over")
    print(utils::head(over_display[, show_cols, drop = FALSE], n),
          row.names = FALSE)
    cat("\n")
  }

  if (type %in% c("all", "under") && object$n_under > 0L) {
    cat("  Under-represented (top", min(n, object$n_under), "):\n")
    under_display <- display_df(under_sorted, "under")
    print(utils::head(under_display[, show_cols, drop = FALSE], n),
          row.names = FALSE)
    cat("\n")
  }

  if (object$n_anomalous == 0L) {
    cat("  No anomalous paths detected.\n")
  } else if (type == "over" && object$n_over == 0L) {
    cat("  No over-represented paths.\n")
  } else if (type == "under" && object$n_under == 0L) {
    cat("  No under-represented paths.\n")
  }

  pick <- function(df, dir) {
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    out <- display_df(df, dir)
    cols <- c("order", "path", "observed", "expected", "ratio", "p_tail")
    out <- out[, intersect(cols, names(out)), drop = FALSE]
    out$direction <- dir
    utils::head(out, n)
  }
  parts <- switch(type,
    all   = list(pick(over_sorted, "over"), pick(under_sorted, "under")),
    over  = list(pick(over_sorted, "over")),
    under = list(pick(under_sorted, "under"))
  )
  combined <- do.call(rbind, Filter(Negate(is.null), parts))
  if (is.null(combined)) {
    combined <- data.frame(path = character(0L),
                           observed = numeric(0L),
                           expected = numeric(0L),
                           ratio = numeric(0L),
                           p_tail = numeric(0L),
                           direction = character(0L),
                           stringsAsFactors = FALSE)
  } else {
    if (type == "all" && nrow(combined) > 0L) {
      idx <- switch(order_by,
        sig   = order(combined$p_tail),
        ratio = order(-abs(log(combined$ratio))),
        freq  = order(-combined$observed),
        path  = order(combined$path)
      )
      combined <- combined[idx, , drop = FALSE]
    }
    row.names(combined) <- NULL
  }
  # The tables above are already printed; return invisibly so they are not
  # printed a second time at top level.
  invisible(combined)
}
