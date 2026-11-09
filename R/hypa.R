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
  # HYPA >= 1 - alpha (LaRock et al. 2020, Sec. 3) <=> P(X > f) <= alpha.
  # As in the original implementation (tlarock/hypa), a path whose order has
  # a single possible path has P(X > f) = 0 and is "over" even when its
  # observed count equals its expected one; pinned by test-hypa-oracles.R.
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
#' @description
#' `build_hypa()` tests which paths in a set of sequences occur more or less
#' often than expected, with HYPA (LaRock et al., 2020). For an order `k` it
#' builds the order-`k` De Bruijn graph, whose nodes are runs of `k` states
#' and whose edges are paths of `k + 1` states, and compares each path's
#' count with a hypergeometric null model that keeps every node's observed
#' in- and out-strength. Paths that are rarer than the null model allows are
#' labelled under-represented, paths that are more common over-represented,
#' and the rest normal. Unobserved paths that the observed transitions make
#' possible are tested too, so a path that never occurs can be
#' under-represented. Several orders can be tested in one call.
#'
#' A path is over-represented when its last state follows its `k`-state
#' context more often than the null model predicts from the strengths of the
#' De Bruijn nodes alone (at order 2, from the first-order transitions).
#' Each path is tested on its own at `alpha`, so with many possible paths
#' some are flagged by chance: `p_adjust` corrects for that, and with small
#' counts the upper-tail test rejects more often than `alpha` (see
#' `p_adjust`). Whether the sequences need a higher order overall is a
#' separate question, tested by [markov_order_test()]; [build_mogen()]
#' selects an order and [build_hon()] builds a network of variable order.
#'
#' @details
#' Internal missing states split a trajectory, so no path spans a gap, and
#' pieces shorter than two states are dropped. State labels must not
#' contain the control characters `"\x01"` or `"\x02"`, which are used
#' internally as separators; such labels raise an error of class
#' `nestimate_hypa_bad_label`.
#'
#' The possible paths of order `k` link each node with an observed outgoing
#' path to every node with an observed incoming path that continues it (its
#' last `k - 1` states are the other's first) through an observed
#' first-order transition; at order 1 they are the observed transitions.
#' The propensity matrix of the null model is fitted
#' on them by continuous iterative proportional fitting; when the fit does
#' not reach its tolerance, a warning of class `nestimate_hypa_no_converge`
#' is raised. The original Python implementation rounds the propensities
#' inside every sweep, so its expected counts can differ from these by a
#' few thousandths.
#'
#' As in the original implementation, when an order has only one possible
#' path, that path's upper-tail p-value is 0, so it is labelled
#' over-represented once it is observed `min_count` times, even when its
#' count equals the expected one.
#'
#' @param data The sequences: a wide data frame (one row per trajectory,
#'   position columns only, named as a numbered run), a character matrix
#'   with the same layout, a list of character vectors, a `tna` model, or a
#'   `netobject` built from sequences. For a
#'   `tna` model or a `netobject`, numeric state codes are converted to
#'   their labels. A frame with other columns beside the positions, a frame
#'   or matrix of only 0 and 1, a `netobject_group`, and a one-hot or
#'   association network raise an error of class `nestimate_bad_input`.
#' @param order The order or orders of the De Bruijn graph, whole numbers of
#'   at least 1 without duplicates. Default `2`. Order `k` tests paths of
#'   `k + 1` states (`k` transitions), each conditioned on the `k`-state
#'   context before its last state. With several orders, one layer is built
#'   per order, in ascending order whatever the order given; an order with
#'   no paths (longer than the sequences allow) is dropped without a
#'   message, and an error is raised when none is left. The network fields
#'   (`weights`, `edges`, `adjacency`, `xi`, `nodes`, `meta`) describe the
#'   lowest order built, while `scores` holds every order.
#' @param alpha The significance threshold for the labels, a number in
#'   (0, 0.5). Default `0.05`. A path is under-represented when its adjusted
#'   lower-tail p-value (`p_adjusted_under`) is below `alpha`, otherwise
#'   over-represented when its adjusted upper-tail p-value
#'   (`p_adjusted_over`) is at most `alpha`, and normal otherwise.
#' @param min_count The smallest count a path needs to be labelled, a whole
#'   number of at least 1. Default `5`. An over-represented path must be
#'   observed at least `min_count` times, an under-represented path expected
#'   at least `min_count` times under the null model, so a path that never
#'   occurs can still be under-represented. Paths below it are normal.
#' @param p_adjust The correction of the p-values for multiple testing: any
#'   of [stats::p.adjust.methods], such as `"BH"` for Benjamini-Hochberg
#'   false discovery rate control. Default `"none"`, which compares the raw
#'   p-values with `alpha`, as LaRock et al. (2020) do. Under- and
#'   over-representation p-values are adjusted separately, within each
#'   order, over every possible path. With small counts the paper's
#'   upper-tail score `P(X > f)` rejects more often than `alpha` under the
#'   null, so neither the raw nor the corrected labels carry an exact type-I
#'   error or false discovery rate guarantee; `min_count` limits the
#'   problem.
#' @param k The former name of `order`, deprecated: when given, it replaces
#'   `order`, with a warning of class `deprecatedWarning`. Default `NULL`.
#' @return A `net_hypa` object (class `c("net_hypa", "cograph_network")`).
#'   `print()` reports the orders, the number of paths tested and how many
#'   are over- and under-represented (per order when there are several);
#'   [summary()] prints and returns the strongest anomalies; [pathways()]
#'   gives the anomalous paths as pathway strings for
#'   `cograph::plot_simplicial()`. The object's fields:
#'   \describe{
#'     \item{scores}{A data frame with one row per possible path of each
#'       order: every observed path and every unobserved one the observed
#'       transitions permit (`observed = 0`). Within each order the
#'       over-represented paths come first (by decreasing ratio), then the
#'       under-represented (by increasing ratio), then the rest. Columns:
#'       `path`, the full state sequence (`"A -> B -> C"`); `from`, the
#'       context (the conditioning states); `to`, the next state;
#'       `observed`; `expected`, the count under the null model, whose
#'       propensity matrix keeps the observed in- and out-strength of every
#'       De Bruijn node (at order 2 this equals the first-order Markov
#'       expectation \eqn{n(a b \cdot) n(\cdot b c) /
#'       n(\cdot b \cdot)}{n(a,b,.) n(.,b,c) / n(.,b,.)});
#'       `ratio`, observed / expected; `p_under`, the hypergeometric
#'       lower-tail probability `P(X <= observed)`, and `p_value`, the same
#'       value under its older name; `p_over`, the strict upper tail
#'       `P(X > observed) = 1 - p_under`, so that `p_over <= alpha` is the
#'       paper's HYPA score of at least `1 - alpha`; `p_adjusted_under` and
#'       `p_adjusted_over`, the same after `p_adjust`; `anomaly` (`"over"`,
#'       `"under"` or `"normal"`); and `order`.}
#'     \item{ho_edges}{The same table as `scores`, under the name the other
#'       higher-order builders use.}
#'     \item{over, under}{The rows of `scores` labelled over- or
#'       under-represented.}
#'     \item{adjacency}{The weighted adjacency matrix (path counts) of the
#'       lowest-order De Bruijn graph.}
#'     \item{weights}{The cograph weight matrix of the lowest-order graph.}
#'     \item{xi}{The fitted propensity matrix of the lowest-order graph.}
#'     \item{edges}{The cograph edge data frame of the lowest-order graph.}
#'     \item{by_order}{A list with the results of each order built, named
#'       by order.}
#'     \item{order}{The orders built, ascending.}
#'     \item{k}{The same as `order`, kept for backward compatibility.}
#'     \item{alpha}{The significance threshold used.}
#'     \item{p_adjust}{The correction used.}
#'     \item{n_anomalous, n_over, n_under}{The number of anomalous, over- and
#'       under-represented paths, over all orders.}
#'     \item{n_edges}{The number of possible paths tested, over all orders.}
#'     \item{nodes}{A data frame (`id`, `label`, `name`) of the nodes of the
#'       lowest-order De Bruijn graph; `label` and `name` are the node's `k`
#'       states joined by the internal separator `"\x01"`, as are the
#'       dimnames of `adjacency`, `weights` and `xi`.}
#'     \item{directed}{Always `TRUE`.}
#'     \item{meta}{The cograph meta list of the lowest-order graph.}
#'     \item{node_groups}{Always `NULL`.}
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
  .refuse_onehot_network(data, "build_hypa()")
  .refuse_binary_matrix(data, "build_hypa()")
  if (is.data.frame(data)) {
    .refuse_undeclared_roles(data, "build_hypa()", via_network = TRUE,
                             whole_frame = TRUE)
  }
  data <- .coerce_sequence_input(data)

  valid_methods <- unique(c(stats::p.adjust.methods, "none"))
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% valid_methods) {
    stop(errorCondition(sprintf("'p_adjust' must be one of: %s",
                 paste(valid_methods, collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }
  .assert(
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
    stop(errorCondition(
      "No valid trajectories (each must have at least 2 states)",
      class = "nestimate_bad_input", call = NULL))
  }

  per_order <- lapply(order, function(o) {
    .hypa_one_order(trajectories, o, alpha, min_count, p_adjust)
  })
  names(per_order) <- as.character(order)
  per_order <- Filter(Negate(is.null), per_order)
  if (length(per_order) == 0L) {
    stop(errorCondition(
      "No edges at any requested order (paths too short or too few)",
      class = "nestimate_bad_input", call = NULL))
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
#' @param x For the `print()` method: a `net_hypa`.
#' @param ... For the `print()` and `summary()` methods: ignored.
#' @param object For the `summary()` method: a `net_hypa`.
#' @param n For the `summary()` method: the most paths to show and return
#'   per direction. Default `10`.
#' @param type For the `summary()` method: which anomalies to show, `"all"`
#'   (default), `"over"` or `"under"`.
#' @param order_by For the `summary()` method: how the paths are ranked
#'   within each direction before the top `n` are kept. `"sig"` (default)
#'   ranks by the tail probability of that direction, `"freq"` (or
#'   `"frequency"`) by observed count, `"ratio"` by observed / expected
#'   (highest first for over-represented paths, lowest first for
#'   under-represented ones), and `"path"` alphabetically. With `type =
#'   "all"` the returned rows of both directions are ranked together the
#'   same way, by how far the ratio is from 1 on the log scale
#'   (`abs(log(ratio))`) for `"ratio"`.
#' @return `print()` returns its input invisibly. `summary()` prints the
#'   settings, the anomaly counts and the top `n` paths of each direction,
#'   and returns invisibly a data frame of those paths, one row per path:
#'   `order` (the De Bruijn order the path was tested at), `path`,
#'   `observed`, `expected`, `ratio`, `p_tail` (the raw tail probability of
#'   its direction, `p_over` for over-represented paths and `p_under` for
#'   under-represented ones) and `direction` (`"over"` or `"under"`). When
#'   nothing is anomalous it is a zero-row data frame with the same columns
#'   except `order`.
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
