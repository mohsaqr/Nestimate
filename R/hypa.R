# ---- HYPA: Hypothesis Testing for Path Anomalies ----
#
# Implements HYPA (LaRock et al. 2020) for detecting anomalous paths in
# sequential data. Uses a multi-hypergeometric null model on k-th order
# De Bruijn graphs to identify over/under-represented paths.

# ---------------------------------------------------------------------------
# Internal: Fit Xi matrix (iterative proportional fitting)
# ---------------------------------------------------------------------------

#' Compute propensity matrix Xi from node strengths
#'
#' Xi_{vw} = s_out(v) * s_in(w) for edges present in the De Bruijn graph.
#' The product of strengths gives N = sum(Xi) >> m = sum(adj), ensuring a
#' non-degenerate hypergeometric null model. This follows the HYPA null
#' model where edge propensity is proportional to the product of endpoint
#' weighted degrees.
#'
#' @param adj Square adjacency matrix of the De Bruijn graph.
#' @return Matrix Xi with same dimensions as adj (sum(Xi) >> sum(adj)).
#' @noRd
.hypa_fit_xi <- function(adj) {
  out_strength <- rowSums(adj)
  in_strength <- colSums(adj)
  mask <- adj > 0

  # Xi_{vw} = s_out(v) * s_in(w) where edges exist in the De Bruijn graph
  outer(out_strength, in_strength) * mask
}

# ---------------------------------------------------------------------------
# Internal: Compute HYPA scores
# ---------------------------------------------------------------------------

#' Compute hypergeometric p-values for each edge
#'
#' For each edge (v,w) with observed weight f, computes:
#'   \code{HYPA(v,w) = P(X <= f)} where \code{X ~ Hypergeometric(N, K, n)},
#'   \code{N = round(sum(Xi))}, \code{K = round(Xi[v,w])}, \code{n = sum(adj)}
#'
#' @param adj Adjacency matrix (edge weights = path frequencies).
#' @param xi Fitted propensity matrix.
#' @return Data frame with from, to, observed, expected, hypa_score, anomaly.
#' @noRd
.hypa_compute_scores <- function(adj, xi) {
  n <- nrow(adj)
  nodes <- rownames(adj)

  # Total pool
  N_total <- round(sum(xi))
  n_draws <- sum(adj)

  # Collect edges
  edge_idx <- which(adj > 0, arr.ind = TRUE)
  if (nrow(edge_idx) == 0L) {
    return(data.frame(path = character(0L), from = character(0L),
                      to = character(0L), observed = integer(0L),
                      expected = numeric(0L), ratio = numeric(0L),
                      hypa_score = numeric(0L), anomaly = character(0L),
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

    # Hypergeometric CDF: P(X <= f_obs)
    hypa_score <- stats::phyper(f_obs, K, N_total - K, n_clamp)

    # Expected value: n * K / N
    expected <- if (N_total > 0) n_draws * K / N_total else 0

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
      hypa_score = hypa_score,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, results)
  # Classify anomalies
  result$anomaly <- vapply(result$hypa_score, function(s) {
    if (s < 0.05) "under" else if (s > 0.95) "over" else "normal"
  }, character(1L))

  result
}

# ---------------------------------------------------------------------------
# Main function: build_hypa
# ---------------------------------------------------------------------------

#' Detect Path Anomalies via HYPA
#'
#' Constructs a k-th order De Bruijn graph from sequential trajectory data and
#' uses a hypergeometric null model to detect paths with anomalous frequencies.
#' Paths occurring more or less often than expected under the null model are
#' flagged as over- or under-represented.
#'
#' @param data A data.frame (rows = trajectories) or list of character vectors.
#' @param k Integer. Order of the De Bruijn graph (default 2). Detects
#'   anomalies in paths of length k.
#' @param alpha Numeric. Significance threshold for anomaly classification
#'   (default 0.05). Paths with HYPA score < alpha are under-represented;
#'   paths with score > 1-alpha are over-represented.
#' @param min_count Integer. Minimum observed count for a path to be
#'   classified as anomalous (default 2). Paths with fewer observations
#'   are always classified as \code{"normal"} regardless of their
#'   HYPA score, since single occurrences are unreliable.
#' @return An object of class \code{net_hypa} with components:
#'   \describe{
#'     \item{scores}{Data frame with path, from, to, observed, expected,
#'       ratio, hypa_score, anomaly columns. The \code{path} column shows
#'       the full state sequence (e.g., "A -> B -> C"); \code{from} is the
#'       context (conditioning states); \code{to} is the next state;
#'       \code{ratio} is observed / expected.}
#'     \item{adjacency}{Weighted adjacency matrix of the De Bruijn graph.}
#'     \item{xi}{Fitted propensity matrix.}
#'     \item{k}{Order of the De Bruijn graph.}
#'     \item{alpha}{Significance threshold used.}
#'     \item{n_anomalous}{Number of anomalous paths detected.}
#'     \item{n_over}{Number of over-represented paths.}
#'     \item{n_under}{Number of under-represented paths.}
#'     \item{n_edges}{Total number of edges.}
#'     \item{nodes}{Node names in the De Bruijn graph.}
#'   }
#'
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' & Schweitzer, F. (2020). HYPA: Efficient Detection of Path Anomalies in
#' Time Series Data on Networks. \emph{SDM 2020}, 460–468.
#'
#' @examples
#' \dontrun{
#' trajs <- list(c("A","B","C"), c("A","B","C"), c("A","B","C"),
#'               c("A","B","D"), c("C","B","D"), c("C","B","A"))
#' h <- build_hypa(trajs, k = 2)
#' print(h)
#' }
#'
#' @export
build_hypa <- function(data, k = 3L, alpha = 0.05, min_count = 5L) {
  k <- as.integer(k)
  min_count <- as.integer(min_count)
  stopifnot(
    "'data' must be a data.frame or list" =
      is.data.frame(data) || is.list(data),
    "'k' must be >= 1" = k >= 1L,
    "'alpha' must be in (0, 0.5)" = alpha > 0 && alpha < 0.5,
    "'min_count' must be >= 1" = min_count >= 1L
  )

  trajectories <- .hon_parse_input(data, collapse_repeats = FALSE)
  if (length(trajectories) == 0L) {
    stop("No valid trajectories (each must have at least 2 states)")
  }

  # Build k-th order De Bruijn graph (reuse MOGen infrastructure)
  kg <- .mogen_count_kgrams(trajectories, k)

  if (nrow(kg$edges) == 0L) {
    stop(sprintf("No edges at order %d (paths too short or too few)", k))
  }

  # Build adjacency matrix (weighted, NOT normalized)
  nodes <- kg$nodes
  n <- length(nodes)
  adj <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
  idx <- cbind(match(kg$edges$from, nodes), match(kg$edges$to, nodes))
  adj[idx] <- kg$edges$weight

  # Compute Xi (propensity matrix)
  xi <- .hypa_fit_xi(adj)

  # Compute HYPA scores
  scores <- .hypa_compute_scores(adj, xi)

  # Classify anomalies: must exceed alpha AND min_count
  scores$anomaly <- vapply(seq_len(nrow(scores)), function(i) {
    if (scores$observed[i] < min_count) return("normal")
    s <- scores$hypa_score[i]
    if (s < alpha) "under" else if (s > (1 - alpha)) "over" else "normal"
  }, character(1L))

  n_over <- sum(scores$anomaly == "over")
  n_under <- sum(scores$anomaly == "under")

  result <- list(
    scores = scores,
    edges = scores,
    adjacency = adj,
    xi = xi,
    k = k,
    alpha = alpha,
    n_anomalous = n_over + n_under,
    n_over = n_over,
    n_under = n_under,
    n_edges = nrow(scores),
    nodes = nodes
  )

  class(result) <- "net_hypa"
  result
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_hypa
#'
#' @param x A \code{net_hypa} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.net_hypa <- function(x, ...) {
  cat("HYPA: Path Anomaly Detection\n")
  cat(sprintf("  Order k:      %d\n", x$k))
  cat(sprintf("  Edges:        %d\n", x$n_edges))
  cat(sprintf("  Anomalous:    %d (alpha=%.2f)\n", x$n_anomalous, x$alpha))
  cat(sprintf("    Over-repr:  %d\n", x$n_over))
  cat(sprintf("    Under-repr: %d\n", x$n_under))
  invisible(x)
}

#' Summary Method for net_hypa
#'
#' @param object A \code{net_hypa} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
summary.net_hypa <- function(object, ...) {
  cat("HYPA Summary\n\n")
  cat(sprintf("  Order: %d | Nodes: %d | Edges: %d\n",
              object$k, length(object$nodes), object$n_edges))
  cat(sprintf("  Alpha: %.2f\n\n", object$alpha))

  if (object$n_anomalous > 0L) {
    anom <- object$scores[object$scores$anomaly != "normal", ]
    anom <- anom[order(anom$hypa_score), ]
    cat("  Anomalous paths:\n")
    print(anom, row.names = FALSE)
  } else {
    cat("  No anomalous paths detected.\n")
  }

  invisible(object)
}

