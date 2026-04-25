# ---- Link Prediction for Network Objects ----

#' Predict Missing or Future Links in a Network
#'
#' @description
#' Computes link prediction scores for all node pairs using one or more
#' structural similarity methods. Accepts \code{netobject}, \code{mcml},
#' \code{cograph_network}, or a raw weight matrix.
#'
#' All methods are fully vectorized using matrix operations — no loops.
#' Supports both weighted and binary adjacency, directed and undirected
#' networks.
#'
#' @param x A \code{netobject}, \code{mcml}, \code{cograph_network}, or
#'   numeric square matrix.
#' @param methods Character vector. One or more of:
#'   \code{"common_neighbors"}, \code{"resource_allocation"},
#'   \code{"adamic_adar"}, \code{"jaccard"}, \code{"preferential_attachment"},
#'   \code{"katz"}. Default: all six methods.
#' @param weighted Logical. If \code{TRUE}, use the weight matrix directly
#'   instead of binarizing. Default: \code{TRUE}.
#' @param top_n Integer or NULL. Return only the top N predictions per method.
#'   Default: \code{NULL} (all pairs).
#' @param exclude_existing Logical. If \code{TRUE}, exclude node pairs that
#'   already have an edge. Default: \code{TRUE}.
#' @param include_self Logical. If \code{TRUE}, include self-loop predictions.
#'   Default: \code{FALSE}.
#' @param katz_damping Numeric or NULL. Attenuation factor for Katz index.
#'   If NULL, auto-computed as \code{0.9 / spectral_radius(A)}.
#'   Default: \code{NULL}.
#'
#' @return An object of class \code{"net_link_prediction"} containing:
#' \describe{
#'   \item{predictions}{Data frame with columns: from, to, method, score, rank.
#'     Sorted by score (descending) within each method.}
#'   \item{scores}{Named list of score matrices (one per method).}
#'   \item{methods}{Character vector of methods used.}
#'   \item{nodes}{Character vector of node names.}
#'   \item{directed}{Logical.}
#'   \item{weighted}{Logical.}
#'   \item{n_nodes}{Integer.}
#'   \item{n_existing}{Integer. Number of existing edges.}
#' }
#'
#' @details
#' ## Methods
#'
#' \describe{
#'   \item{common_neighbors}{Number of shared neighbors. For directed graphs,
#'     sums shared out-neighbors and shared in-neighbors.
#'     Vectorized as \code{A \%*\% t(A) + t(A) \%*\% A}.}
#'   \item{resource_allocation}{Zhou et al. (2009). Like common neighbors but
#'     weights each shared neighbor z by \code{1/degree(z)}.
#'     Penalizes hubs, rewards rare shared connections.}
#'   \item{adamic_adar}{Adamic & Adar (2003). Like resource allocation but
#'     weights by \code{1/log(degree(z))}. Less aggressive penalty than RA.}
#'   \item{jaccard}{Ratio of shared neighbors to total neighbors.
#'     For directed graphs, computed on combined (out+in) neighbor sets.}
#'   \item{preferential_attachment}{Product of source out-degree and target
#'     in-degree. Captures the "rich-get-richer" effect.}
#'   \item{katz}{Katz (1953). Weighted sum of all paths between nodes,
#'     exponentially damped by path length. Computed via matrix inversion:
#'     \code{(I - beta * A)^{-1} - I}. Captures global structure.}
#' }
#'
#' @references
#' Liben-Nowell, D. & Kleinberg, J. (2007). The link-prediction problem for
#' social networks. \emph{JASIST}, 58(7), 1019--1031.
#'
#' Zhou, T., Lu, L. & Zhang, Y.-C. (2009). Network topology and link
#' prediction. \emph{European Physical Journal B}, 71, 623--630.
#'
#' Adamic, L. A. & Adar, E. (2003). Friends and neighbors on the Web.
#' \emph{Social Networks}, 25(3), 211--230.
#'
#' Katz, L. (1953). A new status index derived from sociometric analysis.
#' \emph{Psychometrika}, 18(1), 39--43.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:5], 50, TRUE),
#'   V2 = sample(LETTERS[1:5], 50, TRUE),
#'   V3 = sample(LETTERS[1:5], 50, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' pred <- predict_links(net)
#' print(pred)
#' summary(pred)
#'
#' @seealso \code{\link{evaluate_links}} for prediction evaluation,
#'   \code{\link{build_network}} for network estimation.
#'
#' @export
predict_links <- function(x,
                          methods = c("common_neighbors", "resource_allocation",
                                      "adamic_adar", "jaccard",
                                      "preferential_attachment", "katz"),
                          weighted = TRUE,
                          top_n = NULL,
                          exclude_existing = TRUE,
                          include_self = FALSE,
                          katz_damping = NULL) {

  # ---- Extract weight matrix ----
  parsed <- .lp_parse_input(x)
  W <- parsed$weights
  directed <- parsed$directed
  nodes <- parsed$nodes
  n <- length(nodes)

  stopifnot(is.numeric(W), nrow(W) == n, ncol(W) == n)
  methods <- match.arg(methods, several.ok = TRUE)

  # Binary adjacency
  A <- (W != 0) * 1L
  # Working matrix: weighted or binary
  M <- if (isTRUE(weighted)) W else A

  n_existing <- sum(A != 0)
  if (!include_self) {
    n_existing <- n_existing - sum(diag(A) != 0)
  }

  # ---- Compute scores ----
  score_mats <- list()

  if ("common_neighbors" %in% methods) {
    score_mats$common_neighbors <- .lp_common_neighbors(M, directed)
  }
  if ("resource_allocation" %in% methods) {
    score_mats$resource_allocation <- .lp_resource_allocation(M, A, directed)
  }
  if ("adamic_adar" %in% methods) {
    score_mats$adamic_adar <- .lp_adamic_adar(M, A, directed)
  }
  if ("jaccard" %in% methods) {
    score_mats$jaccard <- .lp_jaccard(M, directed)
  }
  if ("preferential_attachment" %in% methods) {
    score_mats$preferential_attachment <- .lp_preferential_attachment(A)
  }
  if ("katz" %in% methods) {
    score_mats$katz <- .lp_katz(A, katz_damping)
  }

  # Name dimensions
  score_mats <- lapply(score_mats, function(s) {
    dimnames(s) <- list(nodes, nodes)
    s
  })

  # ---- Build predictions data frame ----
  predictions <- .lp_build_predictions(
    score_mats, nodes, A, directed,
    exclude_existing = exclude_existing,
    include_self = include_self,
    top_n = top_n
  )

  # ---- Consensus ranking across methods ----
  consensus <- .lp_build_consensus(predictions, methods)

  structure(list(
    predictions = predictions,
    consensus   = consensus,
    scores      = score_mats,
    adjacency   = A,
    methods     = methods,
    nodes       = nodes,
    directed    = directed,
    weighted    = weighted,
    n_nodes     = n,
    n_existing  = n_existing
  ), class = "net_link_prediction")
}


#' Build consensus ranking across link prediction methods
#' @noRd
.lp_build_consensus <- function(predictions, methods) {
  if (length(methods) <= 1L || nrow(predictions) == 0L) return(NULL)

  # Create pair key for aggregation
  pair_key <- paste0(predictions$from, "\t", predictions$to)

  # Average rank and method count per pair
  avg_rank <- tapply(predictions$rank, pair_key, mean)
  n_methods <- tapply(predictions$rank, pair_key, length)

  keys <- names(avg_rank)
  parts <- strsplit(keys, "\t", fixed = TRUE)

  result <- data.frame(
    from      = vapply(parts, `[`, character(1), 1L),
    to        = vapply(parts, `[`, character(1), 2L),
    avg_rank  = as.numeric(avg_rank),
    n_methods = as.integer(n_methods),
    stringsAsFactors = FALSE
  )
  result <- result[order(result$avg_rank), , drop = FALSE]
  result$consensus_rank <- seq_len(nrow(result))
  rownames(result) <- NULL
  result
}


# ---- Vectorized Score Functions ----

#' @noRd
.lp_common_neighbors <- function(M, directed) {
  if (directed) {
    # Shared out-neighbors + shared in-neighbors
    s <- tcrossprod(M) + crossprod(M)
  } else {
    s <- tcrossprod(M)
  }
  diag(s) <- 0
  s
}

#' @noRd
.lp_resource_allocation <- function(M, A, directed) {
  # RA(i,j) = sum_z M(i,z)*M(z,j) / degree(z)
  deg <- if (directed) rowSums(A) + colSums(A) else rowSums(A)
  inv_deg <- ifelse(deg > 0, 1 / deg, 0)
  D <- diag(inv_deg)
  if (directed) {
    s <- M %*% D %*% t(M) + t(M) %*% D %*% M
  } else {
    s <- M %*% D %*% M
  }
  diag(s) <- 0
  s
}

#' @noRd
.lp_adamic_adar <- function(M, A, directed) {
  deg <- if (directed) rowSums(A) + colSums(A) else rowSums(A)
  inv_log_deg <- ifelse(deg > 1, 1 / log(deg), 0)
  D <- diag(inv_log_deg)
  if (directed) {
    s <- M %*% D %*% t(M) + t(M) %*% D %*% M
  } else {
    s <- M %*% D %*% M
  }
  diag(s) <- 0
  s
}

#' @noRd
.lp_jaccard <- function(M, directed) {
  # Jaccard on combined neighbor sets (out + in for directed)
  if (directed) {
    B <- pmin(M + t(M), 1)  # symmetrized binary presence
  } else {
    B <- M
  }
  intersection <- tcrossprod(B)
  # |N(i) union N(j)| = |N(i)| + |N(j)| - |N(i) intersect N(j)|
  deg <- rowSums(B > 0)
  union_mat <- outer(deg, deg, "+") - intersection
  s <- ifelse(union_mat > 0, intersection / union_mat, 0)
  diag(s) <- 0
  s
}

#' @noRd
.lp_preferential_attachment <- function(A) {
  s <- outer(rowSums(A), colSums(A), "*")
  diag(s) <- 0
  s
}

#' @noRd
.lp_katz <- function(A, damping) {
  n <- nrow(A)
  # Auto-compute damping from spectral radius
  if (is.null(damping)) {
    sr <- max(abs(eigen(A, only.values = TRUE)$values))
    damping <- if (sr > 0) 0.9 / sr else 0.1
  } else {
    sr <- max(abs(eigen(A, only.values = TRUE)$values))
    if (sr > 0 && damping >= 1 / sr) {
      damping <- 0.9 / sr
      warning("katz_damping exceeded convergence bound; auto-adjusted to ",
              round(damping, 4), call. = FALSE)
    }
  }
  S <- tryCatch(
    solve(diag(n) - damping * A) - diag(n),
    error = function(e) {
      # Fallback: power series truncated at 6 terms
      result <- matrix(0, n, n)
      Ak <- A
      for (k in seq_len(6L)) {
        result <- result + damping^k * Ak
        Ak <- Ak %*% A
      }
      result
    }
  )
  diag(S) <- 0
  S
}


# ---- Input Parsing ----

#' @noRd
.lp_parse_input <- function(x) {
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "netobject_group")) {
    stop("predict_links() requires a single network, not a group. ",
         "Use lapply(group, predict_links) for per-group predictions.",
         call. = FALSE)
  }
  if (inherits(x, "tna")) {
    return(list(
      weights  = x$weights,
      directed = !isSymmetric(x$weights),
      nodes    = rownames(x$weights)
    ))
  }
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
  }
  if (inherits(x, "netobject")) {
    return(list(
      weights  = x$weights,
      directed = x$directed,
      nodes    = x$nodes$label
    ))
  }
  if (is.matrix(x) && is.numeric(x)) {
    nodes <- rownames(x) %||% paste0("V", seq_len(nrow(x)))
    return(list(
      weights  = x,
      directed = !isSymmetric(x),
      nodes    = nodes
    ))
  }
  stop("x must be a netobject, cograph_network, or numeric matrix.",
       call. = FALSE)
}


# ---- Predictions Data Frame ----

#' @noRd
.lp_build_predictions <- function(score_mats, nodes, A, directed,
                                   exclude_existing, include_self, top_n) {
  n <- length(nodes)
  # Pre-build row/col indices
  idx <- expand.grid(row = seq_len(n), col = seq_len(n),
                     KEEP.OUT.ATTRS = FALSE)
  idx$from <- nodes[idx$row]
  idx$to <- nodes[idx$col]

  # Filter self-loops
  if (!include_self) {
    idx <- idx[idx$row != idx$col, , drop = FALSE]
  }

  # For undirected, keep only upper triangle
  if (!directed) {
    idx <- idx[idx$row < idx$col, , drop = FALSE]
  }

  # Existing edge mask
  existing <- A[cbind(idx$row, idx$col)] != 0

  dfs <- lapply(names(score_mats), function(m) {
    s <- score_mats[[m]]
    scores <- s[cbind(idx$row, idx$col)]
    df <- data.frame(
      from   = idx$from,
      to     = idx$to,
      method = m,
      score  = scores,
      existing = existing,
      stringsAsFactors = FALSE
    )
    if (exclude_existing) {
      df <- df[!df$existing, , drop = FALSE]
    }
    df <- df[order(-df$score), , drop = FALSE]
    df$rank <- seq_len(nrow(df))
    if (!is.null(top_n) && nrow(df) > top_n) {
      df <- df[seq_len(top_n), , drop = FALSE]
    }
    df
  })
  result <- do.call(rbind, dfs)
  rownames(result) <- NULL
  result
}


# ---- Evaluation ----

#' Evaluate Link Predictions Against Known Edges
#'
#' @description
#' Computes AUC-ROC, precision\eqn{@}k, and average precision for link
#' predictions against a set of known true edges.
#'
#' @param pred A \code{net_link_prediction} object.
#' @param true_edges A data frame with columns \code{from} and \code{to},
#'   or a binary matrix where 1 indicates a true edge.
#' @param k Integer vector. Values of k for precision\eqn{@}k.
#'   Default: \code{c(5, 10, 20)}.
#'
#' @return A data frame with columns: method, auc, average_precision,
#'   and one precision_at_k column per k value.
#'
#' @examples
#' set.seed(42)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:5], 50, TRUE),
#'   V2 = sample(LETTERS[1:5], 50, TRUE),
#'   V3 = sample(LETTERS[1:5], 50, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' pred <- predict_links(net, exclude_existing = FALSE)
#'
#' # Evaluate: predict the network's own edges
#' true <- data.frame(from = pred$predictions$from[1:5],
#'                    to = pred$predictions$to[1:5])
#' evaluate_links(pred, true)
#'
#' @export
evaluate_links <- function(pred, true_edges, k = c(5L, 10L, 20L)) {
  stopifnot(inherits(pred, "net_link_prediction"))

  # Parse true edges to a set of "from->to" keys
  if (is.matrix(true_edges)) {
    nodes <- rownames(true_edges) %||% pred$nodes
    idx <- which(true_edges != 0, arr.ind = TRUE)
    true_keys <- paste0(nodes[idx[, 1]], "\t", nodes[idx[, 2]])
  } else {
    stopifnot(is.data.frame(true_edges),
              all(c("from", "to") %in% names(true_edges)))
    true_keys <- paste0(true_edges$from, "\t", true_edges$to)
  }

  methods <- unique(pred$predictions$method)
  results <- lapply(methods, function(m) {
    df <- pred$predictions[pred$predictions$method == m, , drop = FALSE]
    df <- df[order(-df$score), , drop = FALSE]
    pred_keys <- paste0(df$from, "\t", df$to)
    is_true <- pred_keys %in% true_keys

    # AUC via Wilcoxon-Mann-Whitney
    auc <- .lp_compute_auc(df$score, is_true)

    # Average precision
    ap <- .lp_compute_ap(is_true)

    # Precision@k
    pk <- vapply(k, function(ki) {
      if (ki > length(is_true)) return(NA_real_)
      sum(is_true[seq_len(ki)]) / ki
    }, numeric(1))
    names(pk) <- paste0("precision_at_", k)

    c(list(method = m, auc = auc, average_precision = ap), as.list(pk))
  })
  do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
}

#' @noRd
.lp_compute_auc <- function(scores, labels) {
  # Wilcoxon-Mann-Whitney AUC
  pos <- scores[labels]
  neg <- scores[!labels]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  # Count pairs where pos > neg, plus 0.5 for ties
  comparisons <- outer(pos, neg, "-")
  (sum(comparisons > 0) + 0.5 * sum(comparisons == 0)) /
    (length(pos) * length(neg))
}

#' @noRd
.lp_compute_ap <- function(is_true) {
  if (sum(is_true) == 0) return(0)
  cum_true <- cumsum(is_true)
  precision_at_i <- cum_true / seq_along(is_true)
  sum(precision_at_i * is_true) / sum(is_true)
}


# ---- S3 Methods ----

#' Print Method for net_link_prediction
#'
#' @param x A \code{net_link_prediction} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE),
#'   V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' pred <- predict_links(net)
#' print(pred)
#'
#' @export
print.net_link_prediction <- function(x, ...) {
  dir_lbl <- if (x$directed) "directed" else "undirected"
  wt_lbl <- if (x$weighted) "weighted" else "binary"
  cat(sprintf("Link Prediction  [%s | %s | %d nodes | %d existing edges]\n",
              dir_lbl, wt_lbl, x$n_nodes, x$n_existing))
  cat(sprintf("  Methods: %s\n", paste(x$methods, collapse = ", ")))

  arrow <- if (x$directed) "->" else "--"

  if (!is.null(x$consensus)) {
    # Multi-method: show consensus ranking
    top <- utils::head(x$consensus, 10L)
    n_total <- nrow(x$consensus)
    cat(sprintf("\n  Top predicted links (consensus across %d methods):\n",
                length(x$methods)))
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      cat(sprintf("    %d. %s %s %s  (avg rank: %.1f, agreed: %d/%d)\n",
                  i, r$from, arrow, r$to, r$avg_rank,
                  r$n_methods, length(x$methods)))
    }
    if (n_total > 10L) {
      cat(sprintf("    ... and %d more predictions\n", n_total - 10L))
    }
  } else {
    # Single method: show directly
    df <- x$predictions
    top <- utils::head(df, 10L)
    if (nrow(top) > 0L) {
      cat(sprintf("\n  Top predicted links (%s):\n", x$methods[1]))
      for (i in seq_len(nrow(top))) {
        r <- top[i, ]
        cat(sprintf("    %d. %s %s %s  (score: %.4f)\n",
                    i, r$from, arrow, r$to, r$score))
      }
      if (nrow(df) > 10L) {
        cat(sprintf("    ... and %d more predictions\n", nrow(df) - 10L))
      }
    }
  }
  invisible(x)
}


#' Summary Method for net_link_prediction
#'
#' @param object A \code{net_link_prediction} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with per-method summary statistics, invisibly.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE),
#'   V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' pred <- predict_links(net)
#' summary(pred)
#'
#' @export
summary.net_link_prediction <- function(object, ...) {
  df <- object$predictions
  do.call(rbind, lapply(object$methods, function(m) {
    sub <- df[df$method == m, , drop = FALSE]
    data.frame(
      method        = m,
      n_predictions = nrow(sub),
      score_mean    = mean(sub$score, na.rm = TRUE),
      score_sd      = sd(sub$score, na.rm = TRUE),
      score_max     = max(sub$score, na.rm = TRUE),
      score_min     = min(sub$score, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names     = NULL
    )
  }))
}
