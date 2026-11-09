# ---- Link Prediction for Network Objects ----

#' Predict Missing or Future Links in a Network
#'
#' @description
#' `predict_links()` scores every pair of nodes for how strongly the
#' network's structure suggests a link between them, with one or more
#' similarity indices: common neighbours, resource allocation, Adamic-Adar,
#' Jaccard, preferential attachment and Katz (Liben-Nowell & Kleinberg,
#' 2007). It returns one row per pair and index, ranked within each index,
#' and a consensus ranking across the indices. It takes a `netobject`, a
#' `cograph_network`, a `tna` model or a square weight matrix, weighted or
#' binary, directed or undirected.
#'
#' A high score says that a pair resembles linked pairs in the network's
#' structure (shared neighbours, busy endpoints, many connecting paths), not
#' that the link exists in the population or will form. The indices are
#' heuristics and often disagree; the consensus rank shows the pairs on which
#' they agree. How well an index recovers links is measured by scoring
#' held-out edges with [evaluate_links()]; [pathways()] writes the top
#' predictions as pathway strings.
#'
#' @param x The network: a `netobject`, a `cograph_network`, a `tna` model,
#'   or a numeric square weight matrix. A network keeps its own direction; a
#'   matrix or `tna` model is directed unless it is symmetric. A
#'   `netobject_group` or an `mcml` raises an error of class
#'   `nestimate_bad_input`: predict each network separately. For the `print()`
#'   method: a `predict_links()` result.
#' @param methods The indices to compute, one or more of
#'   `"common_neighbors"`, `"resource_allocation"`, `"adamic_adar"`,
#'   `"jaccard"`, `"preferential_attachment"` and `"katz"` (defined under
#'   *Details*). By default all six.
#' @param weighted Whether common neighbours, resource allocation,
#'   Adamic-Adar and Jaccard use the edge weights (`TRUE`, default) or the
#'   binary adjacency (`FALSE`). Preferential attachment, Katz and the degrees
#'   that resource allocation and Adamic-Adar divide by always use the binary
#'   adjacency.
#' @param top_n The number of predictions kept per index, the highest
#'   scores. Default `NULL` (every pair). The ranks are computed before the
#'   cut, so they stay the same.
#' @param exclude_existing Whether pairs that already have an edge are left
#'   out of the predictions. Default `TRUE`.
#' @param include_self Whether self-pairs (a node with itself) are listed for
#'   a directed network; an undirected network never lists them. All indices
#'   score self-pairs 0. It also decides whether self-loops count in the
#'   number of existing edges. Default `FALSE`.
#' @param katz_damping The attenuation factor \eqn{\beta} of the Katz index.
#'   Default `NULL`, which uses \eqn{0.9 / \rho(A)} for the spectral radius
#'   \eqn{\rho(A)} of the adjacency (`0.1` when it is 0). A value at or above
#'   \eqn{1 / \rho(A)}, where the series diverges, is replaced by
#'   \eqn{0.9 / \rho(A)} with a warning.
#'
#' @return A `net_link_prediction` object. `print()` lists the top ten
#'   predicted links (by consensus rank with several indices, by score with
#'   one) and returns its input invisibly. `summary()` returns a data frame
#'   with one row per index: `method`, `n_predictions`, and the `score_mean`,
#'   `score_sd`, `score_max` and `score_min` of its predictions; an index
#'   with no predictions (every possible link already exists) has
#'   `n_predictions = 0` and `NA` scores. [evaluate_links()] scores the
#'   predictions against known edges, and [pathways()] writes them as pathway
#'   strings. The object's fields:
#' \describe{
#'   \item{predictions}{Data frame, one row per node pair and index:
#'     `from`, `to`, `method`, `score`, `existing` (whether the pair already
#'     has an edge) and `rank`, sorted by score (highest first) within each
#'     index. An undirected network lists each pair once.}
#'   \item{consensus}{Data frame, one row per node pair: `from`, `to`,
#'     `avg_rank` (the pair's mean rank over the indices that list it),
#'     `n_methods` (how many do) and `consensus_rank`, ordered by `avg_rank`.
#'     `NULL` when only one index was computed.}
#'   \item{scores}{Named list of score matrices, one per index.}
#'   \item{adjacency}{The 0/1 adjacency matrix of the input network.}
#'   \item{methods}{The indices computed.}
#'   \item{nodes}{The node names.}
#'   \item{directed, weighted}{Whether the network was read as directed, and
#'     whether weights were used.}
#'   \item{n_nodes}{The number of nodes.}
#'   \item{n_existing}{The number of existing edges: non-zero cells of the
#'     adjacency matrix, each undirected edge counted once, self-loops
#'     counted only with `include_self = TRUE`.}
#' }
#'
#' @details
#' ## Methods
#'
#' \describe{
#'   \item{common_neighbors}{The number of shared neighbours. For a directed
#'     network, shared out-neighbours plus shared in-neighbours, computed as
#'     \code{A \%*\% t(A) + t(A) \%*\% A}.}
#'   \item{resource_allocation}{Zhou et al. (2009). Like common neighbours,
#'     but each shared neighbour z is weighted by \code{1/degree(z)} (in- plus
#'     out-degree for a directed network). It penalises hubs and rewards rare
#'     shared connections.}
#'   \item{adamic_adar}{Adamic & Adar (2003). Like resource allocation, but
#'     weighted by \code{1/log(degree(z))}, a milder penalty; a neighbour of
#'     degree 1 adds nothing.}
#'   \item{jaccard}{Shared neighbours divided by all neighbours of the two
#'     nodes. With `weighted = TRUE` it is the weighted Jaccard (Ruzicka)
#'     index, \eqn{\sum_z \min(|w_{iz}|, |w_{jz}|) / \sum_z \max(|w_{iz}|,
#'     |w_{jz}|)}, which equals the binary index on 0/1 weights and lies
#'     between 0 and 1. For a directed network, computed on the combined (out +
#'     in) neighbourhoods, each neighbour weighted by the larger of its two
#'     directions.}
#'   \item{preferential_attachment}{The source's out-degree times the
#'     target's in-degree (Barabasi & Albert, 1999): the "rich get richer"
#'     effect.}
#'   \item{katz}{Katz (1953). A sum over all paths between the two nodes,
#'     damped exponentially by path length, computed by matrix inversion as
#'     \code{(I - beta * A)^{-1} - I}. It reads the global structure. If the
#'     inversion fails, the series is truncated after six terms.}
#' }
#'
#' Every index is computed with matrix operations on the whole network at
#' once.
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
#' Jaccard, P. (1901). Etude comparative de la distribution florale dans une
#' portion des Alpes et des Jura. \emph{Bulletin de la Societe Vaudoise des
#' Sciences Naturelles}, 37, 547--579.
#'
#' Ioffe, S. (2010). Improved consistent sampling, weighted minhash and L1
#' sketching. \emph{IEEE International Conference on Data Mining},
#' 246--255.
#'
#' Barabasi, A.-L. & Albert, R. (1999). Emergence of scaling in random
#' networks. \emph{Science}, 286(5439), 509--512.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A", "B", "C", "D", "A", "C", "E", "B"),
#'   V2 = c("B", "C", "D", "E", "C", "E", "A", "D"),
#'   V3 = c("C", "D", "E", "A", "D", "A", "B", "E")
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

  .assert(is.numeric(W), nrow(W) == n, ncol(W) == n)
  methods <- match.arg(methods, several.ok = TRUE)

  # Binary adjacency
  A <- (W != 0) * 1L
  # Working matrix: weighted or binary
  M <- if (isTRUE(weighted)) W else A

  # An undirected edge is one pair: count the upper triangle only, so the
  # count matches the network's own n_edges (it used to count both cells).
  counted <- if (directed) A != 0 else upper.tri(A, diag = TRUE) & A != 0
  n_existing <- sum(counted)
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
  # Weighted Jaccard (Ruzicka similarity): sum_z min(w_iz, w_jz) over
  # sum_z max(w_iz, w_jz), on absolute weights (a negative edge is still a
  # neighbour, as in the adjacency). With 0/1 weights it is exactly
  # |N(i) & N(j)| / |N(i) | N(j)| (Jaccard 1901; Ioffe 2010 for the weighted
  # form). The old code divided a weighted intersection by a binary union,
  # which went negative for weights above 1 and was clamped to 0.
  # Directed: one weight per neighbour over both directions, the larger of
  # the two, which reduces to the combined (out + in) neighbour set.
  B <- abs(M)
  if (directed) B <- pmax(B, t(B))
  n <- nrow(B)
  overlap <- function(f) {
    vapply(seq_len(n), function(i) {
      rowSums(f(B, matrix(B[i, ], n, n, byrow = TRUE)))
    }, numeric(n))
  }
  intersection <- overlap(pmin)
  union_mat <- overlap(pmax)
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
  if (inherits(x, "mcml")) {
    .stop_mcml_input("predict_links()", "scores one network")
  }
  if (inherits(x, "netobject_group")) {
    stop(errorCondition(paste(
      "predict_links() scores one network, not a group: pass each group's",
      "network on its own."),
      class = "nestimate_bad_input", call = NULL))
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
  stop(errorCondition(
    "x must be a netobject, cograph_network, or numeric matrix.",
    class = "nestimate_bad_input", call = NULL))
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
#' `evaluate_links()` measures how well each index of a [predict_links()]
#' result recovers a set of known edges. For each index it returns the area
#' under the ROC curve (AUC), the average precision, and the precision among
#' the top `k` predictions.
#'
#' The AUC is the probability that a known edge scores higher than a pair
#' that is not one (ties count half), so 0.5 is chance level. Average
#' precision is the mean of the precisions at the rank of each known edge,
#' and precision\eqn{@}k the share of known edges among the `k`
#' highest-scoring pairs. Only the pairs listed in the predictions are
#' evaluated: with `exclude_existing = TRUE`, the default of
#' [predict_links()], the network's current edges are not scored, so the
#' known edges should be ones held out of the network (a later period, a
#' removed subset), and `top_n` limits the evaluation to the pairs it kept.
#' The AUC is `NA` when the evaluated pairs include no known edge, or
#' nothing but known edges.
#'
#' @param pred A [predict_links()] result.
#' @param true_edges The known edges: a data frame with columns `from` and
#'   `to` (every row is an edge, whatever its other columns, such as a
#'   weight, hold), or a matrix whose non-zero cells are edges. Its row
#'   names name both its rows and its columns (column names are not read),
#'   so its columns must be in the same order as its rows; without row
#'   names, the nodes of `pred` are used. For a directed prediction, pairs
#'   are matched by `from` and `to`, direction included; for an undirected
#'   one, in either order, so `(C, A)` matches the listed pair `(A, C)`.
#' @param k The cut-offs for precision\eqn{@}k. Default `c(5, 10, 20)`. A
#'   `k` larger than an index's number of predictions gives `NA`.
#'
#' @return A data frame with one row per index: `method`, `auc`,
#'   `average_precision` (`0` when no evaluated pair is a known edge), and
#'   one `precision_at_<k>` column per value of `k`. A prediction with no
#'   listed pairs (every possible link already exists) gives `NULL`.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A", "B", "C", "D", "A", "C", "E", "B"),
#'   V2 = c("B", "C", "D", "E", "C", "E", "A", "D"),
#'   V3 = c("C", "D", "E", "A", "D", "A", "B", "E")
#' )
#' net <- build_network(seqs, method = "relative")
#' pred <- predict_links(net, exclude_existing = FALSE)
#'
#' # Evaluate against the network's own edges as the known truth
#' evaluate_links(pred, extract_edges(net, threshold = 0.001))
#'
#' @export
evaluate_links <- function(pred, true_edges, k = c(5L, 10L, 20L)) {
  .assert(
    inherits(pred, "net_link_prediction"),
    "`k` must be positive whole numbers" =
      is.numeric(k) && length(k) >= 1L && all(is.finite(k)) &&
        all(k >= 1) && all(k == round(k))
  )

  # Parse true edges to a set of "from->to" keys
  if (is.matrix(true_edges)) {
    nodes <- rownames(true_edges) %||% pred$nodes
    idx <- which(true_edges != 0, arr.ind = TRUE)
    true_from <- nodes[idx[, 1]]
    true_to <- nodes[idx[, 2]]
  } else {
    .assert(is.data.frame(true_edges),
              all(c("from", "to") %in% names(true_edges)))
    true_from <- as.character(true_edges$from)
    true_to <- as.character(true_edges$to)
  }
  # An undirected pair has no direction: match (C, A) to the listed (A, C).
  pair_key <- function(from, to) {
    if (isTRUE(pred$directed)) return(paste0(from, "\t", to))
    paste0(pmin(from, to), "\t", pmax(from, to))
  }
  true_keys <- pair_key(true_from, true_to)

  methods <- unique(pred$predictions$method)
  results <- lapply(methods, function(m) {
    df <- pred$predictions[pred$predictions$method == m, , drop = FALSE]
    df <- df[order(-df$score), , drop = FALSE]
    pred_keys <- pair_key(df$from, df$to)
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
  if (!length(results)) {
    # no candidate pairs (every pair is already an edge): an empty table with
    # the same columns, not NULL
    empty <- c(list(method = character(0), auc = numeric(0),
                    average_precision = numeric(0)),
               stats::setNames(rep(list(numeric(0)), length(k)),
                               paste0("precision_at_", k)))
    return(as.data.frame(empty, stringsAsFactors = FALSE))
  }
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

#' @rdname predict_links
#' @param ... For the `print()` and `summary()` methods: not used.
#' @param object For `summary()`: a `predict_links()` result.
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


#' @rdname predict_links
#' @export
summary.net_link_prediction <- function(object, ...) {
  df <- object$predictions
  do.call(rbind, lapply(object$methods, function(m) {
    score <- df$score[df$method == m]
    # A complete network has no missing links to score: report NA, not the
    # NaN / -Inf / Inf that mean() / max() / min() return on an empty vector.
    has_scores <- length(score) > 0L
    data.frame(
      method        = m,
      n_predictions = length(score),
      score_mean    = if (has_scores) mean(score) else NA_real_,
      score_sd      = if (has_scores) stats::sd(score) else NA_real_,
      score_max     = if (has_scores) max(score) else NA_real_,
      score_min     = if (has_scores) min(score) else NA_real_,
      stringsAsFactors = FALSE,
      row.names     = NULL
    )
  }))
}
