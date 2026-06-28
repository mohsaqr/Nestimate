# ---- Built-in centrality measures (no external dependencies) ----
#
# All path-based measures (Betweenness, Closeness) are derived from
# all-pairs shortest paths via Floyd-Warshall. For transition/frequency
# networks, weights are inverted so higher weight = shorter distance.

#' All-pairs shortest paths via Floyd-Warshall (vectorized)
#'
#' @param W Square numeric weight matrix (zeros = no edge).
#' @param invert Logical. Convert weights to distances by 1/w? Default TRUE.
#' @return List with \code{D} (distance matrix) and \code{sigma} (shortest
#'   path count matrix).
#' @noRd
.floyd_warshall_sp <- function(W, invert = TRUE) {
  n   <- nrow(W)
  pos <- W > 0

  D              <- matrix(Inf, n, n)
  D[pos]         <- if (invert) 1 / W[pos] else W[pos]
  diag(D)        <- 0

  sigma          <- matrix(0L, n, n)
  sigma[pos]     <- 1L
  diag(sigma)    <- 1L

  Reduce(function(s, k) {
    D     <- s$D
    sigma <- s$sigma
    new_d <- outer(D[, k], D[k, ], "+")
    new_s <- outer(sigma[, k], sigma[k, ], "*")
    # Standard Floyd-Warshall only relaxes (s, t) when both s and t differ
    # from k. With diag(D) = 0, new_d[k, j] would otherwise match D[k, j]
    # under the equal-path test and double-count sigma.
    shorter        <- new_d < D & is.finite(new_d)
    equal          <- (new_d == D) & is.finite(new_d) & new_d > 0
    shorter[k, ]   <- FALSE; shorter[, k] <- FALSE
    equal[k, ]     <- FALSE; equal[, k]   <- FALSE
    sigma[shorter] <- new_s[shorter]
    sigma[equal]   <- sigma[equal] + new_s[equal]
    new_D          <- D
    new_D[shorter] <- new_d[shorter]
    list(D = new_D, sigma = sigma)
  }, seq_len(n), list(D = D, sigma = sigma))
}


#' Betweenness centrality (directed or undirected, weighted)
#'
#' Fraction of all shortest paths passing through each node.
#' Normalized by \code{(n-1)(n-2)} for directed, \code{(n-1)(n-2)/2}
#' for undirected.
#'
#' @param W Square numeric weight matrix.
#' @param directed Logical.
#' @param invert Logical. Invert weights to distances? Default TRUE.
#' @return Named numeric vector of betweenness values.
#' @noRd
.betweenness <- function(W, directed = TRUE, invert = TRUE) {
  n  <- nrow(W)
  if (n < 3L) return(setNames(numeric(n), rownames(W)))
  sp <- .floyd_warshall_sp(W, invert)
  D  <- sp$D
  sg <- sp$sigma

  btw <- vapply(seq_len(n), function(v) {
    idx    <- seq_len(n)[-v]
    d_sv   <- D[idx, v]
    d_vt   <- D[v, idx]
    d_st   <- D[idx, idx]
    sg_sv  <- sg[idx, v]
    sg_vt  <- sg[v, idx]
    sg_st  <- sg[idx, idx]
    d_svt  <- outer(d_sv, d_vt, "+")
    on_path <- is.finite(d_st) & sg_st > 0L &
               abs(d_svt - d_st) < 1e-10
    diag(on_path) <- FALSE
    num  <- outer(sg_sv, sg_vt, "*")
    sum((num / sg_st)[on_path], na.rm = TRUE)
  }, numeric(1))

  norm <- if (directed) (n - 1) * (n - 2) else (n - 1) * (n - 2) / 2
  if (norm > 0) btw <- btw / norm
  setNames(btw, rownames(W))
}


#' Edge betweenness (directed or undirected, weighted)
#'
#' Replaces each edge's weight with the number of shortest paths that
#' traverse it (fractional when shortest paths tie), matching
#' \code{tna::betweenness_network()} / \code{igraph::edge_betweenness()}.
#' For a probability/transition network, weights are inverted to distances
#' (\code{invert = TRUE}) so the geodesic between two states is the most
#' probable route, not the fewest hops.
#'
#' @param W Square numeric weight matrix (zeros = no edge).
#' @param invert Logical. Invert weights to distances by \code{1/w}?
#'   Default \code{TRUE}.
#' @return A numeric matrix the same shape as \code{W}: edge-betweenness on
#'   edge positions, zero elsewhere. Symmetric when \code{W} is symmetric.
#' @noRd
.edge_betweenness <- function(W, invert = TRUE) {
  n   <- nrow(W)
  EB  <- matrix(0, n, n, dimnames = dimnames(W))
  if (n < 2L) return(EB)

  sp  <- .floyd_warshall_sp(W, invert)
  D   <- sp$D
  sg  <- sp$sigma
  pos <- W > 0
  len <- matrix(Inf, n, n)
  len[pos] <- if (invert) 1 / W[pos] else W[pos]

  edges <- which(pos, arr.ind = TRUE)
  vals  <- vapply(seq_len(nrow(edges)), function(e) {
    a <- edges[e, 1L]
    b <- edges[e, 2L]
    # A shortest s->t path uses edge a->b iff
    #   d(s, a) + len(a, b) + d(b, t) == d(s, t).
    # Such paths number sigma(s, a) * sigma(b, t); divide by sigma(s, t).
    through <- outer(D[, a], D[b, ], "+") + len[a, b]
    on_path <- is.finite(D) & sg > 0L & abs(through - D) < 1e-9
    diag(on_path) <- FALSE                       # exclude s == t
    sum((outer(sg[, a], sg[b, ]) / sg)[on_path])
  }, numeric(1))

  EB[edges] <- vals
  EB
}


#' Closeness centrality (directed or undirected, weighted)
#'
#' For directed networks returns both InCloseness and OutCloseness.
#' Closeness = (reachable - 1) / sum(distances to/from reachable nodes).
#' Isolated nodes (no reachable peers) get 0.
#'
#' @param W Square numeric weight matrix.
#' @param directed Logical.
#' @param invert Logical. Invert weights to distances? Default TRUE.
#' @return For directed: named list with \code{InCloseness} and
#'   \code{OutCloseness} vectors. For undirected: named list with
#'   \code{Closeness} vector.
#' @noRd
.closeness <- function(W, directed = TRUE, invert = TRUE) {
  n   <- nrow(W)
  nms <- rownames(W)
  D   <- .floyd_warshall_sp(W, invert)$D

  .cl <- function(d_vec) {
    finite_d <- d_vec[is.finite(d_vec) & d_vec > 0]
    r <- length(finite_d)
    if (r == 0L) 0 else r / sum(finite_d)
  }

  if (directed) {
    list(
      InCloseness  = setNames(vapply(seq_len(n), function(v) .cl(D[, v]),
                                     numeric(1)), nms),
      OutCloseness = setNames(vapply(seq_len(n), function(v) .cl(D[v, ]),
                                     numeric(1)), nms)
    )
  } else {
    list(
      Closeness = setNames(vapply(seq_len(n), function(v) .cl(D[v, ]),
                                  numeric(1)), nms)
    )
  }
}


# ---- Internal centrality() generic (S3 dispatch) ----

#' @noRd
centrality <- function(x, ...) {
  UseMethod("centrality")
}


# ---- Exported net_centrality() ----

#' Compute Centrality Measures for a Network
#'
#' Computes centrality measures from a \code{netobject},
#' \code{netobject_group}, or \code{cograph_network}. For directed networks
#' the default measures are InStrength, OutStrength, and Betweenness. For
#' undirected networks the defaults are Closeness and Betweenness.
#'
#' @param x A \code{netobject}, \code{netobject_group}, or
#'   \code{cograph_network}.
#' @param measures Character vector. Centrality measures to compute.
#'   Built-in: \code{"InStrength"}, \code{"OutStrength"},
#'   \code{"Betweenness"}, \code{"InCloseness"}, \code{"OutCloseness"},
#'   \code{"Closeness"}. \code{"Closeness"} is defined only for
#'   undirected networks; \code{"InCloseness"}/\code{"OutCloseness"}
#'   only for directed networks (requesting the wrong one for the
#'   network's directedness is an error). Default depends on
#'   directedness.
#' @param loops Logical. Include self-loops (diagonal) in computation?
#'   Default: \code{FALSE}.
#' @param centrality_fn Optional function. Custom centrality function that
#'   takes a weight matrix and returns a named list of centrality vectors.
#' @param ... Additional arguments (ignored).
#'
#' @return For a \code{netobject}: a data frame with node names as rows and
#'   centrality measures as columns. For a \code{netobject_group}: a named
#'   list of such data frames (one per group).
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' net_centrality(net)
#'
#' @export
net_centrality <- function(x, measures = NULL, loops = FALSE,
                            centrality_fn = NULL, ...) {
  if (isFALSE(loops)) {
    message("centralities computed excluding loops (diagonal). ",
            "Pass `loops = TRUE` to include self-transitions.")
  }
  centrality(x, measures = measures, loops = loops,
             centrality_fn = centrality_fn, ...)
}


#' @noRd
centrality.netobject <- function(x, measures = NULL, loops = FALSE,
                                  centrality_fn = NULL, ...) {
  mat      <- x$weights
  states   <- x$nodes$label
  directed <- x$directed

  if (is.null(measures)) {
    measures <- if (directed) {
      c("InStrength", "OutStrength", "Betweenness")
    } else {
      c("Closeness", "Betweenness")
    }
  }

  res <- .compute_centralities(mat, states, directed, measures,
                                loops = loops, centrality_fn = centrality_fn)
  as.data.frame(res, row.names = states)
}


#' @noRd
centrality.netobject_group <- function(x, measures = NULL, loops = FALSE,
                                        centrality_fn = NULL, ...) {
  lapply(x, function(net) {
    centrality.netobject(net, measures = measures, loops = loops,
                         centrality_fn = centrality_fn)
  })
}


#' @noRd
centrality.cograph_network <- function(x, measures = NULL, loops = FALSE,
                                        centrality_fn = NULL, ...) {
  centrality.netobject(.as_netobject(x), measures = measures, loops = loops,
                        centrality_fn = centrality_fn)
}


#' @noRd
centrality.mcml <- function(x, measures = NULL, loops = FALSE,
                             centrality_fn = NULL, ...) {
  centrality.netobject_group(as_tna(x), measures = measures, loops = loops,
                              centrality_fn = centrality_fn)
}


# ---- Exported net_edge_betweenness() ----

#' Edge Betweenness Network
#'
#' Builds a network in which each edge's weight is replaced by its
#' betweenness: the number of shortest paths between all node pairs that
#' traverse that edge (fractional when shortest paths tie). This is the
#' Nestimate counterpart of \code{tna::betweenness_network()} and produces
#' identical values for transition networks; the name differs to avoid a
#' clash with \code{tna::betweenness_network()} and
#' \code{igraph::edge_betweenness()}.
#'
#' For a probability/transition network the edge weights are transition
#' probabilities, so they are inverted to distances (\code{invert = TRUE})
#' before path-finding: the geodesic between two states is then the most
#' probable route rather than the one with the fewest hops. Pass
#' \code{invert = FALSE} when the weights already represent distances.
#'
#' Directedness is taken from the network itself. A directed network yields
#' an asymmetric betweenness matrix; an undirected (symmetric) network
#' yields a symmetric one.
#'
#' @param x A \code{netobject} or \code{netobject_group}.
#' @param invert Logical. Invert weights to distances by \code{1/w} before
#'   computing shortest paths? Default \code{TRUE} (correct for probability
#'   and frequency networks).
#' @param ... Additional arguments (ignored).
#'
#' @return For a \code{netobject}: a new \code{netobject} (class
#'   \code{c("netobject", "cograph_network")}) whose \code{$weights} are the
#'   edge-betweenness scores, with \code{method = "edge_betweenness"}. Call
#'   \code{extract_edges()} on it for a tidy per-edge table, or \code{plot()}
#'   to render it. For a \code{netobject_group}: a \code{netobject_group} of
#'   such networks, one per group.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' eb  <- net_edge_betweenness(net)
#' extract_edges(eb)
#'
#' @export
net_edge_betweenness <- function(x, invert = TRUE, ...) {
  UseMethod("net_edge_betweenness")
}

#' @rdname net_edge_betweenness
#' @export
net_edge_betweenness.netobject <- function(x, invert = TRUE, ...) {
  mat <- x$weights
  if (is.null(dimnames(mat)) && !is.null(x$nodes$name)) {
    dimnames(mat) <- list(x$nodes$name, x$nodes$name)
  }
  directed <- if (!is.null(x$directed)) isTRUE(x$directed) else TRUE
  eb <- .edge_betweenness(mat, invert = invert)
  .wrap_netobject(eb, data = x$data, method = "edge_betweenness",
                  directed = directed, inits = x$inits)
}

#' @rdname net_edge_betweenness
#' @export
net_edge_betweenness.netobject_group <- function(x, invert = TRUE, ...) {
  result <- lapply(x, net_edge_betweenness.netobject, invert = invert)
  names(result) <- names(x)
  class(result) <- "netobject_group"
  result
}

#' @rdname net_edge_betweenness
#' @export
net_edge_betweenness.default <- function(x, invert = TRUE, ...) {
  stop("net_edge_betweenness() requires a 'netobject' or 'netobject_group'; ",
       "got '", class(x)[1L], "'.", call. = FALSE)
}
