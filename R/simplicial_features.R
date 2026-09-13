# ==============================================================================
# simplicial_features() -- tidy topological features, one row per network.
#
# build_simplicial() is deliberately one-network-at-a-time. Using topology as
# regression predictors needs the opposite shape: one row per unit, one column
# per feature. This verb owns that reshaping so callers never loop over a
# netobject_group and rbind the pieces by hand.
# ==============================================================================

utils::globalVariables(c("network", "feature", "value"))

#' Tidy Topological Features for One or Many Networks
#'
#' Builds a simplicial complex per network and returns its topological
#' summaries as a tidy \code{data.frame} -- one row per network, one column
#' per feature -- ready to use as regression predictors or to join onto
#' unit-level outcomes.
#'
#' Higher-order structure is reported as \code{d2}, \code{d3}, ... : the
#' number of simplices of that dimension. A 2-simplex is a triangle of three
#' mutually connected states, a 3-simplex a tetrahedron of four. These count
#' \emph{co-participation in a dense region}, not statistical interaction.
#'
#' @param x A \code{netobject}, \code{netobject_group}, \code{mcml}, a square
#'   weight matrix, or a named list of any of these. A group or list yields
#'   one row per member, an \code{mcml} one row per cluster, and a single
#'   network one row.
#' @param threshold Minimum absolute edge weight for an edge to exist
#'   (passed to \code{\link{build_simplicial}}). Topology is a step function
#'   of this value, so a single threshold is a choice, not a result -- pass a
#'   vector to sweep it and get one row per network per threshold.
#' @param max_dim Maximum simplex dimension retained. Default \code{4}.
#' @param normalize Divide simplex counts by the number of nodes, so networks
#'   of different size are comparable. Default \code{FALSE}.
#' @param type Complex type passed to \code{\link{build_simplicial}}.
#'   Default \code{"clique"}.
#'
#' @return A \code{data.frame} with one row per network (per threshold), and
#'   columns \code{network}, \code{threshold}, \code{n_nodes}, \code{n_edges},
#'   \code{b0}, \code{b1} (Betti numbers), \code{euler}, \code{max_q},
#'   \code{d1} ... \code{d<max_dim>} (simplex counts by dimension), and
#'   \code{higher_order} (the total of \code{d2} upward).
#' @seealso \code{\link{build_simplicial}}, \code{\link{betti_numbers}},
#'   \code{\link{q_analysis}}, \code{\link{outcome_model}}
#' @examples
#' m1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3,
#'              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' m2 <- matrix(c(0, .2, 0, .2, 0, .1, 0, .1, 0), 3, 3,
#'              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' simplicial_features(list(dense = m1, sparse = m2), threshold = 0.3)
#'
#' # Sweep the threshold rather than committing to one.
#' simplicial_features(list(dense = m1), threshold = c(0.1, 0.3, 0.5))
#' @export
simplicial_features <- function(x, threshold = 0, max_dim = 4L,
                                normalize = FALSE, type = "clique") {
  stopifnot(
    "`threshold` must be a finite numeric vector" =
      is.numeric(threshold) && length(threshold) >= 1L && all(is.finite(threshold)),
    "`max_dim` must be a single non-negative integer" =
      length(max_dim) == 1L && is.finite(max_dim) && max_dim >= 0,
    "`normalize` must be TRUE or FALSE" =
      isTRUE(normalize) || isFALSE(normalize)
  )
  nets <- .sf_as_list(x)
  grid <- expand.grid(net = seq_along(nets), thr = seq_along(threshold))
  rows <- Map(function(i, j) {
    .sf_one(nets[[i]], names(nets)[i], threshold[j], max_dim, normalize, type)
  }, grid$net, grid$thr)
  out <- do.call(rbind, rows)
  out <- out[order(out$network, out$threshold), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Normalise every accepted input to a named list of networks.
.sf_as_list <- function(x) {
  if (inherits(x, "mcml")) {
    layers <- x$clusters
    if (is.null(layers) || !length(layers)) {
      stop("mcml has no cluster layers to summarise.", call. = FALSE)
    }
    return(lapply(layers, function(z) z$weights))
  }
  if (inherits(x, "netobject_group")) {
    nms <- names(x) %||% paste0("network_", seq_along(x))
    return(stats::setNames(lapply(seq_along(x), function(i) x[[i]]), nms))
  }
  if (inherits(x, "netobject") || inherits(x, "cograph_network") ||
      inherits(x, "tna") || is.matrix(x)) {
    return(list(network_1 = x))
  }
  if (is.list(x)) {
    nms <- names(x) %||% paste0("network_", seq_along(x))
    return(stats::setNames(x, nms))
  }
  stop("Cannot extract networks from '", class(x)[1L], "'. Supply a ",
       "netobject, netobject_group, mcml, matrix, or a named list of these.",
       call. = FALSE)
}

# Features for a single network at a single threshold.
.sf_one <- function(net, label, threshold, max_dim, normalize, type) {
  sc <- build_simplicial(net, type = type, threshold = threshold,
                         max_dim = as.integer(max_dim))
  bn   <- betti_numbers(sc)
  dims <- vapply(sc$simplices, function(s) length(s) - 1L, integer(1L))
  n    <- sc$n_nodes
  scale <- if (isTRUE(normalize) && n > 0) n else 1

  counts <- vapply(seq_len(max_dim), function(d) sum(dims == d), integer(1L))
  names(counts) <- paste0("d", seq_len(max_dim))

  out <- data.frame(
    network   = label,
    threshold = threshold,
    n_nodes   = n,
    n_edges   = sum(dims == 1L),
    b0        = .sf_betti(bn, "b0"),
    b1        = .sf_betti(bn, "b1"),
    euler     = euler_characteristic(sc),
    max_q     = q_analysis(sc)$max_q,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(as.list(counts / scale)))
  out$higher_order <- sum(counts[-1L]) / scale
  out
}

# Betti numbers above the complex's dimension are absent, not missing -- an
# absent cycle count is a zero, so report it as one rather than an NA.
.sf_betti <- function(bn, which) {
  if (which %in% names(bn)) as.numeric(bn[[which]]) else 0
}
