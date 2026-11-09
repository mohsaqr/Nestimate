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
#' @description
#' `simplicial_features()` builds a simplicial complex for each network with
#' [build_simplicial()] and returns its topological summaries as a tidy data
#' frame, one row per network and one column per feature, ready to use as
#' regression predictors or to join onto unit-level outcomes. It takes one
#' network, a group of networks, the within-cluster networks of an `mcml`,
#' or a list of networks, at one threshold or a sweep of them.
#'
#' Higher-order structure is reported as `d2`, `d3`, ...: the number of
#' simplices of that dimension. A 2-simplex is a triangle of three mutually
#' connected states, a 3-simplex a tetrahedron of four. These count
#' co-participation in a dense region, not statistical interaction. The
#' counts are a step function of `threshold`, so one threshold is a choice,
#' not a result: sweep several and check that a conclusion holds across
#' them.
#'
#' @param x The networks: a `netobject`, `cograph_network`, `tna` model or
#'   square weight matrix (one row); a `netobject_group`, or a list of
#'   single networks or matrices (one row per member, named by the list, or
#'   `network_1`, `network_2`, ... when unnamed); or an `mcml` object (one row per
#'   cluster, from its within-cluster networks; one built without them
#'   raises an error).
#' @param threshold The smallest absolute edge weight that makes an edge,
#'   passed to [build_simplicial()]. Default `0` (every non-zero edge). A
#'   vector sweeps it, giving one row per network per threshold. Values must
#'   be finite and not negative. `type = "vr"` and `type = "pathway"` ignore
#'   it, although the `threshold` column still records it.
#' @param max_dim The highest simplex dimension built, a whole number of at
#'   least 1 (a fractional value is truncated). Default `4`. When a network holds a larger clique, the complex
#'   is truncated and [betti_numbers()] warns with class
#'   `nestimate_truncated_complex`; the warning concerns the Betti number of
#'   dimension `max_dim`, so it affects `b1` only when `max_dim = 1`.
#' @param normalize Whether to divide the simplex counts (`d1` onwards and
#'   `higher_order`) by the number of nodes, so networks of different size
#'   are comparable; `n_edges` stays a count. Default `FALSE`.
#' @param type The complex construction passed to [build_simplicial()]:
#'   `"clique"` (default), `"pathway"`, `"vr"` or its alias `"rips"`.
#'   `"pathway"` needs networks of class `netobject` or `tna`, from which a
#'   higher-order network is built; a matrix, a bare `cograph_network` or an
#'   `mcml` (whose cluster layers are matrices) raises an error.
#'
#' @return A data frame with one row per network and threshold, ordered by
#'   network name and then threshold, and columns `network`, `threshold`,
#'   `n_nodes`, `n_edges`, `b0` and `b1` (the Betti numbers: connected
#'   components and independent cycles, 0 when the complex has no such
#'   dimension), `euler` (the Euler characteristic of the complex as
#'   built), `max_q` (the top q-level of [q_analysis()], the dimension of the
#'   largest maximal simplex), `d1`, `d2`, ... up to dimension `max_dim`
#'   (simplex counts by dimension) and `higher_order` (the total from `d2` upward).
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
  .assert(
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
      stop(errorCondition("mcml has no cluster layers to summarise.",
                          class = "nestimate_bad_input", call = NULL))
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
  stop(errorCondition(
    paste0(
      "Cannot extract networks from '",
      class(x)[1L],
      "'. Supply a ",
      "netobject, netobject_group, mcml, matrix, or a named list of these."),
    class = "nestimate_bad_input", call = NULL))
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
