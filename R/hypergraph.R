# ---- Hypergraph construction (EG-4) -------------------------------------
# Build a hypergraph from a network's clique structure, optionally retaining
# pairwise edges. Foundation for higher-order analyses (centrality, walks,
# contagion, public-goods games).

#' Higher-order hypergraph from a network's clique structure
#'
#' @description
#' `build_hypergraph()` turns a network into a hypergraph by promoting its
#' cliques to hyperedges: every set of 3 to `max_size` nodes that are all
#' tied to each other becomes a hyperedge, kept with probability `p`, and
#' the pairwise edges can be kept as hyperedges of size 2. It takes a
#' `netobject`, a bare `cograph_network`, a `simplicial_complex` or an
#' adjacency or weight matrix, and returns a `net_hypergraph`, the input of
#' the hypergraph verbs such as [hypergraph_measures()],
#' [hypergraph_centrality()] and [clique_expansion()].
#'
#' A hyperedge built this way records that every pair of its nodes is tied
#' in the thresholded network; a pairwise network cannot show whether those
#' nodes also act together as a group. When group memberships are observed
#' directly, [bipartite_groups()] builds the hypergraph from them instead.
#' Each clique also contributes its smaller cliques: with `max_size = 4`, a
#' 4-clique gives one hyperedge of size 4 and four of size 3.
#'
#' @param net The network: a `netobject`, bare `cograph_network`,
#'   `simplicial_complex` (its edges are read as ties of weight 1 and the
#'   cliques are found again up to `max_size`), or a numeric or logical adjacency or weight
#'   matrix (nodes are named `V1`, `V2`, ... when it has no row names). A
#'   tie counts in either direction, so a directed network is symmetrised,
#'   and self-loops are ignored.
#' @param p The probability that each clique of 3 or more nodes becomes a
#'   hyperedge, a number in `[0, 1]`; each clique is kept or dropped
#'   independently. Default `1` (every clique is kept). Pairwise edges are
#'   not sampled.
#' @param method How hyperedges are found. `"clique"` (default) promotes
#'   the cliques of the thresholded network. `"vr"` and `"rips"` (a metric
#'   Vietoris-Rips construction) are not implemented and raise an error
#'   instead of falling back to `"clique"`.
#' @param include_pairwise Whether the pairwise edges of the thresholded
#'   network become hyperedges of size 2. `FALSE` gives a hypergraph of
#'   hyperedges of size 3 or more only. Default `TRUE`.
#' @param max_size The largest hyperedge size, at least 2. Default `3`
#'   (triangles); `4` adds the 4-cliques, and so on. With `2`, only the
#'   pairwise edges remain (none, with `include_pairwise = FALSE`).
#' @param threshold The cut-off that turns weights into ties: an edge is
#'   kept when its absolute weight is non-zero and at least `threshold`, so
#'   negative weights count by their size. A negative value raises an
#'   error. Default `0` (every non-zero weight is a tie).
#' @param seed An integer that makes the sampling reproducible when
#'   `0 < p < 1`. The session's own random stream is left as it was.
#'   Default `NULL`.
#'
#' @return A `net_hypergraph` object. `print()` shows the numbers of nodes
#'   and hyperedges, the size distribution and the settings, and returns
#'   its input invisibly. `summary()` prints the counts with the mean and
#'   largest hyperedge size, and returns a data frame with one row per node:
#'   `node` and `degree`, the number of hyperedges that contain it. The
#'   object holds:
#' \describe{
#'   \item{`hyperedges`}{A list of integer vectors, one per hyperedge,
#'     each the sorted indices of the nodes it spans; pairwise edges, when
#'     kept, come first.}
#'   \item{`incidence`}{An integer matrix of nodes by hyperedges, 1 where a
#'     node belongs to a hyperedge and 0 elsewhere. Its row names are the
#'     node names and its column names `h1`, `h2`, ...}
#'   \item{`nodes`}{The node names.}
#'   \item{`n_nodes`, `n_hyperedges`}{The numbers of nodes and
#'     hyperedges.}
#'   \item{`size_distribution`}{A named integer vector with the number of
#'     hyperedges of each size, named `size_2`, `size_3`, ...}
#'   \item{`params`}{The settings used: `method`, `p`,
#'     `include_pairwise`, `max_size`, `threshold` and `seed`.}
#' }
#'
#' @details
#' The construction follows Burgio, Matamalas, Gomez and Arenas (2020) on
#' cooperation in networks and hypergraphs. The cliques are found with
#' [build_simplicial()], and each clique of 3 or more nodes then goes
#' through an independent Bernoulli(`p`) trial that decides whether it
#' becomes a hyperedge. The pairwise edges are not sampled: with
#' `include_pairwise = TRUE` all of them are kept, so the hypergraph holds
#' the pairwise structure and the sampled higher-order interactions
#' together.
#'
#' At the limits:
#' \itemize{
#'   \item `p = 0` with `include_pairwise = TRUE` gives the thresholded
#'     pairwise network as a hypergraph of size-2 edges.
#'   \item `p = 1` with `include_pairwise = FALSE` gives only the
#'     hyperedges of size 3 or more, every clique of the network's clique
#'     complex up to `max_size`.
#' }
#'
#' @seealso [build_simplicial()] (the clique enumeration),
#'   [bipartite_groups()], [build_network()].
#'
#' @examples
#' set.seed(1)
#' n <- 8
#' adj <- matrix(stats::rbinom(n * n, 1, 0.5), n, n)
#' diag(adj) <- 0
#' adj <- (adj + t(adj)) > 0
#' rownames(adj) <- colnames(adj) <- LETTERS[seq_len(n)]
#' hg <- build_hypergraph(adj, p = 1, max_size = 3L)
#' print(hg)
#' summary(hg)
#'
#' @references
#' Burgio, G., Matamalas, J. T., Gomez, S., & Arenas, A. (2020). Evolution
#' of cooperation in the presence of higher-order interactions: from
#' networks to hypergraphs. \emph{Entropy} 22(7), 744.
#' \doi{10.3390/e22070744}
#'
#' @export
build_hypergraph <- function(net,
                              p = 1,
                              method = c("clique", "vr", "rips"),
                              include_pairwise = TRUE,
                              max_size = 3L,
                              threshold = 0,
                              seed = NULL) {
  method <- match.arg(method)
  if (method == "vr" || method == "rips") {
    stop(errorCondition(
      paste0(
        "method = \"vr\" (Vietoris-Rips) is not implemented. A genuine ",
        "metric Vietoris-Rips construction requires a distance-ball ",
        "graph, which this package does not build. Use ",
        "method = \"clique\" (cliques in the binarised adjacency)."),
      class = "nestimate_bad_input", call = NULL))
  }
  .assert(
    is.numeric(p), length(p) == 1L, p >= 0, p <= 1,
    is.logical(include_pairwise), length(include_pairwise) == 1L,
    is.numeric(max_size), length(max_size) == 1L, max_size >= 2L,
    is.numeric(threshold), length(threshold) == 1L
  )
  max_size <- as.integer(max_size)
  if (!is.null(seed)) .local_seed(as.integer(seed))

  # ---- Get adjacency + node names from any supported input ------------
  parsed <- .hg_extract_adj(net)
  adj    <- parsed$adj
  nodes  <- parsed$nodes
  n      <- length(nodes)

  # ---- Find all simplices (cliques) up to max_size --------------------
  sc <- build_simplicial(adj, type = method, threshold = threshold,
                         max_dim = max_size - 1L)
  simplices <- sc$simplices

  sizes   <- vapply(simplices, length, integer(1L))
  edges_2 <- simplices[sizes == 2L]

  # k>=3 simplices: sampled with prob p
  hi_simp <- simplices[sizes >= 3L & sizes <= max_size]
  sampled_hi <- if (length(hi_simp) == 0L || p <= 0) {
    list()
  } else if (p >= 1) {
    hi_simp
  } else {
    keep <- as.logical(stats::rbinom(length(hi_simp), 1L, p))
    hi_simp[keep]
  }

  hyperedges <- if (include_pairwise) c(edges_2, sampled_hi) else sampled_hi
  m <- length(hyperedges)

  # ---- Build incidence matrix [n_nodes x n_hyperedges] ----------------
  if (m == 0L) {
    incidence <- matrix(0L, nrow = n, ncol = 0L,
                        dimnames = list(nodes, NULL))
  } else {
    incidence <- matrix(0L, nrow = n, ncol = m,
                        dimnames = list(nodes, paste0("h", seq_len(m))))
    for (j in seq_len(m)) {
      incidence[hyperedges[[j]], j] <- 1L
    }
  }

  # ---- Size distribution ----------------------------------------------
  he_sizes <- vapply(hyperedges, length, integer(1L))
  size_dist <- if (length(he_sizes)) {
    tab <- table(he_sizes)
    out <- as.integer(tab)
    names(out) <- paste0("size_", names(tab))
    out
  } else {
    integer(0L)
  }

  structure(
    list(
      hyperedges        = hyperedges,
      incidence         = incidence,
      nodes             = nodes,
      n_nodes           = n,
      n_hyperedges      = m,
      size_distribution = size_dist,
      params = list(
        method           = method,
        p                = p,
        include_pairwise = include_pairwise,
        max_size         = max_size,
        threshold        = threshold,
        seed             = seed
      )
    ),
    class = "net_hypergraph"
  )
}

# ---- Internal helpers --------------------------------------------------

#' Extract symmetric numeric adjacency + node names from any supported input
#' @noRd
.hg_extract_adj <- function(net) {
  if (inherits(net, "simplicial_complex")) {
    nodes <- net$nodes
    n     <- length(nodes)
    adj   <- matrix(0, n, n, dimnames = list(nodes, nodes))
    edge_simp <- Filter(function(s) length(s) == 2L, net$simplices)
    for (e in edge_simp) {
      adj[e[1L], e[2L]] <- 1
      adj[e[2L], e[1L]] <- 1
    }
    return(list(adj = adj, nodes = nodes))
  }

  if (inherits(net, "netobject") || inherits(net, "cograph_network")) {
    w <- net$weights
    nm <- if (!is.null(net$nodes$name)) {
      net$nodes$name
    } else if (!is.null(net$nodes$label)) {
      net$nodes$label
    } else if (!is.null(net$nodes$id)) {
      as.character(net$nodes$id)
    } else if (!is.null(rownames(w))) {
      rownames(w)
    } else {
      paste0("V", seq_len(nrow(w)))
    }
    rownames(w) <- colnames(w) <- nm
    return(list(adj = w, nodes = nm))
  }

  if (is.matrix(net) && (is.numeric(net) || is.logical(net))) {
    nm <- rownames(net) %||% paste0("V", seq_len(nrow(net)))
    storage.mode(net) <- "double"
    rownames(net) <- colnames(net) <- nm
    return(list(adj = net, nodes = nm))
  }

  stop(errorCondition(
    paste0(
      "`net` must be a netobject, cograph_network, simplicial_complex, ",
      "or numeric matrix."),
    class = "nestimate_bad_input", call = NULL))
}

# ---- S3 methods ---------------------------------------------------------

#' @rdname build_hypergraph
#' @param x For the `print()` method: a `net_hypergraph`.
#' @param object For the `summary()` method: a `net_hypergraph`.
#' @param ... Not used.
#' @export
print.net_hypergraph <- function(x, ...) {
  cat(sprintf("Hypergraph: %d nodes, %d hyperedges\n",
              x$n_nodes, x$n_hyperedges))
  if (length(x$size_distribution)) {
    cat("Size distribution:\n")
    for (nm in names(x$size_distribution)) {
      cat(sprintf("  %-8s : %d\n", nm, x$size_distribution[[nm]]))
    }
  }
  cat(sprintf("Method: %s, p = %.2f, include_pairwise = %s, max_size = %d\n",
              x$params$method, x$params$p,
              x$params$include_pairwise, x$params$max_size))
  invisible(x)
}

#' @rdname build_hypergraph
#' @export
summary.net_hypergraph <- function(object, ...) {
  he_sizes <- vapply(object$hyperedges, length, integer(1L))
  cat("Hypergraph summary\n")
  cat(sprintf("  Nodes:         %d\n", object$n_nodes))
  cat(sprintf("  Hyperedges:    %d\n", object$n_hyperedges))
  if (length(he_sizes)) {
    cat(sprintf("  Mean size:     %.2f\n", mean(he_sizes)))
    cat(sprintf("  Max size:      %d\n", max(he_sizes)))
  }

  nodes <- if (!is.null(object$nodes)) object$nodes else
    rownames(object$incidence)
  if (is.null(nodes)) nodes <- paste0("n", seq_len(object$n_nodes))
  if (object$n_hyperedges > 0L) {
    deg <- as.integer(rowSums(object$incidence > 0))
  } else {
    deg <- rep(0L, object$n_nodes)
  }
  data.frame(node = as.character(nodes), degree = deg,
             stringsAsFactors = FALSE, row.names = NULL)
}
