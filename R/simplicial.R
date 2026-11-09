# ---- Simplicial Complex Analysis ----
#
# Construction, homology, centrality, and Q-analysis for simplicial
# complexes built from networks and higher-order pathway data.
#
# Clique finding verified against igraph::cliques().
# Betti numbers verified against known topological invariants.

# =========================================================================
# Core constructor
# =========================================================================

#' Build a Simplicial Complex
#'
#' @description
#' `build_simplicial()` builds a simplicial complex, the higher-order
#' structure of a network: besides nodes (0-simplices) and edges
#' (1-simplices) it holds filled triangles (2-simplices), tetrahedra and
#' larger groups of nodes. Three constructions are available through `type`:
#'
#' \itemize{
#'   \item \strong{Clique complex} (`"clique"`): every clique of the network
#'     whose edges reach `threshold` becomes a simplex, the standard bridge
#'     from a graph to algebraic topology.
#'   \item \strong{Pathway complex} (`"pathway"`): the states on each
#'     higher-order pathway of a `net_hon`, `net_hypa` or `net_mogen` (or of
#'     one built from a `tna` model or `netobject`) become a simplex.
#'   \item \strong{Vietoris-Rips complex} (`"vr"`, alias `"rips"`): the input
#'     is a non-negative distance matrix, every group of nodes whose pairwise
#'     distances are all at most `max_scale` becomes a simplex, and each
#'     simplex records the scale at which it enters (its largest pairwise
#'     distance), so [persistent_homology()] can read the filtration.
#' }
#'
#' A clique complex is fixed by the thresholded graph, since every clique is
#' filled in; its Betti numbers count the holes that remain, such as cycles
#' of edges that no triangles fill. Its shape depends on `threshold`, so
#' read [betti_numbers()] at more than one threshold, or follow the features
#' across all of them with [persistent_homology()]. A pathway complex is
#' built from paths observed in the sequences, so its simplices say which
#' states occur together on a path.
#' [simplicial_degree()] and [q_analysis()] describe how the nodes take part
#' in the complex.
#'
#' @details
#' For the clique and Vietoris-Rips complexes the weight matrix is made
#' symmetric first (the larger absolute weight of the two directions for a
#' clique complex, the larger distance for a Vietoris-Rips complex) and its
#' diagonal is not read. Cliques are found with `igraph::cliques()` when
#' \pkg{igraph} is installed, and by Bron-Kerbosch otherwise.
#'
#' @param x The network: a square weight matrix, a `tna` model, a
#'   `netobject` or `cograph_network`, or a `net_hon`, `net_hypa` or
#'   `net_mogen`. With `type = "vr"`, the matrix holds non-negative
#'   distances. The clique and Vietoris-Rips constructions read the weight
#'   matrix of any of these networks; with `type = "pathway"`, a bare
#'   `cograph_network` or a matrix raises an error. For the `print()` and
#'   `plot()` methods: a `simplicial_complex` object.
#' @param type The construction: `"clique"` (default), `"pathway"` or
#'   `"vr"` (alias `"rips"`, stored as `"vr"`).
#' @param threshold For `type = "clique"`: the smallest absolute edge weight
#'   an edge needs to enter (edges at exactly `threshold` are included). Zero
#'   weights are never edges. Default `0`. Ignored for `"vr"`, which uses
#'   `max_scale`, and for `"pathway"`.
#' @param max_dim The highest simplex dimension built (a k-simplex has k + 1
#'   nodes), a non-negative whole number. Default `10`. When a larger clique
#'   is cut at this dimension, the clique or Vietoris-Rips complex is marked
#'   truncated; [betti_numbers()] then warns with class
#'   `nestimate_truncated_complex`, as does [persistent_homology()] on a
#'   Vietoris-Rips complex. A pathway with more than `max_dim + 1` distinct
#'   states is cut to its faces of that size, without the truncated mark.
#' @param max_pathways For `type = "pathway"`, how many pathways to keep:
#'   the most frequent higher-order paths for a `net_hon`, the most frequent
#'   transitions of [mogen_transitions()] at the optimal order for a
#'   `net_mogen` (the highest order fitted when the optimal one is 0), and
#'   for a `net_hypa` the anomalous paths ranked by their observed /
#'   expected ratio (see `anomaly`). Default `NULL` (all).
#' @param anomaly For `type = "pathway"` on a `net_hypa`: which anomalous
#'   paths become simplices: `"all"` (default), `"over"` or `"under"`.
#'   Over-represented and `"all"` paths are ranked by the largest ratio,
#'   under-represented ones by the smallest. Possible paths that never occur
#'   are left out, since a pathway simplex is a set of states observed
#'   together on a path.
#' @param max_scale For `type = "vr"`: the largest distance at which an edge
#'   enters; longer edges are left out. Default `NULL`, the largest finite
#'   distance.
#' @param ... For `type = "pathway"` on a `tna` model or `netobject`: `method`
#'   (`"hon"`, the default, `"hypa"` or `"mogen"`) chooses the higher-order
#'   model built from it, and the other arguments go to [build_hon()],
#'   [build_hypa()] or [build_mogen()]. Any other construction or input does
#'   not read `...` and warns with class `nestimate_ignored_argument` when it
#'   is given. For the `print()` and `plot()` methods: ignored.
#' @return A `simplicial_complex` object. `print()` shows the counts, the
#'   density, the f-vector (simplices per dimension), the Euler
#'   characteristic and the non-zero Betti numbers; `plot()` draws the
#'   f-vector, the Betti numbers and the simplicial degree of each node;
#'   [betti_numbers()], [euler_characteristic()], [simplicial_degree()],
#'   [q_analysis()] and [persistent_homology()] read it. The object's
#'   fields:
#' \describe{
#'   \item{simplices}{List of integer vectors, one per simplex, each holding
#'     the sorted indices of the nodes it spans. Every face of every simplex
#'     is present, as are all 0-simplices (isolated nodes included).}
#'   \item{nodes}{Character vector of node labels, which `simplices` index.}
#'   \item{n_nodes, n_simplices}{Counts of nodes and simplices.}
#'   \item{dimension}{Highest simplex dimension present.}
#'   \item{f_vector}{Named integer vector `dim_0`, `dim_1`, ...: the number
#'     of simplices of each dimension.}
#'   \item{density}{`n_simplices` divided by the number of simplices a
#'     complete complex on `n_nodes` nodes would have up to `dimension`.}
#'   \item{mean_dim}{Mean simplex dimension.}
#'   \item{type}{`"clique"`, `"pathway"` or `"vr"`.}
#'   \item{truncated}{For the clique and Vietoris-Rips complexes: whether
#'     `max_dim` cut a larger clique.}
#' }
#'   A Vietoris-Rips complex also holds `filtration` (the scale at which each
#'   simplex enters, in the order of `simplices`) and `max_scale` (the cap
#'   used), which [persistent_homology()] reads.
#'
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' print(sc)
#' betti_numbers(sc)
#'
#' # Vietoris-Rips on a distance matrix:
#' d <- 1 - mat
#' diag(d) <- 0
#' sc_vr <- build_simplicial(d, type = "vr", max_scale = 0.6)
#'
#' @seealso \code{\link{betti_numbers}}, \code{\link{persistent_homology}},
#'   \code{\link{simplicial_degree}}, \code{\link{q_analysis}}
#'
#' @export
build_simplicial <- function(x, type = "clique", threshold = 0,
                              max_dim = 10L, max_pathways = NULL,
                              anomaly = c("all", "over", "under"),
                              max_scale = NULL, ...) {
  type <- match.arg(type, c("clique", "pathway", "vr", "rips"))
  anomaly <- match.arg(anomaly)
  .assert(
    "`threshold` must be a single non-negative number" =
      is.numeric(threshold) && length(threshold) == 1L &&
        !is.na(threshold) && threshold >= 0,
    "`max_dim` must be a single non-negative whole number" =
      is.numeric(max_dim) && length(max_dim) == 1L &&
        !is.na(max_dim) && max_dim >= 0 && max_dim == as.integer(max_dim)
  )
  if (type != "pathway") .sc_warn_ignored_dots("build_simplicial()", ...)

  if (type == "vr" || type == "rips") {
    d <- .sc_extract_matrix(x)
    fc <- .filter_vr_complex(d, max_dim = max_dim, max_scale = max_scale,
                             detect_truncation = TRUE)
    sc <- .make_simplicial_complex(fc$simplices, fc$nodes, "vr")
    # Reorder filtration to match the simplex order in sc$simplices.
    sc$filtration <- .align_filtration(fc, sc)
    sc$max_scale <- fc$max_w
    sc$truncated <- isTRUE(fc$truncated)
    return(sc)
  }

  if (type == "pathway") {
    return(.build_simplicial_pathway(x, max_dim, max_pathways,
                                     anomaly = anomaly, ...))
  }

  mat <- .sc_extract_matrix(x)
  .build_simplicial_clique(mat, threshold, max_dim)
}

# `...` reaches a higher-order builder only for a pathway complex built from a
# tna model or netobject; anywhere else it used to be dropped without a word,
# so a misspelled argument had no effect.
#' @noRd
.sc_warn_ignored_dots <- function(fn, ...) {
  if (...length() == 0L) return(invisible(NULL))
  given <- ...names() %||% rep("", ...length())
  given[is.na(given) | !nzchar(given)] <- "<unnamed>"
  warning(warningCondition(
    sprintf(paste0(
      "%s: `...` is read only for type = \"pathway\" on a tna model or ",
      "netobject; ignored here: %s."),
      fn, paste(given, collapse = ", ")),
    class = "nestimate_ignored_argument", call = NULL
  ))
  invisible(NULL)
}

#' @noRd
.align_filtration <- function(fc, sc) {
  # .make_simplicial_complex may reorder simplices and add isolated vertices;
  # re-key the filtration vector to match sc$simplices order.
  fc_keys <- fc$key
  sc_keys <- vapply(sc$simplices, function(s) paste(sort(s), collapse = ","),
                    character(1))
  out <- numeric(length(sc_keys))
  m <- match(sc_keys, fc_keys)
  out[!is.na(m)] <- fc$filt_asc[m[!is.na(m)]]
  out[is.na(m)] <- 0  # isolated vertices added by .make_simplicial_complex
  out
}

# =========================================================================
# Clique complex - verified against igraph::cliques()
# =========================================================================

#' @noRd
.build_simplicial_clique <- function(mat, threshold = 0, max_dim = 10L,
                                     inclusive = TRUE) {
  .assert(
    "the network must be a square numeric weight matrix" =
      is.matrix(mat) && (is.numeric(mat) || is.logical(mat)) &&
        nrow(mat) == ncol(mat)
  )
  nodes <- rownames(mat) %||% paste0("V", seq_len(nrow(mat)))
  n <- nrow(mat)

  adj <- .sc_threshold_adjacency(mat, threshold, inclusive = inclusive)
  diag(adj) <- FALSE
  adj <- adj | t(adj)

  simplices <- .find_all_cliques(adj, max_dim)

  sc <- .make_simplicial_complex(simplices, nodes, "clique")
  # max_dim cut a larger clique: the top-dimension cycles have no filling
  # simplices in the complex, so betti_numbers() warns about them
  sc$truncated <- .clique_exceeds(adj, simplices, max_dim)
  sc
}

#' Does the graph hold a clique larger than max_dim + 1 vertices?
#'
#' True when some (max_dim + 1)-vertex clique already enumerated extends by
#' one vertex adjacent to all of its members. Reads only the simplices that
#' were built, so it costs (top simplices x vertices) rather than an exact
#' maximum-clique search, which took seconds on dense graphs.
#' @param adj Symmetric logical adjacency, zero diagonal.
#' @param simplices Integer vertex sets of the complex.
#' @param max_dim The construction's max_dim.
#' @noRd
.clique_exceeds <- function(adj, simplices, max_dim) {
  tops <- Filter(function(s) length(s) == max_dim + 1L, simplices)
  any(vapply(tops, function(s) {
    # a member is never adjacent to itself, so only outsiders can reach |s|
    any(colSums(adj[s, , drop = FALSE]) == length(s))
  }, logical(1L)))
}

#' Is a filtered complex still truncated below a scale cap?
#'
#' A complex built at a larger scale may be truncated there but not within
#' the cap persistent_homology() is asked to use: recheck with the edges and
#' top simplices whose filtration value is at most \code{cap}.
#' @noRd
.sc_truncated_within <- function(x, cap) {
  if (is.null(cap)) return(TRUE)
  inside <- x$filtration <= cap
  simp <- x$simplices[inside]
  adj <- matrix(FALSE, x$n_nodes, x$n_nodes)
  edges <- Filter(function(s) length(s) == 2L, simp)
  if (length(edges)) {
    ij <- do.call(rbind, edges)
    adj[ij] <- TRUE
    adj[ij[, 2:1, drop = FALSE]] <- TRUE
  }
  .clique_exceeds(adj, simp, x$dimension)
}

#' @noRd
.sc_threshold_adjacency <- function(mat, threshold, inclusive = TRUE) {
  weights <- abs(mat)
  diag(weights) <- 0
  weights <- pmax(weights, t(weights))

  if (isTRUE(inclusive)) {
    weights > 0 & weights >= threshold
  } else {
    weights > 0 & weights > threshold
  }
}

#' Find all cliques (all sizes) via igraph or fallback Bron-Kerbosch
#'
#' When igraph is available, delegates to igraph::cliques() for
#' correctness and speed. Results are verified to match on package tests.
#' @noRd
.find_all_cliques <- function(adj, max_dim = 10L) {
  n <- nrow(adj)

  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                              diag = FALSE)
    raw <- igraph::cliques(g, min = 1, max = max_dim + 1L)
    simplices <- lapply(raw, function(cl) sort(as.integer(cl)))
  } else {
    # Fallback: Bron-Kerbosch + expand to all faces
    maximal <- .bron_kerbosch_all(adj) # nocov
    simplices <- .expand_to_faces(maximal, max_dim) # nocov
  }

  simplices
}

# =========================================================================
# Pathway complex (from HON/HYPA)
# =========================================================================

#' @noRd
.build_simplicial_pathway <- function(x, max_dim = 10L,
                                       max_pathways = NULL,
                                       anomaly = c("all", "over", "under"),
                                       ...) {
  anomaly <- match.arg(anomaly)
  if (inherits(x, c("net_hon", "net_hypa", "net_mogen"))) {
    .sc_warn_ignored_dots("build_simplicial()", ...)
  }

  if (inherits(x, "net_hon")) {
    edges <- x$ho_edges
    ho <- edges[edges$from_order > 1L, , drop = FALSE]
    ho <- ho[order(-ho$count), , drop = FALSE]
    if (!is.null(max_pathways) && nrow(ho) > max_pathways) {
      ho <- ho[seq_len(max_pathways), , drop = FALSE]
    }
    nodes <- x$first_order_states
    raw_paths <- ho$path
  } else if (inherits(x, "net_hypa")) {
    # a pathway simplex is a set of states observed together on a path, so a
    # possible path that never occurs (scored "under" since 0.9.26) is not one
    scores <- x$scores[x$scores$observed > 0, , drop = FALSE]
    if (anomaly == "all") {
      anom <- scores[scores$anomaly != "normal", , drop = FALSE]
    } else {
      anom <- scores[scores$anomaly == anomaly, , drop = FALSE]
    }
    if ("ratio" %in% names(anom)) {
      if (anomaly == "under") {
        anom <- anom[order(anom$ratio), , drop = FALSE]
      } else {
        anom <- anom[order(-anom$ratio), , drop = FALSE]
      }
    }
    if (!is.null(max_pathways) && nrow(anom) > max_pathways) { # nocov
      anom <- anom[seq_len(max_pathways), , drop = FALSE] # nocov
    }
    parts <- strsplit(
      gsub("\x01", " -> ", x$nodes$label, fixed = TRUE), " -> ", fixed = TRUE
    )
    nodes <- sort(unique(unlist(parts)))
    raw_paths <- anom$path
  } else if (inherits(x, "net_mogen")) {
    # Use mogen_transitions() at optimal (or highest available) order
    # Its $path column is already in "A -> B -> C" format
    order_used <- x$optimal_order
    if (order_used < 1L) order_used <- max(x$orders[x$orders >= 1L], 0L)
    if (order_used < 1L) {
      stop(errorCondition(
        "MOGen model has no higher-order transitions (optimal_order = 0)",
        class = "nestimate_bad_input", call = NULL))
    }
    trans <- mogen_transitions(x, order = order_used)
    if (nrow(trans) == 0L) {
      return(.make_simplicial_complex(list(), x$states, "pathway"))
    }
    if (!is.null(max_pathways) && nrow(trans) > max_pathways) {
      trans <- trans[seq_len(max_pathways), , drop = FALSE]
    }
    nodes <- x$states
    raw_paths <- trans$path
  } else if (inherits(x, c("tna", "netobject"))) {
    dots <- list(...)
    method <- match.arg(dots$method %||% "hon", c("hon", "hypa", "mogen"))
    dots$method <- NULL
    ho_obj <- switch(method,
      hon   = do.call(build_hon,   c(list(x), dots)),
      hypa  = do.call(build_hypa,  c(list(x), dots)),
      mogen = do.call(build_mogen, c(list(x), dots))
    )
    return(.build_simplicial_pathway(ho_obj, max_dim, max_pathways,
                                     anomaly = anomaly))
  } else {
    stop(errorCondition(
      paste0(
        "For type='pathway', x must be a net_hon, net_hypa, net_mogen, ",
        "tna, or netobject."),
      class = "nestimate_bad_input", call = NULL))
  }

  if (length(raw_paths) == 0L) {
    return(.make_simplicial_complex(list(), nodes, "pathway"))
  }

  node_idx <- setNames(seq_along(nodes), nodes)
  simplices_raw <- lapply(raw_paths, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    unique(parts)
  })

  simplices <- lapply(simplices_raw, function(s) {
    idx <- node_idx[s]
    sort(idx[!is.na(idx)])
  })
  simplices <- simplices[vapply(simplices, length, integer(1)) >= 2L]
  simplices <- .expand_to_faces(simplices, max_dim)

  .make_simplicial_complex(simplices, nodes, "pathway")
}

# =========================================================================
# Matrix extraction
# =========================================================================

#' @noRd
.sc_extract_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (inherits(x, "tna")) return(x$weights)
  if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    return(x$weights)
  }
  if (inherits(x, "net_hon")) return(x$matrix)
  stop(errorCondition(
    paste0("Cannot extract adjacency matrix from '", class(x)[1], "'."),
    class = "nestimate_bad_input", call = NULL))
}

# =========================================================================
# Bron-Kerbosch (fallback when igraph unavailable)
# =========================================================================

#' @noRd
.bron_kerbosch_all <- function(adj) { # nocov start
  n <- nrow(adj)
  neighbors <- lapply(seq_len(n), function(i) which(adj[i, ]))
  cliques <- list()

  .bk <- function(R, P, X) {
    if (length(P) == 0L && length(X) == 0L) {
      cliques[[length(cliques) + 1L]] <<- sort(R)
      return(invisible(NULL))
    }
    union_px <- c(P, X)
    pivot <- union_px[which.max(
      vapply(union_px, function(v) sum(P %in% neighbors[[v]]), integer(1))
    )]
    # A loop: each branch shrinks P and grows X for the branches after it.
    for (v in setdiff(P, neighbors[[pivot]])) {
      nbrs <- neighbors[[v]]
      .bk(c(R, v), intersect(P, nbrs), intersect(X, nbrs))
      P <- setdiff(P, v)
      X <- c(X, v)
    }
  }

  .bk(integer(0), seq_len(n), integer(0))
  cliques
} # nocov end

#' Expand maximal cliques to all sub-simplices
#' @noRd
.expand_to_faces <- function(simplices, max_dim = 10L) {
  # every face of every simplex, smallest first within a simplex; a face shared
  # by several simplices is kept where it first appears
  faces <- unlist(lapply(simplices, function(simplex) {
    simplex <- sort(as.integer(simplex))
    max_size <- min(length(simplex), max_dim + 1L)
    unlist(lapply(seq_len(max_size), function(size) {
      utils::combn(simplex, size, simplify = FALSE)
    }), recursive = FALSE)
  }), recursive = FALSE)
  if (is.null(faces)) return(list())
  keys <- vapply(faces, paste, character(1), collapse = ",")
  faces[!duplicated(keys)]
}

# =========================================================================
# Constructor
# =========================================================================

#' @noRd
.make_simplicial_complex <- function(simplices, nodes, type) {
  # Ensure all 0-simplices are present
  keys <- vapply(simplices, function(s) paste(sort(s), collapse = ","),
                 character(1))
  absent <- setdiff(as.character(seq_along(nodes)), keys)
  simplices <- c(simplices, as.list(as.integer(absent)))

  dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))
  max_d <- if (length(dims) > 0L) max(dims) else 0L
  f_vec <- vapply(0:max_d, function(d) sum(dims == d), integer(1))
  names(f_vec) <- paste0("dim_", 0:max_d)

  n_simplices <- length(simplices)
  n_nodes <- length(nodes)

  # Density: fraction of possible simplices that exist
  # Max possible k-simplices = C(n, k+1)
  max_possible <- sum(vapply(0:max_d, function(d) {
    choose(n_nodes, d + 1L)
  }, numeric(1)))
  density <- if (max_possible > 0) n_simplices / max_possible else 0

  # Mean simplex dimension
  mean_dim <- if (n_simplices > 0L) mean(dims) else 0

  structure(list(
    simplices = simplices,
    nodes = nodes,
    n_nodes = n_nodes,
    n_simplices = n_simplices,
    dimension = max_d,
    f_vector = f_vec,
    density = density,
    mean_dim = mean_dim,
    type = type
  ), class = "simplicial_complex")
}

# =========================================================================
# Homology
# =========================================================================

#' Betti Numbers
#'
#' @description
#' `betti_numbers()` counts the holes of a simplicial complex in each
#' dimension: \eqn{\beta_0} is the number of connected components,
#' \eqn{\beta_1} the number of independent loops (cycles not filled in by
#' triangles), \eqn{\beta_2} the number of enclosed voids, and so on up to
#' the dimension of the complex.
#'
#' Betti numbers describe the complex as it was built; a clique complex
#' built at a different threshold can give different numbers.
#' [persistent_homology()] follows each hole across all thresholds and shows
#' how long it lasts.
#'
#' @details
#' \eqn{\beta_k} is the dimension of the k-th homology group, the number of
#' k-simplices minus the rank of the k-th boundary matrix, minus the rank of
#' the (k + 1)-th. The ranks are computed with `qr()` on the signed boundary
#' matrices.
#'
#' @param sc A `simplicial_complex` from [build_simplicial()]. When it was
#'   built with a `max_dim` below its largest clique, the top Betti number
#'   counts cycles whose filling simplices were not built, so it is an upper
#'   bound, and a warning of class `nestimate_truncated_complex` says so.
#' @return A named integer vector `c(b0 = , b1 = , ...)`, one entry per
#'   dimension from 0 to the dimension of the complex.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' betti_numbers(sc)
#'
#' @export
betti_numbers <- function(sc) {
  .assert("`sc` must be a simplicial_complex (see build_simplicial())" =
            inherits(sc, "simplicial_complex"))
  if (isTRUE(sc$truncated)) {
    warning(warningCondition(
      sprintf(paste0(
        "the complex was built with max_dim = %d below its largest clique; ",
        "b%d counts cycles whose filling simplices were not built and is an ",
        "upper bound. Rebuild with a larger `max_dim`."),
        sc$dimension, sc$dimension),
      class = "nestimate_truncated_complex", call = NULL
    ))
  }
  .compute_betti(sc)
}

#' @noRd
.compute_betti <- function(sc) {
  max_d <- sc$dimension
  dims <- vapply(sc$simplices, function(s) length(s) - 1L, integer(1))
  by_dim <- lapply(0:max_d, function(d) sc$simplices[dims == d])

  # rank of the boundary map from the d-simplices to their (d - 1)-faces
  boundary_rank <- function(d) {
    k_simplices <- by_dim[[d + 1L]]
    km1_simplices <- by_dim[[d]]
    if (length(k_simplices) == 0L || length(km1_simplices) == 0L) return(0L) # nocov

    km1_keys <- vapply(km1_simplices, function(s) {
      paste(sort(s), collapse = ",")
    }, character(1))

    # one row per (simplex j, omitted vertex i): the face it leaves and its sign
    cells <- do.call(rbind, lapply(seq_along(k_simplices), function(j) {
      simplex <- sort(k_simplices[[j]])
      omitted <- seq_along(simplex)
      face_keys <- vapply(omitted, function(i) {
        paste(simplex[-i], collapse = ",")
      }, character(1))
      cbind(row = match(face_keys, km1_keys), col = j,
            sign = (-1)^(omitted + 1L))
    }))
    cells <- cells[!is.na(cells[, "row"]), , drop = FALSE]

    bmat <- matrix(0, nrow = length(km1_simplices),
                   ncol = length(k_simplices))
    bmat[cells[, c("row", "col"), drop = FALSE]] <- cells[, "sign"]
    qr(bmat)$rank
  }
  boundary_ranks <- c(0L, vapply(seq_len(max_d), boundary_rank, integer(1)),
                      0L)

  betti <- vapply(0:max_d, function(d) {
    n_k <- length(by_dim[[d + 1L]])
    nullity_k <- n_k - boundary_ranks[d + 1L]
    rank_kp1 <- if (d < max_d) boundary_ranks[d + 2L] else 0L
    as.integer(max(nullity_k - rank_kp1, 0L))
  }, integer(1))

  names(betti) <- paste0("b", 0:max_d)
  betti
}

#' Euler Characteristic
#'
#' @description
#' `euler_characteristic()` returns the Euler characteristic of a simplicial
#' complex, \eqn{\chi = \sum_{k=0}^{d} (-1)^k f_k}, where \eqn{f_k} is the
#' number of k-simplices and \eqn{d} the dimension of the complex. By the
#' Euler-Poincare theorem it also equals \eqn{\sum_{k} (-1)^k \beta_k}, the
#' alternating sum of the [betti_numbers()].
#'
#' @param sc A `simplicial_complex` from [build_simplicial()].
#' @return A single integer.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' euler_characteristic(sc)
#'
#' @export
euler_characteristic <- function(sc) {
  .assert("`sc` must be a simplicial_complex (see build_simplicial())" =
            inherits(sc, "simplicial_complex"))
  signs <- (-1L)^(seq_along(sc$f_vector) - 1L)
  as.integer(sum(signs * sc$f_vector))
}

# =========================================================================
# Persistent homology
#
# Algorithm: full boundary-matrix reduction over Z/2 (Edelsbrunner, Letscher
# & Zomorodian 2000). Filtered complex is built once at full graph;
# filtration values are assigned per simplex (vertices at 0, k-simplex at
# max_w - min edge weight in sigma for the clique filtration; max pairwise
# distance in sigma for the VR filtration). Persistence pairs are read off the
# reduction directly - the previous Betti-difference heuristic that mispaired
# features born/dying between adjacent grid steps is gone.
# =========================================================================

#' @noRd
.filter_clique_complex <- function(mat, max_dim = 3L) {
  # Clique filtration in similarity (descending) semantics.
  # filt_asc(sigma) = max_w - min(edge weights in sigma), vertex = 0.
  # So a high-weight simplex has a small filt_asc (enters early in ascending).
  n <- nrow(mat)
  nodes <- rownames(mat) %||% paste0("V", seq_len(n))
  max_w <- max(mat)

  adj <- mat > 0
  diag(adj) <- FALSE

  if (!any(adj)) {
    simps <- as.list(seq_len(n))
    return(list(
      simplices = simps, dim = integer(n), filt_asc = numeric(n),
      key = as.character(seq_len(n)), nodes = nodes,
      max_filt = 0, max_w = max_w, mode = "clique"
    ))
  }

  all_simp <- .find_all_cliques(adj, max_dim)
  dims <- vapply(all_simp, function(s) length(s) - 1L, integer(1))
  filt <- vapply(seq_along(all_simp), function(j) {
    s <- all_simp[[j]]
    if (length(s) == 1L) return(0)
    pairs <- utils::combn(s, 2L)
    max_w - min(mat[cbind(pairs[1L, ], pairs[2L, ])])
  }, numeric(1))

  ord <- order(filt, dims)
  simplices <- all_simp[ord]
  dims <- dims[ord]
  filt <- filt[ord]
  keys <- vapply(simplices, function(s) paste(s, collapse = ","), character(1))

  list(
    simplices = simplices, dim = dims, filt_asc = filt,
    key = keys, nodes = nodes,
    max_filt = if (length(filt) > 0L) max(filt) else 0,
    max_w = max_w, mode = "clique"
  )
}

#' @noRd
.filter_vr_complex <- function(d, max_dim = 3L, max_scale = NULL,
                               detect_truncation = FALSE) {
  # Vietoris-Rips filtration on a non-negative distance matrix.
  # filt(sigma) = max pairwise distance in sigma; vertex = 0.
  .assert(is.matrix(d), nrow(d) == ncol(d))
  n <- nrow(d)
  nodes <- rownames(d) %||% paste0("V", seq_len(n))
  d <- pmax(d, t(d))
  diag(d) <- 0
  if (any(d < 0, na.rm = TRUE)) {
    stop(errorCondition(
      "VR filtration requires a non-negative distance matrix.",
      class = "nestimate_bad_input", call = NULL))
  }
  finite_d <- d
  finite_d[!is.finite(finite_d)] <- NA_real_
  cap <- if (is.null(max_scale)) {
    if (all(is.na(finite_d))) 0 else max(finite_d, na.rm = TRUE)
  } else {
    .assert(is.numeric(max_scale), length(max_scale) == 1L,
              !is.na(max_scale), max_scale >= 0)
    max_scale
  }

  # Edges within cap. d(i,j) == 0 for i != j is a valid pseudometric case
  # (duplicate points / equivalence classes), so include zero-distance
  # off-diagonal edges and rely on diag(adj) <- FALSE to exclude self-loops.
  adj <- !is.na(finite_d) & finite_d >= 0 & finite_d <= cap
  diag(adj) <- FALSE
  all_simp <- .find_all_cliques(adj, max_dim)
  # only build_simplicial() reports truncation; persistence does not need it
  truncated <- isTRUE(detect_truncation) &&
    .clique_exceeds(adj, all_simp, max_dim)
  dims <- vapply(all_simp, function(s) length(s) - 1L, integer(1))
  filt <- vapply(seq_along(all_simp), function(j) {
    s <- all_simp[[j]]
    if (length(s) == 1L) return(0)
    pairs <- utils::combn(s, 2L)
    max(d[cbind(pairs[1L, ], pairs[2L, ])])
  }, numeric(1))

  ord <- order(filt, dims)
  simplices <- all_simp[ord]
  dims <- dims[ord]
  filt <- filt[ord]
  keys <- vapply(simplices, function(s) paste(s, collapse = ","), character(1))

  list(
    simplices = simplices, dim = dims, filt_asc = filt,
    key = keys, nodes = nodes,
    max_filt = if (length(filt) > 0L) max(filt) else 0,
    max_w = cap, mode = "vr", truncated = truncated
  )
}

#' @noRd
.fc_from_filtered_complex <- function(sc, max_dim = 3L, max_scale = NULL) {
  # Convert a simplicial_complex with attached $filtration into the internal
  # filtered-complex shape consumed by .persistence_pairs_z2(). Honors max_dim
  # by dropping simplices above that dimension and re-orders by (filt, dim) so
  # boundary reduction is well-defined.
  .assert(
    "`x` must be a simplicial_complex with one filtration value per simplex" =
      inherits(sc, "simplicial_complex") && !is.null(sc$filtration) &&
        length(sc$filtration) == length(sc$simplices)
  )
  simplices <- sc$simplices
  filt <- as.numeric(sc$filtration)
  dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))

  # Drop above max_dim
  keep <- dims <= max_dim
  simplices <- simplices[keep]
  filt <- filt[keep]
  dims <- dims[keep]

  # Apply max_scale cap if requested (and the complex is a VR build)
  mode <- if (identical(sc$type, "vr")) "vr" else "clique"
  if (!is.null(max_scale)) {
    .assert(is.numeric(max_scale), length(max_scale) == 1L,
              !is.na(max_scale), max_scale >= 0)
    keep <- filt <= max_scale
    simplices <- simplices[keep]
    filt <- filt[keep]
    dims <- dims[keep]
  }

  # Order by (filt asc, dim asc) so faces precede cofaces at the same filt
  ord <- order(filt, dims)
  simplices <- simplices[ord]
  filt <- filt[ord]
  dims <- dims[ord]
  keys <- vapply(simplices, function(s) paste(sort(s), collapse = ","),
                 character(1))

  max_w <- if (mode == "vr") {
    if (!is.null(sc$max_scale)) sc$max_scale
    else if (length(filt) > 0L) max(filt) else 0
  } else {
    if (length(filt) > 0L) max(filt) else 0
  }

  list(
    simplices = simplices, dim = dims, filt_asc = filt,
    key = keys, nodes = sc$nodes,
    max_filt = if (length(filt) > 0L) max(filt) else 0,
    max_w = max_w, mode = mode
  )
}

#' @noRd
.persistence_pairs_z2 <- function(fc) {
  # Standard left-to-right boundary-matrix reduction over Z/2.
  # The j-loop is sequential by construction (column j depends on reduced
  # columns 1..j-1) - this is the package's second documented for-loop
  # exception alongside the permutation loop in sequence_compare.R.
  simplices <- fc$simplices
  dims <- fc$dim
  keys <- fc$key
  filt <- fc$filt_asc
  N <- length(simplices)

  if (N == 0L) {
    empty <- data.frame(dimension = integer(0), birth = numeric(0),
                        death = numeric(0), persistence = numeric(0),
                        stringsAsFactors = FALSE)
    return(list(pairs = empty, essential = empty))
  }

  key_to_idx <- setNames(seq_len(N), keys)

  # Boundary columns: k-simplex (k>=1) has (k+1) (k-1)-faces.
  D <- lapply(seq_len(N), function(j) {
    if (dims[j] < 1L) return(integer(0))
    s <- simplices[[j]]
    face_keys <- vapply(seq_along(s), function(i) {
      paste(s[-i], collapse = ",")
    }, character(1))
    sort.int(as.integer(key_to_idx[face_keys]))
  })

  low_to_col <- integer(N) # low_to_col[r] = column j with low(j)=r, 0 if none
  paired_b <- integer(0L)
  paired_d <- integer(0L)

  for (j in seq_len(N)) {
    col <- D[[j]]
    while (length(col) > 0L) {
      l <- col[length(col)]
      i <- low_to_col[l]
      if (i == 0L) break
      # XOR with D[[i]]: symmetric difference of sorted integer vectors
      col <- sort.int(c(setdiff(col, D[[i]]), setdiff(D[[i]], col)))
    }
    D[[j]] <- col
    if (length(col) > 0L) {
      l <- col[length(col)]
      low_to_col[l] <- j
      paired_b <- c(paired_b, l)
      paired_d <- c(paired_d, j)
    }
  }

  essential_idx <- setdiff(seq_len(N), c(paired_b, paired_d))

  pairs_df <- if (length(paired_b) == 0L) {
    data.frame(dimension = integer(0), birth = numeric(0),
               death = numeric(0), persistence = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    data.frame(
      dimension = dims[paired_b],
      birth = filt[paired_b],
      death = filt[paired_d],
      persistence = filt[paired_d] - filt[paired_b],
      stringsAsFactors = FALSE
    )
  }
  pairs_df <- pairs_df[pairs_df$persistence > 0, , drop = FALSE]

  essential_df <- if (length(essential_idx) == 0L) {
    data.frame(dimension = integer(0), birth = numeric(0),
               death = numeric(0), persistence = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    data.frame(
      dimension = dims[essential_idx],
      birth = filt[essential_idx],
      death = Inf, persistence = Inf,
      stringsAsFactors = FALSE
    )
  }
  list(pairs = pairs_df, essential = essential_df)
}

#' @noRd
.betti_curve_from_pairs <- function(pers, thresholds, max_dim, mode) {
  grid <- expand.grid(threshold = thresholds, dimension = 0:max_dim,
                      KEEP.OUT.ATTRS = FALSE)
  if (nrow(pers) == 0L) {
    grid$betti <- 0L
    return(grid[, c("threshold", "dimension", "betti")])
  }
  # Clique mode: descending thresholds; alive at t iff birth >= t AND death < t.
  # VR mode: ascending thresholds; alive at t iff birth <= t AND death > t.
  alive_per_row <- if (mode == "clique") {
    vapply(seq_len(nrow(grid)), function(k) {
      t <- grid$threshold[k]
      sub <- pers[pers$dimension == grid$dimension[k], , drop = FALSE]
      sum(sub$birth >= t & sub$death < t)
    }, integer(1))
  } else {
    vapply(seq_len(nrow(grid)), function(k) {
      t <- grid$threshold[k]
      sub <- pers[pers$dimension == grid$dimension[k], , drop = FALSE]
      sum(sub$birth <= t & sub$death > t)
    }, integer(1))
  }
  grid$betti <- alive_per_row
  grid[, c("threshold", "dimension", "betti")]
}

#' Persistent Homology
#'
#' @description
#' `persistent_homology()` follows the holes of a network's simplicial
#' complex as the complex is grown edge by edge: each component, loop or
#' void is born when it first appears and dies when it is filled in, and the
#' result lists every such feature with its birth, death and persistence (how
#' long it lasts). Two filtrations are available:
#' \describe{
#'   \item{`type = "clique"`}{A weighted clique filtration on similarities:
#'     edges enter from the strongest to the weakest, and a simplex enters
#'     at the smallest absolute weight among its edges, so the thresholds run
#'     from high to low.}
#'   \item{`type = "vr"`}{A Vietoris-Rips filtration on a non-negative
#'     distance matrix: a simplex enters at the largest distance among its
#'     nodes, so the thresholds run from low to high. `max_scale` caps the
#'     distances considered.}
#' }
#' A `simplicial_complex` from `build_simplicial(type = "vr")`, which stores
#' the scale at which each simplex enters, is also read directly.
#'
#' Features with long persistence last over a wide range of thresholds;
#' features that die soon after they are born are commonly read as noise.
#' Essential features never die: the component that remains at the end, and
#' any hole that no simplex in the filtration fills. [bottleneck_distance()]
#' compares the diagrams of two networks, and [persistence_landscape()]
#' turns a diagram into functions that can be averaged.
#'
#' @details
#' The persistence pairs come from the full reduction of the boundary
#' matrix over \eqn{\mathbb{Z}/2} (Edelsbrunner, Letscher & Zomorodian,
#' 2000): each feature is paired with the simplex whose entry creates it
#' (birth) and the simplex whose entry destroys it (death). Pairs born and
#' killed at the same value are dropped. The diagram is exact; `n_steps`
#' only sets the grid on which the Betti curve is reported.
#'
#' In the clique filtration the weight matrix is made symmetric with the
#' larger absolute weight of the two directions, its diagonal is not read,
#' and every node enters at the largest weight. An essential feature has
#' `death = 0` and `persistence = birth`; the Betti curve is read on
#' `n_steps` thresholds from the largest weight down to 1% of it. In the
#' Vietoris-Rips filtration the matrix is made symmetric with the larger
#' distance, every node enters at 0, and an essential feature has
#' `death = Inf` and `persistence = Inf`; the curve runs from 0 to the
#' largest entry value. The filtration holds simplices one dimension above
#' `max_dim`, so a feature of dimension `max_dim` can die (the Ripser
#' convention).
#'
#' @param x The network: a square weight matrix (similarities for
#'   `type = "clique"`, non-negative distances for `type = "vr"`), a `tna`
#'   model, a `netobject` or `cograph_network` (a `net_hon` among them), or
#'   a `simplicial_complex` from `build_simplicial(type = "vr")`. The
#'   complex's stored entry values are reduced as they are, so its own type
#'   is used and `type` is ignored;
#'   a complex without them (a clique or pathway complex) raises an error.
#'   A complex built with a `max_dim` below its largest clique warns with
#'   class `nestimate_truncated_complex`, since its top-dimension features
#'   cannot die, when the complex's dimension is at most `max_dim` and the
#'   truncated part lies within `max_scale`. All-zero weights raise an error
#'   in the clique filtration, and negative distances in the Vietoris-Rips
#'   one. So does a Vietoris-Rips input in which no edge can enter: every
#'   off-diagonal distance missing or infinite, or `max_scale = 0` with all
#'   distances positive. An all-zero distance matrix gives a single
#'   component, and distances that all exceed `max_scale` give one essential
#'   component per node. For the `print()` and `plot()` methods: a
#'   `persistent_homology` object.
#' @param n_steps The number of thresholds on which the Betti curve is
#'   reported, a positive whole number. Default `20`. The persistence
#'   diagram does not depend on it.
#' @param max_dim The highest feature dimension reported, a non-negative
#'   whole number. Default `3`. A `simplicial_complex` input supplies only
#'   the simplices it holds.
#' @param type The filtration: `"clique"` (default), on similarity weights,
#'   or `"vr"` (Vietoris-Rips), on distances.
#' @param max_scale For `type = "vr"`, and for a Vietoris-Rips
#'   `simplicial_complex`: the largest distance at which an edge enters;
#'   longer edges and the simplices that hold them are left out. Default
#'   `NULL`, the largest finite distance (for a complex, every simplex it
#'   holds).
#'
#' @return A `persistent_homology` object. `print()` shows the threshold
#'   range, the number of features per dimension and the three
#'   longest-lived; `plot()` draws the Betti curve and the persistence
#'   diagram; [bottleneck_distance()] and [persistence_landscape()] read it.
#'   The object's fields:
#' \describe{
#'   \item{betti_curve}{Data frame with one row per threshold and dimension:
#'     `threshold`, `dimension` and `betti`, the number of features alive at
#'     that threshold.}
#'   \item{persistence}{Data frame with one row per feature, sorted by
#'     decreasing persistence: `dimension`, `birth`, `death` and
#'     `persistence`. Essential features are included, with `death = 0` in
#'     the clique filtration and `death = Inf` in the Vietoris-Rips one.}
#'   \item{thresholds}{The thresholds of the Betti curve.}
#'   \item{mode}{`"clique"` or `"vr"`.}
#' }
#'
#' @references
#' Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2000). Topological
#' persistence and simplification. In \emph{Proceedings of the 41st Annual
#' Symposium on Foundations of Computer Science}, 454-463. Journal version:
#' \emph{Discrete & Computational Geometry} (2002) \strong{28}, 511-533.
#'
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' ph <- persistent_homology(mat, n_steps = 10)
#' print(ph)
#'
#' @export
persistent_homology <- function(x, n_steps = 20L, max_dim = 3L,
                                type = "clique", max_scale = NULL) {
  .assert(
    "`n_steps` must be a single positive whole number" =
      is.numeric(n_steps) && length(n_steps) == 1L &&
        !is.na(n_steps) && n_steps >= 1L && n_steps == as.integer(n_steps),
    "`max_dim` must be a single non-negative whole number" =
      is.numeric(max_dim) && length(max_dim) == 1L &&
        !is.na(max_dim) && max_dim >= 0 && max_dim == as.integer(max_dim)
  )
  type <- match.arg(type, c("clique", "vr"))
  if (inherits(x, "simplicial_complex") && is.null(x$filtration)) {
    # a clique or pathway complex records no scale at which a simplex enters
    stop(errorCondition(paste0(
      "This ", x$type, " complex carries no filtration, so its persistence ",
      "cannot be read. Pass the network or weight matrix itself, or a ",
      "complex built with build_simplicial(type = \"vr\")."),
      class = "nestimate_bad_input", call = NULL))
  }

  # Filtered-complex handoff: build_simplicial(type = "vr") attaches a
  # $filtration vector. Consume it directly instead of rebuilding from a
  # matrix - this is the workflow advertised by the build_simplicial docs.
  if (inherits(x, "simplicial_complex") && !is.null(x$filtration)) {
    if (isTRUE(x$truncated) && x$dimension <= max_dim &&
        .sc_truncated_within(x, max_scale)) {
      warning(warningCondition(
        sprintf(paste0(
          "the complex was built with max_dim = %d below its largest clique, ",
          "so its %d-dimensional classes cannot die and are reported as ",
          "essential. Rebuild it with `max_dim = %d`."),
          x$dimension, x$dimension, max_dim + 1L),
        class = "nestimate_truncated_complex", call = NULL
      ))
    }
    fc <- .fc_from_filtered_complex(x, max_dim = max_dim + 1L,
                                    max_scale = max_scale)
  } else if (type == "vr") {
    d <- if (is.matrix(x)) x else .sc_extract_matrix(x)
    fc <- .filter_vr_complex(d, max_dim = max_dim + 1L, max_scale = max_scale)
    if (fc$max_w == 0 && fc$max_filt == 0 &&
        all(fc$dim == 0L)) {
      stop(errorCondition(paste0(
        "No edge enters the Vietoris-Rips filtration: every off-diagonal ",
        "distance is missing, infinite or larger than `max_scale`."),
        class = "nestimate_bad_input", call = NULL))
    }
  } else {
    mat <- .sc_extract_matrix(x)
    .assert(
      "the network must be a square numeric weight matrix" =
        (is.numeric(mat) || is.logical(mat)) && nrow(mat) == ncol(mat)
    )
    mat <- abs(mat)
    mat <- pmax(mat, t(mat))
    diag(mat) <- 0
    if (max(mat) == 0) {
      stop(errorCondition("All weights are zero; cannot build filtration.",
                        class = "nestimate_bad_input", call = NULL))
    }
    fc <- .filter_clique_complex(mat, max_dim = max_dim + 1L)
  }

  # Every filtration above holds simplices up to max_dim + 1, so a
  # max_dim-dimensional class can die (the Ripser convention); building only
  # to max_dim left every top-dimension cycle unfilled and essential.
  red <- .persistence_pairs_z2(fc)

  # Translate to user-facing scale and assemble persistence table
  if (fc$mode == "clique") {
    max_w <- fc$max_w
    pairs <- red$pairs
    if (nrow(pairs) > 0L) {
      bd_asc <- pairs$birth
      dd_asc <- pairs$death
      pairs$birth <- max_w - bd_asc
      pairs$death <- max_w - dd_asc
      pairs$persistence <- pairs$birth - pairs$death
    }
    ess <- red$essential
    if (nrow(ess) > 0L) {
      ess$birth <- max_w - ess$birth
      ess$death <- 0
      ess$persistence <- ess$birth
    }
    thresholds <- seq(max_w, max_w * 0.01, length.out = n_steps)
  } else {
    pairs <- red$pairs
    ess <- red$essential
    # Keep essential death = Inf so the Betti curve correctly counts them as
    # alive at the final threshold (death > t holds for any finite t). The
    # plot path caps Inf at max_filt for display.
    if (nrow(ess) > 0L) {
      ess$persistence <- Inf
    }
    thresholds <- seq(0, max(fc$max_filt, .Machine$double.eps),
                      length.out = n_steps)
  }

  persistence <- rbind(pairs, ess)
  persistence <- persistence[order(-persistence$persistence), , drop = FALSE]

  # Report homology up to max_dim; the (max_dim + 1)-simplices only kill
  persistence <- persistence[persistence$dimension <= max_dim, , drop = FALSE]
  rownames(persistence) <- NULL

  betti_curve <- .betti_curve_from_pairs(persistence, thresholds, max_dim,
                                         fc$mode)

  structure(list(
    betti_curve = betti_curve,
    persistence = persistence,
    thresholds = thresholds,
    mode = fc$mode
  ), class = "persistent_homology")
}

# =========================================================================
# Simplicial centrality
# =========================================================================

#' Simplicial Degree
#'
#' @description
#' `simplicial_degree()` counts, for each node of a simplicial complex, how
#' many simplices of each dimension contain it: edges (dimension 1),
#' triangles (dimension 2) and so on. It is the higher-order counterpart of
#' a node's degree, and says which nodes take part in many group-level
#' structures rather than only in pairs.
#'
#' The counts describe the complex as built, so they move with the
#' threshold or pathways [build_simplicial()] was given. A node in a clique
#' of m nodes is in `choose(m - 1, k)` of its k-simplices, so one large
#' clique raises the node's counts in every dimension up to the clique's
#' own;
#' compare nodes within a dimension.
#'
#' @param sc A `simplicial_complex` from [build_simplicial()].
#' @param normalized Whether to divide each count of dimension k by the
#'   largest possible one, `choose(n_nodes - 1, k)`, the number of
#'   k-simplices a node can be in. Default `FALSE`.
#'
#' @return A data frame with one row per node, sorted by decreasing `total`:
#'   `node`, then `d0`, `d1`, ... up to the dimension of the complex (the
#'   number of simplices of that dimension containing the node; `d0` is
#'   always 1), and `total`, the sum of `d1` onwards.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' simplicial_degree(sc)
#'
#' @export
simplicial_degree <- function(sc, normalized = FALSE) {
  .assert(
    "`sc` must be a simplicial_complex (see build_simplicial())" =
      inherits(sc, "simplicial_complex"),
    "`normalized` must be TRUE or FALSE" =
      (is.logical(normalized) || is.numeric(normalized)) &&
        length(normalized) == 1L && !is.na(normalized)
  )
  n <- sc$n_nodes
  max_d <- sc$dimension

  dims <- vapply(sc$simplices, function(s) length(s) - 1L, integer(1))

  # cell (v, d + 1) counts the d-simplices that contain node v
  members <- unlist(sc$simplices, use.names = FALSE)
  member_dim <- rep(dims, times = dims + 1L)
  mat <- matrix(
    tabulate((member_dim * n) + members, nbins = n * (max_d + 1L)),
    nrow = n, ncol = max_d + 1L
  )

  if (normalized && n > 1L) {
    denom <- choose(n - 1L, 0:max_d)
    denom[!(denom > 0)] <- 1
    mat <- sweep(mat, 2L, denom, "/")
  }

  df <- as.data.frame(mat)
  names(df) <- paste0("d", 0:max_d)
  df <- cbind(data.frame(node = sc$nodes, stringsAsFactors = FALSE), df)
  df$total <- rowSums(mat[, -1L, drop = FALSE])
  df[order(-df$total), ]
}

# =========================================================================
# Q-analysis
# =========================================================================

#' Q-Analysis
#'
#' @description
#' `q_analysis()` describes how the maximal simplices of a simplicial complex
#' (those not contained in a larger one) hang together, after Atkin (1974).
#' Two maximal simplices are q-near when they share a face of dimension at
#' least q (q + 1 shared nodes), and q-connected when a chain of q-near
#' simplices joins them. The result has two parts:
#' \itemize{
#'   \item \strong{Q-vector}: \eqn{Q_q}, the number of q-connected
#'     components among the maximal simplices of dimension at least q, for
#'     each q from the largest simplex dimension down to 0. Simplices of lower
#'     dimension do not exist at level q and are not counted.
#'   \item \strong{Structure vector}: here the dimension of the largest
#'     simplex each node belongs to. (Atkin uses "structure vector" for the
#'     Q-vector itself; the name is kept for compatibility.)
#' }
#'
#' A Q-vector of 1 at every level means the maximal simplices stay in one
#' piece however many shared nodes are required. Requiring more shared
#' nodes splits the complex more, but fewer simplices take part at higher
#' q, so the count can rise or fall with q; `print()` reports the lowest
#' level at which there is more than one. Like
#' the complex itself, the vectors depend on the threshold or pathways
#' [build_simplicial()] was given.
#'
#' @param sc A `simplicial_complex` from [build_simplicial()].
#'
#' @return A `q_analysis` object. `print()` shows the components at each
#'   level, the level at which the complex first splits and the structure
#'   vector; `plot()` draws the two vectors. The object's fields:
#'   \describe{
#'     \item{q_vector}{Named integer vector, one entry per level from the
#'       highest down (`q_2`, `q_1`, `q_0` for a complex of dimension 2):
#'       the number of q-connected components at that level.}
#'     \item{max_q}{The largest simplex dimension, the highest level.}
#'     \item{structure_vector}{Named integer vector, one entry per node: the
#'       dimension of the largest simplex containing it.}
#'   }
#'
#' @references
#' Atkin, R. H. (1974). \emph{Mathematical Structure in Human Affairs}.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' q_analysis(sc)
#'
#' @export
q_analysis <- function(sc) {
  .assert("`sc` must be a simplicial_complex (see build_simplicial())" =
            inherits(sc, "simplicial_complex"))

  simplices <- sc$simplices
  dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))

  # Structure vector: max simplex dimension per node
  sv <- vapply(seq_len(sc$n_nodes), function(v) {
    d <- vapply(simplices, function(s) {
      if (v %in% s) length(s) - 1L else -1L
    }, integer(1))
    max(d)
  }, integer(1))
  names(sv) <- sc$nodes

  # Find maximal simplices
  n_s <- length(simplices)
  is_maximal <- vapply(seq_len(n_s), function(i) {
    si <- simplices[[i]]
    !any(vapply(seq_len(n_s), function(j) {
      if (j == i || dims[j] <= dims[i]) return(FALSE)
      all(si %in% simplices[[j]])
    }, logical(1)))
  }, logical(1))

  maximal <- simplices[is_maximal]
  n_max <- length(maximal)
  max_q <- if (n_max > 0L) max(vapply(maximal, length, integer(1))) - 1L else 0L

  q_levels <- max_q:0

  if (n_max <= 1L) {
    q_vec <- setNames(rep(1L, length(q_levels)), paste0("q_", q_levels))
    return(structure(list(
      q_vector = q_vec, max_q = max_q,
      structure_vector = sv
    ), class = "q_analysis"))
  }

  # Shared face dimension between pairs
  shared_dim <- vapply(seq_len(n_max), function(i) {
    vapply(seq_len(n_max), function(j) {
      if (i == j) return(-1L)
      length(intersect(maximal[[i]], maximal[[j]])) - 1L
    }, integer(1))
  }, integer(n_max))

  # Q_q counts components among the maximal simplices of dimension >= q only
  max_dims <- vapply(maximal, length, integer(1)) - 1L
  q_vec <- vapply(q_levels, function(q) {
    active <- max_dims >= q
    adj_q <- shared_dim[active, active, drop = FALSE] >= q
    diag(adj_q) <- FALSE
    .count_components(adj_q)
  }, integer(1))
  names(q_vec) <- paste0("q_", q_levels)

  structure(list(
    q_vector = q_vec,
    max_q = max_q,
    structure_vector = sv
  ), class = "q_analysis")
}

#' @noRd
.count_components <- function(adj) {
  n <- nrow(adj)
  visited <- logical(n)
  n_comp <- 0L
  # A loop: a breadth-first search, where each start depends on the nodes the
  # earlier searches already reached.
  for (start in seq_len(n)) {
    if (visited[start]) next
    n_comp <- n_comp + 1L
    queue <- start
    visited[start] <- TRUE
    while (length(queue) > 0L) {
      v <- queue[1L]
      queue <- queue[-1L]
      nbrs <- which(adj[v, ] & !visited)
      visited[nbrs] <- TRUE
      queue <- c(queue, nbrs)
    }
  }
  n_comp
}

# =========================================================================
# Verification helper
# =========================================================================

#' Verify Simplicial Complex Against igraph
#'
#' @description
#' `verify_simplicial()` builds the clique complex of `mat` with
#' [build_simplicial()] and checks it against two independent references:
#' its simplices
#' against the cliques `igraph::cliques()` finds in the same thresholded
#' graph, and its Euler characteristic against the alternating sum of its
#' Betti numbers (the Euler-Poincare theorem). It prints the result of each
#' check and the Betti numbers, and is meant for testing. It needs
#' \pkg{igraph}, and raises an error without it.
#'
#' @details
#' The complex is built with the default `max_dim = 10`, while
#' `igraph::cliques()` has no size limit, so a graph with a clique of more
#' than 11 nodes reports a mismatch, with a warning of class
#' `nestimate_truncated_complex` from [betti_numbers()].
#'
#' @param mat A square weight matrix, or a network object whose weight
#'   matrix [build_simplicial()] reads (a `netobject`, `tna` model,
#'   `cograph_network` or `net_hon`).
#' @param threshold The smallest absolute edge weight an edge needs to
#'   enter, as in [build_simplicial()]. Default `0`.
#'
#' @return A list, invisibly: `cliques_match` (whether the simplices equal
#'   the cliques from `igraph::cliques()` exactly), `n_simplices_ours`,
#'   `n_simplices_igraph`, `betti`, `euler` and `f_vector`. The two checks
#'   are printed as a side effect.
#' @examplesIf requireNamespace("igraph", quietly = TRUE)
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' verify_simplicial(mat, threshold = 0.3)
#' @export
verify_simplicial <- function(mat, threshold = 0) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required for verification.", call. = FALSE) # nocov
  }

  mat <- .sc_extract_matrix(mat)
  sc <- build_simplicial(mat, threshold = threshold)

  # igraph clique count
  adj <- .sc_threshold_adjacency(mat, threshold, inclusive = TRUE)
  diag(adj) <- FALSE
  adj <- adj | t(adj)
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                            diag = FALSE)
  ig_cliques <- igraph::cliques(g)

  # Compare sorted simplex sets
  our_keys <- sort(vapply(sc$simplices, function(s) {
    paste(sort(s), collapse = ",")
  }, character(1)))
  ig_keys <- sort(vapply(ig_cliques, function(cl) {
    paste(sort(as.integer(cl)), collapse = ",")
  }, character(1)))

  betti <- betti_numbers(sc)
  euler <- euler_characteristic(sc)

  result <- list(
    cliques_match = identical(our_keys, ig_keys),
    n_simplices_ours = length(sc$simplices),
    n_simplices_igraph = length(ig_cliques),
    betti = betti,
    euler = euler,
    f_vector = sc$f_vector
  )

  clique_ok <- result$cliques_match
  dims <- seq_along(betti) - 1L
  euler_from_betti <- as.integer(sum((-1L)^dims * betti))
  euler_ok <- euler == euler_from_betti

  b_str <- paste(sprintf("%s=%d", names(betti), betti), collapse = " ")
  cat(sprintf("  Cliques:  %s (%d simplices)\n",
              if (clique_ok) "MATCH" else "MISMATCH", result$n_simplices_ours))
  cat(sprintf("  Betti:    %s\n", b_str))
  cat(sprintf("  Euler:    %d (Euler-Poincare: %s)\n",
              euler, if (euler_ok) "VERIFIED" else "FAILED"))

  invisible(result)
}

# =========================================================================
# Print methods
# =========================================================================

#' @rdname build_simplicial
#' @param combined For the `plot()` method: `TRUE` (default) draws the four
#'   panels (f-vector, Betti numbers, simplicial degree per node and degree
#'   by dimension) as one 2 x 2 figure, which needs \pkg{gridExtra};
#'   `FALSE` returns the four ggplots in a named list (`f_vector`, `betti`,
#'   `degree`, `degree_heatmap`) to print, save or arrange separately. A
#'   complex without edges (dimension 0) has no simplicial degree to draw
#'   and raises an error of class `nestimate_bad_input`.
#' @return The `print()` method returns its input invisibly. `plot()`
#'   returns the drawn grid grob invisibly with `combined = TRUE`, and the
#'   named list of four ggplots invisibly with `combined = FALSE`.
#' @export
print.simplicial_complex <- function(x, ...) {
  labels <- c("clique" = "Clique Complex",
              "pathway" = "Pathway Complex",
              "vr"      = "Vietoris-Rips Complex")
  cat(labels[x$type] %||% "Simplicial Complex", "\n")

  betti <- .compute_betti(x)
  chi <- euler_characteristic(x)

  cat(sprintf("  %d nodes, %d simplices, dimension %d\n",
              x$n_nodes, x$n_simplices, x$dimension))
  cat(sprintf("  Density: %.1f%%  |  Mean dim: %.2f  |  Euler: %d\n",
              x$density * 100, x$mean_dim, chi))

  # f-vector: compact
  f_str <- paste(sprintf("f%d=%d", seq_along(x$f_vector) - 1L,
                          x$f_vector), collapse = " ")
  cat(sprintf("  f-vector: (%s)\n", f_str))

  # Betti: only non-zero
  nz <- which(betti > 0)
  if (length(nz) == 0L) {
    cat("  Betti: all zero (contractible)\n") # nocov
  } else {
    b_str <- paste(sprintf("%s=%d", names(betti)[nz], betti[nz]),
                   collapse = " ")
    cat(sprintf("  Betti: %s\n", b_str))
  }

  if (x$n_nodes <= 15L) {
    cat("  Nodes:", paste(x$nodes, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @rdname persistent_homology
#' @param ... For the `print()` and `plot()` methods: ignored.
#' @param combined For the `plot()` method: `TRUE` (default) draws the Betti
#'   curve and the persistence diagram side by side, which needs
#'   \pkg{gridExtra}; `FALSE` returns them as a named list of ggplots
#'   (`betti_curve`, `persistence`).
#' @return The `print()` method returns its input invisibly. `plot()`
#'   returns the drawn grid grob invisibly with `combined = TRUE`, and the
#'   named list of two ggplots invisibly with `combined = FALSE`. The
#'   diagram plots each feature's birth against its death; essential
#'   features are drawn at the end of the filtration, at `death = 0` in the
#'   clique filtration and, for the Vietoris-Rips one, at the largest finite
#'   value in the diagram or on the threshold grid.
#' @export
print.persistent_homology <- function(x, ...) {
  cat("Persistent Homology\n")
  # the clique filtration descends through weights, the Vietoris-Rips one
  # ascends through distances; essential classes never die (death 0 in the
  # clique convention, Inf in the Vietoris-Rips one)
  vr <- identical(x$mode, "vr")
  ends <- if (vr) range(x$thresholds) else rev(range(x$thresholds))
  cat(sprintf("  %d filtration steps [%.4f \u2192 %.4f]\n",
              length(x$thresholds), ends[1L], ends[2L]))

  if (nrow(x$persistence) > 0L) {
    dims <- sort(unique(x$persistence$dimension))
    parts <- vapply(dims, function(d) {
      sub <- x$persistence[x$persistence$dimension == d, ]
      n_p <- if (vr) sum(is.infinite(sub$death)) else sum(sub$death == 0)
      sprintf("b%d: %d (%d persistent)", d, nrow(sub), n_p)
    }, character(1))
    cat("  Features:", paste(parts, collapse = "  |  "), "\n")

    # Top 3 only
    top <- head(x$persistence[x$persistence$persistence > 0, ], 3)
    if (nrow(top) > 0L) {
      cat("  Longest-lived:\n")
      cat(sprintf("    b%d: %.4f \u2192 %.4f (life: %.4f)\n",
                  top$dimension, top$birth, top$death, top$persistence),
          sep = "")
    }
  }
  invisible(x)
}

#' @rdname q_analysis
#' @param x For the `print()` and `plot()` methods: a `q_analysis` object.
#' @param ... For the `print()` and `plot()` methods: ignored.
#' @param combined For the `plot()` method: `TRUE` (default) draws the
#'   Q-vector and the structure vector side by side, which needs
#'   \pkg{gridExtra}; `FALSE` returns them as a named list of ggplots
#'   (`q_vector`, `structure_vector`).
#' @return The `print()` method returns its input invisibly. `plot()`
#'   returns the drawn grid grob invisibly with `combined = TRUE`, and the
#'   named list of two ggplots invisibly with `combined = FALSE`.
#' @export
print.q_analysis <- function(x, ...) {
  cat(sprintf("Q-Analysis (max q = %d)\n", x$max_q))

  # Q-vector as compact line
  q_levels <- as.integer(sub("^q_", "", names(x$q_vector)))
  q_str <- paste(sprintf("q%d:%d", q_levels, x$q_vector), collapse = " ")
  cat(sprintf("  Components: %s\n", q_str))

  # Fragmentation: lowest q where components > 1. q_vector is stored from
  # the highest q down, so scan it in ascending q.
  ascending <- order(q_levels)
  asc_levels <- q_levels[ascending]
  asc_counts <- x$q_vector[ascending]
  frag_idx <- which(asc_counts > 1L)[1]
  if (!is.na(frag_idx)) {
    if (frag_idx > 1L) {
      cat(sprintf("  Fragments at q = %d (%d \u2192 %d components)\n",
                  asc_levels[frag_idx], asc_counts[frag_idx - 1L],
                  asc_counts[frag_idx]))
    } else {
      cat(sprintf("  Fragments at q = %d (%d components)\n",
                  asc_levels[frag_idx], asc_counts[frag_idx]))
    }
  } else {
    cat("  Fully connected at all q levels\n")
  }

  # Structure vector: compact sorted
  sv <- sort(x$structure_vector, decreasing = TRUE)
  sv_str <- paste(sprintf("%s:%d", names(sv), sv), collapse = " ")
  cat(sprintf("  Structure: %s\n", sv_str))
  invisible(x)
}

# =========================================================================
# Plot methods
# =========================================================================

.sc_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "grey40",
                                             size = base_size - 2),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' @rdname build_simplicial
#' @export
plot.simplicial_complex <- function(x, combined = TRUE, ...) {
  .assert("`combined` must be TRUE or FALSE" =
            is.logical(combined) && length(combined) == 1L && !is.na(combined))
  if (x$dimension < 1L) {
    # the degree panels show participation in simplices of dimension 1 and up
    stop(errorCondition(paste0(
      "This complex has no edges (dimension 0), so there is no simplicial ",
      "degree to plot. Build it with a lower `threshold`."),
      class = "nestimate_bad_input", call = NULL))
  }

  deg <- simplicial_degree(x)
  betti <- .compute_betti(x)

  # --- Panel 1: f-vector ---
  fdf <- data.frame(dim = factor(seq_along(x$f_vector) - 1L),
                     count = as.integer(x$f_vector))
  p1 <- ggplot2::ggplot(fdf, ggplot2::aes(x = dim, y = count)) +
    ggplot2::geom_col(fill = "#4A7FB5", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3,
                        size = 3.5) +
    ggplot2::labs(title = "f-vector",
                  subtitle = "Simplices per dimension",
                  x = "Dimension", y = "Count") +
    .sc_theme()

  # --- Panel 2: degree ranking ---
  deg$node <- factor(deg$node, levels = deg$node)
  p2 <- ggplot2::ggplot(deg, ggplot2::aes(x = node, y = total,
                                            fill = total)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_gradient(low = "#81B1D3", high = "#E8734A") +
    ggplot2::geom_text(ggplot2::aes(label = total), vjust = -0.3,
                        size = 3.2) +
    ggplot2::labs(title = "Simplicial Degree",
                  subtitle = "Higher-order participation (dim 1+)",
                  x = NULL, y = "Total") +
    .sc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
                                                         hjust = 1))

  # --- Panel 3: degree heatmap ---
  d_cols <- paste0("d", seq_len(x$dimension))
  deg_long <- stats::reshape(deg[, c("node", d_cols)],
                              direction = "long",
                              varying = d_cols,
                              v.names = "count", timevar = "dim",
                              times = seq_len(x$dimension))
  deg_long$node <- factor(deg_long$node, levels = rev(deg$node))

  p3 <- ggplot2::ggplot(deg_long, ggplot2::aes(x = factor(dim), y = node,
                                                  fill = count)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = count), size = 3.2) +
    ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#E8734A",
                                  guide = "none") +
    ggplot2::labs(title = "Degree by Dimension",
                  subtitle = "Simplex participation per node",
                  x = "Dimension", y = NULL) +
    .sc_theme()

  # --- Panel 4: Betti numbers ---
  bdf <- data.frame(dim = factor(seq_along(betti) - 1L),
                     betti = as.integer(betti))
  b_subtitle <- sprintf("Euler characteristic: %d", euler_characteristic(x))
  p4 <- ggplot2::ggplot(bdf, ggplot2::aes(x = dim, y = betti)) +
    ggplot2::geom_col(fill = "#6AAB9C", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = betti), vjust = -0.3,
                        size = 3.5) +
    ggplot2::labs(title = "Betti Numbers",
                  subtitle = b_subtitle,
                  x = "Dimension", y = expression(beta[k])) +
    .sc_theme()

  panels <- list(f_vector = p1, betti = p4, degree = p2,
                 degree_heatmap = p3)
  if (!combined) return(invisible(panels))
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("combined = TRUE requires the gridExtra package.", call. = FALSE) # nocov
  }
  combined_grob <- gridExtra::arrangeGrob(p1, p4, p2, p3, ncol = 2,
                                          padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined_grob)
  invisible(combined_grob)
}

#' @rdname persistent_homology
#' @export
plot.persistent_homology <- function(x, combined = TRUE, ...) {
  .assert("`combined` must be TRUE or FALSE" =
            is.logical(combined) && length(combined) == 1L && !is.na(combined))

  filt <- x$betti_curve
  filt$dim_label <- factor(paste0("B", filt$dimension))
  # the Vietoris-Rips filtration runs over distances, the clique one over weights
  vr <- identical(x$mode, "vr")
  curve_subtitle <- if (vr) {
    "Betti numbers across distance thresholds"
  } else {
    "Betti numbers across inclusive weight thresholds"
  }
  curve_x <- if (vr) "Distance Threshold" else "Weight Threshold"

  # --- Panel 1: Betti curve ---
  p1 <- ggplot2::ggplot(filt, ggplot2::aes(x = threshold, y = betti,
                                              color = dim_label)) +
    ggplot2::geom_step(linewidth = 1.1, direction = "vh") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(title = "Betti Curve",
                  subtitle = curve_subtitle,
                  x = curve_x, y = "Betti Number",
                  color = NULL) +
    .sc_theme() +
    ggplot2::theme(legend.position = "top")

  # --- Panel 2: Persistence diagram ---
  pers <- x$persistence[x$persistence$persistence > 0, ]
  # Cap both Inf death and Inf persistence for display (essential VR classes).
  # ggplot's continuous size scale drops Inf-valued rows, which would silently
  # erase essential features (e.g., the surviving H_0 in VR mode). Cap to the
  # finite max so every row renders.
  inf_death <- !is.finite(pers$death)
  inf_pers  <- !is.finite(pers$persistence)
  if (any(inf_death) || any(inf_pers)) {
    finite_max <- max(c(pers$birth,
                        pers$death[!inf_death],
                        pers$persistence[!inf_pers],
                        x$thresholds), na.rm = TRUE)
    pers$death[inf_death]      <- finite_max
    pers$persistence[inf_pers] <- finite_max
  }

  if (nrow(pers) > 0L) {
    pers$dim_label <- factor(paste0("B", pers$dimension))
    lim <- max(c(pers$birth, pers$death)) * 1.15

    p2 <- ggplot2::ggplot(pers, ggplot2::aes(x = birth, y = death,
                                                color = dim_label,
                                                size = persistence)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                            color = "grey60") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_size_continuous(range = c(1.5, 6), guide = "none") +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
      ggplot2::labs(title = "Persistence Diagram",
                    subtitle = "Boundary-matrix reduction over Z/2",
                    x = "Birth", y = "Death", color = NULL) +
      .sc_theme() +
      ggplot2::theme(legend.position = "top")
  } else { # nocov start
    p2 <- ggplot2::ggplot() +
      ggplot2::labs(title = "Persistence Diagram",
                    subtitle = "No features detected") +
      .sc_theme() # nocov end
  }

  panels <- list(betti_curve = p1, persistence = p2)
  if (!combined) return(invisible(panels))
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("combined = TRUE requires the gridExtra package.", call. = FALSE) # nocov
  }
  combined_grob <- gridExtra::arrangeGrob(p1, p2, ncol = 2,
                                          padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined_grob)
  invisible(combined_grob)
}

#' @rdname q_analysis
#' @export
plot.q_analysis <- function(x, combined = TRUE, ...) {
  .assert("`combined` must be TRUE or FALSE" =
            is.logical(combined) && length(combined) == 1L && !is.na(combined))

  # --- Panel 1: Q-vector ---
  qdf <- data.frame(q = as.integer(sub("^q_", "", names(x$q_vector))),
                    components = as.integer(x$q_vector))

  p1 <- ggplot2::ggplot(qdf, ggplot2::aes(x = q, y = components)) +
    ggplot2::geom_step(linewidth = 1.2, color = "#4A7FB5",
                        direction = "vh") +
    ggplot2::geom_point(size = 3, color = "#E8734A") +
    ggplot2::geom_text(ggplot2::aes(label = components), vjust = -1,
                        size = 3.5) +
    ggplot2::scale_x_continuous(breaks = qdf$q) +
    ggplot2::labs(title = "Q-Vector",
                  subtitle = "Connected components at each q-level",
                  x = "q (shared face dimension)", y = "Components") +
    .sc_theme()

  # --- Panel 2: Structure vector ---
  sv <- x$structure_vector
  svdf <- data.frame(node = names(sv), dim = as.integer(sv),
                      stringsAsFactors = FALSE)
  svdf <- svdf[order(-svdf$dim, svdf$node), ]
  svdf$node <- factor(svdf$node, levels = svdf$node)

  p2 <- ggplot2::ggplot(svdf, ggplot2::aes(x = node, y = dim, fill = dim)) +
    ggplot2::geom_col(show.legend = FALSE, width = 0.7) +
    ggplot2::scale_fill_gradient(low = "#81B1D3", high = "#E8734A") +
    ggplot2::geom_text(ggplot2::aes(label = dim), vjust = -0.3,
                        size = 3.5) +
    ggplot2::labs(title = "Structure Vector",
                  subtitle = "Highest simplex dimension per node",
                  x = NULL, y = "Max Dimension") +
    .sc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
                                                         hjust = 1))

  panels <- list(q_vector = p1, structure_vector = p2)
  if (!combined) return(invisible(panels))
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("combined = TRUE requires the gridExtra package.", call. = FALSE) # nocov
  }
  combined_grob <- gridExtra::arrangeGrob(p1, p2, ncol = 2,
                                          padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined_grob)
  invisible(combined_grob)
}
