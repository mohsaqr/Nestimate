# ---- Hypergraph structural measures (HON-9) ------------------------------
# Node-level, hyperedge-level, and global statistics for a net_hypergraph.
# All measures are fast matrix operations on the incidence matrix B.

#' Structural measures for a hypergraph
#'
#' @description
#' `hypergraph_measures()` computes the structural statistics of a
#' hypergraph at three levels: for each node (how many hyperedges it
#' belongs to, how large they are, and how many it shares with each other
#' node), for each hyperedge and pair of hyperedges (their sizes and how
#' much they overlap), and for the hypergraph as a whole (density, mean
#' hyperedge size, the intersection profile). It reads the `net_hypergraph`
#' objects that [build_hypergraph()] and [bipartite_groups()] return.
#'
#' The measures describe the hypergraph as observed: they carry no
#' uncertainty and test nothing. The intersection profile shows whether
#' hyperedges mostly share no or few nodes or share substantial cores (Do
#' et al., 2020; Lee, Choe & Shin, 2021). To rank nodes by their position
#' in the hypergraph, use [hypergraph_centrality()]; to read it as an
#' ordinary pairwise network, project it with [clique_expansion()].
#'
#' @param hg A `net_hypergraph`, from [build_hypergraph()] or
#'   [bipartite_groups()]. Any other object raises an error.
#'
#' @return A `hypergraph_measures` object, a named list; `print()` shows the
#'   density, mean hyperedge size and pairwise participation, the minimum,
#'   median and maximum of hyperdegree, node strength and hyperedge size,
#'   and the intersection profile, and returns its input invisibly.
#'   Its components, node-level first:
#'
#' \describe{
#'   \item{**Node-level (one value per node)**}{}
#'   \item{`hyperdegree`}{The number of hyperedges containing each node.}
#'   \item{`node_strength`}{Total participation: for node \eqn{i},
#'     \eqn{\sum_{e \ni i} |e|}{sum of |e| over the hyperedges e
#'     containing i}. A node in many large hyperedges has
#'     high strength.}
#'   \item{`max_edge_size`}{The size of the largest hyperedge containing
#'     each node (0 for a node in none).}
#'   \item{`co_degree`}{An `n_nodes` x `n_nodes` matrix:
#'     `co_degree[i, j] = |{e : i, j in e}|`, the number of hyperedges
#'     containing both nodes \eqn{i} and \eqn{j}. The diagonal is zero.}
#'   \item{**Hyperedge-level (one value per hyperedge, or m x m)**}{}
#'   \item{`edge_sizes`}{The hyperedge sizes \eqn{|e|}.}
#'   \item{`edge_pairwise_overlap`}{An `m` x `m` matrix of
#'     `|e_i intersect e_j|`. The diagonal is zero.}
#'   \item{`overlap_coefficient`}{An `m` x `m` matrix of
#'     `|e_i intersect e_j| / min(|e_i|, |e_j|)`: how much of the smaller
#'     hyperedge lies in the larger. The diagonal is zero.}
#'   \item{`jaccard`}{An `m` x `m` matrix of the symmetric overlap index
#'     `|e_i intersect e_j| / |e_i union e_j|`. The diagonal is zero.}
#'   \item{**Global (scalars)**}{}
#'   \item{`density`}{For a hypergraph whose hyperedges all have size `k`:
#'     `m / choose(n, k)`. For mixed sizes:
#'     `sum(|e|) / (n * m)`, the mean fraction of nodes per hyperedge.}
#'   \item{`avg_edge_size`}{The mean of `edge_sizes`.}
#'   \item{`size_distribution`}{The count of hyperedges of each size,
#'     passed through from `hg`.}
#'   \item{`intersection_profile`}{The number of hyperedge pairs at each
#'     intersection size that occurs, named `overlap_<size>` (sizes with no
#'     pairs are left out); empty with fewer than two hyperedges.}
#'   \item{`pairwise_participation`}{The fraction of node pairs that
#'     appear together in at least one hyperedge (0 with fewer than two
#'     nodes).}
#'   \item{`n_nodes`, `n_hyperedges`}{The numbers of nodes and
#'     hyperedges.}
#' }
#'
#' @details
#' All measures are matrix operations on the binary incidence matrix
#' \eqn{B = (b_{ij})}, where \eqn{b_{ij} = 1} if node \eqn{i} is in
#' hyperedge \eqn{j}; a weighted incidence is binarised first:
#' \itemize{
#'   \item `hyperdegree = rowSums(B)`,
#'         `edge_sizes = colSums(B)`
#'   \item `co_degree = tcrossprod(B)` (with zero diagonal)
#'   \item `edge_pairwise_overlap = crossprod(B)` (with zero diagonal)
#'   \item `overlap_coefficient[i, j] = overlap[i, j] /
#'         min(edge_sizes[i], edge_sizes[j])`
#'   \item `jaccard[i, j] = overlap[i, j] /
#'         (edge_sizes[i] + edge_sizes[j] - overlap[i, j])`
#' }
#' A ratio whose denominator is 0 is set to 0.
#'
#' A hypergraph without hyperedges (`n_hyperedges == 0`) gives zeros for
#' the node measures, empty hyperedge matrices, a density and pairwise
#' participation of 0, and an `avg_edge_size` of `NA`.
#'
#' @seealso [build_hypergraph()], [bipartite_groups()],
#'   [hypergraph_centrality()], [clique_expansion()].
#'
#' @examples
#' df <- data.frame(
#'   player  = c("A", "B", "C", "A", "B", "D", "C", "D"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3")
#' )
#' hg <- bipartite_groups(df, "player", "session")
#' m  <- hypergraph_measures(hg)
#' print(m)   # node-level, hyperedge-level and global measures
#'
#' @references
#' Lee, G., Choe, M., & Shin, K. (2021). How do hyperedges overlap in
#' real-world hypergraphs? Patterns, measures, and generators.
#' \emph{Proceedings of the Web Conference 2021}, 3396-3407.
#'
#' Do, M. T., Yoon, S., Hooi, B., & Shin, K. (2020). Structural patterns
#' and generative models of real-world hypergraphs. arXiv:2006.07060.
#'
#' @export
hypergraph_measures <- function(hg) {
  .assert(inherits(hg, "net_hypergraph"))

  B <- hg$incidence
  n <- hg$n_nodes
  m <- hg$n_hyperedges

  if (m == 0L) {
    zeros_n <- stats::setNames(numeric(n), hg$nodes)
    empty_mm <- matrix(0, 0L, 0L)
    return(structure(
      list(
        hyperdegree            = stats::setNames(integer(n), hg$nodes),
        node_strength          = zeros_n,
        max_edge_size          = stats::setNames(integer(n), hg$nodes),
        co_degree              = matrix(0, n, n,
                                        dimnames = list(hg$nodes, hg$nodes)),
        edge_sizes             = integer(0),
        edge_pairwise_overlap  = empty_mm,
        overlap_coefficient    = empty_mm,
        jaccard                = empty_mm,
        density                = 0,
        avg_edge_size          = NA_real_,
        size_distribution      = hg$size_distribution,
        intersection_profile   = integer(0),
        pairwise_participation = 0,
        n_nodes                = n,
        n_hyperedges           = 0L
      ),
      class = "hypergraph_measures"
    ))
  }

  B_bin <- (B > 0) * 1.0
  edge_sizes <- as.integer(colSums(B_bin))
  edge_names <- colnames(B) %||% paste0("h", seq_len(m))

  # ---- Node-level ----
  hyperdegree <- stats::setNames(as.integer(rowSums(B_bin)), hg$nodes)
  node_strength <- stats::setNames(as.numeric(B_bin %*% edge_sizes), hg$nodes)
  max_edge_size <- stats::setNames(
    vapply(seq_len(n), function(i) {
      mem <- B_bin[i, ] > 0
      if (any(mem)) max(edge_sizes[mem]) else 0L
    }, integer(1L)),
    hg$nodes
  )

  # ---- Pairwise co-degree (n x n) ----
  co_degree <- tcrossprod(B_bin)
  diag(co_degree) <- 0
  dimnames(co_degree) <- list(hg$nodes, hg$nodes)

  # ---- Edge-level pairwise overlap (m x m) ----
  edge_pairwise_overlap <- crossprod(B_bin)
  diag(edge_pairwise_overlap) <- 0
  dimnames(edge_pairwise_overlap) <- list(edge_names, edge_names)

  # ---- Overlap coefficient & Jaccard (m x m) ----
  size_min   <- outer(edge_sizes, edge_sizes, pmin)
  size_union <- outer(edge_sizes, edge_sizes, `+`) - edge_pairwise_overlap
  overlap_coefficient <- ifelse(size_min > 0,
                                edge_pairwise_overlap / size_min, 0)
  jaccard <- ifelse(size_union > 0,
                    edge_pairwise_overlap / size_union, 0)
  diag(overlap_coefficient) <- 0
  diag(jaccard) <- 0
  dimnames(overlap_coefficient) <- list(edge_names, edge_names)
  dimnames(jaccard) <- list(edge_names, edge_names)

  # ---- Global ----
  is_uniform <- length(unique(edge_sizes)) == 1L
  density <- if (is_uniform) {
    k <- edge_sizes[1L]
    if (k > n) 0 else m / choose(n, k)
  } else {
    sum(edge_sizes) / (n * m)
  }
  avg_edge_size <- mean(edge_sizes)

  intersection_profile <- if (m >= 2L) {
    upper <- edge_pairwise_overlap[upper.tri(edge_pairwise_overlap)]
    tab <- table(upper)
    out <- as.integer(tab)
    names(out) <- paste0("overlap_", names(tab))
    out
  } else {
    integer(0L)
  }

  pairwise_participation <- if (n >= 2L) {
    sum(co_degree[upper.tri(co_degree)] > 0) / choose(n, 2L)
  } else {
    0
  }

  structure(
    list(
      hyperdegree            = hyperdegree,
      node_strength          = node_strength,
      max_edge_size          = max_edge_size,
      co_degree              = co_degree,
      edge_sizes             = edge_sizes,
      edge_pairwise_overlap  = edge_pairwise_overlap,
      overlap_coefficient    = overlap_coefficient,
      jaccard                = jaccard,
      density                = density,
      avg_edge_size          = avg_edge_size,
      size_distribution      = hg$size_distribution,
      intersection_profile   = intersection_profile,
      pairwise_participation = pairwise_participation,
      n_nodes                = n,
      n_hyperedges           = m
    ),
    class = "hypergraph_measures"
  )
}

#' @rdname hypergraph_measures
#' @param x For the `print()` method: a `hypergraph_measures` object.
#' @param ... For the `print()` method: not used.
#' @export
print.hypergraph_measures <- function(x, ...) {
  cat(sprintf("Hypergraph measures: %d nodes, %d hyperedges\n",
              x$n_nodes, x$n_hyperedges))
  cat(sprintf("  Density:                %.4f\n", x$density))
  cat(sprintf("  Mean edge size:         %.2f\n",
              if (length(x$avg_edge_size)) x$avg_edge_size else NA))
  cat(sprintf("  Pairwise participation: %.4f\n",
              x$pairwise_participation))
  if (length(x$hyperdegree)) {
    cat(sprintf("\n  Hyperdegree:    min=%d  med=%.1f  max=%d\n",
                min(x$hyperdegree),
                stats::median(x$hyperdegree),
                max(x$hyperdegree)))
    cat(sprintf("  Node strength:  min=%.1f  med=%.1f  max=%.1f\n",
                min(x$node_strength),
                stats::median(x$node_strength),
                max(x$node_strength)))
  }
  if (length(x$edge_sizes)) {
    cat(sprintf("  Edge sizes:     min=%d  med=%.1f  max=%d\n",
                min(x$edge_sizes),
                stats::median(x$edge_sizes),
                max(x$edge_sizes)))
  }
  if (length(x$intersection_profile)) {
    cat("\nIntersection profile (pair counts by overlap size):\n")
    for (nm in names(x$intersection_profile)) {
      cat(sprintf("  %-12s: %d\n", nm, x$intersection_profile[[nm]]))
    }
  }
  invisible(x)
}
