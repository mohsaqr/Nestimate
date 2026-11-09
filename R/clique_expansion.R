# ---- Hypergraph -> pairwise network projection (HON-7) -------------------
# Standard clique expansion: for each hyperedge, every pair of members
# contributes a unit (or weighted) tally to the resulting pairwise network.
# Closes the I/O cycle with build_hypergraph() / bipartite_groups().

#' Clique expansion of a hypergraph
#'
#' @description
#' `clique_expansion()` projects a hypergraph onto an ordinary pairwise
#' network, its *clique expansion* (also called the "downgrade" of a
#' hypergraph to a dyadic graph): each hyperedge links every pair of its
#' members, and the weight between two nodes sums these links over the
#' hyperedges. It takes a `net_hypergraph` from [build_hypergraph()] or
#' [bipartite_groups()] and returns an undirected
#' [netobject][build_network]. With binary incidence, the weight `W[i, j]`
#' is the number of hyperedges containing both `i` and `j`; with weighted
#' incidence, it is the sum over hyperedges of the product of the two
#' members' incidence values.
#'
#' The projection keeps which pairs co-occur and how often, and drops the
#' groups themselves: three hyperedges `{A, B}`, `{B, C}` and `{A, C}` give
#' the same network as one hyperedge `{A, B, C}`. Read its edges as
#' pairwise co-occurrence; [hypergraph_measures()] and
#' [hypergraph_centrality()] work on the groups.
#'
#' @param hg The hypergraph: a `net_hypergraph` from [build_hypergraph()] or
#'   [bipartite_groups()].
#' @param weighted Whether to use the hypergraph's incidence values. With
#'   `TRUE` (default), a weighted hypergraph (from [bipartite_groups()] with
#'   `weight`) gives a weighted projection; with `FALSE`, the incidence is
#'   binarised first, so each weight counts the hyperedges the pair shares.
#'   For an unweighted hypergraph the two agree.
#'
#' @return A `netobject` (class `c("netobject", "cograph_network")`) with
#'   method `"clique_expansion"`, undirected. `print()`, `summary()`,
#'   [extract_edges()], [net_centrality()] and `cograph::splot()` read it as
#'   any network from [build_network()]. It carries no data, so the
#'   resampling verbs, such as [bootstrap_network()], and sequence
#'   clustering refuse it. Its weight matrix is
#'   `W = incidence %*% t(incidence)` with a zero diagonal, and it holds the
#'   standard `netobject` fields (`weights`, `nodes`, `edges`, `n_nodes`,
#'   `n_edges`, `meta`). `edges` has one row per non-zero upper-triangle
#'   cell, with integer `from`/`to` node indices and `weight`. `params`
#'   records `source`, `weighted`, `n_hyperedges` and
#'   `hypergraph_size_distribution`.
#'
#' @details
#' The weight matrix is computed in one BLAS call, `tcrossprod(incidence)`,
#' in `O(n_nodes^2 * n_hyperedges)` time. With binary incidence (an
#' unweighted hypergraph, or `weighted = FALSE`), its off-diagonal counts
#' are those of [cooccurrence()] on the same groups, which also keeps the
#' diagonal.
#'
#' @seealso [build_hypergraph()], [bipartite_groups()], [build_network()].
#'
#' @examples
#' df <- data.frame(
#'   player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
#' )
#' hg  <- bipartite_groups(df, player = "player", group = "session")
#' net <- clique_expansion(hg)
#' extract_edges(net, threshold = 1)
#'
#' @references
#' Tian, H., & Zafarani, R. (2024). Higher-order networks representation and
#' learning: A survey. \emph{ACM SIGKDD Explorations Newsletter} 26(1), 1-18.
#'
#' @note Experimental. Validated against `tcrossprod(incidence)` with a zero
#'   diagonal. No external R package exposes clique expansion as a primitive;
#'   the implementation restates the definition in one line.
#'
#' @export
clique_expansion <- function(hg, weighted = TRUE) {
  .assert(
    "`hg` must be a net_hypergraph (build_hypergraph/bipartite_groups)" =
      inherits(hg, "net_hypergraph"),
    "`weighted` must be TRUE or FALSE" =
      is.logical(weighted) && length(weighted) == 1L
  )

  if (hg$n_hyperedges == 0L) {
    weights <- matrix(0, hg$n_nodes, hg$n_nodes,
                      dimnames = list(hg$nodes, hg$nodes))
  } else {
    inc <- hg$incidence
    if (!weighted) {
      inc <- (inc > 0) * 1L
    }
    storage.mode(inc) <- "double"
    weights <- tcrossprod(inc)
    diag(weights) <- 0
    rownames(weights) <- colnames(weights) <- hg$nodes
  }

  net <- .wrap_netobject(weights,
                         method   = "clique_expansion",
                         directed = FALSE)
  net$params <- list(
    source                       = "clique_expansion",
    weighted                     = weighted,
    n_hyperedges                 = hg$n_hyperedges,
    hypergraph_size_distribution = hg$size_distribution
  )
  net
}
