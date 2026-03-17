#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @import tna
#' @importFrom stats cor var sd median quantile setNames aggregate
#'   approx spline cmdscale coef confint cov fitted lm na.omit
#'   p.adjust pnorm predict pt qnorm qt residuals rnorm runif
#'   rbinom rpois reshape symnum ecdf complete.cases
#' @importFrom utils head tail combn modifyList
#' @importFrom igraph graph_from_adjacency_matrix graph_from_edgelist
#'   V E vcount ecount degree strength distances get.edgelist
#'   as_edgelist as_adjacency_matrix neighbors shortest_paths
#'   all_shortest_paths is_directed set_vertex_attr set_edge_attr
#'   graph_from_data_frame delete_edges delete_vertices
#'   transitivity reciprocity closeness betweenness eigen_centrality
#'   page_rank hub_score authority_score edge_density dyad_census
#' @importFrom glasso glasso
#' @importFrom data.table data.table setDT rbindlist
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators()
}

# Global variable declarations for R CMD check
utils::globalVariables(c(
  "from", "to", "weight", "value", "name", "variable",
  "edge", "node", "state", "time", "onset", "terminus",
  "lower", "upper", "significant", "p_value", "ci_lower", "ci_upper",
  "metric", "centrality", "vertex", "group"
))
