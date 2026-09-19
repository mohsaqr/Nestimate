# ==============================================================================
# macro_network() -- the cluster-level network, optionally with one cluster
# expanded back into its own states.
#
# An mcml's macro layer collapses every cluster to a single node, which is the
# point of it -- but it also hides the one cluster you may care about. This
# verb rebuilds the network at a mixed resolution: the chosen cluster keeps its
# member states as separate nodes, every other cluster stays collapsed.
#
# It re-counts from the sequence data rather than splitting the macro matrix,
# because a k x k aggregate cannot be disaggregated after the fact.
# ==============================================================================

#' Cluster-Level Network, With One Cluster Expanded
#'
#' Returns the macro (cluster-level) network of an \code{mcml}, optionally
#' with one or more clusters expanded back into their member states. Every
#' other cluster stays collapsed to a single node, so the result is a network
#' at mixed resolution: the cluster of interest in detail, its context in
#' summary.
#'
#' @param x An \code{mcml} built from sequence data, or an
#'   \code{mcml_pc} from \code{\link{build_mcml_pc}}. A matrix-derived
#'   \code{mcml} carries no node-level data and cannot be expanded.
#' @param expand Names of clusters to expand into their member states.
#'   \code{NULL} (default) collapses every cluster, reproducing the macro
#'   layer; \code{"all"} or \code{TRUE} expands every cluster, so each state
#'   is its own node.
#' @param method Estimator passed to \code{\link{build_network}}. Default
#'   \code{"relative"} (row-normalised transitions).
#' @param ... Further arguments passed to \code{\link{build_network}}.
#'
#' @return A \code{netobject} (also a \code{cograph_network}) whose nodes are
#'   the collapsed clusters plus the member states of any expanded cluster,
#'   with weights re-counted from the sequence data by
#'   \code{\link{build_network}}. \code{$node_groups} is a two-column data
#'   frame (\code{node}, \code{group}) mapping every node to its cluster, and
#'   the same labels are a factor in \code{$nodes$groups}, so the result
#'   plots grouped; an expanded cluster's states each map to that cluster, a
#'   collapsed cluster maps to itself. \code{$expanded} records the cluster
#'   names that were expanded (\code{NULL} when none were).
#' @seealso \code{\link{build_mcml}}, \code{\link{as_tna}}
#' @examples
#' seqs <- data.frame(
#'   t1 = c("A", "C", "A", "B"), t2 = c("B", "D", "C", "A"),
#'   t3 = c("C", "A", "D", "C"), stringsAsFactors = FALSE
#' )
#' mc <- build_mcml(seqs, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))
#' macro_network(mc)                    # every cluster collapsed
#' macro_network(mc, expand = "G2")     # G2 shown as C and D
#' @export
macro_network <- function(x, expand = NULL, method = "relative", ...) {
  stopifnot(
    "`x` must be an mcml or an mcml_pc" =
      inherits(x, "mcml") || inherits(x, "mcml_pc"),
    "`expand` must be a character vector, TRUE, \"all\" or NULL" =
      is.null(expand) || is.character(expand) || isTRUE(expand)
  )
  # A psychometric fit already holds its cluster-level network as a full
  # netobject -- there is nothing to re-count, so the verb hands it back.
  # Mixed-resolution expansion is a transition-network operation: it re-counts
  # sequences, which an association network has none of.
  if (inherits(x, "mcml_pc")) {
    if (!is.null(expand)) {
      stop(errorCondition(
        paste0("`expand` is not available for a psychometric (mcml_pc) fit: ",
               "there are no sequences to re-count. Use as_networks() to get ",
               "the macro network together with every within-cluster network."),
        class = "nestimate_no_expand", call = NULL))
    }
    return(x$macro)
  }
  members <- x$cluster_members
  # "all" / TRUE expands every cluster, so the macro carries each state as its
  # own node. Spelling it out as names(x$cluster_members) would make the
  # caller reach into the object for something the verb can name itself.
  if (isTRUE(expand) || identical(expand, "all")) expand <- names(members)
  unknown <- setdiff(expand, names(members))
  if (length(unknown) > 0L) {
    stop("Unknown cluster(s) in `expand`: ",
         paste(utils::head(unknown, 5L), collapse = ", "),
         ". Available: ", paste(names(members), collapse = ", "),
         call. = FALSE)
  }
  data <- attr(x, "htna_source")
  if (is.null(data)) {
    stop(errorCondition(
      paste0("This mcml carries no node-level sequence data (it was built ",
             "from a matrix or aggregate), so its clusters cannot be ",
             "expanded. Rebuild it from sequences."),
      class = "nestimate_no_expand_source", call = NULL))
  }

  # A state maps to its own label when its cluster is expanded, to the cluster
  # name otherwise. Recoding the sequences and re-counting gives exact
  # transitions at the mixed resolution.
  cluster_of <- stats::setNames(rep(names(members), lengths(members)),
                                unlist(members, use.names = FALSE))
  map  <- cluster_of
  keep <- unlist(members[names(members) %in% expand], use.names = FALSE)
  map[keep] <- keep

  recoded <- as.data.frame(
    lapply(data, function(col) unname(map[as.character(col)])),
    stringsAsFactors = FALSE
  )
  net <- build_network(recoded, method = method, ...)

  # Label every node with the cluster it belongs to, so the result plots
  # grouped: an expanded state carries its cluster, a collapsed cluster itself.
  # `cluster_of` still maps every state to its cluster, so an expanded state
  # is grouped with its siblings rather than with itself.
  node_cluster <- ifelse(net$nodes$label %in% names(members),
                         net$nodes$label,
                         unname(cluster_of[net$nodes$label]))
  node_cluster[is.na(node_cluster)] <- net$nodes$label[is.na(node_cluster)]
  net$node_groups <- data.frame(node = net$nodes$label, group = node_cluster,
                                stringsAsFactors = FALSE)
  net$nodes$groups <- factor(node_cluster, levels = unique(node_cluster))
  net$expanded <- expand
  net
}
