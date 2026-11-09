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
#' @description
#' `macro_network()` returns the cluster-level (macro) network of an `mcml`,
#' optionally with one or more clusters expanded back into their member
#' states. Every other cluster stays collapsed to a single node, so the
#' result is a network at mixed resolution: the cluster of interest in
#' detail, its context in summary, with weights counted from the sequences.
#' For an `mcml_pc` from [build_mcml_pc()], it returns the cluster-level
#' network that fit already holds.
#'
#' @details
#' The weights are counted again from the sequences, after each state is
#' recoded to its cluster (or kept, for an expanded cluster), because a
#' cluster-level matrix cannot be split back into its states.
#'
#' @param x An `mcml` built from sequence data, or an `mcml_pc` from
#'   [build_mcml_pc()]. An `mcml` built from a matrix or an edge list holds
#'   no sequences to count, so it raises an error of class
#'   `nestimate_no_expand_source`, even with `expand = NULL`; a name in
#'   `expand` that is not a cluster is reported first.
#' @param expand The clusters to expand into their member states. `NULL`
#'   (default) collapses every cluster, which reproduces the macro layer;
#'   `"all"` or `TRUE` expands every cluster, so each state is its own
#'   node. A name that is not a cluster raises an error. For an `mcml_pc`,
#'   any value other than `NULL` raises an error of class
#'   `nestimate_no_expand`: there are no sequences to count again.
#' @param method The estimator passed to [build_network()]. Default
#'   `"relative"` (row-normalised transitions). For an `mcml_pc`, giving it
#'   raises an error: the estimator was set by [build_mcml_pc()].
#' @param scaling,threshold,params Passed to [build_network()] for the
#'   re-counted network: `scaling` (default `NULL`, none), `threshold`
#'   (default `0`) and `params` (default `list()`), as there. For an
#'   `mcml_pc`, giving any raises an error, as for `method`.
#' @param ... Further arguments for [build_network()], such as `start` or
#'   `end`. For an `mcml_pc`, giving any raises an error, as for `method`.
#'
#' @return For an `mcml`: a `netobject` (also a `cograph_network`) whose
#'   nodes are the collapsed clusters plus the member states of any
#'   expanded cluster, with weights counted by [build_network()]. It is
#'   printed, resampled (for example with [bootstrap_network()]) and drawn
#'   like any other `netobject`; its node grouping is read by cograph's
#'   group-aware plots, such as `cograph::plot_mtna()`. Beside the usual
#'   `netobject` fields it holds `node_groups`, a data frame with the
#'   columns `node` and `group` that maps each node to its cluster (an
#'   expanded cluster's states map to that cluster, a collapsed cluster maps
#'   to itself), the same labels as a factor in the `groups` column of
#'   `nodes`, and
#'   `expanded`, the names of the clusters expanded (`NULL` when none
#'   were).
#'
#'   For an `mcml_pc`: its cluster-level network, the `netobject` estimated
#'   by [build_mcml_pc()], unchanged.
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
macro_network <- function(x, expand = NULL, method = "relative",
                          scaling = NULL, threshold = 0, params = list(),
                          ...) {
  .assert(
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
    if (!missing(method) || !missing(scaling) || !missing(threshold) ||
        !missing(params) || ...length() > 0L) {
      stop(errorCondition(
        paste0(
          "`method`, `scaling`, `threshold`, `params` and `...` do not apply ",
          "to an mcml_pc: its cluster-level ",
          "network was estimated by build_mcml_pc(); set the estimator there."),
        class = "nestimate_bad_input", call = NULL))
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
    stop(errorCondition(paste0("Unknown cluster(s) in `expand`: ",
         paste(utils::head(unknown, 5L), collapse = ", "),
         ". Available: ", paste(names(members), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  data <- attr(x, "htna_source")
  if (is.null(data)) {
    stop(errorCondition(
      paste0("This mcml carries no node-level sequence data (only an mcml ",
             "built from sequence data keeps it, not one built from a ",
             "matrix, an edge list or a network), so its clusters cannot be ",
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
  net <- build_network(.mark_positions(recoded), method = method,
                       scaling = scaling, threshold = threshold,
                       params = params, ...)

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
