# ---- MCML for Partial-Correlation / Psychometric Networks (EXPERIMENTAL) ----

#' Multi-Cluster Multi-Level Aggregation for Psychometric Networks
#'
#' @description
#' \strong{Experimental.} Aggregates a node-level psychometric network
#' (correlation, partial correlation, or EBICglasso) into a cluster-level
#' macro network plus per-cluster within networks — the MCML view that
#' \code{\link{build_mcml}} provides for transition networks, adapted to
#' the statistics of undirected association networks. The API and the
#' exact aggregation formulas may change between releases.
#'
#' Unlike transition counts, partial correlations do not aggregate by
#' arithmetic: the submatrix of a pcor matrix is \emph{not} the pcor
#' network of the subsystem (the conditioning set changes), and averaging
#' pcor entries across blocks is descriptive only. The five
#' \code{aggregation} methods therefore have explicitly different
#' statuses:
#'
#' \describe{
#'   \item{\code{"average"} (descriptive)}{Macro edge A--B = mean of the
#'     signed node-level weights between members of A and members of B;
#'     macro diagonal = mean within-block off-diagonal weight. Needs only
#'     the weight matrix (works on data-less netobjects). Caveat: signed
#'     averaging can cancel opposite-sign edges; interpret as net average
#'     association, not connectivity strength.}
#'   \item{\code{"composite"} (re-estimated)}{Per-observation cluster
#'     scores = (optionally standardized) mean of member variables; the
#'     chosen \code{method} is then re-fit on the k composite columns,
#'     so the macro network is a genuine correlation / pcor / EBICglasso
#'     network among clusters. Requires raw data.}
#'   \item{\code{"loadings"} (re-estimated, connectivity-weighted)}{As
#'     \code{"composite"}, but member variables are weighted by their
#'     within-cluster connection strength in the node-level network (the
#'     absolute summed weight to the other members of their own cluster,
#'     normalized to sum to 1 per cluster). Nodes that anchor their
#'     cluster contribute more to its composite. This is Nestimate's own
#'     weighting — related in spirit to network loadings
#'     (Christensen & Golino 2021) but not a reimplementation of any
#'     EGA-family estimator. Requires raw data.}
#'   \item{\code{"rv"} (descriptive, multivariate)}{Macro edge A--B =
#'     Escoufier's RV coefficient between the member blocks — a matrix
#'     correlation in [0, 1] computed from the block covariance
#'     structure. No composites are formed and no estimator is re-fit,
#'     so nothing is lost to averaging; signs are not represented.
#'     Requires raw data.}
#'   \item{\code{"canonical"} (descriptive, multivariate)}{Macro edge
#'     A--B = the first canonical correlation between the member blocks:
#'     the strongest linear relationship any weighting of A's items can
#'     have with any weighting of B's items — an upper bound on what
#'     composite methods can recover. Requires raw data.}
#' }
#'
#' @details
#' \strong{Item diagnostics.} Whenever raw data or a node-level network
#' is available, every item's connection strength to \emph{every} cluster
#' is computed. The \code{$loadings} table reports, per item: its signed
#' own-cluster loading, its composite weight, its strongest cross-cluster
#' loading, and a \code{misfit} flag set when the cross-cluster loading
#' exceeds the own-cluster loading — evidence the item is assigned to the
#' wrong cluster. Misfit items trigger a warning; every aggregation
#' silently inherits a bad membership, so fix the assignment rather than
#' ignoring the flag.
#'
#' \strong{Reverse-keyed items.} With \code{signed = TRUE} (default),
#' items whose summed within-cluster association is negative are flipped
#' (their standardized values enter composites with weight sign -1), so a
#' reverse-keyed item reinforces its cluster composite instead of
#' cancelling it. Flips are reported in the \code{sign} column and via a
#' warning.
#'
#' \strong{Missing data.} Composites are per-row weighted means over the
#' \emph{observed} members (weights renormalized per row); rows with no
#' observed member yield \code{NA} and are dropped by the estimator with
#' a message. Node-level estimation applies the estimators' own
#' complete-case handling.
#'
#' \strong{Ordinal items.} \code{cor_method = "polychoric"} (requires the
#' \pkg{lavaan} package) estimates the node-level and within-cluster
#' networks from polychoric correlations — appropriate for Likert items.
#' Composites are continuous sums, so the macro re-estimation uses
#' Pearson correlations of the composites regardless.
#'
#' \strong{Within-cluster networks} follow \code{within}:
#' \code{"reestimate"} (default) re-fits the estimator on the member
#' columns alone — the honest conditional structure of the subsystem;
#' \code{"subnetwork"} slices the node-level weight matrix and is
#' descriptive (for pcor/glasso it retains conditioning on
#' out-of-cluster nodes). Modes without raw data force
#' \code{"subnetwork"}. Singleton clusters get no within network
#' (\code{NULL}).
#'
#' \strong{Uncertainty.} The composite/loadings macro network is a full
#' netobject carrying its composite data, so
#' \code{\link{bootstrap_network}(fit$macro)} (edge-weight CIs) and
#' \code{\link{vertex_bootstrap}(fit$macro)} (network-level CIs) work
#' directly; \code{\link{vertex_compare}(fit1$macro, fit2$macro)}
#' compares two groups. \code{\link{loading_stability}} quantifies how
#' stable the composite weights themselves are under case resampling.
#'
#' All constituent networks are undirected
#' (\code{meta$directed = FALSE}), so renderers that auto-detect
#' directedness (e.g. \code{cograph::plot_mcml()}) draw the result
#' without arrowheads.
#'
#' @param x A \code{netobject} estimated with an undirected association
#'   method (\code{cor}, \code{pcor}, \code{glasso}; aliases accepted) —
#'   its \code{$data} and \code{$method} are reused — or a numeric
#'   data.frame of raw observations (then \code{method} decides the
#'   estimator).
#' @param clusters Cluster membership: a named list of node-label vectors
#'   (names = cluster labels), or a vector of cluster labels named by
#'   node. Every node must be assigned to exactly one cluster.
#' @param aggregation Character. \code{"average"}, \code{"composite"},
#'   \code{"loadings"}, \code{"rv"}, or \code{"canonical"} (see
#'   Description). Default \code{"average"}.
#' @param method Character. Network estimator for the re-estimation
#'   paths and for data.frame input: \code{"pcor"} (default),
#'   \code{"glasso"}, or \code{"cor"} — the same vocabulary as
#'   \code{\link{build_network}}. Ignored (with the netobject's own
#'   method used instead) when \code{x} is a netobject and
#'   \code{method} is not given.
#' @param within Character. \code{"reestimate"} (default) or
#'   \code{"subnetwork"} (see Details).
#' @param weighting Character. How items are weighted inside their
#'   cluster composite (\code{aggregation = "composite"} only):
#'   \describe{
#'     \item{\code{"equal"}}{(default) 1/m per item — the scale as
#'       scored.}
#'     \item{\code{"strength"}}{Mean absolute connection to the other
#'       own-cluster members in the node-level network (what
#'       \code{aggregation = "loadings"} selects).}
#'     \item{\code{"eigen"}}{Leading-eigenvector weights of the
#'       within-cluster block of the node-level network — like strength
#'       but giving extra weight to items connected to other
#'       well-connected items.}
#'     \item{\code{"pca"}}{First principal component of the member
#'       items' correlation matrix — a data-statistical weighting,
#'       blind to the estimated network.}
#'     \item{\code{"factor"}}{Standardized loadings of a one-factor
#'       model per cluster — the classical latent-variable weighting.
#'       The extraction method is chosen by \code{fa_method} and the
#'       correlation input respects \code{cor_method} (so polychoric
#'       factor analysis of ordinal items is one call). Clusters with
#'       fewer than 3 items (or non-converging fits) fall back to
#'       \code{"pca"} with a warning.}
#'     \item{\code{"closeness"}, \code{"betweenness"}}{Within-cluster
#'       closeness / betweenness centrality of the item in the node-level
#'       network block (absolute weights). Betweenness can be all-zero in
#'       densely connected clusters; equal weights are used then, with a
#'       warning.}
#'     \item{\code{"expected_influence"}}{Mean \emph{signed} connection
#'       to the other own-cluster members (Robinaugh et al. 2016) — like
#'       strength, but opposite-sign connections subtract, and negative
#'       expected influence marks reverse-keyed items.}
#'     \item{\code{"specificity"}}{The misfit margin: own-cluster
#'       strength minus the strongest cross-cluster strength, floored at
#'       0. Items that belong as much to another cluster contribute
#'       nothing to their composite — the weighting twin of the
#'       \code{misfit} diagnostic.}
#'     \item{\code{"item_total"}}{Corrected item-total correlation:
#'       each item against the mean of the other (standardized) members —
#'       the classical scale-construction weighting.}
#'   }
#'   \strong{Custom weightings:} \code{weighting} also accepts a
#'   \emph{named numeric vector} (one entry per item; absolute values
#'   are normalized within each cluster, signs flip items when
#'   \code{signed = TRUE}) or a \emph{function}
#'   \code{function(W_block, data_block, nodes)} returning one numeric
#'   weight per item, evaluated per cluster.
#'
#'   Network-based and data-based weightings answer different questions;
#'   comparing them is informative — divergence means the network's view
#'   of the cluster differs from its latent-variable view. Schemes with
#'   inherently non-negative weights (equal, closeness, betweenness,
#'   specificity) keep the eigenvector-based item signs; sign-carrying
#'   schemes (eigen, pca, factor, expected_influence, item_total, custom)
#'   use their own.
#' @param scale Logical. Standardize member variables before forming
#'   composites (default \code{TRUE}). Ignored for \code{"average"},
#'   \code{"rv"}, and \code{"canonical"}.
#' @param cor_method Character. Correlation type for estimation from raw
#'   data: \code{"pearson"} (default), \code{"spearman"}, or
#'   \code{"polychoric"} (needs \pkg{lavaan}; see Details).
#' @param signed Logical. Flip reverse-keyed items in composites
#'   (default \code{TRUE}; see Details).
#' @param fa_method Character. Extraction method for
#'   \code{weighting = "factor"}:
#'   \describe{
#'     \item{\code{"ml"}}{(default) Maximum likelihood
#'       (\code{stats::factanal} on the \code{cor_method}-consistent
#'       correlation matrix).}
#'     \item{\code{"paf"}}{Iterated principal axis factoring (SMC
#'       start, communalities iterated on the reduced correlation
#'       matrix).}
#'     \item{\code{"minres"}}{Minimum residual / unweighted least
#'       squares (uniquenesses optimized to minimize squared
#'       off-diagonal residuals).}
#'     \item{\code{"cfa"}}{One-factor confirmatory model in
#'       \pkg{lavaan} (standardized loadings); with
#'       \code{cor_method = "polychoric"} the items are declared
#'       ordered, giving the categorical (DWLS) factor model.}
#'   }
#'   On well-behaved unidimensional clusters the four agree closely;
#'   divergence indicates Heywood-prone or non-unidimensional clusters.
#' @param ... Further arguments forwarded directly to the
#'   \code{fa_method} backend, exactly as that backend spells them
#'   (requires \code{weighting = "factor"}): lavaan arguments for
#'   \code{"cfa"} (\code{estimator = "WLSMV"}, \code{missing = "fiml"},
#'   \code{se = "robust"}, ...), \code{stats::factanal()} arguments for
#'   \code{"ml"}, and
#'   \code{max_iter} / \code{tol} for \code{"paf"}. Arguments managed
#'   internally (\code{model}, \code{data}, \code{covmat},
#'   \code{n.obs}, \code{factors}) are ignored with a warning — the
#'   model is always the one-factor model per cluster, because the
#'   composite needs exactly one weight per item.
#' @param id_col Character vector or NULL. Identifier column(s) to drop
#'   from data.frame input before analysis (e.g., the \code{rid}/actor
#'   columns produced by
#'   \code{\link{convert_sequence_format}(format = "frequency")}, whose
#'   output otherwise feeds this function directly as per-actor behavior
#'   profiles). Same convention as the association estimators.
#'
#' @return An object of class \code{"mcml_pc"} containing:
#' \describe{
#'   \item{macro}{Cluster-level netobject (k x k, undirected).}
#'   \item{clusters}{Named list of within-cluster netobjects
#'     (\code{NULL} for singleton clusters).}
#'   \item{cluster_members}{Named list of member node labels.}
#'   \item{loadings}{Tidy item-diagnostic data frame (one row per node):
#'     \code{node}, \code{cluster}, \code{loading} (signed own-cluster),
#'     \code{weight}, \code{sign}, \code{max_cross},
#'     \code{cross_cluster}, \code{misfit}. \code{NULL} only when no
#'     node-level network is available.}
#'   \item{node_network}{The node-level netobject the aggregation was
#'     based on (kept for diagnostics and \code{loading_stability}).}
#'   \item{data}{The raw data (data.frame) when available, else
#'     \code{NULL}.}
#'   \item{meta}{List: \code{aggregation}, \code{method},
#'     \code{within}, \code{scale}, \code{cor_method}, \code{signed},
#'     \code{n_nodes}, \code{n_clusters}, \code{cluster_sizes},
#'     \code{n_misfit}, \code{n_flipped}, \code{directed = FALSE},
#'     \code{source = "pc"}, \code{experimental = TRUE}.}
#' }
#'
#' @references
#' Escoufier, Y. (1973). Le traitement des variables vectorielles.
#' \emph{Biometrics}, 29(4), 751-760.
#'
#' Hotelling, H. (1936). Relations between two sets of variates.
#' \emph{Biometrika}, 28(3/4), 321-377.
#'
#' Robinaugh, D. J., Millner, A. J., & McNally, R. J. (2016). Identifying
#' highly influential nodes in the complicated grief network.
#' \emph{Journal of Abnormal Psychology}, 125(6), 747-757.
#'
#' Christensen, A. P., & Golino, H. (2021). On the equivalency of factor
#' and network loadings. \emph{Behavior Research Methods}, 53, 1563-1580.
#'
#' Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
#' correlation networks. \emph{Psychological Methods}, 23(4), 617-634.
#'
#' @seealso \code{\link{build_mcml}} for transition networks,
#'   \code{\link{loading_stability}} for composite-weight stability,
#'   \code{\link{bootstrap_network}} and \code{\link{vertex_bootstrap}}
#'   for uncertainty on the macro network.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' sigma <- matrix(0.15, 6, 6)
#' sigma[1:3, 1:3] <- 0.5
#' sigma[4:6, 4:6] <- 0.5
#' diag(sigma) <- 1
#' z <- matrix(rnorm(n * 6), n, 6) %*% chol(sigma)
#' df <- as.data.frame(z)
#' names(df) <- c("a1", "a2", "a3", "b1", "b2", "b3")
#' clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#'
#' fit <- build_mcml_pc(df, clusters, aggregation = "composite",
#'                      method = "cor")
#' fit$macro$weights
#' fit$loadings
#'
#' @export
build_mcml_pc <- function(x,
                          clusters,
                          aggregation = c("average", "composite",
                                          "loadings", "rv", "canonical"),
                          method = c("pcor", "glasso", "cor"),
                          within = c("reestimate", "subnetwork"),
                          weighting = c("equal", "strength", "eigen",
                                        "closeness", "betweenness",
                                        "expected_influence", "specificity",
                                        "pca", "factor", "item_total"),
                          scale = TRUE,
                          cor_method = c("pearson", "spearman",
                                         "polychoric"),
                          signed = TRUE,
                          id_col = NULL,
                          fa_method = c("ml", "paf", "minres", "cfa"),
                          ...) {
  aggregation <- match.arg(aggregation)
  weighting_given <- !missing(weighting)
  custom_weighting <- !is.character(weighting)
  if (custom_weighting) {
    if (!(is.function(weighting) ||
          (is.numeric(weighting) && !is.null(names(weighting))))) {
      stop("'weighting' must be a scheme name, a named numeric vector ",
           "of item weights, or a function(W_block, data_block, nodes).",
           call. = FALSE)
    }
  } else {
    weighting <- match.arg(weighting)
  }
  # "loadings" is the composite path with network-strength weighting
  if (aggregation == "loadings") {
    aggregation <- "composite"
    if (!weighting_given) weighting <- "strength"
  }
  within <- match.arg(within)
  cor_method <- match.arg(cor_method)
  fa_method <- match.arg(fa_method)
  fa_args <- list(...)
  if (length(fa_args) > 0 && (is.null(names(fa_args)) ||
                              any(!nzchar(names(fa_args))))) {
    stop("Arguments passed via ... must be named.", call. = FALSE)
  }
  if ("estimator" %in% names(fa_args) && fa_method != "cfa") {
    stop("'estimator' in ... is the lavaan estimator (e.g. \"WLSMV\") ",
         "and applies only to fa_method = \"cfa\". The network ",
         "estimator is selected by 'method'.", call. = FALSE)
  }
  if (length(fa_args) > 0 && !identical(weighting, "factor")) {
    stop("Unused argument(s): ", paste(names(fa_args), collapse = ", "),
         ". Arguments in ... are forwarded to the factor-analysis ",
         "backend and require weighting = \"factor\".", call. = FALSE)
  }
  if (length(fa_args) > 0) {
    protected <- switch(fa_method,
      cfa = c("model", "data"),
      ml = c("covmat", "n.obs", "factors", "x", "data"),
      paf = "R",
      minres = "R")
    bad <- intersect(names(fa_args), protected)
    if (length(bad) > 0) {
      warning("fa_args element(s) ignored (managed internally): ",
              paste(bad, collapse = ", "), call. = FALSE)
      fa_args[bad] <- NULL
    }
  }
  stopifnot(
    is.logical(scale), length(scale) == 1L, !is.na(scale),
    is.logical(signed), length(signed) == 1L, !is.na(signed)
  )

  pc_methods <- c("pcor", "glasso", "cor")
  method_given <- !missing(method)
  method <- match.arg(method)

  # ---- Resolve input: node-level network, raw data, estimator ----
  if (inherits(x, c("netobject", "cograph_network"))) {
    node_network <- x
    data <- x$data
    # Association netobjects store $data as a numeric matrix
    if (is.matrix(data)) data <- as.data.frame(data)
    net_method <- .resolve_method_alias(x$method)
    if (!method_given && net_method %in% pc_methods) {
      method <- net_method
    }
    if (isTRUE(x$directed)) {
      stop("'x' is a directed network (method = \"", x$method, "\"). ",
           "build_mcml_pc() is for undirected association networks; ",
           "use build_mcml() for transition networks.", call. = FALSE)
    }
  } else if (is.data.frame(x)) {
    if (!is.null(id_col)) {
      stopifnot(is.character(id_col))
      missing_ids <- setdiff(id_col, names(x))
      if (length(missing_ids) > 0) {
        stop("id_col column(s) not found in data: ",
             paste(missing_ids, collapse = ", "), call. = FALSE)
      }
      x <- x[, setdiff(names(x), id_col), drop = FALSE]
    }
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("Data frame input must be all-numeric (observations x ",
           "variables). Non-numeric columns: ",
           paste(names(x)[!numeric_cols], collapse = ", "), call. = FALSE)
    }
    data <- x
    node_network <- .pc_estimate(data, method, cor_method)
  } else {
    stop("'x' must be an undirected netobject or a numeric data.frame.",
         call. = FALSE)
  }
  W <- node_network$weights
  if (!is.matrix(W) || is.null(rownames(W))) {
    stop("Could not extract a labeled weight matrix from 'x'.",
         call. = FALSE)
  }
  node_labels <- rownames(W)

  has_data <- is.data.frame(data) && nrow(data) > 0L &&
    all(node_labels %in% colnames(data))
  needs_data <- c("composite", "loadings", "rv", "canonical")
  if (aggregation %in% needs_data && !has_data) {
    stop("aggregation = \"", aggregation, "\" requires raw data with one ",
         "column per node. Pass a netobject that carries $data, or the ",
         "raw data.frame itself. For weight-only input use ",
         "aggregation = \"average\".", call. = FALSE)
  }
  if (within == "reestimate" && !has_data) {
    within <- "subnetwork"
  }

  # ---- Normalize cluster membership ----
  members <- .pc_normalize_clusters(clusters, node_labels)
  k <- length(members)
  if (k < 2) {
    stop("At least 2 clusters are required.", call. = FALSE)
  }
  cluster_names <- names(members)

  # ---- Item diagnostics + composite weights under the chosen scheme ----
  loadings_df <- .pc_network_loadings(W, members, signed = signed)
  if (aggregation == "composite" && !identical(weighting, "strength")) {
    loadings_df <- .pc_apply_weighting(loadings_df, W, data, members,
                                       weighting, signed,
                                       fa_method = fa_method,
                                       cor_method = cor_method,
                                       fa_args = fa_args)
  }
  misfit_items <- loadings_df$node[loadings_df$misfit]
  if (length(misfit_items) > 0) {
    warning("Item(s) more strongly connected to another cluster than ",
            "their own (possible misassignment): ",
            paste(misfit_items, collapse = ", "),
            ". See $loadings (misfit, cross_cluster).", call. = FALSE)
  }
  flipped_items <- loadings_df$node[loadings_df$sign < 0]
  if (signed && length(flipped_items) > 0 && aggregation == "composite") {
    warning("Reverse-keyed item(s) flipped in composites: ",
            paste(flipped_items, collapse = ", "),
            ". See $loadings (sign).", call. = FALSE)
  }

  # ---- Macro network ----
  macro <- switch(aggregation,
    average = .wrap_netobject(.pc_average_macro(W, members), data = NULL,
                              method = "pc_average", directed = FALSE),
    composite = {
      composites <- .pc_composites(data, members, loadings_df, scale,
                                   signed = signed)
      .pc_estimate(composites, method, "pearson")
    },
    rv = .wrap_netobject(.pc_block_macro(data, members, "rv"),
                         data = NULL, method = "pc_rv", directed = FALSE),
    canonical = .wrap_netobject(.pc_block_macro(data, members, "canonical"),
                                data = NULL, method = "pc_canonical",
                                directed = FALSE)
  )

  # ---- Within-cluster networks ----
  within_nets <- lapply(cluster_names, function(cl) {
    nodes <- members[[cl]]
    if (length(nodes) < 2) return(NULL)
    if (within == "reestimate" && has_data) {
      .pc_estimate(data[, nodes, drop = FALSE], method, cor_method)
    } else {
      .wrap_netobject(W[nodes, nodes, drop = FALSE], data = NULL,
                      method = "pc_subnetwork", directed = FALSE)
    }
  })
  names(within_nets) <- cluster_names

  result <- list(
    macro = macro,
    clusters = within_nets,
    cluster_members = members,
    loadings = loadings_df,
    node_network = node_network,
    data = if (has_data) data else NULL,
    meta = list(
      aggregation = aggregation,
      method = if (aggregation == "composite") method
        else NA_character_,
      within = within,
      weighting = if (aggregation != "composite") NA_character_
        else if (is.function(weighting)) "custom (function)"
        else if (is.numeric(weighting)) "custom (vector)"
        else weighting,
      fa_method = if (identical(weighting, "factor")) fa_method
        else NA_character_,
      fa_args = if (identical(weighting, "factor") && length(fa_args) > 0)
        fa_args else NULL,
      scale = scale,
      cor_method = cor_method,
      signed = signed,
      n_nodes = length(node_labels),
      n_clusters = k,
      cluster_sizes = vapply(members, length, integer(1)),
      n_misfit = length(misfit_items),
      n_flipped = length(flipped_items),
      directed = FALSE,
      source = "pc",
      experimental = TRUE
    )
  )
  class(result) <- "mcml_pc"
  result
}


#' Estimate an association network, honoring cor_method incl. polychoric
#' @noRd
.pc_estimate <- function(data, estimator, cor_method) {
  if (identical(cor_method, "polychoric")) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      stop("cor_method = \"polychoric\" requires the 'lavaan' package.",
           call. = FALSE)
    }
    df <- as.data.frame(data)
    S <- lavaan::lavCor(df, ordered = names(df))
    S <- unclass(S)
    attributes(S)[setdiff(names(attributes(S)),
                          c("dim", "dimnames"))] <- NULL
    n_obs <- sum(stats::complete.cases(df))
    return(build_network(S, method = estimator, n = n_obs))
  }
  build_network(as.data.frame(data), method = estimator,
                cor_method = cor_method)
}


#' Normalize cluster input to a named list of member labels
#' @noRd
.pc_normalize_clusters <- function(clusters, node_labels) {
  if (is.list(clusters)) {
    if (is.null(names(clusters)) || any(!nzchar(names(clusters)))) {
      stop("'clusters' list must have non-empty names.", call. = FALSE)
    }
    members <- lapply(clusters, as.character)
  } else if (is.atomic(clusters) && !is.null(names(clusters))) {
    members <- split(names(clusters), as.character(clusters))
  } else {
    stop("'clusters' must be a named list of node labels or a vector of ",
         "cluster labels named by node.", call. = FALSE)
  }
  assigned <- unlist(members, use.names = FALSE)
  if (anyDuplicated(assigned)) {
    stop("Each node may belong to exactly one cluster. Duplicated: ",
         paste(unique(assigned[duplicated(assigned)]), collapse = ", "),
         call. = FALSE)
  }
  unknown <- setdiff(assigned, node_labels)
  if (length(unknown) > 0) {
    stop("Cluster members not present in the network: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  missing_nodes <- setdiff(node_labels, assigned)
  if (length(missing_nodes) > 0) {
    stop("Nodes without a cluster assignment: ",
         paste(missing_nodes, collapse = ", "), call. = FALSE)
  }
  members
}


#' Block-average macro matrix (descriptive aggregation)
#' @noRd
.pc_average_macro <- function(W, members) {
  k <- length(members)
  cluster_names <- names(members)
  macro <- matrix(0, k, k, dimnames = list(cluster_names, cluster_names))
  pairs <- expand.grid(a = seq_len(k), b = seq_len(k))
  vals <- vapply(seq_len(nrow(pairs)), function(r) {
    a <- members[[pairs$a[r]]]
    b <- members[[pairs$b[r]]]
    block <- W[a, b, drop = FALSE]
    if (pairs$a[r] == pairs$b[r]) {
      off <- block[row(block) != col(block)]
      if (length(off) == 0) 0 else mean(off)
    } else {
      mean(block)
    }
  }, numeric(1))
  macro[cbind(pairs$a, pairs$b)] <- vals
  (macro + t(macro)) / 2
}


#' Block-level macro via RV coefficient or first canonical correlation
#' @noRd
.pc_block_macro <- function(data, members, type) {
  k <- length(members)
  cluster_names <- names(members)
  X <- lapply(members, function(nodes) {
    as.matrix(data[, nodes, drop = FALSE])
  })
  keep <- stats::complete.cases(do.call(cbind, X))
  X <- lapply(X, function(m) m[keep, , drop = FALSE])

  macro <- diag(1, k)
  dimnames(macro) <- list(cluster_names, cluster_names)
  idx <- which(upper.tri(macro), arr.ind = TRUE)
  vals <- vapply(seq_len(nrow(idx)), function(r) {
    A <- X[[idx[r, 1]]]
    B <- X[[idx[r, 2]]]
    if (type == "rv") {
      .pc_rv_coefficient(A, B)
    } else {
      suppressWarnings(stats::cancor(A, B)$cor[1])
    }
  }, numeric(1))
  macro[idx] <- vals
  macro[idx[, c(2, 1), drop = FALSE]] <- vals
  macro
}


#' Escoufier's RV coefficient between two observation blocks
#' @noRd
.pc_rv_coefficient <- function(A, B) {
  A <- base::scale(A, center = TRUE, scale = FALSE)
  B <- base::scale(B, center = TRUE, scale = FALSE)
  sab <- crossprod(A, B)
  saa <- crossprod(A)
  sbb <- crossprod(B)
  num <- sum(sab * sab)
  den <- sqrt(sum(saa * saa) * sum(sbb * sbb))
  if (den == 0) return(NA_real_)
  num / den
}


#' Item diagnostics: signed own-cluster loadings, cross-loadings, misfit
#'
#' For every node, the signed sum of its node-level weights to the other
#' members of each cluster. Own-cluster loading drives the composite
#' weight (absolute share, sums to 1 per cluster); the strongest
#' off-cluster absolute loading drives the misfit flag.
#' @noRd
.pc_network_loadings <- function(W, members, signed = TRUE) {
  cluster_names <- names(members)
  rows <- lapply(cluster_names, function(cl) {
    nodes <- members[[cl]]
    singleton <- length(nodes) == 1L

    # Item signs from the leading eigenvector of the within-block weight
    # matrix, oriented so the majority of items are positive. This is
    # robust to a single reverse-keyed item: pairwise signed sums are
    # not (a reversed neighbor drags a good item's signed sum to ~0).
    item_sign <- rep(1, length(nodes))
    names(item_sign) <- nodes
    if (signed && !singleton) {
      block <- W[nodes, nodes, drop = FALSE]
      if (any(block[row(block) != col(block)] != 0)) {
        v <- eigen(block, symmetric = TRUE)$vectors[, 1]
        if (sum(v) < 0) v <- -v
        item_sign <- ifelse(v < 0, -1, 1)
        names(item_sign) <- nodes
      }
    }

    per_node <- lapply(nodes, function(i) {
      # mean absolute connection strength to every cluster -- means, not
      # sums, so the comparison is not biased toward larger clusters
      to_cluster <- vapply(cluster_names, function(target) {
        targets <- setdiff(members[[target]], i)
        if (length(targets) == 0) return(0)
        mean(abs(W[i, targets]))
      }, numeric(1))
      own <- to_cluster[[cl]]
      cross <- to_cluster[setdiff(cluster_names, cl)]
      if (length(cross) == 0) {
        max_cross <- 0
        cross_cluster <- NA_character_
      } else {
        max_cross <- max(cross)
        cross_cluster <- names(cross)[which.max(cross)]
      }
      data.frame(
        node = i, cluster = cl,
        loading = unname(item_sign[i]) * own,
        sign = unname(item_sign[i]),
        max_cross = max_cross, cross_cluster = cross_cluster,
        # A singleton has no within-cluster connections by construction;
        # misfit is undefined for it, not evidence of misassignment.
        misfit = if (singleton) FALSE else max_cross > own,
        stringsAsFactors = FALSE
      )
    })
    cl_df <- do.call(rbind, per_node)
    if (singleton) {
      cl_df$weight <- 1
      return(cl_df)
    }
    total <- sum(abs(cl_df$loading))
    if (total == 0) {
      warning("Cluster '", cl, "' has all-zero within-cluster loadings; ",
              "using equal composite weights.", call. = FALSE)
      cl_df$weight <- 1 / nrow(cl_df)
    } else {
      cl_df$weight <- abs(cl_df$loading) / total
    }
    cl_df
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("node", "cluster", "loading", "weight", "sign",
          "max_cross", "cross_cluster", "misfit")]
}


#' Per-observation cluster composites (NA-tolerant, sign-corrected)
#' @noRd
.pc_composites <- function(data, members, loadings_df, scale, signed) {
  cols <- lapply(names(members), function(cl) {
    nodes <- members[[cl]]
    mat <- as.matrix(data[, nodes, drop = FALSE])
    if (scale) mat <- base::scale(mat)
    idx <- match(nodes, loadings_df$node)
    w <- loadings_df$weight[idx]
    if (signed) w <- w * loadings_df$sign[idx]
    # Row-wise weighted mean over observed members: renormalize the
    # absolute weights across non-missing cells per row.
    wmat <- matrix(w, nrow(mat), length(w), byrow = TRUE)
    wmat[is.na(mat)] <- NA
    num <- rowSums(mat * wmat, na.rm = TRUE)
    den <- rowSums(abs(wmat), na.rm = TRUE)
    out <- num / den
    out[den == 0] <- NA_real_
    out
  })
  names(cols) <- names(members)
  as.data.frame(cols, stringsAsFactors = FALSE)
}


#' Composite-Weight Stability Under Case Resampling
#'
#' @description
#' \strong{Experimental.} Bootstraps the item weights of a
#' \code{\link{build_mcml_pc}} fit: rows of the raw data are resampled,
#' the node-level network is re-estimated each time, and the
#' connectivity-based composite weights are recomputed. Wide intervals
#' mean the weighting (and therefore the \code{"loadings"} macro
#' network) should not be over-interpreted.
#'
#' @param x An \code{mcml_pc} object that carries raw data.
#' @param iter Integer. Bootstrap replicates (default 200; node-level
#'   re-estimation makes this heavier than a plain bootstrap).
#' @param ci_level Numeric. Significance level for percentile CIs
#'   (default 0.05).
#' @param seed Integer or NULL. RNG seed.
#'
#' @return An object of class \code{"pc_loading_stability"}: a list with
#'   \code{summary} (tidy data frame: \code{node}, \code{cluster},
#'   \code{weight}, \code{boot_mean}, \code{boot_sd}, \code{ci_lower},
#'   \code{ci_upper}, \code{sign_flips} — the proportion of replicates
#'   in which the item's sign differed from the observed one),
#'   \code{boot_weights} (iter x n_nodes matrix), \code{iter}, and
#'   \code{ci_level}. Has print and plot methods.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df <- as.data.frame(matrix(rnorm(600), 100, 6))
#' names(df) <- c("a1", "a2", "a3", "b1", "b2", "b3")
#' cl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#' fit <- build_mcml_pc(df, cl, aggregation = "loadings",
#'                      estimator = "cor")
#' ls <- loading_stability(fit, iter = 50, seed = 1)
#' ls$summary
#' }
#'
#' @export
loading_stability <- function(x, iter = 200L, ci_level = 0.05,
                              seed = NULL) {
  if (!inherits(x, "mcml_pc")) {
    stop("'x' must be an mcml_pc object from build_mcml_pc().",
         call. = FALSE)
  }
  if (is.null(x$data)) {
    stop("This mcml_pc object carries no raw data; loading_stability() ",
         "needs the observations to resample.", call. = FALSE)
  }
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(ci_level), length(ci_level) == 1,
    ci_level > 0, ci_level < 1
  )
  iter <- as.integer(iter)
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  data <- x$data[, x$loadings$node, drop = FALSE]
  members <- x$cluster_members
  estimator <- .resolve_method_alias(x$node_network$method)
  if (!estimator %in% c("pcor", "glasso", "cor")) {
    estimator <- "cor"
  }
  cor_method <- x$meta$cor_method
  signed <- x$meta$signed
  node_order <- x$loadings$node
  n_rows <- nrow(data)

  boot <- vapply(seq_len(iter), function(b) {
    rows <- sample.int(n_rows, n_rows, replace = TRUE)
    net_b <- .pc_estimate(data[rows, , drop = FALSE], estimator,
                          cor_method)
    ld_b <- suppressWarnings(
      .pc_network_loadings(net_b$weights, members, signed = signed)
    )
    idx <- match(node_order, ld_b$node)
    ld_b$weight[idx] * ld_b$sign[idx]
  }, numeric(length(node_order)))
  boot <- t(boot)
  colnames(boot) <- node_order

  observed <- x$loadings$weight * x$loadings$sign
  probs <- c(ci_level / 2, 1 - ci_level / 2)
  ci <- apply(boot, 2, quantile, probs = probs, na.rm = TRUE)
  sign_flips <- colMeans(sweep(sign(boot), 2, sign(observed), "!="),
                         na.rm = TRUE)

  summary_df <- data.frame(
    node = node_order,
    cluster = x$loadings$cluster,
    weight = observed,
    boot_mean = colMeans(boot, na.rm = TRUE),
    boot_sd = apply(boot, 2, sd, na.rm = TRUE),
    ci_lower = ci[1, ],
    ci_upper = ci[2, ],
    sign_flips = unname(sign_flips),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  result <- list(summary = summary_df, boot_weights = boot,
                 iter = iter, ci_level = ci_level)
  class(result) <- "pc_loading_stability"
  result
}


#' Print Composite-Weight Stability
#'
#' @param x A \code{pc_loading_stability} object.
#' @param digits Number of digits to display (default 3).
#' @param ... Additional arguments (ignored).
#' @return \code{x}, invisibly.
#' @export
print.pc_loading_stability <- function(x, digits = 3, ...) {
  cat("Composite-Weight Stability (case bootstrap, experimental)\n")
  cat(sprintf("  %d replicates | %.0f%% percentile CIs\n\n",
              x$iter, 100 * (1 - x$ci_level)))
  df <- x$summary
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  print(df, row.names = FALSE)
  invisible(x)
}


#' Plot Composite-Weight Stability
#'
#' Signed composite weights with bootstrap percentile intervals, faceted
#' by cluster.
#'
#' @param x A \code{pc_loading_stability} object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.pc_loading_stability <- function(x, ...) {
  df <- x$summary
  ggplot2::ggplot(df, ggplot2::aes(x = .data$weight,
                                   y = stats::reorder(.data$node,
                                                      .data$weight))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_lower,
                                        xmax = .data$ci_upper),
                           width = 0.3, color = "#4A6FE3",
                           linewidth = 0.5) +
    ggplot2::geom_point(size = 2, color = "#D33F6A") +
    ggplot2::facet_wrap(~cluster, scales = "free_y") +
    ggplot2::labs(
      x = "Signed composite weight", y = NULL,
      title = "Composite-weight stability under case resampling",
      subtitle = sprintf("%d replicates | %.0f%% percentile CIs",
                         x$iter, 100 * (1 - x$ci_level))
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Print an MCML-PC Object
#'
#' @param x An \code{mcml_pc} object.
#' @param digits Number of digits to display (default 3).
#' @param ... Additional arguments (ignored).
#' @return \code{x}, invisibly.
#' @export
print.mcml_pc <- function(x, digits = 3, ...) {
  m <- x$meta
  cat("MCML for Psychometric Networks (experimental)\n")
  est <- if (is.na(m$method)) "none (no re-estimation)" else m$method
  wt <- if (is.na(m$weighting)) "" else paste0(" | Weights: ", m$weighting)
  cat(sprintf("  Aggregation: %s | Estimator: %s | Within: %s | cor: %s%s\n",
              m$aggregation, est, m$within, m$cor_method, wt))
  cat(sprintf("  %d nodes in %d clusters (%s)\n",
              m$n_nodes, m$n_clusters,
              paste(names(m$cluster_sizes), m$cluster_sizes,
                    sep = ": ", collapse = ", ")))
  if (m$n_misfit > 0) {
    cat(sprintf("  ! %d misfit item(s) - see $loadings\n", m$n_misfit))
  }
  if (m$n_flipped > 0) {
    cat(sprintf("  ! %d reverse-keyed item(s) flipped - see $loadings\n",
                m$n_flipped))
  }
  cat("\nMacro weights:\n")
  print(round(x$macro$weights, digits))
  invisible(x)
}


#' Summarize an MCML-PC Object
#'
#' @param object An \code{mcml_pc} object.
#' @param ... Additional arguments (ignored).
#' @return Tidy data frame with one row per macro edge (upper triangle),
#'   columns \code{from}, \code{to}, \code{weight}.
#' @export
summary.mcml_pc <- function(object, ...) {
  W <- object$macro$weights
  idx <- which(upper.tri(W), arr.ind = TRUE)
  data.frame(
    from = rownames(W)[idx[, 1]],
    to = colnames(W)[idx[, 2]],
    weight = W[idx],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Plot an MCML-PC Object
#'
#' Heatmap of the macro (cluster-level) weights with the package's
#' diverging palette. For the two-layer network rendering use
#' \code{cograph::plot_mcml()}, which accepts \code{mcml_pc} objects and
#' draws them undirected.
#'
#' @param x An \code{mcml_pc} object.
#' @param digits Number of digits for tile labels (default 2).
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.mcml_pc <- function(x, digits = 2, ...) {
  W <- x$macro$weights
  df <- data.frame(
    from = rep(rownames(W), times = ncol(W)),
    to = rep(colnames(W), each = nrow(W)),
    weight = as.vector(W),
    stringsAsFactors = FALSE
  )
  lim <- max(abs(df$weight))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$to, y = .data$from,
                                   fill = .data$weight)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = round(.data$weight, digits)),
                       size = 3.5) +
    ggplot2::scale_fill_gradient2(low = "#D33F6A", mid = "white",
                                  high = "#4A6FE3", midpoint = 0,
                                  limits = c(-lim, lim)) +
    ggplot2::scale_y_discrete(limits = rev(rownames(W))) +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "Weight",
      title = "Cluster-level (macro) network weights",
      subtitle = sprintf("aggregation = \"%s\" (experimental)",
                         x$meta$aggregation)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Promote a psychometric MCML result to a network group
#'
#' \code{as_networks()} is the psychometric-network counterpart of
#' \code{\link{as_tna}}. It promotes the cluster-level (macro) and
#' within-cluster networks produced by \code{\link{build_mcml_pc}} into a
#' single \code{netobject_group}, so the result flows into the same
#' downstream verbs as any other group of networks (\code{print()},
#' \code{summary()}, \code{plot()}, \code{\link{net_centrality}}).
#'
#' Where \code{as_tna()} promotes \emph{transition} networks (directed,
#' row-normalised, with initial probabilities) and re-wraps raw matrices,
#' \code{as_networks()} promotes \emph{psychometric} networks (undirected;
#' correlation / partial-correlation / glasso). The macro and within-cluster
#' components of an \code{mcml_pc} object are already full netobjects carrying
#' their estimator, directedness and data, so this function assembles them
#' into a group rather than re-wrapping matrices.
#'
#' @param x An object to convert. The \code{mcml_pc} method (from
#'   \code{\link{build_mcml_pc}}) is the primary path.
#'
#' @return A \code{netobject_group}: a named list whose first element is
#'   \code{macro} (the cluster-level network), followed by one netobject per
#'   non-singleton cluster.
#'
#' @seealso \code{\link{build_mcml_pc}} to create the input,
#'   \code{\link{as_tna}} for the transition-network counterpart.
#'
#' @examples
#' set.seed(1)
#' df <- as.data.frame(matrix(stats::rnorm(200 * 6), 200, 6))
#' names(df) <- c("a1", "a2", "a3", "b1", "b2", "b3")
#' clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#' fit <- build_mcml_pc(df, clusters, aggregation = "composite", method = "cor")
#' nets <- as_networks(fit)
#' nets
#' nets$macro$weights
#' @export
as_networks <- function(x) {
  UseMethod("as_networks")
}

#' @rdname as_networks
#' @return The \code{mcml_pc} method returns a \code{netobject_group};
#'   singleton clusters (no within-network) are dropped with a
#'   \code{warning()}.
#' @export
as_networks.mcml_pc <- function(x) {
  clusters <- x$clusters
  is_singleton <- vapply(clusters, is.null, logical(1))
  if (any(is_singleton)) {
    warning("Dropped singleton clusters with no within-network: ",
            paste(names(clusters)[is_singleton], collapse = ", "),
            call. = FALSE)
  }
  result <- c(list(macro = x$macro), clusters[!is_singleton])
  class(result) <- "netobject_group"
  result
}

#' @rdname as_networks
#' @return The default method returns the input unchanged if it is already a
#'   \code{netobject_group}, otherwise it errors.
#' @export
as_networks.default <- function(x) {
  if (inherits(x, "netobject_group")) {
    return(x)
  }
  stop("Cannot convert object of class '", class(x)[1],
       "' to a network group.", call. = FALSE)
}


#' Override loading weights/signs with the chosen weighting scheme
#'
#' "strength" weights come from .pc_network_loadings directly; this
#' handles "equal", "eigen", "pca", "factor". Signs for "equal" keep the
#' eigenvector-based signs already in loadings_df; vector-based schemes
#' carry their own signs.
#' @noRd
.pc_apply_weighting <- function(loadings_df, W, data, members, weighting,
                                signed, fa_method = "ml",
                                cor_method = "pearson",
                                fa_args = list()) {
  # Magnitude-only schemes have non-negative weight vectors; they keep the
  # eigenvector-based item signs already in loadings_df. Sign-carrying
  # schemes derive both weight and sign from their own vector.
  magnitude_only <- c("equal", "closeness", "betweenness", "specificity")

  per_cluster <- lapply(names(members), function(cl) {
    nodes <- members[[cl]]
    m <- length(nodes)
    idx <- match(nodes, loadings_df$node)
    out <- loadings_df[idx, , drop = FALSE]
    if (m == 1L) {
      out$weight <- 1
      return(out)
    }
    block <- W[nodes, nodes, drop = FALSE]

    fallback_equal <- function(reason) {
      warning("Cluster '", cl, "': ", reason,
              "; using equal weights.", call. = FALSE)
      rep(1, m)
    }

    v <- if (is.function(weighting)) {
      val <- weighting(block, data[, nodes, drop = FALSE], nodes)
      if (!is.numeric(val) || length(val) != m || any(!is.finite(val))) {
        stop("Custom weighting function must return ", m, " finite ",
             "numeric values for cluster '", cl, "'.", call. = FALSE)
      }
      val
    } else if (is.numeric(weighting)) {
      missing_w <- setdiff(nodes, names(weighting))
      if (length(missing_w) > 0) {
        stop("Custom weighting vector lacks entries for: ",
             paste(missing_w, collapse = ", "), call. = FALSE)
      }
      unname(weighting[nodes])
    } else {
      switch(weighting,
        equal = rep(1, m),
        eigen = eigen(block, symmetric = TRUE)$vectors[, 1],
        closeness = {
          cv <- .closeness(abs(block), directed = FALSE)$Closeness
          if (all(cv == 0)) fallback_equal("all-zero closeness") else cv
        },
        betweenness = {
          bv <- .betweenness(abs(block), directed = FALSE)
          if (all(bv == 0)) {
            fallback_equal(
              "all-zero betweenness (no item lies between others)")
          } else {
            bv
          }
        },
        expected_influence = {
          # mean SIGNED connection to the other own-cluster members
          # (Robinaugh et al. 2016); negative EI marks reverse-keyed items
          vapply(seq_len(m), function(i) {
            mean(block[i, -i])
          }, numeric(1))
        },
        specificity = {
          # the misfit margin: own-cluster strength minus the strongest
          # cross-cluster strength; items that belong elsewhere get ~0
          margin <- pmax(abs(out$loading) - out$max_cross, 0)
          if (all(margin == 0)) {
            fallback_equal("no item is cluster-specific (all margins <= 0)")
          } else {
            margin
          }
        },
        pca = .pc_first_pc(data[, nodes, drop = FALSE], cl),
        factor = {
          if (m < 3) {
            warning("Cluster '", cl, "' has fewer than 3 items; ",
                    "weighting = \"factor\" falls back to \"pca\".",
                    call. = FALSE)
            .pc_first_pc(data[, nodes, drop = FALSE], cl)
          } else {
            tryCatch(
              .pc_fa_loadings(data[, nodes, drop = FALSE], fa_method,
                              cor_method, fa_args),
              error = function(e) {
                warning("Factor extraction (", fa_method, ") failed for ",
                        "cluster '", cl, "' (", conditionMessage(e),
                        "); falling back to \"pca\".", call. = FALSE)
                .pc_first_pc(data[, nodes, drop = FALSE], cl)
              }
            )
          }
        },
        item_total = {
          # corrected item-total correlation: item vs the mean of the
          # OTHER (standardized) members. Columns are pre-oriented by the
          # eigenvector-based signs already in loadings_df so a reversed
          # member cannot contaminate the rest-mean (in small clusters
          # that collapses the good items' item-totals to ~0); the eigen
          # signs are then reapplied so reversed items keep sign -1.
          z <- base::scale(as.matrix(data[, nodes, drop = FALSE]))
          z <- sweep(z, 2, out$sign, "*")
          r_corrected <- vapply(seq_len(m), function(i) {
            rest <- rowMeans(z[, -i, drop = FALSE], na.rm = TRUE)
            r <- suppressWarnings(
              stats::cor(z[, i], rest, use = "pairwise.complete.obs"))
            if (is.na(r)) 0 else r
          }, numeric(1))
          r_corrected * out$sign
        }
      )
    }

    if (sum(v) < 0) v <- -v
    total <- sum(abs(v))
    out$weight <- if (total == 0) rep(1 / m, m) else abs(v) / total
    sign_carrying <- is.function(weighting) || is.numeric(weighting) ||
      !(weighting %in% magnitude_only)
    if (signed && sign_carrying) {
      out$sign <- ifelse(v < 0, -1, 1)
    }
    out
  })
  out <- do.call(rbind, per_cluster)
  rownames(out) <- NULL
  out
}


#' First principal component of a data block's correlation matrix
#' @noRd
.pc_first_pc <- function(block, cl) {
  R <- stats::cor(as.matrix(block), use = "pairwise.complete.obs")
  if (any(!is.finite(R))) {
    warning("Cluster '", cl, "' correlation matrix has non-finite ",
            "entries; using equal PCA weights.", call. = FALSE)
    return(rep(1, ncol(R)))
  }
  eigen(R, symmetric = TRUE)$vectors[, 1]
}


#' Correlation matrix of a data block per cor_method
#' @noRd
.pc_cor_matrix <- function(block, cor_method) {
  df <- as.data.frame(block)
  if (identical(cor_method, "polychoric")) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      stop("cor_method = \"polychoric\" requires the 'lavaan' package.",
           call. = FALSE)
    }
    S <- lavaan::lavCor(df, ordered = names(df))
    S <- unclass(S)
    attributes(S)[setdiff(names(attributes(S)),
                          c("dim", "dimnames"))] <- NULL
    list(R = S, n = sum(stats::complete.cases(df)))
  } else {
    list(R = stats::cor(as.matrix(df), use = "pairwise.complete.obs",
                        method = cor_method),
         n = sum(stats::complete.cases(df)))
  }
}


#' One-factor loadings of a data block under the chosen extraction method
#'
#' "ml" = maximum likelihood (stats::factanal on the correlation matrix);
#' "paf" = iterated principal axis factoring; "minres" = minimum residual
#' (ULS) via L-BFGS-B on the uniquenesses; "cfa" = one-factor CFA in
#' lavaan (DWLS on declared-ordered items when cor_method =
#' "polychoric", ML otherwise), standardized loadings. All operate on
#' the cor_method-consistent correlation structure.
#' @noRd
.pc_fa_loadings <- function(block, fa_method, cor_method,
                            fa_args = list()) {
  # Arguments that would break the one-weight-per-item contract or the
  # internal model construction are protected per backend.
  drop_protected <- function(args, protected) {
    bad <- intersect(names(args), protected)
    if (length(bad) > 0) args[bad] <- NULL
    args
  }

  if (fa_method == "cfa") {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      stop("fa_method = \"cfa\" requires the 'lavaan' package.",
           call. = FALSE)
    }
    df <- as.data.frame(block)
    # syntactic-safe internal names for the lavaan model string
    orig <- names(df)
    names(df) <- paste0(".it", seq_along(orig))
    model <- paste("F =~", paste(names(df), collapse = " + "))
    args <- list(
      model = model, data = df, std.lv = TRUE,
      ordered = if (identical(cor_method, "polychoric")) names(df) else NULL
    )
    extra <- drop_protected(fa_args, c("model", "data"))
    args <- utils::modifyList(args, extra)
    fit <- do.call(lavaan::cfa, args)
    std <- lavaan::standardizedSolution(fit)
    lam <- std$est.std[std$op == "=~"]
    return(lam[match(names(df), std$rhs[std$op == "=~"])])
  }

  cm <- .pc_cor_matrix(block, cor_method)
  R <- cm$R
  switch(fa_method,
    ml = {
      args <- list(covmat = R, n.obs = cm$n, factors = 1)
      extra <- drop_protected(fa_args,
                              c("covmat", "n.obs", "factors", "x", "data"))
      args <- utils::modifyList(args, extra)
      as.numeric(do.call(stats::factanal, args)$loadings)
    },
    paf = do.call(.pc_fa_paf,
                  c(list(R = R),
                    drop_protected(fa_args, "R")[
                      intersect(names(fa_args), c("max_iter", "tol"))])),
    minres = .pc_fa_minres(R)
  )
}


#' Iterated principal axis factoring (one factor)
#' @noRd
.pc_fa_paf <- function(R, max_iter = 100L, tol = 1e-6) {
  h2 <- 1 - 1 / diag(solve(R))  # squared multiple correlations start
  for (it in seq_len(max_iter)) {
    Rr <- R
    diag(Rr) <- h2
    e <- eigen(Rr, symmetric = TRUE)
    lambda <- e$vectors[, 1] * sqrt(max(e$values[1], 0))
    h2_new <- pmin(lambda^2, 0.998)
    if (max(abs(h2_new - h2)) < tol) break
    h2 <- h2_new
  }
  lambda
}


#' Minimum-residual (ULS) one-factor extraction
#' @noRd
.pc_fa_minres <- function(R) {
  m <- ncol(R)
  lambda_for <- function(psi) {
    Rr <- R
    diag(Rr) <- 1 - psi
    e <- eigen(Rr, symmetric = TRUE)
    e$vectors[, 1] * sqrt(max(e$values[1], 0))
  }
  obj <- function(psi) {
    resid <- R - tcrossprod(lambda_for(psi))
    diag(resid) <- 0
    sum(resid^2)
  }
  start <- pmin(pmax(1 / diag(solve(R)), 0.005), 0.995)
  fit <- stats::optim(start, obj, method = "L-BFGS-B",
                      lower = rep(0.005, m), upper = rep(0.995, m))
  lambda_for(fit$par)
}
