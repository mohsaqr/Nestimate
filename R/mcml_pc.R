# ---- MCML for Partial-Correlation / Psychometric Networks (EXPERIMENTAL) ----

#' Multi-Cluster Multi-Level Aggregation for Psychometric Networks
#'
#' @description
#' `build_mcml_pc()` collapses a psychometric network whose nodes are assigned
#' to clusters into a cluster-level (macro) network and one within-cluster
#' network per cluster, and returns an `mcml_pc` object. It is the counterpart
#' of [build_mcml()] for undirected association networks: it reads a
#' `netobject` estimated with `cor`, `pcor` or `glasso`, or a numeric data
#' frame of observations, from which it estimates the node-level network
#' first. \strong{Experimental}: the API and the aggregation formulas may
#' change between releases.
#'
#' Partial correlations do not aggregate by arithmetic. The submatrix of a
#' partial-correlation matrix is \emph{not} the partial-correlation network of
#' that subsystem, because the conditioning set changes, and averaging entries
#' across blocks is descriptive only. The aggregations (see `aggregation`)
#' therefore differ in status:
#'
#' \describe{
#'   \item{Score paths, re-estimated (`"scaled"`, `"composite"`, `"mean"`,
#'     `"median"`, `"loadings"`)}{Each observation gets one score per
#'     cluster, formed from the cluster's member items (see `weighting`), and
#'     the estimator in `method` is re-fitted on these scores, so the macro
#'     network is a genuine correlation, partial-correlation or EBICglasso
#'     network among clusters. Needs raw data.}
#'   \item{`"average"` (descriptive)}{The macro edge between clusters A and B
#'     is the mean of the signed node-level weights between members of A and
#'     members of B; the macro diagonal is the mean within-block
#'     off-diagonal weight. Needs only the weight matrix, so it works on a
#'     netobject without data. Signed averaging can cancel opposite-sign
#'     edges: read it as a net average association, not a connection
#'     strength.}
#'   \item{`"escoufier"` (descriptive, multivariate)}{The macro edge is
#'     Escoufier's RV coefficient between the member blocks, a matrix
#'     correlation in `[0, 1]` computed from the block covariance structure.
#'     No scores are formed and nothing is re-estimated, so nothing is lost to
#'     averaging; signs are not represented. Needs raw data.}
#'   \item{`"cancor"` (descriptive, multivariate)}{The macro edge is the first
#'     canonical correlation between the member blocks: the strongest linear
#'     relationship any weighting of A's items can have with any weighting of
#'     B's items, an upper bound on what the score paths can recover. Needs
#'     raw data.}
#' }
#'
#' Every aggregation inherits the cluster assignment. Before reading the macro
#' network, check the item table with [item_loadings()]: an item more strongly
#' connected to another cluster than to its own is flagged as misfit and
#' named in a warning. A flagged item means every aggregation carries a
#' doubtful membership, so fix the assignment rather than ignore the flag.
#' How far the composite weights move under case resampling is checked with
#' [loading_stability()]. The macro network of a score path carries its
#' scores as data, so [bootstrap_network()] (edge-weight intervals) and,
#' with at least three clusters, [vertex_bootstrap()] (network-level
#' intervals) work on it directly, as in
#' `bootstrap_network(macro_network(fit))`, and [vertex_compare()] (also with
#' at least three clusters) compares the macro networks of two fits;
#' [macro_network()] extracts it.
#'
#' @details
#' \strong{Item diagnostics.} Every item's mean absolute connection to the
#' other members of \emph{every} cluster is computed from the node-level
#' network. [item_loadings()] reports, per item, its own-cluster loading
#' (signed by the eigenvector orientation; a sign-carrying weighting is
#' recorded in the `sign` column instead), its composite weight, its strongest cross-cluster loading and the
#' cluster that loading points to, and a `misfit` flag set when the cross-cluster
#' loading exceeds the own-cluster loading. A singleton cluster has no
#' within-cluster connections, so its item is never flagged.
#'
#' \strong{Reverse-keyed items.} With `signed = TRUE` (default), the item
#' signs come from the leading eigenvector of the within-cluster block of the
#' node-level network, oriented so most items are positive. An item with a
#' negative sign enters its cluster score with sign -1, so a reverse-keyed
#' item reinforces its cluster score instead of cancelling it. Flips are
#' reported in the `sign` column of [item_loadings()] and, on the score
#' paths, in a warning. On the score paths a sign-carrying weighting sets the
#' signs from its own values instead (see `weighting`).
#'
#' \strong{Missing data.} Cluster scores are row-wise weighted means over the
#' \emph{observed} members, with the weights renormalised per row; a row with
#' no observed member is `NA` and the estimator drops it with a message.
#' Node-level estimation applies the estimators' own complete-case handling,
#' and `"escoufier"` and `"cancor"` use the rows complete on every item.
#'
#' \strong{Ordinal items.} `cor_method = "polychoric"` (needs the
#' \pkg{lavaan} package) estimates the node-level and within-cluster networks
#' from polychoric correlations, which suits Likert items. Cluster scores are
#' continuous, so the macro network is always re-estimated from their Pearson
#' correlations.
#'
#' \strong{Within-cluster networks} follow `within`. `"reestimate"` re-fits
#' the estimator on the member columns alone, giving the conditional
#' structure of the subsystem; `"subnetwork"` slices the node-level weight
#' matrix and is descriptive (for `pcor` and `glasso` it keeps the
#' conditioning on items outside the cluster). Input without raw data always
#' uses `"subnetwork"`. A singleton cluster gets no within-cluster network
#' (`NULL`).
#'
#' All constituent networks are undirected (`meta$directed = FALSE`), so
#' renderers that detect directedness, such as `cograph::plot_mcml()`, draw
#' them without arrowheads.
#'
#' @param x The network or the data. A `netobject` (or a `cograph_network`)
#'   estimated with an undirected association method (`cor`, `pcor`, `glasso`, or their
#'   aliases), whose weights, stored data and method are reused; a directed
#'   network raises an error and points to [build_mcml()]. Or a numeric data
#'   frame of observations (one row per observation, one column per item),
#'   from which the node-level network is estimated with `method` and
#'   `cor_method`; a non-numeric column raises an error. Any other input, a
#'   matrix included, raises an error. For the `print()` and `plot()`
#'   methods: an `mcml_pc` object.
#' @param clusters The cluster membership, in one of three forms: a named
#'   list of node-label vectors (the names are the cluster labels), a data
#'   frame whose first column holds the node names and second column the
#'   cluster labels (read by position, column names are not used), or a
#'   vector of cluster labels named by node. Every node must belong to
#'   exactly one cluster, every member must be a node, and at least two
#'   clusters are needed; anything else raises an error.
#' @param aggregation How the clusters are collapsed into the macro network.
#'   Default `"scaled"`. The score paths re-estimate the network on one score
#'   per cluster and observation:
#'   \describe{
#'     \item{`"scaled"`}{The weighted mean of the \emph{standardised} member
#'       items.}
#'     \item{`"composite"`}{The former name of `"scaled"`, kept as an
#'       alias.}
#'     \item{`"mean"`}{As `"scaled"`, on the raw (unstandardised) items.}
#'     \item{`"median"`}{The row-wise median of the raw, sign-corrected
#'       items. The median is unweighted, so the weight magnitudes from
#'       `weighting` are not used, but a sign-carrying weighting still sets
#'       the item signs.}
#'     \item{`"loadings"`}{As `"scaled"`, with the items weighted by their
#'       within-cluster network strength (`weighting = "strength"`) unless
#'       `weighting` is given.}
#'   }
#'   The descriptive ones relate the blocks without forming scores:
#'   `"average"` (the mean item-pair edge of each between-cluster block),
#'   `"escoufier"` (Escoufier's RV coefficient) and `"cancor"` (the first
#'   canonical correlation); the Description says how each can be read.
#'   Every path except `"average"` needs raw data; without it, an error
#'   suggests `"average"`.
#' @param method The network estimator for the node-level network of a data
#'   frame, the macro network of the score paths and the re-estimated
#'   within-cluster networks: `"pcor"` (default), `"glasso"` or `"cor"`, as
#'   in [build_network()]. When `x` is a `cor`, `pcor` or `glasso` netobject
#'   and `method` is not given, the netobject's own method is used; any other
#'   undirected netobject (an `ising` one, say) uses the default `"pcor"`.
#' @param within How the within-cluster networks are formed: `"reestimate"`
#'   (default) re-fits the estimator on each cluster's items, `"subnetwork"`
#'   slices the node-level weight matrix (see Details). Input without raw data
#'   uses `"subnetwork"` whatever is set.
#' @param weighting How the items are weighted inside their cluster score on
#'   the score paths; on `"median"` only the signs of a sign-carrying scheme
#'   are used. Default `"equal"` (`"strength"` under
#'   `aggregation = "loadings"`); the descriptive aggregations ignore it, and
#'   their item table keeps the `"strength"` weights. The weights are
#'   absolute values normalised to sum to 1 within each cluster.
#'   \describe{
#'     \item{`"equal"`}{1/m per item: the scale as scored.}
#'     \item{`"strength"`}{The item's mean absolute connection to the other
#'       members of its cluster in the node-level network (what
#'       `aggregation = "loadings"` selects), so items that anchor their
#'       cluster count more in its score. This is Nestimate's own weighting,
#'       related in spirit to network loadings (Christensen & Golino 2021)
#'       but not a reimplementation of any EGA-family estimator. A cluster
#'       whose within-cluster connections are all zero gets equal weights,
#'       with a warning.}
#'     \item{`"eigen"`}{The leading eigenvector of the cluster's block of the
#'       node-level network: like strength, with extra weight for items
#'       connected to other well-connected items.}
#'     \item{`"pca"`}{The first principal component of the member items'
#'       correlation matrix: a weighting from the data, blind to the
#'       estimated network. A correlation block with non-finite entries gives
#'       equal weights, with a warning.}
#'     \item{`"factor"`}{The standardised loadings of a one-factor model per
#'       cluster, the classical latent-variable weighting. `fa_method` sets
#'       the extraction and the correlations follow `cor_method`, so
#'       polychoric factor analysis of ordinal items is one call. A cluster
#'       with fewer than 3 items, or whose factor fit fails, falls back to
#'       `"pca"` with a warning.}
#'     \item{`"closeness"`, `"betweenness"`}{The item's closeness or
#'       betweenness centrality within its cluster's block of the node-level
#'       network (absolute weights). Betweenness can be zero for every item of
#'       a densely connected cluster; equal weights are then used, with a
#'       warning, and the same holds for all-zero closeness.}
#'     \item{`"expected_influence"`}{The item's mean \emph{signed} connection
#'       to the other members of its cluster (Robinaugh et al. 2016): like
#'       strength, but opposite-sign connections subtract, and a negative
#'       value marks a reverse-keyed item.}
#'     \item{`"specificity"`}{The misfit margin: own-cluster strength minus
#'       the strongest cross-cluster strength, floored at 0, so an item that
#'       belongs as much to another cluster adds nothing to its score; the
#'       weighting twin of the `misfit` flag. When no item of a cluster has a
#'       positive margin, equal weights are used, with a warning.}
#'     \item{`"item_total"`}{The corrected item-total correlation of each
#'       item with the mean of the other (standardised, sign-corrected)
#'       members: the classical scale-construction weighting.}
#'   }
#'   `weighting` also takes a \emph{named numeric vector} with one entry per
#'   item (its absolute values are normalised within each cluster and, when
#'   `signed = TRUE`, its signs set the item signs after the orientation
#'   described below), or a \emph{function}
#'   `function(W_block, data_block, nodes)` that returns one finite weight per
#'   item, evaluated per cluster; a vector missing an item, or a function
#'   returning anything else, raises an error.
#'
#'   Network-based and data-based weightings answer different questions, and
#'   comparing them is informative: when they diverge, the network's view of
#'   the cluster differs from its latent-variable view. The schemes whose
#'   weights cannot be negative (`"equal"`, `"strength"`, `"closeness"`,
#'   `"betweenness"`, `"specificity"`) keep the eigenvector-based item
#'   signs; the sign-carrying ones (`"eigen"`, `"pca"`, `"factor"`,
#'   `"expected_influence"`, `"item_total"`, a custom vector or function)
#'   set the signs from their own values. Within each cluster those values
#'   are oriented first, negated when they sum to a negative number, so the
#'   sign carrying the larger total counts as positive and only the items of
#'   the other sign are flipped: marking the reverse-keyed items of a
#'   cluster negative works only when they do not outweigh the rest in
#'   total.
#' @param cor_method The correlation used when estimating from raw data:
#'   `"pearson"` (default), `"spearman"` or `"polychoric"` (needs
#'   \pkg{lavaan}; see Details).
#' @param signed Whether reverse-keyed items are flipped in the cluster
#'   scores (see Details). Default `TRUE`.
#' @param fa_method The extraction method for `weighting = "factor"`; ignored
#'   otherwise.
#'   \describe{
#'     \item{`"ml"`}{(default) Maximum likelihood, `stats::factanal()` on the
#'       correlation matrix that `cor_method` gives.}
#'     \item{`"paf"`}{Iterated principal axis factoring, starting from the
#'       squared multiple correlations and iterating the communalities on the
#'       reduced correlation matrix.}
#'     \item{`"minres"`}{Minimum residual (unweighted least squares): the
#'       uniquenesses are chosen to minimise the squared off-diagonal
#'       residuals.}
#'     \item{`"cfa"`}{A one-factor confirmatory model fitted with
#'       `lavaan::cfa()` (standardised loadings); with
#'       `cor_method = "polychoric"` the items are declared ordered, which
#'       gives the categorical (DWLS) factor model.}
#'   }
#'   On well-behaved unidimensional clusters the four agree closely; a
#'   divergence points to a cluster prone to Heywood cases or not
#'   unidimensional.
#' @param ... Further arguments passed to the `fa_method` backend, spelled as
#'   that backend spells them, and allowed only with `weighting = "factor"`
#'   (otherwise they raise an error, as does an unnamed one): `lavaan::cfa()`
#'   arguments for `"cfa"` (`estimator = "WLSMV"`, `missing = "fiml"`,
#'   `se = "robust"`, ...), `stats::factanal()` arguments for `"ml"`, and
#'   `max_iter` / `tol` for `"paf"`; `"minres"` uses none. A name the
#'   backend does not read is dropped with a warning of class
#'   `nestimate_ignored_argument`. A `"paf"` or `"minres"` extraction that
#'   stops before converging keeps its last loadings and raises a warning of
#'   class `nestimate_no_converge`. `estimator` is lavaan's and raises an error for
#'   any other `fa_method`, since the network estimator is set by `method`.
#'   Arguments the backend call sets itself (`model` and `data` for `"cfa"`;
#'   `covmat`, `n.obs`, `factors`, `x` and `data` for `"ml"`; `R` for
#'   `"paf"` and `"minres"`) are dropped with a warning: the model is always
#'   the one-factor model per cluster, because the score needs exactly one
#'   weight per item. For the `print()`, `summary()` and `plot()` methods:
#'   ignored.
#' @param id_col Identifier columns to drop from data-frame input before the
#'   analysis. Default `NULL`. Use it for the id columns that
#'   [convert_sequence_format()] adds with `format = "frequency"`, whose
#'   per-sequence counts otherwise feed this function directly as behaviour
#'   profiles. A `rid` column is left out of the node-level network without
#'   being listed, as in the association estimators of [build_network()],
#'   except with `cor_method = "polychoric"`, where it must be listed too;
#'   any other id column, including the `.rid` a wide table gets, is read
#'   as an item unless listed. A column not in the data raises an error; a
#'   netobject ignores the argument.
#'
#' @return An `mcml_pc` object. `print()` shows the aggregation, estimator,
#'   cluster sizes, the misfit and flipped-item counts (when there are any)
#'   and the macro weights;
#'   `summary()` returns one row per macro edge; `plot()` draws the macro
#'   weights as a heatmap. [item_loadings()] gives the item table,
#'   [composites()] the cluster scores of a score path, [macro_network()] the
#'   macro network, [as_networks()] the macro and within-cluster networks as
#'   a `netobject_group`, [loading_stability()] the stability of the
#'   composite weights, and `cograph::plot_mcml()` draws the two-layer
#'   network. The object's fields:
#' \describe{
#'   \item{macro}{The cluster-level netobject (k x k, undirected).}
#'   \item{clusters}{Named list of within-cluster netobjects (`NULL` for a
#'     singleton cluster).}
#'   \item{cluster_members}{Named list of each cluster's node labels.}
#'   \item{loadings}{The item table that [item_loadings()] returns, one row
#'     per node.}
#'   \item{node_network}{The node-level netobject the aggregation was based
#'     on, kept for the diagnostics and [loading_stability()].}
#'   \item{data}{The raw data as a data frame when available, else `NULL`.}
#'   \item{composites}{The cluster scores of a score path, else `NULL`.}
#'   \item{meta}{List: `aggregation` (the name as the caller spelled it),
#'     `method`, `within`, `weighting`, `weighting_spec` (the scheme itself,
#'     a custom vector or function included), `fa_method`, `fa_args`,
#'     `scale`, `cor_method`, `signed`, `n_nodes`, `n_clusters`,
#'     `cluster_sizes`, `n_misfit`, `n_flipped`, `directed = FALSE`,
#'     `source = "pc"` and `experimental = TRUE`. `method` and `weighting`
#'     are `NA` on the descriptive paths, which re-estimate no macro
#'     network, and `fa_method` / `fa_args` are set only for
#'     `weighting = "factor"`.}
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
#' fit            # macro (cluster-level) weights
#' summary(fit)   # one row per macro edge
#' plot(fit)      # heatmap of the same weights
#'
#' @export
build_mcml_pc <- function(x,
                          clusters,
                          aggregation = c("scaled", "composite", "mean",
                                          "median",
                                          "loadings", "average", "escoufier",
                                          "cancor"),
                          method = c("pcor", "glasso", "cor"),
                          within = c("reestimate", "subnetwork"),
                          weighting = c("equal", "strength", "eigen",
                                        "closeness", "betweenness",
                                        "expected_influence", "specificity",
                                        "pca", "factor", "item_total"),
                          cor_method = c("pearson", "spearman",
                                         "polychoric"),
                          signed = TRUE,
                          id_col = NULL,
                          fa_method = c("ml", "paf", "minres", "cfa"),
                          ...) {
  aggregation <- match.arg(aggregation)
  # Public label kept for printing / metadata; the strategy logic below runs
  # on the internal "composite" token.
  aggregation_label <- aggregation
  # "scaled" / "composite" / "mean" / "median" are the re-estimated score path:
  # they differ only in how each cluster score is formed (z-scored mean / raw
  # mean / raw median). Resolve them to the internal composite strategy.
  score_fn <- "mean"
  scale <- TRUE
  if (aggregation %in% c("scaled", "composite", "mean", "median")) {
    score_fn    <- if (aggregation == "median") "median" else "mean"
    scale       <- aggregation %in% c("scaled", "composite")
    aggregation <- "composite"
  }
  # Public names for the two descriptive similarity coefficients map onto the
  # internal strategy tokens.
  if (aggregation == "escoufier") aggregation <- "rv"
  if (aggregation == "cancor")    aggregation <- "canonical"
  weighting_given <- !missing(weighting)
  custom_weighting <- !is.character(weighting)
  if (custom_weighting) {
    if (!(is.function(weighting) ||
          (is.numeric(weighting) && !is.null(names(weighting))))) {
      stop(errorCondition(
        paste0(
          "'weighting' must be a scheme name, a named numeric vector ",
          "of item weights, or a function(W_block, data_block, nodes)."),
        class = "nestimate_bad_input", call = NULL))
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
    stop(errorCondition(
      "Arguments passed via ... must be named.",
      class = "nestimate_bad_input", call = NULL))
  }
  if ("estimator" %in% names(fa_args) && fa_method != "cfa") {
    stop(errorCondition(
      paste0(
        "'estimator' in ... is the lavaan estimator (e.g. \"WLSMV\") ",
        "and applies only to fa_method = \"cfa\". The network ",
        "estimator is selected by 'method'."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (length(fa_args) > 0 && !identical(weighting, "factor")) {
    stop(errorCondition(
      paste0(
        "Unused argument(s): ",
        paste(names(fa_args), collapse = ", "),
        ". Arguments in ... are forwarded to the factor-analysis ",
        "backend and require weighting = \"factor\"."),
      class = "nestimate_bad_input", call = NULL))
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
    # "paf" reads max_iter and tol, "minres" takes no argument: say which of
    # the rest are dropped instead of dropping them without a word.
    unused <- switch(fa_method,
      paf = setdiff(names(fa_args), c("max_iter", "tol")),
      minres = names(fa_args),
      character())
    if (length(unused) > 0) {
      warning(warningCondition(
        paste0("Argument(s) in ... ignored by fa_method = \"", fa_method,
               "\": ", paste(unused, collapse = ", "), "."),
        class = "nestimate_ignored_argument", call = NULL))
    }
  }
  .assert(
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
      stop(errorCondition(
        paste0(
          "'x' is a directed network (method = \"",
          x$method,
          "\"). ",
          "build_mcml_pc() is for undirected association networks; ",
          "use build_mcml() for transition networks."),
        class = "nestimate_bad_input", call = NULL))
    }
  } else if (is.data.frame(x)) {
    if (!is.null(id_col)) {
      .assert(is.character(id_col))
      missing_ids <- setdiff(id_col, names(x))
      if (length(missing_ids) > 0) {
        stop(errorCondition(paste0("id_col column(s) not found in data: ",
             paste(missing_ids, collapse = ", ")),
             class = "nestimate_bad_input", call = NULL))
      }
      x <- x[, setdiff(names(x), id_col), drop = FALSE]
    }
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop(errorCondition(
        paste0(
          "Data frame input must be all-numeric (observations x ",
          "variables). Non-numeric columns: ",
          paste(names(x)[!numeric_cols], collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
    data <- x
    node_network <- .pc_estimate(data, method, cor_method)
  } else {
    stop(errorCondition(
      "'x' must be an undirected netobject or a numeric data.frame.",
      class = "nestimate_bad_input", call = NULL))
  }
  W <- node_network$weights
  if (!is.matrix(W) || is.null(rownames(W))) {
    stop(errorCondition(
      "Could not extract a labeled weight matrix from 'x'.",
      class = "nestimate_bad_input", call = NULL))
  }
  node_labels <- rownames(W)

  has_data <- is.data.frame(data) && nrow(data) > 0L &&
    all(node_labels %in% colnames(data))
  needs_data <- c("composite", "loadings", "rv", "canonical")
  if (aggregation %in% needs_data && !has_data) {
    stop(errorCondition(
      paste0(
        "aggregation = \"",
        aggregation_label,
        "\" requires raw data with ",
        "one column per node. Pass a netobject that carries $data, or the ",
        "raw data.frame itself. For weight-only input use ",
        "aggregation = \"average\"."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (within == "reestimate" && !has_data) {
    within <- "subnetwork"
  }

  # ---- Normalize cluster membership ----
  members <- .pc_normalize_clusters(clusters, node_labels)
  k <- length(members)
  if (k < 2) {
    stop(errorCondition(
      "At least 2 clusters are required.",
      class = "nestimate_bad_input", call = NULL))
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
            ". See item_loadings(fit, misfit = TRUE).", call. = FALSE)
  }
  flipped_items <- loadings_df$node[loadings_df$sign < 0]
  if (signed && length(flipped_items) > 0 && aggregation == "composite") {
    warning("Reverse-keyed item(s) flipped in composites: ",
            paste(flipped_items, collapse = ", "),
            ". See item_loadings() (column sign).", call. = FALSE)
  }

  # ---- Macro network ----
  composites <- NULL
  macro <- switch(aggregation,
    average = .wrap_netobject(.pc_average_macro(W, members), data = NULL,
                              method = "pc_average", directed = FALSE),
    composite = {
      composites <- .pc_composites(data, members, loadings_df, scale,
                                   signed = signed, score_fn = score_fn)
      # keep the input's row identities: the estimator drops incomplete rows
      # and renumbers, so composites() must not read the estimator's data
      rownames(composites) <- rownames(data)
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
    composites = composites,
    meta = list(
      aggregation = aggregation_label,
      method = if (aggregation == "composite") method
        else NA_character_,
      within = within,
      weighting = if (aggregation != "composite") NA_character_
        else if (is.function(weighting)) "custom (function)"
        else if (is.numeric(weighting)) "custom (vector)"
        else weighting,
      # Retain custom vectors/functions as well as the printable scheme name
      # so case-resampling can replay the fitted weighting exactly.
      weighting_spec = if (aggregation == "composite") weighting else "strength",
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
      stop(errorCondition(
        "cor_method = \"polychoric\" requires the 'lavaan' package.",
        class = "nestimate_missing_package", call = NULL))
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
  if (is.data.frame(clusters)) {
    # Two-column data.frame read by position: first column = node names
    # (the item / variable columns), second column = group labels. Column
    # names are irrelevant, matching build_mcml()'s data.frame handling.
    if (ncol(clusters) < 2L) {
      stop(errorCondition(
        paste0(
          "'clusters' data.frame must have at least two columns ",
          "(node names, then group labels)."),
        class = "nestimate_bad_input", call = NULL))
    }
    members <- split(as.character(clusters[[1L]]),
                     as.character(clusters[[2L]]))
  } else if (is.list(clusters)) {
    if (is.null(names(clusters)) || any(!nzchar(names(clusters)))) {
      stop(errorCondition(
        "'clusters' list must have non-empty names.",
        class = "nestimate_bad_input", call = NULL))
    }
    members <- lapply(clusters, as.character)
  } else if (is.atomic(clusters) && !is.null(names(clusters))) {
    members <- split(names(clusters), as.character(clusters))
  } else {
    stop(errorCondition(
      paste0(
        "'clusters' must be a named list of node labels, a two-column ",
        "data.frame (node, group), or a vector of cluster labels named ",
        "by node."),
      class = "nestimate_bad_input", call = NULL))
  }
  assigned <- unlist(members, use.names = FALSE)
  if (anyDuplicated(assigned)) {
    stop(errorCondition(
      paste0(
        "Each node may belong to exactly one cluster. Duplicated: ",
        paste(unique(assigned[duplicated(assigned)]), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  unknown <- setdiff(assigned, node_labels)
  if (length(unknown) > 0) {
    stop(errorCondition(paste0("Cluster members not present in the network: ",
         paste(unknown, collapse = ", ")),
         class = "nestimate_bad_input", call = NULL))
  }
  missing_nodes <- setdiff(node_labels, assigned)
  if (length(missing_nodes) > 0) {
    stop(errorCondition(paste0("Nodes without a cluster assignment: ",
         paste(missing_nodes, collapse = ", ")),
         class = "nestimate_bad_input", call = NULL))
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
.pc_composites <- function(data, members, loadings_df, scale, signed,
                           score_fn = "mean") {
  cols <- lapply(names(members), function(cl) {
    nodes <- members[[cl]]
    mat <- as.matrix(data[, nodes, drop = FALSE])
    if (scale) mat <- base::scale(mat)
    idx <- match(nodes, loadings_df$node)
    sgn <- if (signed) loadings_df$sign[idx] else rep(1, length(nodes))
    if (identical(score_fn, "median")) {
      # Row-wise median over observed members (sign-corrected). Median is
      # unweighted by definition, so item weights do not apply here.
      smat <- sweep(mat, 2L, sgn, "*")
      out <- apply(smat, 1L, stats::median, na.rm = TRUE)
      out[rowSums(!is.na(smat)) == 0L] <- NA_real_
      return(out)
    }
    # Row-wise weighted mean over observed members: renormalize the
    # absolute weights across non-missing cells per row.
    w <- loadings_df$weight[idx] * sgn
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
#' `loading_stability()` bootstraps the item weights of a [build_mcml_pc()]
#' fit and returns a `pc_loading_stability` object with each item's observed
#' signed weight, its bootstrap mean and standard deviation, a percentile
#' interval and how often its sign flipped. Each replicate resamples the rows
#' of the raw data with replacement, re-estimates the node-level network
#' with the fit's estimator (`cor`, `pcor` or `glasso`; any other node-level
#' method is re-estimated as `cor`) and correlation type, and reapplies the
#' fitted weighting scheme, data-dependent and custom weights and the
#' factor-analysis options included. A fit with a descriptive aggregation
#' (`"average"`, `"escoufier"`, `"cancor"`) has its item-table weights, the
#' `"strength"` weights, resampled. \strong{Experimental.}
#'
#' Wide intervals mean the weighting should not be over-interpreted. Equal
#' weights keep their magnitudes in every replicate, so for them only the
#' sign flips carry information; a sign that flips often marks an item whose
#' reverse-keying is not settled by the data. Supplied numeric weights are
#' fixed in magnitude and sign, so their intervals have zero width and
#' resampling them tells nothing. The intervals describe the weights, not
#' the macro network: edge-weight intervals for it come from
#' [bootstrap_network()] on its macro network.
#'
#' @param x An `mcml_pc` object that carries its raw data; one built from a
#'   network without data raises an error, and so does an older saved fit
#'   with custom weights, which must be rebuilt with [build_mcml_pc()] to
#'   keep its weighting. For the `print()` and `plot()` methods: a
#'   `pc_loading_stability` object.
#' @param iter The number of bootstrap replicates, at least 2. Default `200`;
#'   each replicate re-estimates the node-level network, so this is slower
#'   than a plain bootstrap.
#' @param ci_level The tail level of the percentile intervals: `0.05`
#'   (default) gives 95% intervals.
#' @param seed An integer that makes the resampling reproducible. Default
#'   `NULL`. The session's own random stream is left as it was.
#'
#' @return A `pc_loading_stability` object. `print()` shows the item table
#'   and `summary()` returns it;
#'   `plot()` draws each item's signed weight with its interval, one panel
#'   per cluster. The object's fields:
#'   \describe{
#'     \item{summary}{A data frame with one row per item: `node`, `cluster`,
#'       `weight` (the observed signed weight), `boot_mean`, `boot_sd`,
#'       `ci_lower`, `ci_upper` and `sign_flips` (the share of replicates in
#'       which the item's sign differed from the observed one).}
#'     \item{boot_weights}{The `iter` x items matrix of resampled signed
#'       weights.}
#'     \item{iter, ci_level}{The settings used.}
#'   }
#'
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(100)
#' g <- stats::rnorm(100)
#' df <- data.frame(a1 = f + stats::rnorm(100), a2 = f + stats::rnorm(100),
#'                  a3 = f + stats::rnorm(100), b1 = g + stats::rnorm(100),
#'                  b2 = g + stats::rnorm(100), b3 = g + stats::rnorm(100))
#' cl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#' fit <- build_mcml_pc(df, cl, aggregation = "loadings",
#'                      method = "cor")
#' stability <- loading_stability(fit, iter = 50, seed = 1)
#' stability
#'
#' @export
loading_stability <- function(x, iter = 200L, ci_level = 0.05,
                              seed = NULL) {
  if (!inherits(x, "mcml_pc")) {
    stop(errorCondition(
      "'x' must be an mcml_pc object from build_mcml_pc().",
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.null(x$data)) {
    stop(errorCondition(
      paste0(
        "This mcml_pc object carries no raw data; loading_stability() ",
        "needs the observations to resample."),
      class = "nestimate_bad_input", call = NULL))
  }
  .assert(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(ci_level), length(ci_level) == 1,
    ci_level > 0, ci_level < 1
  )
  iter <- as.integer(iter)
  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
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
  weighting <- x$meta$weighting_spec %||% x$meta$weighting
  # Older descriptive fits stored NA, and their diagnostics used strength.
  if (is.null(weighting) ||
      (is.character(weighting) && length(weighting) == 1L && is.na(weighting))) {
    weighting <- "strength"
  }
  if (is.character(weighting) && any(grepl("^custom", weighting))) {
    stop(errorCondition(
      paste0(
        "This saved fit does not retain its custom weighting. Rebuild it ",
        "with build_mcml_pc() before calling loading_stability()."),
      class = "nestimate_bad_input", call = NULL))
  }

  # A factor extraction that stops unconverged warns once per cluster and
  # replicate; they are counted here and reported in one warning below.
  n_unconverged <- 0L
  boot <- vapply(seq_len(iter), function(b) {
    rows <- sample.int(n_rows, n_rows, replace = TRUE)
    data_b <- data[rows, , drop = FALSE]
    net_b <- .pc_estimate(data_b, estimator, cor_method)
    ld_b <- suppressWarnings(
      .pc_network_loadings(net_b$weights, members, signed = signed)
    )
    if (!identical(weighting, "strength")) {
      ld_b <- withCallingHandlers(
        .pc_apply_weighting(ld_b, net_b$weights, data_b, members,
                            weighting, signed,
                            fa_method = x$meta$fa_method %||% "ml",
                            cor_method = cor_method,
                            fa_args = x$meta$fa_args %||% list()),
        nestimate_no_converge = function(w) {
          n_unconverged <<- n_unconverged + 1L
          invokeRestart("muffleWarning")
        })
    }
    idx <- match(node_order, ld_b$node)
    ld_b$weight[idx] * ld_b$sign[idx]
  }, numeric(length(node_order)))
  boot <- t(boot)
  colnames(boot) <- node_order
  if (n_unconverged > 0L) {
    warning(warningCondition(
      sprintf(paste0("loading_stability(): the factor extraction stopped ",
                     "unconverged in %d cluster fit(s) over %d replicates; ",
                     "their last-iteration loadings are used."),
              n_unconverged, iter),
      class = "nestimate_no_converge", call = NULL))
  }

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


#' @rdname loading_stability
#' @param digits For the `print()` method: the decimal places of the
#'   numeric columns. Default `3`.
#' @param ... For the `print()`, `summary()` and `plot()` methods: ignored.
#' @return For the methods: `print()` returns `x` invisibly, `plot()` a
#'   `ggplot`.
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


#' @rdname loading_stability
#' @param object For the `summary()` method: a `pc_loading_stability` object.
#' @return For the `summary()` method: the item table, a data frame with one
#'   row per item and the columns listed under `summary` above.
#' @export
summary.pc_loading_stability <- function(object, ...) {
  object$summary
}


#' @rdname loading_stability
#' @section Methods:
#' * `plot()` draws each item's observed signed composite weight with its
#'   bootstrap percentile interval, one panel per cluster, with a dashed line
#'   at zero.
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


#' @rdname build_mcml_pc
#' @param digits For the `print()` and `plot()` methods: the decimal places
#'   of the printed macro weights (default `3`) or of the tile labels
#'   (default `2`).
#' @param object For the `summary()` method: an `mcml_pc` object.
#' @return For the methods: `print()` returns `x` invisibly; `summary()`
#'   returns a data frame with one row per macro edge (the upper triangle)
#'   and columns `from`, `to` and `weight`; `plot()` returns a `ggplot`.
#' @export
print.mcml_pc <- function(x, digits = 3, ...) {
  m <- x$meta
  cat("MCML for Psychometric Networks (experimental)\n")
  # "Estimator" is the macro network's. A descriptive aggregation re-estimates
  # no macro network, but within = "reestimate" still re-estimates the
  # within-cluster ones: name their estimator beside it.
  est <- if (is.na(m$method)) "none (no re-estimation)" else m$method
  wt <- if (is.na(m$weighting)) "" else paste0(" | Weights: ", m$weighting)
  within_method <- unique(unlist(lapply(x$clusters, function(net) net$method)))
  within_lab <- if (is.na(m$method) && identical(m$within, "reestimate") &&
                    length(within_method) == 1L) {
    sprintf("%s (%s)", m$within, within_method)
  } else {
    m$within
  }
  cat(sprintf("  Aggregation: %s | Estimator: %s | Within: %s | cor: %s%s\n",
              m$aggregation, est, within_lab, m$cor_method, wt))
  cat(sprintf("  %d nodes in %d clusters (%s)\n",
              m$n_nodes, m$n_clusters,
              paste(names(m$cluster_sizes), m$cluster_sizes,
                    sep = ": ", collapse = ", ")))
  if (m$n_misfit > 0) {
    cat(sprintf("  ! %d misfit item(s) - see item_loadings()\n", m$n_misfit))
  }
  # Items are flipped only where scores are formed (the composite path with
  # signed = TRUE); the descriptive aggregations flip nothing.
  if (m$n_flipped > 0 && !is.na(m$method) && !isFALSE(m$signed)) {
    cat(sprintf("  ! %d reverse-keyed item(s) flipped - see item_loadings()\n",
                m$n_flipped))
  }
  cat("\nMacro weights:\n")
  print(round(x$macro$weights, digits))
  invisible(x)
}


#' @rdname build_mcml_pc
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


#' @rdname build_mcml_pc
#' @section Methods:
#' * `plot()` draws the macro (cluster-level) weights as a heatmap in the
#'   package's diverging palette, with each weight printed in its tile. For
#'   the two-layer network drawing use `cograph::plot_mcml()`, which accepts
#'   `mcml_pc` objects and draws them undirected.
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


#' Item Diagnostics From a Psychometric MCML Fit
#'
#' @description
#' `item_loadings()` returns the item table behind a [build_mcml_pc()] fit:
#' one row per item, with how strongly the item connects to its own cluster,
#' the weight it carries into that cluster's score, and whether it connects
#' more strongly to some other cluster. The loadings are network loadings,
#' read from the node-level network (an item's mean absolute connection to
#' the members of a cluster), not factor loadings.
#'
#' An item flagged as `misfit` connects more strongly to another cluster than
#' to its own, which is evidence that it is assigned to the wrong cluster;
#' every aggregation of the fit inherits that membership. A `sign` of -1
#' marks an item treated as reverse-keyed. How stable the weights are under
#' resampling is checked with [loading_stability()].
#'
#' @param x An `mcml_pc` object from [build_mcml_pc()]. A fit that carries no
#'   item table raises an error of class `nestimate_no_loadings`.
#' @param misfit Which items to return: `NULL` (default) returns every item,
#'   `TRUE` only those whose strongest cross-cluster connection exceeds their
#'   own-cluster loading, and `FALSE` only those that fit where they were
#'   assigned.
#' @param ... Ignored.
#'
#' @return A data frame with one row per item (only the misfitting or the
#'   fitting items when `misfit` is set) and columns `node`, `cluster`,
#'   `loading` (the mean absolute connection to the other members of its own
#'   cluster, signed by the leading eigenvector of the cluster's block; under
#'   a sign-carrying weighting the `sign` column holds the sign actually
#'   used, which can differ), `weight` (its weight in the cluster
#'   score, or its strength weight for a descriptive aggregation; the weights
#'   of a cluster sum to 1), `sign` (+1, or -1 for a
#'   reverse-keyed item), `max_cross` (its strongest mean absolute
#'   connection to any other cluster), `cross_cluster` (which cluster that
#'   is) and `misfit` (logical). An item alone in its cluster has loading 0,
#'   weight 1 and is never flagged.
#'
#' @seealso \code{\link{build_mcml_pc}} to create the fit,
#'   \code{\link{loading_stability}} for the weights' sampling
#'   uncertainty, \code{\link{composites}} for the scores these weights
#'   produce.
#'
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(200)
#' g <- stats::rnorm(200)
#' df <- data.frame(a1 = f + stats::rnorm(200), a2 = f + stats::rnorm(200),
#'                  a3 = f + stats::rnorm(200), b1 = g + stats::rnorm(200),
#'                  b2 = g + stats::rnorm(200), b3 = g + stats::rnorm(200))
#' clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#' fit <- build_mcml_pc(df, clusters, aggregation = "loadings",
#'                      method = "cor")
#' item_loadings(fit)
#' item_loadings(fit, misfit = TRUE)
#'
#' @export
item_loadings <- function(x, ...) {
  UseMethod("item_loadings")
}

#' @rdname item_loadings
#' @export
item_loadings.mcml_pc <- function(x, misfit = NULL, ...) {
  .assert(
    "`misfit` must be TRUE, FALSE or NULL" =
      is.null(misfit) ||
      (is.logical(misfit) && length(misfit) == 1L && !is.na(misfit))
  )
  table <- x$loadings
  if (is.null(table)) {
    stop(errorCondition(
      paste0("This fit carries no item diagnostics: it was built without a ",
             "node-level network."),
      class = "nestimate_no_loadings", call = NULL))
  }
  if (is.null(misfit)) {
    return(table)
  }
  out <- table[table$misfit == misfit, , drop = FALSE]
  row.names(out) <- NULL
  out
}


#' Cluster Scores From a Psychometric MCML Fit
#'
#' @description
#' `composites()` returns the per-observation cluster scores that the macro
#' network of a re-estimated [build_mcml_pc()] fit was estimated from: each
#' respondent's weighted, sign-corrected score on every cluster. These are
#' the scores to carry into a profile analysis, a regression or any other
#' model that needs one number per cluster per respondent.
#'
#' A score is only as good as the items behind it: check the item weights
#' and misfit flags with [item_loadings()] and the stability of the weights
#' with [loading_stability()] before using the scores elsewhere.
#'
#' @param x An `mcml_pc` fit with a score aggregation (`"scaled"`,
#'   `"composite"`, `"mean"`, `"median"` or `"loadings"`). The descriptive
#'   aggregations (`"average"`, `"escoufier"`, `"cancor"`) relate the clusters
#'   without ever forming a score and raise an error of class
#'   `nestimate_no_composites`.
#' @param ... Ignored.
#'
#' @return A data frame with one row per row of the input data, in input
#'   order and with the input's row names, and one numeric column per
#'   cluster, named by the cluster (as a syntactic name, so `"Clu A"` becomes
#'   `Clu.A`). A row whose members of a cluster are all missing is `NA` in
#'   that column; items missing only in part are averaged over the observed
#'   ones. The macro network is estimated on the complete rows, so
#'   `build_network(composites(fit), method = ...)` with the fit's estimator
#'   reproduces it.
#'
#' @seealso \code{\link{build_mcml_pc}} to create the fit,
#'   \code{\link{item_loadings}} for the item weights behind these scores.
#'
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(200)
#' g <- stats::rnorm(200)
#' df <- data.frame(a1 = f + stats::rnorm(200), a2 = f + stats::rnorm(200),
#'                  a3 = f + stats::rnorm(200), b1 = g + stats::rnorm(200),
#'                  b2 = g + stats::rnorm(200), b3 = g + stats::rnorm(200))
#' clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#' fit <- build_mcml_pc(df, clusters, aggregation = "loadings",
#'                      method = "cor")
#' head(composites(fit))
#'
#' @export
composites <- function(x, ...) {
  UseMethod("composites")
}

#' @rdname composites
#' @export
composites.mcml_pc <- function(x, ...) {
  scores <- x$composites
  if (is.null(scores)) {
    stop(errorCondition(
      sprintf(paste0("aggregation = \"%s\" relates the item blocks directly ",
                     "and forms no cluster scores. Use a re-estimated ",
                     "aggregation (\"scaled\", \"composite\", \"mean\", ",
                     "\"median\" or \"loadings\") to get composites."),
              x$meta$aggregation),
      class = "nestimate_no_composites", call = NULL))
  }
  as.data.frame(scores, stringsAsFactors = FALSE)
}


#' Promote a psychometric MCML result to a network group
#'
#' @description
#' `as_networks()` gathers the cluster-level (macro) network and the
#' within-cluster networks of a [build_mcml_pc()] fit into one
#' `netobject_group`, so they go through the same verbs as any other group of
#' networks: `print()`, `summary()`, [net_centrality()], and
#' `cograph::splot()` for drawing. It is the psychometric-network
#' counterpart of \code{\link[=as_tna]{as_tna()}}.
#'
#' \code{\link[=as_tna]{as_tna()}} converts \emph{transition} networks (directed, row-normalised,
#' with initial probabilities) and re-wraps raw matrices. The macro and
#' within-cluster networks of an `mcml_pc` object are already full
#' netobjects of \emph{psychometric} networks (undirected; correlation,
#' partial correlation or glasso) carrying their estimator and directedness,
#' and their data when they were re-estimated, so `as_networks()` only
#' assembles them into a group.
#'
#' @param x The object to convert: an `mcml_pc` fit from [build_mcml_pc()],
#'   or a `netobject_group`, which is returned unchanged. Any other object
#'   raises an error.
#'
#' @return A `netobject_group`, read by `print()`, `summary()` and
#'   [net_centrality()] and drawn by `cograph::splot()`: a named list whose
#'   first element, `macro`, is the cluster-level network, followed by one
#'   netobject per cluster with more than one member. A singleton cluster has no within-cluster network
#'   and is left out, with a warning naming it.
#'
#' @seealso \code{\link{build_mcml_pc}} to create the input,
#'   \code{\link{as_tna}} for the transition-network counterpart.
#'
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(200)
#' g <- stats::rnorm(200)
#' df <- data.frame(a1 = f + stats::rnorm(200), a2 = f + stats::rnorm(200),
#'                  a3 = f + stats::rnorm(200), b1 = g + stats::rnorm(200),
#'                  b2 = g + stats::rnorm(200), b3 = g + stats::rnorm(200))
#' clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
#' fit <- build_mcml_pc(df, clusters, aggregation = "composite", method = "cor")
#' nets <- as_networks(fit)
#' nets
#' summary(nets)
#' @export
as_networks <- function(x) {
  UseMethod("as_networks")
}

#' @rdname as_networks
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
#' @export
as_networks.default <- function(x) {
  if (inherits(x, "netobject_group")) {
    return(x)
  }
  stop(errorCondition(paste0("Cannot convert object of class '", class(x)[1],
       "' to a network group."), class = "nestimate_bad_input", call = NULL))
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
        stop(errorCondition(
          paste0(
            "Custom weighting function must return ",
            m,
            " finite ",
            "numeric values for cluster '",
            cl,
            "'."),
          class = "nestimate_bad_input", call = NULL))
      }
      val
    } else if (is.numeric(weighting)) {
      missing_w <- setdiff(nodes, names(weighting))
      if (length(missing_w) > 0) {
        stop(errorCondition(
          paste0(
            "Custom weighting vector lacks entries for: ",
            paste(missing_w, collapse = ", ")),
          class = "nestimate_bad_input", call = NULL))
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
      stop(errorCondition(
        "cor_method = \"polychoric\" requires the 'lavaan' package.",
        class = "nestimate_missing_package", call = NULL))
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
      stop(errorCondition(
        "fa_method = \"cfa\" requires the 'lavaan' package.",
        class = "nestimate_missing_package", call = NULL))
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
  converged <- FALSE
  # A loop on purpose: a fixed-point iteration, each step starts from the
  # communalities of the one before and stops at convergence.
  for (it in seq_len(max_iter)) {
    Rr <- R
    diag(Rr) <- h2
    e <- eigen(Rr, symmetric = TRUE)
    lambda <- e$vectors[, 1] * sqrt(max(e$values[1], 0))
    h2_new <- pmin(lambda^2, 0.998)
    if (max(abs(h2_new - h2)) < tol) {
      converged <- TRUE
      break
    }
    h2 <- h2_new
  }
  if (!converged) {
    warning(warningCondition(
      sprintf(paste0("Principal axis factoring did not converge in %d ",
                     "iterations; the loadings of the last iteration are ",
                     "used. Raise max_iter or tol."), as.integer(max_iter)),
      class = "nestimate_no_converge", call = NULL))
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
  if (fit$convergence != 0L) {
    warning(warningCondition(
      sprintf(paste0("Minimum-residual factoring did not converge (optim ",
                     "code %d%s); the loadings at the last iterate are used."),
              fit$convergence,
              if (is.null(fit$message)) "" else paste0(": ", fit$message)),
      class = "nestimate_no_converge", call = NULL))
  }
  lambda_for(fit$par)
}
