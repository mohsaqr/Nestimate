# Cluster Metrics for Network Analysis
# Summary measures for between/within clusters and multilayer networks


# Tagged-list constructor for an mcml layer (macro or one cluster). The
# class buys us a `print.mcml_layer` method so `print(mc$macro)` shows a
# tidy summary instead of dumping every named element of the list.
# Field access (`$weights`, `$inits`, `$labels`, `$data`) is unchanged.
.mcml_layer <- function(weights, inits, labels, data = NULL) {
  if (!is.matrix(weights) || !is.numeric(weights)) {
    stop("'weights' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(weights) != ncol(weights)) {
    stop("'weights' must be a square matrix.", call. = FALSE)
  }
  .validate_mcml_matrix(weights)
  if (!is.character(labels) || length(labels) != nrow(weights) ||
      any(is.na(labels)) || any(!nzchar(labels))) {
    stop("'labels' must be a non-missing character vector with one value per node.",
         call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop("'labels' must be unique.", call. = FALSE)
  }
  if (!is.null(inits)) {
    if (!is.numeric(inits) || length(inits) != length(labels) ||
        any(is.na(inits) | !is.finite(inits))) {
      stop("'inits' must be a finite numeric vector with one value per node.",
           call. = FALSE)
    }
    if (!is.null(names(inits))) {
      missing_inits <- setdiff(labels, names(inits))
      extra_inits <- setdiff(names(inits), labels)
      if (length(missing_inits) > 0L || length(extra_inits) > 0L) {
        stop("'inits' names must match layer labels.", call. = FALSE)
      }
      inits <- inits[labels]
    } else {
      names(inits) <- labels
    }
  }
  obj <- list(weights = weights, inits = inits, labels = labels, data = data)
  class(obj) <- "mcml_layer"
  obj
}


#' @rdname build_mcml
#' @export
print.mcml_layer <- function(x, ...) {
  .mcml_check_unused_dots("print.mcml_layer", ...)
  w  <- x$weights
  if (is.null(w) || nrow(w) == 0L) {
    cat("MCML layer  [empty]\n"); return(invisible(x))
  }

  cat("MCML layer (transition probabilities)  [directed]\n")
  cat(sprintf("  Nodes: %d  |  Non-zero edges: %d\n",
              nrow(w), sum(w != 0)))

  # ---- Weight summary line ----
  nz <- w[w != 0]
  if (length(nz) > 0) {
    cat(sprintf("  Weights: [%.3f, %.3f]  |  mean: %.3f\n",
                min(nz), max(nz), mean(nz)))
  }

  # ---- Weight matrix ----
  cat("\n  Weight matrix:\n")
  digits <- if (all(nz == floor(nz))) 0L else 3L
  mat_r <- round(w, digits)
  if (!is.null(x$labels)) dimnames(mat_r) <- list(x$labels, x$labels)
  formatted <- utils::capture.output(print(mat_r))
  cat(paste0("  ", formatted, collapse = "\n"), "\n")

  # ---- Initial probabilities (bar plot, same style as netobject) ----
  if (!is.null(x$inits) && length(x$inits) > 0L) {
    cat("\n  Initial probabilities:\n")
    init  <- x$inits
    nm    <- names(init)
    if (is.null(nm) && !is.null(x$labels)) nm <- x$labels
    if (is.null(nm)) nm <- as.character(seq_along(init))
    ord   <- order(init, decreasing = TRUE)
    bar_w <- 40L
    max_v <- max(init, na.rm = TRUE)
    for (i in ord) {
      bars <- if (is.finite(max_v) && max_v > 0)
        strrep("#", round(init[i] / max_v * bar_w)) else ""
      cat(sprintf("  %-14s  %.3f  %s\n", nm[i], init[i], bars))
    }
  }

  # ---- Data dimensions ----
  d <- x$data
  if (!is.null(d)) {
    data_dim <- if (is.data.frame(d) || is.matrix(d))
      sprintf("%d x %d", nrow(d), ncol(d))
    else sprintf("length %d", length(d))
    cat(sprintf("\n  Data: %s\n", data_dim))
  }

  invisible(x)
}

# ==============================================================================
# 1. Edge Weight Aggregation
# ==============================================================================

#' Aggregate Edge Weights
#'
#' @description
#' `net_aggregate_weights()` reduces a vector of edge weights to one number:
#' their sum, mean, median, maximum, minimum, product, density or geometric
#' mean. Missing and zero weights are removed first, so every method works on
#' the edges that are present. [cluster_summary()] and [build_mcml()] use it to
#' combine the edges between two clusters into one cluster-level weight.
#'
#' @param w The edge weights, a numeric vector. `NA` is allowed and removed
#'   with the zeros; an infinite value raises an error.
#' @param method How to combine the weights: `"sum"` (default), `"mean"`,
#'   `"median"`, `"max"`, `"min"`, `"prod"`, `"density"` or `"geomean"`. A
#'   single value; any other raises an error. `"density"` divides the sum by
#'   `n_possible`, and `"mean"` by the number of non-zero weights, so the two
#'   agree when every possible edge is present, and `"density"` is the
#'   smaller when some are zero or missing. `"geomean"` is the geometric mean of the
#'   positive weights; negative weights are left out, and the result is `0`
#'   when none is positive.
#' @param n_possible For `method = "density"`: the number of possible edges,
#'   such as the block size `n_i * n_j` of two clusters, a single finite
#'   number. Default `NULL`, under which (as with a value of 0 or less)
#'   `"density"` returns the mean of the non-zero weights. Other methods
#'   ignore it.
#' @return A single number: the chosen aggregate of the non-zero, non-missing
#'   weights, or `0` when none remain.
#' @details
#' `"sum"`, `"prod"`, `"min"`, `"max"`, `"mean"` and `"median"` mean what
#' they mean in igraph's `edge.attr.comb` argument.
#' @export
#' @examples
#' w <- c(0.5, 0.8, 0.3, 0.9)
#' net_aggregate_weights(w, "sum")   # 2.5
#' net_aggregate_weights(w, "mean")  # 0.625
#' net_aggregate_weights(w, "max")   # 0.9
#' net_aggregate_weights(w, "density", n_possible = 9)  # 2.5 / 9
net_aggregate_weights <- function(w, method = "sum", n_possible = NULL) {
  if (!is.numeric(w)) {
    stop(errorCondition("'w' must be a numeric vector.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (any(!is.na(w) & !is.finite(w))) {
    stop(errorCondition("'w' must contain only finite values or NA.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop(errorCondition(
      "'method' must be a single non-missing character value.",
      class = "nestimate_bad_input", call = NULL))
  }
  valid_methods <- c("sum", "mean", "median", "max", "min",
                     "prod", "density", "geomean")
  if (!method %in% valid_methods) {
    stop(errorCondition(paste0("Unknown method: ", method),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(n_possible) &&
      (!is.numeric(n_possible) || length(n_possible) != 1L ||
       is.na(n_possible) || !is.finite(n_possible))) {
    stop(errorCondition(
      "'n_possible' must be a single finite numeric value or NULL.",
      class = "nestimate_bad_input", call = NULL))
  }

  # Remove NA and zero weights
  w <- w[!is.na(w) & w != 0]
  if (length(w) == 0) return(0)

  switch(method,
    "sum"     = sum(w),
    "mean"    = mean(w),
    "median"  = stats::median(w),
    "max"     = max(w),
    "min"     = min(w),
    "prod"    = prod(w),
    "density" = if (!is.null(n_possible) && n_possible > 0) {
      sum(w) / n_possible
    } else {
      sum(w) / length(w)
    },
    "geomean" = {
      pos_w <- w[w > 0]
      if (length(pos_w) == 0) 0 else exp(mean(log(pos_w)))
    },
    stop("Unknown method: ", method, call. = FALSE)
  )
}


# ==============================================================================
# 2. Cluster Summary (Between/Within Aggregates)
# ==============================================================================

#' Cluster Summary Statistics
#'
#' @description
#' `cluster_summary()` aggregates a node-level weight matrix to clusters of
#' nodes and returns an `mcml` object with two levels: a macro network whose
#' entry (i, j) combines every edge from a node of cluster i to a node of
#' cluster j (the diagonal combines the edges inside a cluster), and one
#' within-cluster network per cluster, the submatrix of its nodes. It takes a
#' weight matrix, a `netobject`, a bare `cograph_network` or a `tna` model,
#' and the clusters as a membership vector, a node-to-cluster data frame or a
#' named list. This is the matrix route of Multi-Cluster Multi-Level (MCML)
#' analysis; [build_mcml()] calls it for matrix input.
#'
#' The macro weights are arithmetic aggregates (see `method`) of weights that
#' already exist. They are not row-normalised and are not the transition
#' probabilities of a process over clusters: the transitions between
#' clusters can be counted only from the sequences, which [build_mcml()]
#' does (its default `type = "tna"` gives transition probabilities between
#' clusters). \code{\link[=as_tna]{as_tna()}} keeps the macro weights as
#' they are and does not normalise them.
#'
#' @param x The node-level network: a numeric square matrix whose row and
#'   column names (identical, in the same order) label the nodes, with
#'   finite weights such as transition counts, co-occurrence frequencies or
#'   probabilities; a `netobject` from [build_network()] or a bare
#'   `cograph_network`, whose weight matrix is used; or a `tna` model, whose
#'   weights are used. An `mcml` object is returned unchanged. A matrix
#'   without names labels its nodes `"1"`, `"2"`, ....
#' @param clusters The cluster of each node, required (`NULL`, the default,
#'   raises an error; [build_mcml()] is the verb that reads clusters from a
#'   `netobject`'s node table). Accepted forms:
#'   \describe{
#'     \item{membership vector}{One value per node, in the order of the
#'       matrix rows, such as `c(1, 1, 2, 2, 3, 3)` or `c("A", "A", "B")`. A
#'       named vector is matched to the nodes by name instead. Numeric values
#'       give clusters in sorted order, character values in order of first
#'       appearance.}
#'     \item{data frame}{Node names in the first column, cluster names in the
#'       second, such as
#'       `data.frame(node = c("A", "B", "C"), group = c("G1", "G1", "G2"))`;
#'       the clusters come in alphabetical order.}
#'     \item{named list}{Cluster names mapped to their nodes, such as
#'       `list(Alpha = c("A", "B"), Beta = c("C", "D"))`; the clusters keep
#'       the list's order.}
#'   }
#'   Every node must belong to exactly one cluster: an unknown node, a node
#'   listed twice, a node left out or an empty cluster raises an error.
#' @param method How the edges between (or within) two clusters are combined
#'   into one weight, with [net_aggregate_weights()]; zero edges are left out
#'   first.
#'   \describe{
#'     \item{`"sum"`}{(default) The total weight, which keeps the total flow.
#'       Suits counts such as transition frequencies.}
#'     \item{`"mean"`}{The mean weight, so a cluster of 5 nodes does not get
#'       5 times the weight of a cluster of 1. Suits a matrix of
#'       probabilities (rows summing to 1), where a sum would grow with
#'       cluster size.}
#'     \item{`"median"`}{The median weight, less affected by a few extreme
#'       weights.}
#'     \item{`"max"`}{The strongest connection.}
#'     \item{`"min"`}{The weakest connection.}
#'     \item{`"density"`}{The total weight divided by the number of possible
#'       edges between the two clusters (`n_i * n_j`). It equals `"mean"`
#'       when the block has no zero edges and is smaller than it when it has
#'       some.}
#'     \item{`"geomean"`}{The geometric mean of the positive weights, for
#'       multiplicative processes.}
#'   }
#' @param directed Whether the network is directed, so that A->B and B->A are
#'   separate edges. Default `TRUE`. `FALSE` replaces the matrix by the mean
#'   of it and its transpose before aggregating, which makes every layer
#'   symmetric.
#' @param compute_within Whether to build the within-cluster networks, one
#'   n_i x n_i matrix per cluster. Default `TRUE`; `FALSE` skips them when
#'   only the macro network is needed, and [summary()] then reports their
#'   totals as `NA`.
#'
#' @return An `mcml` object. Print it for the aggregation method, the node
#'   and cluster counts, the clusters and the macro weights; [summary()]
#'   returns a data frame with one row per cluster with its size, total absolute
#'   within-cluster weight and between-cluster weight sent and received
#'   (described under [build_mcml()]);
#'   \code{\link[=as_tna]{as_tna()}} turns the layers into a `netobject_group`;
#'   [simplicial_features()] summarises the topology of each within-cluster
#'   network. [macro_network()], which opens chosen clusters into their
#'   states, needs an `mcml` built from sequences by [build_mcml()]; on this
#'   result it raises an error of class `nestimate_no_expand_source`. The
#'   object's fields:
#'   \describe{
#'     \item{macro}{The cluster-level layer (class `mcml_layer`): `weights`,
#'       the k x k matrix of aggregated weights, diagonal included and rows
#'       not normalised; `inits`, the column sums of that matrix normalised
#'       to sum to 1 (equal shares when every weight is zero); `labels`, the
#'       cluster names; and `data`, `NULL` on this route.}
#'     \item{clusters}{A named list of one `mcml_layer` per cluster: the
#'       n_i x n_i submatrix of its nodes, self-loops included, and its
#'       normalised column sums as `inits` (1 for a single-node cluster).
#'       `NULL` when `compute_within = FALSE`.}
#'     \item{cluster_members}{A named list of each cluster's node labels,
#'       such as `list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))`.}
#'     \item{edges}{`NULL`: a matrix carries no node-level transitions. The
#'       sequence and edge-list routes of [build_mcml()] fill it.}
#'     \item{meta}{A list with `type` (always `"aggregate"` here), `method`,
#'       `directed`, `n_nodes`, `n_clusters`, `cluster_sizes` (a named
#'       integer vector) and `source` (`"matrix"`).}
#'   }
#'
#' @details
#' ## Workflow
#'
#' A typical MCML analysis:
#' \preformatted{
#' # 1. Create the node-level network
#' net <- build_network(data, method = "relative")
#'
#' # 2. Aggregate its edges to cluster level (arithmetic aggregation)
#' cs <- cluster_summary(net, clusters = group_assignments, method = "sum")
#'
#' # 3. Read the result
#' print(cs)      # macro weights
#' summary(cs)    # one row per cluster
#'
#' # 4. Promote the layers to netobjects for downstream verbs
#' nets <- as_tna(cs)
#' }
#'
#' ## The macro matrix
#'
#' The macro matrix has clusters as both rows and columns. An off-diagonal
#' entry (i, j) is the aggregated weight from cluster i to cluster j; a
#' diagonal entry (i, i) aggregates the edges inside cluster i.
#'
#' Rows are not normalised: each entry is the aggregate `method` produces.
#' Row-normalising an arithmetic aggregate here (the old `type = "tna"`
#' combined with `method = "min"`, `"mean"` and so on) gave rows that sum to
#' 1 but are not a probability distribution over any process, which is why
#' `type` was removed from the matrix route. The sequence and edge-list
#' routes of [build_mcml()] keep `type`: there the aggregate is always a
#' count, and `type` chooses between well-defined network constructions.
#'
#' ## Choosing method
#'
#' \tabular{lll}{
#'   \strong{Input data} \tab \strong{Suggested method} \tab \strong{Reason} \cr
#'   Edge counts \tab \code{"sum"} \tab Keeps the total flow between clusters \cr
#'   Transition matrix \tab \code{"mean"} \tab Avoids cluster-size bias \cr
#'   Correlation matrix \tab \code{"mean"} \tab Average correlations \cr
#'   Dense weighted \tab \code{"max"} / \code{"median"} \tab Less sensitive to single weights \cr
#' }
#'
#' @export
#' @seealso
#'   \code{\link{build_mcml}} to build an mcml from raw transitions instead
#'   of a weight matrix,
#'   \code{\link{as_tna}} to promote the layers to netobjects,
#'   \code{\link{macro_network}} for the cluster-level network with one
#'   cluster expanded back into its member states
#'
#' @examples
#' # -----------------------------------------------------
#' # Basic usage with matrix and cluster vector
#' # -----------------------------------------------------
#' set.seed(1)
#' mat <- matrix(runif(100), 10, 10)
#' rownames(mat) <- colnames(mat) <- LETTERS[1:10]
#'
#' cs <- cluster_summary(mat, c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
#' cs            # cluster-level (macro) weights
#' summary(cs)   # one row per cluster
#'
#' # -----------------------------------------------------
#' # Named list clusters (more readable)
#' # -----------------------------------------------------
#' clusters <- list(
#'   Alpha = c("A", "B", "C"),
#'   Beta = c("D", "E", "F"),
#'   Gamma = c("G", "H", "I", "J")
#' )
#' cluster_summary(mat, clusters)
#'
#' # -----------------------------------------------------
#' # A netobject as input: its weight matrix is aggregated
#' # -----------------------------------------------------
#' seqs <- data.frame(
#'   T1 = c("A", "C", "B", "D"), T2 = c("B", "D", "A", "C"),
#'   T3 = c("C", "A", "D", "B")
#' )
#' net <- build_network(seqs, method = "relative")
#' cluster_summary(net, list(G1 = c("A", "B"), G2 = c("C", "D")))
#'
#' # -----------------------------------------------------
#' # Different aggregation methods
#' # -----------------------------------------------------
#' summary(cluster_summary(mat, clusters, method = "sum"))   # total flow
#' summary(cluster_summary(mat, clusters, method = "mean"))  # average
#' summary(cluster_summary(mat, clusters, method = "max"))   # strongest
#'
#' # -----------------------------------------------------
#' # Skip within-cluster computation for speed
#' # -----------------------------------------------------
#' cluster_summary(mat, clusters, compute_within = FALSE)
#'
#' # -----------------------------------------------------
#' # Promote the layers to netobjects
#' # (as_tna() stores the aggregated weights as they are)
#' # -----------------------------------------------------
#' as_tna(cluster_summary(mat, clusters, method = "sum"))
cluster_summary <- function(x,
                            clusters = NULL,
                            method = c("sum", "mean", "median", "max",
                                       "min", "density", "geomean"),
                            directed = TRUE,
                            compute_within = TRUE) {
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop(errorCondition("'directed' must be TRUE or FALSE.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.logical(compute_within) || length(compute_within) != 1L ||
      is.na(compute_within)) {
    stop(errorCondition("'compute_within' must be TRUE or FALSE.",
      class = "nestimate_bad_input", call = NULL))
  }

  # If already an mcml object, return as-is
  if (inherits(x, "mcml")) {
    return(x) # nocov
  }

  method <- match.arg(method)

  # Store original for cluster extraction
  x_orig <- x

  # Extract matrix from various input types
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
    mat <- x$weights
  } else if (inherits(x, "netobject") || inherits(x, "netobject_ml")) {
    mat <- x$weights
  } else if (inherits(x, "tna")) {
    mat <- x$weights
  } else {
    mat <- x
  }

  if (is.null(clusters)) {
    stop(errorCondition("clusters argument is required",
      class = "nestimate_bad_input", call = NULL))
  }

  # Validate input matrix
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop(errorCondition("x must be a netobject, tna object, or numeric matrix",
      class = "nestimate_bad_input", call = NULL)) # nocov
  }
  if (nrow(mat) != ncol(mat)) {
    stop(errorCondition("x must be a square matrix",
      class = "nestimate_bad_input", call = NULL))
  }
  .validate_mcml_matrix(mat)

  # Symmetrize matrix when treating it as undirected, matching the docstring
  # contract. Without this, callers passing `directed = FALSE` got the same
  # output as `directed = TRUE`, which silently violated the contract.
  if (!isTRUE(directed)) mat <- (mat + t(mat)) / 2

  n <- nrow(mat)
  node_names <- rownames(mat)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))

  # Convert clusters to list format
  cluster_list <- .normalize_clusters(clusters, node_names)
  n_clusters <- length(cluster_list)
  cluster_names <- names(cluster_list)
  names(cluster_list) <- cluster_names

  # Get node indices for each cluster
  cluster_indices <- lapply(cluster_list, function(nodes_vec) {
    match(nodes_vec, node_names)
  })

  # ============================================================================
  # Between-cluster computation (always computed)
  # ============================================================================

  # Aggregate between-cluster weights
  between_raw <- matrix(0, n_clusters, n_clusters,
                        dimnames = list(cluster_names, cluster_names))

  for (i in seq_len(n_clusters)) {
    idx_i <- cluster_indices[[i]]
    n_i <- length(idx_i)

    for (j in seq_len(n_clusters)) {
      idx_j <- cluster_indices[[j]]
      n_j <- length(idx_j)

      # Both diagonal (within-cluster) and off-diagonal (between-cluster)
      w_ij <- mat[idx_i, idx_j]
      n_possible <- n_i * n_j
      between_raw[i, j] <- net_aggregate_weights(as.vector(w_ij), method, n_possible)
    }
  }

  # Matrix path is aggregation only -- no post-processing. Caller normalises
  # downstream if they need a Markov chain (e.g. via as_tna()).
  between_weights <- between_raw

  # Compute inits from column sums
  col_sums <- colSums(between_raw)
  total <- sum(col_sums)
  if (total > 0) {
    between_inits <- col_sums / total
  } else {
    between_inits <- rep(1 / n_clusters, n_clusters)
  }
  names(between_inits) <- cluster_names

  # Build $macro
  between <- .mcml_layer(
    weights = between_weights,
    inits   = between_inits,
    labels  = cluster_names,
    data    = NULL
  )

  # ============================================================================
  # Within-cluster computation (optional)
  # ============================================================================

  within_data <- NULL
  if (isTRUE(compute_within)) {
    within_data <- lapply(seq_len(n_clusters), function(i) {
      idx_i <- cluster_indices[[i]]
      n_i <- length(idx_i)
      cl_nodes <- cluster_list[[i]]

      if (n_i <= 1) {
        # Single node: preserve self-loop value
        within_raw <- mat[idx_i, idx_i, drop = FALSE]
        dimnames(within_raw) <- list(cl_nodes, cl_nodes)
        within_weights_i <- within_raw
        within_inits_i <- setNames(1, cl_nodes)
      } else {
        # Extract within-cluster raw weights (preserves self-loops)
        within_raw <- mat[idx_i, idx_i]
        dimnames(within_raw) <- list(cl_nodes, cl_nodes)

        within_weights_i <- within_raw

        # Within-cluster inits (handle NAs)
        col_sums_w <- colSums(within_raw, na.rm = TRUE)
        total_w <- sum(col_sums_w, na.rm = TRUE)
        within_inits_i <- if (!is.na(total_w) && total_w > 0) {
          col_sums_w / total_w
        } else {
          rep(1 / n_i, n_i)
        }
        names(within_inits_i) <- cl_nodes
      }

      .mcml_layer(
        weights = within_weights_i,
        inits   = within_inits_i,
        labels  = cl_nodes,
        data    = NULL
      )
    })
    names(within_data) <- cluster_names
  }

  # ============================================================================
  # Build result
  # ============================================================================

  result <- structure(
    list(
      macro = between,
      clusters = within_data,
      cluster_members = cluster_list,
      edges = NULL,    # NULL placeholder; transitions paths fill this in
      meta = list(
        type = "aggregate",
        method = method,
        directed = directed,
        n_nodes = n,
        n_clusters = n_clusters,
        cluster_sizes = vapply(cluster_list, length, integer(1)),
        source = "matrix"
      )
    ),
    class = "mcml"
  )

  result
}

.validate_mcml_matrix <- function(mat) {
  if (any(is.na(mat) | !is.finite(mat))) {
    stop(errorCondition("x matrix must contain finite non-missing weights.",
      class = "nestimate_bad_input", call = NULL))
  }
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  has_row_names <- !is.null(row_names)
  has_col_names <- !is.null(col_names)
  if (has_row_names && (any(is.na(row_names)) || any(!nzchar(row_names)))) {
    stop(errorCondition(
      "x matrix row names must not contain missing or empty values.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (has_col_names && (any(is.na(col_names)) || any(!nzchar(col_names)))) {
    stop(errorCondition(
      "x matrix column names must not contain missing or empty values.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (has_row_names && anyDuplicated(row_names)) {
    stop(errorCondition("x matrix row names must be unique.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (has_col_names && anyDuplicated(col_names)) {
    stop(errorCondition("x matrix column names must be unique.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (has_row_names && has_col_names && !identical(row_names, col_names)) {
    stop(errorCondition(
      "x matrix row and column names must be identical and in the same order.",
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}


# ==============================================================================
# 2a-bis. Sequence preprocessing owned by build_mcml()
#
# These back the `exclude`, `trim`, `end` and `end_by` arguments so callers
# never hand-assemble the equivalent filter/truncate/append ritual. Order is
# fixed and matters: exclude -> trim -> end. Appending the terminal marker
# last is what keeps a trim from silently eating the markers it is supposed
# to preserve.
# ==============================================================================

# Drop unwanted states from long-format event data (one row per event).
.mcml_drop_states <- function(df, action, exclude) {
  keep <- !(as.character(df[[action]]) %in% exclude)
  df[keep, , drop = FALSE]
}

# Blank unwanted states out of a wide sequence frame (one row per sequence).
.mcml_blank_states <- function(df, exclude) {
  df[] <- lapply(df, function(col) {
    col[as.character(col) %in% exclude] <- NA
    col
  })
  df
}

# Truncate a wide sequence frame to its first `cut` time columns.
.mcml_trim_sequences <- function(df, trim) {
  cut <- .trim_cut(as.matrix(df), trim)
  if (cut >= ncol(df)) return(df)
  df[, seq_len(cut), drop = FALSE]
}

# Resolve `end` (FALSE / TRUE / label) to a label or NULL.
.mcml_end_label <- function(end) {
  if (isFALSE(end)) return(NULL)
  if (isTRUE(end)) return("End")
  as.character(end)[1L]
}

# Which sequence rows carry the terminal marker. With `end_by = NULL` every
# sequence is marked (the build_network `end` contract). With `end_by` naming
# coarser column(s), only the sequence holding that group's LAST EVENT is
# marked.
#
# Note it is the last *event*, not the last-starting session: sessions within
# a group can overlap, so the session that begins last is not always the one
# that ends last. Choosing by start position misplaces the marker on exactly
# those groups (47 of 118,501 skills in the Levebee event log).
.mcml_end_rows <- function(meta, long, end_by, n_rows) {
  if (is.null(end_by)) return(seq_len(n_rows))
  missing_cols <- setdiff(end_by, names(meta))
  if (length(missing_cols) > 0L) {
    stop(errorCondition(paste0("`end_by` column(s) not found in the data: ",
         paste(missing_cols, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  # Position of each session's final event in the ORIGINAL event order.
  # `.order` is prepare()'s record of that; row position in `long_data` is
  # not a substitute, because prepare() orders sessions lexicographically by
  # label -- so try "10" sorts before try "9" and the arithmetically later
  # session lands earlier in the frame.
  sid      <- as.character(long$.session_id)
  ord      <- if (!is.null(long$.order)) long$.order else seq_along(sid)
  last_pos <- tapply(ord, sid, max)
  pos      <- as.numeric(last_pos[as.character(meta$.session_id)])

  key <- do.call(paste, c(lapply(end_by, function(k) meta[[k]]), sep = "\r"))
  idx <- seq_len(n_rows)
  vapply(split(idx, factor(key, levels = unique(key))),
         function(g) g[which.max(pos[g])], integer(1L), USE.NAMES = FALSE)
}

# Append the terminal state into the cell after each marked row's last
# non-NA entry, widening the frame by one column when needed.
.mcml_append_end <- function(df, label, rows) {
  mat <- as.matrix(df)
  storage.mode(mat) <- "character"
  filled <- !is.na(mat)
  last <- max.col(filled * col(filled), ties.method = "last")
  last[!apply(filled, 1L, any)] <- 0L

  if (max(last[rows]) >= ncol(mat)) {
    mat <- cbind(mat, NA_character_)
  }
  mat[cbind(rows, last[rows] + 1L)] <- label
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(out) <- paste0("T", seq_len(ncol(out)))
  out
}

# Validate the four preprocessing arguments up front.
.mcml_check_seq_args <- function(exclude, trim, end, end_by) {
  .assert(
    "`exclude` must be a character vector or NULL" =
      is.null(exclude) || is.character(exclude),
    "`trim` must be a single positive number or NULL" =
      is.null(trim) || (is.numeric(trim) && length(trim) == 1L &&
                        !is.na(trim) && trim > 0),
    "`end` must be TRUE, FALSE, or a single label" =
      isTRUE(end) || isFALSE(end) ||
        (length(end) == 1L && !is.na(end) && nzchar(as.character(end))),
    "`end_by` must be a character vector of column names or NULL" =
      is.null(end_by) || is.character(end_by)
  )
  if (!is.null(end_by) && isFALSE(end)) {
    stop(errorCondition(
      paste0(
        "`end_by` groups the terminal marker, so it needs `end` as well ",
        "(e.g. end = \"Conclude\", end_by = \"SkillID\")."),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

# ==============================================================================
# 2b. Build MCML from Raw Transition Data
# ==============================================================================


#' Build MCML from Raw Transition Data
#'
#' @description
#' `build_mcml()` builds a Multi-Cluster Multi-Level (MCML) model from
#' transitions: it recodes every state to its cluster and counts the
#' transitions that occur, which gives the Markov chain over clusters (the
#' macro network) and one network per cluster of the transitions among its
#' own states. It reads sequences (one row per sequence), event logs (with
#' `action`), edge lists (with `from` and `to`), `tna` models and
#' `netobject`s that carry their data, and it re-partitions an existing
#' `mcml` (see `combine` and `expand`). A weight matrix has no transitions to
#' count, so it is passed to [cluster_summary()], which aggregates the
#' weights instead.
#'
#' With the defaults (`type = "tna"`, `method = "sum"`), a macro edge is the
#' probability that the
#' next state lies in cluster j given that the current one lies in cluster i;
#' its diagonal is the probability of staying in a cluster, whether in the
#' same state or a different one. Both levels depend on the partition, which
#' is the analyst's: compare a few partitions with `combine` and `expand`
#' rather than reading one. The weights are point estimates;
#' \code{\link[=as_tna]{as_tna()}}
#' turns the levels into networks that carry their sequences, so
#' [bootstrap_network()] can show how far each edge moves under resampling.
#'
#' @param x The data. Accepted forms:
#'   \describe{
#'     \item{data frame with `from` and `to` set}{An edge list: one row per
#'       transition, with an optional `weight` column.}
#'     \item{data frame with `action` set}{An event log, turned into
#'       sequences by [prepare()].}
#'     \item{data frame of positions}{Sequences, one row per sequence and one
#'       column per time step: either a numbered run of columns (`T1`, `T2`,
#'       ...) and nothing else, or the columns named in `cols`. Each pair of
#'       adjacent columns (t, t + 1) whose two states are both observed is a
#'       transition; a missing cell breaks the sequence.}
#'     \item{`tna` model}{Its stored sequences are decoded and counted; a
#'       model without data is passed to [cluster_summary()].}
#'     \item{`netobject` or bare `cograph_network`}{A network of a counting
#'       method (`relative`, `frequency`, `co_occurrence` and their aliases)
#'       is re-counted from the data it carries. A network of another method
#'       (such as an association or attention network), a one-hot network, or one
#'       without data is passed to [cluster_summary()], which aggregates its
#'       weights.}
#'     \item{`mcml`}{Returned unchanged, unless `clusters`, `combine` or
#'       `expand` changes its partition: it is then re-estimated from the
#'       sequences it carries (see `combine`).}
#'     \item{square numeric matrix}{Passed to [cluster_summary()].}
#'     \item{non-square or character matrix}{Read as sequences.}
#'   }
#'   A data frame that declares no role and is not a numbered run of
#'   positions, or whose columns hold only 0 and 1, raises an error of class
#'   `nestimate_bad_input`. For the `print()` methods: an `mcml` object or
#'   one of its layers (an `mcml_layer`).
#' @param clusters The cluster of each state. Accepted forms:
#'   \describe{
#'     \item{named list}{Cluster names mapped to their states, such as
#'       `list(A = c("N1", "N2"), B = c("N3", "N4"))`.}
#'     \item{data frame}{State names in the first column, cluster names in
#'       the second, such as
#'       `data.frame(node = c("N1", "N2", "N3"), group = c("A", "A", "B"))`.}
#'     \item{membership vector}{One cluster per state, such as
#'       `c("A", "A", "B", "B")`, matched to the states of the data in
#'       alphabetical order (a named vector is matched by name).}
#'     \item{column name}{For an edge list only: a column holding each row's
#'       cluster, given to both its endpoints. It suits edge lists whose
#'       edges never cross clusters; a node that gets two clusters raises an
#'       error, so pass a list or a data frame for any other edge list.}
#'     \item{`NULL`}{(default) Read from a `netobject`'s node table (a
#'       column `clusters`, `cluster`, `groups` or `group`) or from its
#'       `node_groups` (aligned to the node table by label). The values are
#'       then taken in the order of the node table, which must list the
#'       states in the order [build_network()] gives them (alphabetical for
#'       sequence data). Other input needs an explicit assignment and raises
#'       an error without one.}
#'   }
#'   Every state must belong to exactly one cluster, the `end` state
#'   included. On the sequence and edge-list routes the clusters of the
#'   result come in alphabetical order.
#' @param method How the transitions between (or within) two clusters are
#'   combined into one weight, with [net_aggregate_weights()]: `"sum"`
#'   (default), `"mean"`, `"median"`, `"max"`, `"min"`, `"density"` or
#'   `"geomean"`. On sequences and event logs every transition weighs 1, so
#'   only `"sum"` keeps the counts: `"mean"`, `"median"`, `"max"`, `"min"`
#'   and `"geomean"` give 1 for every pair of clusters with a transition, and
#'   `"density"` divides each count by `n_i * n_j`. The other methods suit
#'   weighted edge lists and weight matrices, where each row is already a
#'   measurement rather than one observation. On sequences and event logs the
#'   within-cluster networks likewise get 1 for each observed pair of states
#'   under every method but `"sum"` (`"density"` too, which is not divided by
#'   `n_i * n_j` there); on an edge list they hold the aggregated weights.
#' @param type How the aggregated counts are turned into weights:
#'   \describe{
#'     \item{`"tna"`}{(default) Each row is normalised to sum to 1, giving
#'       first-order Markov transition probabilities.}
#'     \item{`"raw"`}{The counts as they are.}
#'     \item{`"frequency"`}{The same counts as `"raw"`, under the name
#'       frequency networks use.}
#'     \item{`"cooccurrence"`}{The mean of the count matrix and its
#'       transpose (undirected co-occurrence); the result is recorded as
#'       undirected.}
#'   }
#'   `"semi_markov"` raises an error: the package has no semi-Markov
#'   (holding-time) construction. For input passed to [cluster_summary()]
#'   (a weight matrix, a `tna` model or `netobject` without data, or a
#'   `netobject` of a method that is not a count, such as `"attention"` or
#'   `"glasso"`) `type` has no effect; setting it gives a warning for a
#'   matrix, or for a `tna` model or `netobject` without data.
#' @param directed Whether transitions are directed. Default `TRUE`. `FALSE`
#'   replaces the counts by the mean of the count matrix and its transpose
#'   before `type` is applied.
#' @param compute_within Whether to build the within-cluster networks.
#'   Default `TRUE`. An `mcml` built with `FALSE` keeps no per-cluster
#'   sequences, so it cannot be re-partitioned (see `combine`).
#' @param actor,action,time,order,session,time_threshold For an event log:
#'   the columns and the session gap (in seconds, default `900`) passed to
#'   [prepare()], which forms the sequences that are then counted, as
#'   `prepare(...) |> build_network() |> build_mcml()` would. Setting
#'   `action` declares the data an event log.
#' @param exclude States to drop before anything is counted, such as a
#'   technical marker. In an event log the matching events are removed
#'   before sequences are formed, so a dropped state never takes a position;
#'   in wide sequences the matching cells become missing. Default `NULL`.
#'   Only for wide sequences and event logs: an edge list, a `netobject` or
#'   a `tna` model raises an error. The same holds for `trim` and `end`. On
#'   an existing `mcml` they raise an error together with a re-partition
#'   (`clusters`, `combine` or `expand`); without one the `mcml` is returned
#'   unchanged and they have no effect (see `combine`).
#' @param trim How much of each sequence to keep. `NULL` (default) keeps
#'   every time step. A fraction in (0, 1) keeps the columns up to that
#'   quantile of the sequence lengths (the last observed position of each
#'   row), so `trim = 0.95` cuts only the longest 5% of sequences; a value of
#'   1 or more keeps that many leading columns (`trim = 10` keeps the first
#'   10). It changes the data, so it only ever cuts the end of a sequence;
#'   the fraction in [sequence_plot()]'s `trim` works differently because it
#'   changes only the drawing.
#' @param end A terminal state put in the cell after each sequence's last
#'   observed state. `FALSE` (default) adds none, `TRUE` adds one labelled
#'   `"End"`, and a single string adds one with that label; assign it to a
#'   cluster in `clusters`. It is added after `trim`, so trimming never
#'   removes it.
#' @param end_by Columns grouping the terminal marker at a coarser unit than
#'   the sequence: only the sequence holding each group's last event is
#'   marked. For example,
#'   `session = "AttemptID", end = "Conclude", end_by = "SkillID"`
#'   closes each skill once, not each attempt. It needs
#'   `end` and an event log, and raises an error otherwise. Default `NULL`
#'   (every sequence is marked).
#' @param from,to,weight The columns of an edge list: its source, its target
#'   and an optional non-negative weight (each transition counts 1 without
#'   it). Names give no role: a frame is read as an edge list only when
#'   `from` and `to` are set, and setting them with `action` or `cols` raises
#'   an error of class `nestimate_bad_input`.
#' @param cols For a wide data frame: the columns holding the positions,
#'   needed when the frame is anything but a numbered run (`T1`, `T2`, ...).
#'   Default `NULL`.
#' @param labels New display labels for the states of the within-cluster
#'   networks, applied after estimation: a two-column data frame
#'   (`name`, `label`), a named character vector `c(name = "label")` or a
#'   named list. States not listed keep their names. The new labels also
#'   replace the state names in the cluster membership and the edge table;
#'   the macro network keeps its cluster names.
#' @param combine,expand Change the partition. `combine` merges clusters into
#'   one: a character vector merges one group, and a list merges several,
#'   labelled by its names (default label `"A + B"`); the merged cluster
#'   lists its states in cluster-name order. `expand` then splits the named
#'   clusters (or every cluster, with `"all"` or `TRUE`) into one cluster per
#'   state, named by the state. On new input they apply to `clusters` (in any
#'   of its forms, auto-detection included) before estimation, so
#'   `build_mcml(data, clusters = cl, combine = c("A", "B"))` equals a build
#'   with `A` and `B` merged in `cl`; this works for every input, matrices
#'   included. On an existing `mcml` they re-estimate the model (macro
#'   network, within-cluster networks and stored sequences) from the
#'   sequences it carries, with its original `type`, `method` and `directed`
#'   unless these are passed; a new `clusters` does the same under that
#'   partition, and an unchanged partition returns an equal object. Any
#'   session split, `exclude`, `trim` or `end` applied at build time is
#'   already in those sequences, so passing one of them (or `actor`,
#'   `action`, `time`, `order`, `session`, `time_threshold`, `end_by`,
#'   `labels`, `from`, `to`, `weight`, `cols`) with a re-partition raises an
#'   error, as does passing `clusters` together with `combine` or `expand`.
#'   An `mcml` built from a matrix, from an edge list (it keeps only the
#'   within-cluster edges) or with `compute_within = FALSE` raises an error
#'   of class `nestimate_mcml_no_sequences`. [sequence_plot()] takes the
#'   same two arguments as display options: its `combine` draws what it
#'   draws for `build_mcml(x, combine = )`, while its `expand` opens
#'   clusters in the Summary panel only and keeps one panel per cluster.
#' @param ... For the `print()` and `summary()` methods: not used; any extra
#'   argument raises an error.
#' @param object For `summary()`: an `mcml` object.
#'
#' @return An `mcml` object with the layout described under
#'   [cluster_summary()]. Print it for the type, the method, the node and
#'   cluster counts, the clusters,
#'   the number of transitions between and within clusters and the macro
#'   weights; printing one layer (an `mcml_layer`) shows its node and edge
#'   counts, weight range, weight matrix, initial probabilities and the size
#'   of the data it carries. [summary()] returns a data frame with one row
#'   per cluster, with
#'   columns `cluster`, `size`, `within_total` (the total absolute weight of
#'   the within-cluster network, `NA` when it was not built), `between_out`
#'   and `between_in` (the absolute macro weight leaving and entering the
#'   cluster, diagonal excluded); for an undirected macro network, `between_out` is
#'   the total incident weight and `between_in` is `NA`.
#'   \code{\link[=as_tna]{as_tna()}} turns
#'   the levels into a `netobject_group`, [as_htna()] rebuilds the
#'   node-level network with its clusters, [macro_network()] opens chosen
#'   clusters into their states, and [sequence_plot()] draws the sequences
#'   by cluster. The `print()` methods return their input invisibly.
#'
#'   On the sequence and edge-list routes, the object differs from the
#'   matrix route's in these fields:
#'   \describe{
#'     \item{edges}{A data frame with one row per observed node-level
#'       transition and columns `from`, `to`, `weight`, `cluster_from`,
#'       `cluster_to` and `type` (`"within"` or `"between"`).}
#'     \item{macro}{Its `inits` are the shares of the clusters among the
#'       first-column states (missing cells left out), or the normalised
#'       column sums for an edge list; its `data` holds the sequences recoded to
#'       cluster labels, or one `from`/`to` row per transition for an edge
#'       list.}
#'     \item{clusters}{Each layer counts only the transitions between states
#'       of its own cluster. Its `inits` are the shares of its states among
#'       the first-column states that belong to it (normalised column sums
#'       for an edge list), and its `data` the sequences with the states of other
#'       clusters set missing (for an edge list, its within-cluster
#'       transitions).}
#'     \item{meta}{`source` is `"transitions"`, `type` records the `type`
#'       used, and `directed` whether the stored weights are directed
#'       (`FALSE` for `type = "cooccurrence"`).}
#'   }
#'   Input passed to [cluster_summary()] gives its result: `source`
#'   `"matrix"` in `meta`, and no edge table.
#'
#' @order 1
#' @export
#' @seealso \code{\link{cluster_summary}} for matrix-based aggregation,
#'   \code{\link{as_tna}} to promote the layers to netobjects,
#'   \code{\link{macro_network}} for the cluster-level network with one
#'   cluster expanded
#'
#' @examples
#' # Edge list with clusters: its columns are declared
#' edges <- data.frame(
#'   source = c("A", "A", "B", "C", "C", "D"),
#'   target = c("B", "C", "A", "D", "D", "A"),
#'   count  = c(1, 2, 1, 3, 1, 2)
#' )
#' clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
#' build_mcml(edges, clusters, from = "source", to = "target", weight = "count")
#'
#' # Sequence data with clusters
#' seqs <- data.frame(
#'   T1 = c("A", "C", "B"),
#'   T2 = c("B", "D", "A"),
#'   T3 = c("C", "C", "D"),
#'   T4 = c("D", "A", "C")
#' )
#' cs <- build_mcml(seqs, clusters, type = "raw")
#' cs
#' summary(cs)
#'
#' # Change the partition while building ...
#' three <- build_mcml(seqs, list(G1 = "A", G2 = "B", G3 = c("C", "D")))
#' build_mcml(seqs, list(G1 = "A", G2 = "B", G3 = c("C", "D")),
#'            combine = c("G1", "G2"))
#'
#' # ... or re-partition an existing mcml: the model is re-estimated
#' build_mcml(three, combine = c("G1", "G2"))           # G1 + G2 as one cluster
#' build_mcml(three, combine = list(AB = c("G1", "G2"))) # named merge
#' build_mcml(three, expand = "G3")                      # C and D as clusters
build_mcml <- function(x,
                       clusters = NULL,
                       method = c("sum", "mean", "median", "max",
                                  "min", "density", "geomean"),
                       type = c("tna", "frequency", "cooccurrence",
                                "raw"),
                       directed = TRUE,
                       compute_within = TRUE,
                       actor = NULL,
                       action = NULL,
                       time = NULL,
                       order = NULL,
                       session = NULL,
                       time_threshold = 900,
                       exclude = NULL,
                       trim = NULL,
                       end = FALSE,
                       end_by = NULL,
                       labels = NULL,
                       combine = NULL,
                       expand = NULL,
                       from = NULL,
                       to = NULL,
                       weight = NULL,
                       cols = NULL) {
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop(errorCondition("'directed' must be TRUE or FALSE.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.logical(compute_within) || length(compute_within) != 1L ||
      is.na(compute_within)) {
    stop(errorCondition("'compute_within' must be TRUE or FALSE.",
      class = "nestimate_bad_input", call = NULL))
  }

  # An existing mcml: unchanged unless the partition is changed, in which
  # case it is re-estimated from the sequences it carries.
  if (inherits(x, "mcml")) {
    if (is.null(clusters) && is.null(combine) && is.null(expand)) {
      return(x)
    }
    # Everything that shaped the sequences was applied when x was built and
    # is baked into the sequences it carries; accepting these again would
    # silently do nothing.
    shaping <- c(actor = !missing(actor), action = !missing(action),
                 time = !missing(time), order = !missing(order),
                 session = !missing(session),
                 time_threshold = !missing(time_threshold),
                 exclude = !missing(exclude), trim = !missing(trim),
                 end = !missing(end), end_by = !missing(end_by),
                 labels = !missing(labels), from = !missing(from),
                 to = !missing(to), weight = !missing(weight),
                 cols = !missing(cols))
    if (any(shaping)) {
      stop(errorCondition(
        paste0(
          "`",
          paste(names(shaping)[shaping], collapse = "`, `"),
          "` cannot be applied when re-partitioning an mcml: the sequences ",
          "it carries already reflect how it was built. Build from the ",
          "original data to change them."),
        class = "nestimate_bad_input", call = NULL))
    }
    return(.mcml_rebuild(
      x, clusters = clusters, combine = combine, expand = expand,
      method = if (missing(method)) x$meta$method else match.arg(method),
      type = if (missing(type)) x$meta$type else match.arg(type),
      directed = if (missing(directed)) x$meta$directed %||% directed else directed,
      compute_within = compute_within))
  }
  # combine/expand on fresh input: resolve the partition exactly as the plain
  # build does (every `clusters` form, including auto-detection, ends up as
  # the fit's named `cluster_members`), re-partition that list, and build
  # once more from the same input. Identical to building with the merged
  # list by construction, and works for matrix input too.
  if (!is.null(combine) || !is.null(expand)) {
    args <- list(x = x, clusters = clusters, directed = directed,
                 compute_within = compute_within, actor = actor,
                 action = action, time = time, order = order,
                 session = session, time_threshold = time_threshold,
                 exclude = exclude, trim = trim, end = end, end_by = end_by,
                 from = from, to = to, weight = weight, cols = cols)
    # The first build only reads the partition: without `labels` (its
    # cluster_members must hold the raw state names the second build maps)
    # and without method/type (which do not change the partition, and whose
    # matrix-path warning must fire once, in the real build).
    base <- do.call(build_mcml, args)
    args$clusters <- .mcml_partition(base$cluster_members, combine, expand)
    args$labels   <- labels
    if (!missing(method)) args$method <- method
    if (!missing(type))   args$type   <- type
    return(do.call(build_mcml, args))
  }

  # Remember whether the caller passed `type` explicitly -- used below to
  # warn when type is applied to matrix input, where it has no effect.
  type_explicit <- !missing(type)

  # No semi-Markov / holding-time construction exists in the package. Reject
  # "semi_markov" explicitly instead of silently aliasing "tna" (which the
  # shared row-normalisation branch used to do while print.mcml claimed a
  # semi-Markov model).
  if (is.character(type) && length(type) == 1L &&
      !is.na(type) && type == "semi_markov") {
    stop(errorCondition(
      paste0(
        "type = \"semi_markov\" is not implemented: Nestimate has no ",
        "semi-Markov / holding-time construction. Use type = \"tna\" for a ",
        "first-order Markov transition matrix, or \"raw\"/\"cooccurrence\"."),
      class = "nestimate_bad_input", call = NULL))
  }

  type <- match.arg(type)
  method <- match.arg(method)

  .mcml_check_seq_args(exclude, trim, end, end_by)
  end_label <- .mcml_end_label(end)
  is_edge_list <- .mcml_check_roles(x, action, actor, time, order, session, from,
                                    to, weight, cols)
  if (is.data.frame(x) && !is.null(cols)) x <- .mark_positions(x[cols])

  # Long-format event-log shortcut: if `action` is supplied on a data.frame,
  # delegate to prepare() and feed the resulting wide sequence into the
  # existing sequence path. Behaves identically to
  # prepare(...) |> build_network() |> build_mcml().
  #
  # `exclude` runs before prepare() so dropped states never occupy a sequence
  # position; `trim` and `end` run after, on prepare()'s sessions, because
  # prepare() is what decides where a session begins and ends.
  if (is.data.frame(x) && !is.null(action)) {
    if (!is.null(exclude)) x <- .mcml_drop_states(x, action, exclude)
    prepared <- prepare(x, actor = actor, action = action, time = time,
                        order = order, session = session,
                        time_threshold = time_threshold)
    seq_df <- prepared$sequence_data
    if (!is.null(trim)) seq_df <- .mcml_trim_sequences(seq_df, trim)
    if (!is.null(end_label)) {
      seq_df <- .mcml_append_end(
        seq_df, end_label,
        .mcml_end_rows(prepared$meta_data, prepared$long_data, end_by,
                       nrow(seq_df)))
    }
    x <- seq_df
  } else if (!is.null(exclude) || !is.null(trim) || !is.null(end_label)) {
    # Wide sequence input supports exclude/trim/end; end_by needs the
    # long-format columns prepare() carries, so it has nothing to group by.
    if (!is.null(end_by)) {
      stop(errorCondition(
        paste0(
          "`end_by` needs long-format input (supply `action`), because the ",
          "grouping column it names exists only there."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (!is.data.frame(x) && !is.matrix(x)) {
      stop(errorCondition(
        paste0(
          "`exclude`, `trim` and `end` need sequence data (a wide ",
          "data.frame/matrix, or long-format input via `action`); got '",
          class(x)[1L],
          "'."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (is_edge_list) {
      stop(errorCondition(
        paste0(
          "`exclude`, `trim` and `end` apply to sequence data, not to an ",
          "edge list."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x)) {
      stop(errorCondition(
        paste0(
          "`exclude`, `trim` and `end` apply to sequence data, not to a ",
          "weight matrix: there are no sequences to drop states from."),
        class = "nestimate_bad_input", call = NULL))
    }
    seq_df <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!is.null(exclude)) seq_df <- .mcml_blank_states(seq_df, exclude)
    if (!is.null(trim))    seq_df <- .mcml_trim_sequences(seq_df, trim)
    if (!is.null(end_label)) {
      seq_df <- .mcml_append_end(seq_df, end_label, seq_len(nrow(seq_df)))
    }
    x <- seq_df
  }

  # Coerce cograph_network so downstream branches see a netobject
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)

  input_type <- if (is_edge_list) "edgelist" else .detect_mcml_input(x)

  is_matrix_path <- input_type %in% c("matrix", "tna_matrix",
                                       "netobject_matrix")
  if (type_explicit && is_matrix_path) {
    warning(
      "`type` is ignored for matrix input: the matrix path of build_mcml() ",
      "is aggregation-only. Pass `method` to choose the aggregation and ",
      "normalise downstream if a Markov chain is needed (e.g. via as_tna()).",
      call. = FALSE
    )
  }

  result <- switch(input_type,
"edgelist" = .build_mcml_edgelist(x, clusters, method, type,
                                       directed, compute_within,
                                       from = from, to = to, weight = weight),
    "sequence" = .build_mcml_sequence(x, clusters, method, type,
                                       directed, compute_within),
    # tna stores sequences integer-coded; decode to state labels first.
    "tna_data" = .build_mcml_sequence(.mmm_decode_int_data(x$data, x$labels),
                                       clusters, method, type,
                                       directed, compute_within),
    "tna_matrix" = cluster_summary(x, clusters, method = method,
                                    directed = directed,
                                    compute_within = compute_within),
    "netobject_data" = {
      # Auto-detect clusters from network if not provided
      if (is.null(clusters)) {
        clusters <- .auto_detect_clusters(x) # nocov
      }
      # Respect the netobject's method. Count-based methods (relative,
      # frequency, co_occurrence) are well-defined to re-derive from raw
      # data, so we route through the data path for an exact conditional
      # Markov chain. Model-derived methods (attention, ising, glasso,
      # mgm, pcor, cor, mlvar_*) carry weights that re-counting would
      # silently discard, so we route through the matrix path on $weights.
      net_method <- x$method %||% "relative"
      count_methods <- c("relative", "frequency", "co_occurrence",
                         "tna", "ftna", "cna", "wcna", "wtna",
                         "wtna_cooccurrence", "wtna_transition")
      # a one-hot network's data are 0 and 1 per code, not positions of states
      if (net_method %in% count_methods && !.is_onehot_params(x$params)) {
        data <- x$data
        sub_type <- .detect_mcml_input(data)
        if (sub_type == "edgelist") {
          .build_mcml_edgelist(data, clusters, method, type,
                                directed, compute_within)
        } else {
          .build_mcml_sequence(data, clusters, method, type,
                                directed, compute_within)
        }
      } else {
        cluster_summary(x, clusters, method = method,
                        directed = directed,
                        compute_within = compute_within)
      }
    },
    "netobject_matrix" = {
      if (is.null(clusters)) {
        clusters <- .auto_detect_clusters(x) # nocov
      }
      cluster_summary(x, clusters, method = method,
                       directed = directed, compute_within = compute_within)
    },
    "matrix" = cluster_summary(x, clusters, method = method,
                                directed = directed,
                                compute_within = compute_within),
    stop(errorCondition(
      paste0("Cannot build MCML from input of class '", class(x)[1], "'"),
      class = "nestimate_bad_input", call = NULL))
  )

  if (!is.null(labels)) {
    result <- .apply_node_labels(result, labels)
  }
  # Stash the node-level source so as_htna(mcml) can expand back to the node
  # level without the caller re-passing data. `x` here is already in
  # build_network-ready wide form (long input was widened by prepare()
  # above), so no actor/action/time roles are needed to rebuild. Restricted
  # to sequence input: matrix inputs carry no node-level data, and the
  # edgelist -> build_network(method = "relative") round-trip is unreliable,
  # so those fall back to the explicit "pass data" path in as_htna.mcml().
  if (identical(input_type, "sequence")) {
    attr(result, "htna_source") <- x
  }
  result
}

# Re-partition a named cluster list: `combine` merges clusters (a merged
# group is labelled by its list name or "A + B"; its states keep the order of
# the clusters it merges), then `expand` splits clusters into one singleton
# cluster per member state, named by the state. build_mcml() orders clusters
# by name, so list position is not kept. Shared by build_mcml() on fresh input
# and on an existing mcml.
.mcml_partition <- function(clusters, combine = NULL, expand = NULL) {
  groups <- .mcml_resolve_combine(combine, names(clusters))
  if (length(groups) > 0L) {
    owner <- stats::setNames(names(clusters), names(clusters))
    owner[unlist(groups, use.names = FALSE)] <- rep(names(groups), lengths(groups))
    new_names <- unique(unname(owner))
    clusters <- stats::setNames(lapply(new_names, function(g) {
      unlist(clusters[names(owner)[owner == g]], use.names = FALSE)
    }), new_names)
  }
  opened <- .mcml_resolve_expand(expand, names(clusters))
  if (length(opened) > 0L) {
    kept   <- setdiff(names(clusters), opened)
    states <- unlist(clusters[opened], use.names = FALSE)
    clash  <- intersect(states, kept)
    if (length(clash) > 0L) {
      stop(errorCondition(
        paste0(
          "Expanding would give a state the name of an existing cluster: ",
          paste(clash, collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
    clusters <- do.call(c, lapply(names(clusters), function(k) {
      if (k %in% opened) stats::setNames(as.list(clusters[[k]]), clusters[[k]]) else clusters[k]
    }))
  }
  clusters
}

# Re-estimate an mcml under a changed partition. The stored per-cluster
# matrices are the full sequences with other clusters' states blanked, and the
# clusters partition the states, so coalescing them recovers the sequences the
# mcml was built from (after any session split, exclude, trim or end, which are
# therefore not re-applied). `combine` merges clusters, then `expand` splits
# clusters into singleton clusters, one per member state -- the same order as
# sequence_plot(combine =, expand =).
.mcml_rebuild <- function(x, clusters, combine, expand, method, type,
                          directed, compute_within) {
  datas <- lapply(x$clusters, function(z) z$data)
  if (!length(datas) || any(vapply(datas, is.null, logical(1L)))) {
    stop(errorCondition(
      paste0("This mcml carries no sequences (built from a matrix, or with ",
             "compute_within = FALSE), so it cannot be re-partitioned. ",
             "Rebuild it from the original data with the new `clusters`."),
      class = "nestimate_mcml_no_sequences", call = NULL))
  }
  if (any(vapply(datas, function(d) identical(attr(d, "source"), "edgelist"),
                 logical(1L)))) {
    stop(errorCondition(
      paste0("This mcml was built from an edge list: it stores only the ",
             "within-cluster edges, so the between-cluster edges needed to ",
             "re-partition it are gone. Rebuild from the original edge list, ",
             "e.g. build_mcml(edges, clusters = ..., combine = ...)."),
      class = "nestimate_mcml_no_sequences", call = NULL))
  }
  if (!is.null(clusters) && (!is.null(combine) || !is.null(expand))) {
    stop(errorCondition(
      "Pass either `clusters` or `combine`/`expand`, not both.",
      class = "nestimate_bad_input", call = NULL))
  }
  seqs <- Reduce(function(a, b) { a[is.na(a)] <- b[is.na(a)]; a },
                 lapply(datas, as.data.frame, stringsAsFactors = FALSE))

  if (is.null(clusters)) {
    clusters <- .mcml_partition(x$cluster_members, combine, expand)
  }
  out <- build_mcml(seqs, clusters = clusters, method = method, type = type,
                    directed = directed, compute_within = compute_within,
                    cols = names(seqs))
  # The node-level source (used by as_htna) does not depend on the partition;
  # keep the original's, not the coalesced (possibly relabelled) sequences.
  attr(out, "htna_source") <- attr(x, "htna_source")
  out
}

#' Detect input type for build_mcml
#' @noRd
.detect_mcml_input <- function(x) {
  if (inherits(x, "tna")) {
    if (!is.null(x$data)) return("tna_data")
    return("tna_matrix")
  }

  if (inherits(x, "netobject") || inherits(x, "netobject_ml")) {
    if (!is.null(x$data)) return("netobject_data")
    return("netobject_matrix")
  }

  # Names give no role: an edge list is declared, and stored ones are tagged.
  if (is.data.frame(x)) {
    return(if (identical(attr(x, "source"), "edgelist")) "edgelist" else "sequence")
  }

  if (is.matrix(x)) {
    if (is.numeric(x) && nrow(x) == ncol(x)) return("matrix")
    return("sequence")
  }

  "unknown"
}

#' Validate the declared roles of build_mcml() and refuse an undeclared frame
#'
#' Names give no role. An edge list is declared with `from` and `to`, wide data
#' with `cols` or as nothing but a numbered run of positions, an event log with
#' `action`. Reading any other frame as positions would turn its ids, weights and
#' covariates into states, and reading positions as endpoints would drop columns.
#' @return `TRUE` when `x` is a declared edge list.
#' @noRd
.mcml_check_roles <- function(x, action, actor, time, order, session, from, to,
                              weight, cols) {
  edge_args <- c(from = !is.null(from), to = !is.null(to),
                 weight = !is.null(weight))
  if (!is.data.frame(x)) {
    if (any(edge_args) || !is.null(cols)) {
      stop(errorCondition(
        "`from`, `to`, `weight` and `cols` apply to a data frame.",
        class = "nestimate_bad_input", call = NULL))
    }
    return(FALSE)
  }
  .check_declared_columns(x, "build_mcml()", action = action, actor = actor,
                          time = time, order = order, session = session,
                          from = from, to = to, weight = weight, cols = cols)
  is_edge_list <- any(edge_args)
  if (is_edge_list) {
    single <- all(lengths(list(from, to)) == 1L) && length(weight) <= 1L
    if (!all(edge_args[c("from", "to")]) || !single) {
      stop(errorCondition(
        "An edge list needs one `from` and one `to` column, and at most one `weight`.",
        class = "nestimate_bad_input", call = NULL))
    }
    if (!is.null(action) || !is.null(cols)) {
      stop(errorCondition(
        "Declare an edge list (`from`, `to`), an event log (`action`) or positions (`cols`), not several.",
        class = "nestimate_bad_input", call = NULL))
    }
    return(TRUE)
  }
  if (!is.null(action) || !is.null(cols) || !is.null(.marked_positions(x))) {
    return(FALSE)
  }
  # columns of 0 and 1 could be one-hot codes or positions: no reading is guessed
  if (all(vapply(x, .is_binary_column, logical(1)))) {
    stop(errorCondition(
      sprintf(paste0("build_mcml(): these columns hold only 0 and 1, so they could ",
                     "be one-hot codes or positions: %s. Set `cols` if they are ",
                     "positions; for one-hot data pass ",
                     "build_network(format = \"onehot\")."),
              paste(sprintf("`%s`", utils::head(names(x), 6L)), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  positions <- c(.categorical_positions(x), .numeric_position_run(x))
  if (!setequal(positions, names(x))) {
    stop(errorCondition(
      sprintf(paste0("build_mcml(): no column role is declared, and these columns ",
                     "are not all a numbered run of positions (T1, T2, ...): %s. ",
                     "Set `from` and `to` for an edge list, `action` for an event ",
                     "log, or `cols` for wide sequences."),
              paste(sprintf("`%s`", utils::head(names(x), 6L)), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  FALSE
}

#' Auto-detect clusters from netobject
#' @noRd
.auto_detect_clusters <- function(x) {
  # Prefer x$nodes-embedded cluster column. x$nodes is the canonical node
  # table - if it carries a cluster column, the row order is by definition
  # aligned with the weight matrix, so positional use is safe.
  clusters <- NULL
  if (!is.null(x$nodes)) {
    cluster_cols <- c("clusters", "cluster", "groups", "group")
    for (col in cluster_cols) {
      if (col %in% names(x$nodes)) {
        clusters <- x$nodes[[col]]
        .validate_auto_clusters(clusters, "x$nodes")
        break
      }
    }
  }

  # Fall back to x$node_groups, but require explicit label alignment.
  # audit_mcml #1: node_groups was previously read positionally, which
  # silently mis-assigned nodes whenever the rows were in a different
  # order than x$nodes (a common case for externally constructed
  # netobjects). The accepted shapes are now:
  #   * data.frame with a node identifier column AND a cluster column
  #   * named character/factor vector (names = node labels, values = cluster)
  if (is.null(clusters) && !is.null(x$node_groups)) {
    ng <- x$node_groups

    target_labels <- if (!is.null(x$nodes) && "label" %in% names(x$nodes)) {
      as.character(x$nodes$label)
    } else if (!is.null(rownames(x$weights))) {
      rownames(x$weights)
    } else {
      stop(errorCondition(paste0("Cannot align node_groups: x$nodes$label and ",
           "rownames(x$weights) are both unavailable."),
        class = "nestimate_bad_input", call = NULL))
    }

    if (is.data.frame(ng)) {
      cluster_col <- intersect(c("cluster", "group", "layer"), names(ng))
      if (length(cluster_col) == 0L) {
        stop(errorCondition(
          paste0(
            "'node_groups' data.frame is missing a recognised cluster ",
            "column (one of: cluster, group, layer)."),
          class = "nestimate_bad_input", call = NULL))
      }
      node_col <- intersect(c("node", "name", "label", "id"), names(ng))
      if (length(node_col) == 0L) {
        stop(errorCondition(
          paste0(
            "'node_groups' must include a node identifier column ",
            "(one of: node, name, label, id) so cluster assignments can ",
            "be aligned with x$nodes by label rather than by row order. ",
            "Add the column or pass clusters explicitly as a named list ",
            "or vector."),
          class = "nestimate_bad_input", call = NULL))
      }
      node_values <- as.character(ng[[node_col[1L]]])
      cluster_values <- as.character(ng[[cluster_col[1L]]])
      if (any(is.na(node_values) | !nzchar(node_values))) {
        stop(errorCondition(
          "'node_groups' node column must not contain missing or empty values.",
          class = "nestimate_bad_input", call = NULL))
      }
      .validate_auto_clusters(cluster_values, "node_groups")
      duplicated_nodes <- unique(node_values[duplicated(node_values)])
      if (length(duplicated_nodes) > 0L) {
        stop(errorCondition(
          paste0(
            "'node_groups' assigns duplicate rows to node(s): ",
            paste(utils::head(duplicated_nodes, 5L), collapse = ", ")),
          class = "nestimate_bad_input", call = NULL))
      }
      lookup <- setNames(cluster_values, node_values)
    } else if (is.atomic(ng) && !is.null(names(ng))) {
      node_values <- names(ng)
      cluster_values <- as.character(ng)
      if (any(is.na(node_values) | !nzchar(node_values))) {
        stop(errorCondition(
          "'node_groups' names must not contain missing or empty values.",
          class = "nestimate_bad_input", call = NULL))
      }
      .validate_auto_clusters(cluster_values, "node_groups")
      duplicated_nodes <- unique(node_values[duplicated(node_values)])
      if (length(duplicated_nodes) > 0L) {
        stop(errorCondition(
          paste0(
            "'node_groups' assigns duplicate rows to node(s): ",
            paste(utils::head(duplicated_nodes, 5L), collapse = ", ")),
          class = "nestimate_bad_input", call = NULL))
      }
      lookup <- setNames(cluster_values, node_values)
    } else {
      stop(errorCondition(
        paste0(
          "'node_groups' must be a data.frame with node + cluster ",
          "columns, or a named atomic vector keyed by node label. ",
          "Unnamed vectors cannot be safely aligned to x$nodes."),
        class = "nestimate_bad_input", call = NULL))
    }

    missing_nodes <- setdiff(target_labels, names(lookup))
    if (length(missing_nodes) > 0L) {
      stop(errorCondition(sprintf(
        paste0("node_groups is missing cluster assignments for %d node(s) ",
               "present in x$nodes/weights: %s"),
        length(missing_nodes),
        paste(utils::head(missing_nodes, 5L), collapse = ", ")
      ),
        class = "nestimate_bad_input", call = NULL))
    }
    clusters <- unname(lookup[target_labels])
  }

  if (is.null(clusters)) {
    stop(errorCondition(paste0("No clusters found in netobject. ",
         "Add a 'clusters' column to nodes or provide clusters argument."),
      class = "nestimate_bad_input", call = NULL))
  }
  clusters
}

.validate_auto_clusters <- function(clusters, source) {
  clusters <- as.character(clusters)
  if (any(is.na(clusters) | !nzchar(clusters))) {
    stop(errorCondition(
      paste0(
        source,
        " cluster assignments must not contain missing or empty values."),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

#' Build node-to-cluster lookup from cluster specification
#' @noRd
.build_cluster_lookup <- function(clusters, all_nodes) {
  if (is.data.frame(clusters)) { # nocov start
    # Defensive: .normalize_clusters converts df to list before this is called
    stopifnot(ncol(clusters) >= 2)
    nodes <- as.character(clusters[[1]])
    groups <- as.character(clusters[[2]])
    clusters <- split(nodes, groups)
  } # nocov end

  if (is.list(clusters) && !is.data.frame(clusters)) {
    # Named list: cluster_name -> node vector
    .validate_cluster_partition(clusters, all_nodes)
    lookup <- character(0)
    for (cl_name in names(clusters)) {
      nodes <- clusters[[cl_name]]
      lookup[nodes] <- cl_name
    }
    # Verify all nodes are mapped
    unmapped <- setdiff(all_nodes, names(lookup))
    if (length(unmapped) > 0) {
      stop("Unmapped nodes: ",
           paste(utils::head(unmapped, 5), collapse = ", "),
           call. = FALSE)
    }
    return(lookup)
  }

  if (is.character(clusters) || is.factor(clusters)) {
    clusters <- as.character(clusters)
    if (length(clusters) != length(all_nodes)) {
      stop("Membership vector length (", length(clusters),
           ") must equal number of unique nodes (", length(all_nodes), ")",
           call. = FALSE)
    }
    lookup <- setNames(clusters, all_nodes)
    return(lookup)
  }

  if (is.numeric(clusters) || is.integer(clusters)) {
    if (length(clusters) != length(all_nodes)) {
      stop("Membership vector length (", length(clusters),
           ") must equal number of unique nodes (", length(all_nodes), ")",
           call. = FALSE)
    }
    lookup <- setNames(as.character(clusters), all_nodes)
    return(lookup)
  }

  stop("clusters must be a named list, character/numeric vector, or column name",
       call. = FALSE)
}

#' Empirical first-state distribution from a wide sequence data.frame.
#' Matches `tna::tna()`'s inits: tabulate the values of column 1 (the first
#' time-step), dropping NAs, optionally restricted to a `restrict_to`
#' whitelist (anything outside the whitelist is treated as NA so it is
#' dropped). Normalised to sum to 1. Returns a zero-vector with
#' `length(levels)` entries when no rows contribute.
#' @keywords internal
#' @noRd
.first_state_distribution <- function(wide_df, levels, restrict_to = NULL) {
  empty <- setNames(rep(0, length(levels)), levels)
  if (!is.data.frame(wide_df) || nrow(wide_df) == 0L || ncol(wide_df) == 0L) {
    return(empty)
  }
  first <- as.character(wide_df[[1]])
  if (!is.null(restrict_to)) first[!first %in% restrict_to] <- NA_character_
  first <- first[!is.na(first)]
  if (length(first) == 0L) return(empty)
  fr <- table(factor(first, levels = levels))
  out <- as.numeric(fr) / sum(fr)
  names(out) <- levels
  out
}

#' Build cluster_summary from transition vectors
#' @noRd
.build_from_transitions <- function(from_nodes, to_nodes, weights,
                                     cluster_lookup, cluster_list,
                                     method, type, directed,
                                     compute_within, data = NULL) {
  if (!is.character(from_nodes) || !is.character(to_nodes)) {
    stop("'from_nodes' and 'to_nodes' must be character vectors.",
         call. = FALSE)
  }
  if (length(from_nodes) != length(to_nodes) ||
      length(from_nodes) != length(weights)) {
    stop("'from_nodes', 'to_nodes', and 'weights' must have the same length.",
         call. = FALSE)
  }
  if (any(is.na(from_nodes) | !nzchar(from_nodes) |
          is.na(to_nodes) | !nzchar(to_nodes))) {
    stop("'from_nodes' and 'to_nodes' must not contain missing or empty values.",
         call. = FALSE)
  }
  if (!is.numeric(weights) || any(is.na(weights) | !is.finite(weights))) {
    stop("'weights' must be a finite non-missing numeric vector.",
         call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("'weights' must not contain negative values.", call. = FALSE)
  }
  if (!is.list(cluster_list) || is.data.frame(cluster_list)) {
    stop("'cluster_list' must be a named list.", call. = FALSE)
  }
  cluster_nodes <- sort(unique(unlist(cluster_list, use.names = FALSE)))
  .validate_cluster_partition(cluster_list, cluster_nodes)
  if (!is.character(cluster_lookup) || is.null(names(cluster_lookup))) {
    stop("'cluster_lookup' must be a named character vector.", call. = FALSE)
  }
  if (any(is.na(names(cluster_lookup)) | !nzchar(names(cluster_lookup))) ||
      any(is.na(cluster_lookup) | !nzchar(cluster_lookup))) {
    stop("'cluster_lookup' names and values must not be missing or empty.",
         call. = FALSE)
  }
  if (anyDuplicated(names(cluster_lookup))) {
    stop("'cluster_lookup' names must be unique.", call. = FALSE)
  }
  missing_lookup <- setdiff(unique(c(from_nodes, to_nodes)),
                            names(cluster_lookup))
  if (length(missing_lookup) > 0L) {
    stop("'cluster_lookup' is missing node(s): ",
         paste(utils::head(missing_lookup, 5L), collapse = ", "),
         call. = FALSE)
  }
  unknown_clusters <- setdiff(unique(unname(cluster_lookup)),
                              names(cluster_list))
  if (length(unknown_clusters) > 0L) {
    stop("'cluster_lookup' contains unknown cluster(s): ",
         paste(utils::head(unknown_clusters, 5L), collapse = ", "),
         call. = FALSE)
  }
  if (!is.logical(compute_within) || length(compute_within) != 1L ||
      is.na(compute_within)) {
    stop("'compute_within' must be TRUE or FALSE.", call. = FALSE)
  }

  # Sort clusters alphabetically (TNA convention)
  cluster_list <- cluster_list[order(names(cluster_list))]
  cluster_names <- names(cluster_list)
  n_clusters <- length(cluster_names)

  # Recode to cluster labels
  from_clusters <- cluster_lookup[from_nodes]
  to_clusters <- cluster_lookup[to_nodes]

  # ---- Between-cluster matrix (includes diagonal = within-cluster loops) ----
  between_raw <- matrix(0, n_clusters, n_clusters,
                        dimnames = list(cluster_names, cluster_names))

  # Include ALL transitions -- node-level self-loops (A->A) are valid
  # cluster-level self-loops (e.g. discuss->discuss = Social->Social)
  b_from <- from_clusters
  b_to <- to_clusters
  b_w <- weights

  if (length(b_from) > 0) {
    cluster_sizes <- vapply(cluster_list, length, integer(1L))
    n_possible_mat <- if (method == "density") {
      outer(cluster_sizes, cluster_sizes, "*")
    } else NULL

    fc <- factor(b_from, levels = cluster_names)
    tc <- factor(b_to, levels = cluster_names)
    agg <- tapply(b_w, list(fc, tc), function(w) {
      net_aggregate_weights(w, method)
    }, default = NA_real_)
    agg[is.na(agg)] <- 0
    between_raw[] <- as.numeric(agg)
    if (!is.null(n_possible_mat)) {
      sums <- tapply(b_w, list(fc, tc), sum, default = 0)
      sums[is.na(sums)] <- 0
      between_raw[] <- as.numeric(sums) / as.numeric(n_possible_mat)
      between_raw[!is.finite(between_raw)] <- 0
    }
  }

  # Process between weights
  between_weights <- .process_weights(between_raw, type, directed)

  # Compute initial-state probabilities. When sequence data is available we
  # use the empirical first-state distribution (matches tna::tna()'s inits);
  # otherwise (edgelist or no data) we fall back to the column-sum proxy
  # documented earlier -- there is no first-state to read off an edgelist.
  is_edgelist_data <- identical(attr(data, "source"), "edgelist")
  is_seq_for_inits <- is.data.frame(data) && !is_edgelist_data
  between_inits <- if (is_seq_for_inits) {
    # tna-equivalent macro inits: tabulate column 1 (first time-step) of the
    # wide sequence, drop NAs, map each first state to its cluster, normalise.
    first_state <- as.character(data[[1]])
    first_state <- first_state[!is.na(first_state)]
    if (length(first_state) == 0L) {
      setNames(rep(0, n_clusters), cluster_names)
    } else {
      first_clu <- cluster_lookup[first_state]
      fr <- table(factor(first_clu, levels = cluster_names))
      out <- as.numeric(fr) / sum(fr)
      names(out) <- cluster_names
      out
    }
  } else {
    col_sums <- colSums(between_raw)
    total <- sum(col_sums)
    out <- if (total > 0) col_sums / total else rep(1 / n_clusters, n_clusters)
    names(out) <- cluster_names
    out
  }

  # ---- Build recoded sequence data for bootstrap compatibility ----
  # The bootstrap fast path (.bootstrap_transition) needs a data.frame
  # that it can resample row-wise. Two cases:
  #
  #   * Wide-format sequence input (rows = actors, cols = time-steps):
  #     recode every cell to the cluster label (macro) or keep only
  #     nodes of that cluster (per-cluster).
  #
  #   * Edgelist input (each row is a single from -> to transition):
  #     build 2-column pseudo-sequences, each row representing one
  #     transition. Resampling rows = resampling transitions, which
  #     is the correct bootstrap at the edge level.
  between_seq_data <- NULL
  within_seq_data_list <- NULL
  is_seq <- is.data.frame(data) && !is_edgelist_data

  if (is_seq) {
    # Recode sequences to cluster labels
    between_seq_data <- as.data.frame(
      lapply(data, function(col) {
        recoded <- unname(cluster_lookup[as.character(col)])
        recoded
      }),
      stringsAsFactors = FALSE
    )

    # Filter sequences per cluster (keep only that cluster's nodes, NA others)
    within_seq_data_list <- lapply(cluster_list, function(cl_nodes) {
      as.data.frame(
        lapply(data, function(col) {
          vals <- as.character(col)
          vals[!vals %in% cl_nodes] <- NA_character_
          vals
        }),
        stringsAsFactors = FALSE
      )
    })
  } else if (length(from_nodes) > 0L) {
    # Edgelist branch: every transition -> one 2-row pseudo-sequence.
    # Columns named `from`/`to` so the data is self-describing on inspection.
    # Tagged with source = "edgelist" so bootstrap_network() can warn --
    # edgelist bootstrap treats transitions as independent, ignoring
    # within-actor correlation, so CIs are anti-conservative.
    between_seq_data <- data.frame(
      from = unname(from_clusters),
      to = unname(to_clusters),
      stringsAsFactors = FALSE
    )
    attr(between_seq_data, "source") <- "edgelist"

    within_seq_data_list <- lapply(cluster_list, function(cl_nodes) {
      keep <- from_nodes %in% cl_nodes & to_nodes %in% cl_nodes
      df <- data.frame(
        from = from_nodes[keep],
        to = to_nodes[keep],
        stringsAsFactors = FALSE
      )
      attr(df, "source") <- "edgelist"
      df
    })
  }

  between <- .mcml_layer(
    weights = between_weights,
    inits   = between_inits,
    labels  = cluster_names,
    data    = between_seq_data
  )

  # ---- Within-cluster matrices ----
  within_data <- NULL
  if (isTRUE(compute_within)) {
    # Filter within-cluster transitions (same cluster, including self-loops)
    is_within <- from_clusters == to_clusters
    w_from <- from_nodes[is_within]
    w_to <- to_nodes[is_within]
    w_w <- weights[is_within]

    within_data <- lapply(seq_len(n_clusters), function(i) {
      cl_name <- cluster_names[i]
      cl_nodes <- cluster_list[[cl_name]]
      n_i <- length(cl_nodes)

      if (n_i <= 1) {
        # Singleton cluster: aggregate self-loop weight, then route through
        # .process_weights so type = "tna" row-normalises just like the
        # multi-node branch (1x1 [N] -> [1.0]; 1x1 [0] stays [0]).
        keep_self <- w_from %in% cl_nodes & w_to %in% cl_nodes
        self_weight <- if (any(keep_self)) {
          # one node pair: one possible edge, so "density" is the sum, as
          # the macro layer's sum / (n_i * n_j)
          net_aggregate_weights(w_w[keep_self], method, n_possible = 1)
        } else { 0 }
        within_raw <- matrix(self_weight, 1, 1,
                             dimnames = list(cl_nodes, cl_nodes))
        within_weights_i <- .process_weights(within_raw, type, directed)
        within_inits_i <- if (is_seq_for_inits) {
          .first_state_distribution(data, levels = cl_nodes, restrict_to = cl_nodes)
        } else {
          setNames(1, cl_nodes)
        }
      } else {
        # Filter transitions for this cluster (keep self-loops)
        keep <- w_from %in% cl_nodes & w_to %in% cl_nodes
        cf <- w_from[keep]
        ct <- w_to[keep]
        cw <- w_w[keep]

        within_raw <- matrix(0, n_i, n_i,
                              dimnames = list(cl_nodes, cl_nodes))

        if (length(cf) > 0) {
          fc <- factor(cf, levels = cl_nodes)
          tc <- factor(ct, levels = cl_nodes)
          agg <- tapply(cw, list(fc, tc), function(w) {
            net_aggregate_weights(w, method, n_possible = 1)
          }, default = NA_real_)
          agg[is.na(agg)] <- 0
          within_raw[] <- as.numeric(agg)
        }

        within_weights_i <- .process_weights(within_raw, type, directed)

        # Within-cluster inits: first cluster-k state per row, NA-masked
        # against other clusters (matches tna::tna() on the masked sequence).
        # For edgelist or matrix-only input, fall back to the column-sum proxy.
        within_inits_i <- if (is_seq_for_inits) {
          .first_state_distribution(data, levels = cl_nodes,
                                 restrict_to = cl_nodes)
        } else {
          col_sums_w <- colSums(within_raw, na.rm = TRUE)
          total_w <- sum(col_sums_w, na.rm = TRUE)
          out <- if (!is.na(total_w) && total_w > 0) col_sums_w / total_w
                 else rep(1 / n_i, n_i)
          names(out) <- cl_nodes
          out
        }
      }

      # Attach filtered sequence data for this cluster
      cl_seq_data <- if (!is.null(within_seq_data_list)) {
        within_seq_data_list[[cl_name]]
      } else {
        NULL
      }

      .mcml_layer(
        weights = within_weights_i,
        inits   = within_inits_i,
        labels  = cl_nodes,
        data    = cl_seq_data
      )
    })
    names(within_data) <- cluster_names
  }

  # ---- Edges data.frame ----
  edge_type <- ifelse(from_clusters == to_clusters, "within", "between")
  edges <- data.frame(
    from = from_nodes,
    to = to_nodes,
    weight = weights,
    cluster_from = unname(from_clusters),
    cluster_to = unname(to_clusters),
    type = edge_type,
    stringsAsFactors = FALSE
  )

  # ---- Assemble result ----
  all_nodes <- sort(unique(unlist(cluster_list, use.names = FALSE)))
  n_nodes <- length(all_nodes)

  structure(
    list(
      macro = between,
      clusters = within_data,
      cluster_members = cluster_list,
      edges = edges,
      meta = list(
        type = type,
        method = method,
        # Effective directedness of the stored weights: .process_weights
        # symmetrizes for type = "cooccurrence" (and for directed = FALSE),
        # so record what the weights ARE, not what was asked for.
        directed = isTRUE(directed) && type != "cooccurrence",
        n_nodes = n_nodes,
        n_clusters = n_clusters,
        cluster_sizes = vapply(cluster_list, length, integer(1)),
        source = "transitions"
      )
    ),
    class = "mcml"
  )
}

#' Build MCML from edge list data.frame
#' @noRd
.build_mcml_edgelist <- function(df, clusters, method, type,
                                  directed, compute_within, from = "from",
                                  to = "to", weight = NULL) {
  if (identical(from, to)) {
    stop(errorCondition(
      "Edge-list input must include at least two endpoint columns.",
      class = "nestimate_bad_input", call = NULL))
  }
  has_weight <- !is.null(weight)
  from_vals <- as.character(df[[from]])
  to_vals <- as.character(df[[to]])
  if (any(is.na(from_vals) | !nzchar(from_vals) |
          is.na(to_vals) | !nzchar(to_vals))) {
    stop(errorCondition(
      "Edge-list source and target columns must not contain missing or empty values.",
      class = "nestimate_bad_input", call = NULL))
  }

  if (has_weight) {
    weights <- df[[weight]]
    if (!is.numeric(weights)) {
      stop(errorCondition("Edge-list weight column must be numeric.",
        class = "nestimate_bad_input", call = NULL))
    }
    if (any(is.na(weights) | !is.finite(weights))) {
      stop(errorCondition(
        "Edge-list weight column must contain finite non-missing values.",
        class = "nestimate_bad_input", call = NULL))
    }
    if (any(weights < 0)) {
      stop(errorCondition(
        "Edge-list weight column must not contain negative values.",
        class = "nestimate_bad_input", call = NULL))
    }
  } else {
    weights <- rep(1, nrow(df))
  }

  all_nodes <- sort(unique(c(from_vals, to_vals)))

  # Handle clusters parameter
  if (is.character(clusters) && length(clusters) == 1 &&
      clusters %in% names(df)) {
    # Column name: build lookup from both from+group and to+group
    group_col <- as.character(df[[clusters]])
    if (any(is.na(group_col) | !nzchar(group_col))) {
      stop(errorCondition(
        "Edge-list cluster column must not contain missing or empty values.",
        class = "nestimate_bad_input", call = NULL))
    }

    node_group <- data.frame(
      node = c(from_vals, to_vals),
      group = c(group_col, group_col),
      stringsAsFactors = FALSE
    )
    node_group <- unique(node_group)
    conflicts <- unique(node_group$node[duplicated(node_group$node)])
    if (length(conflicts) > 0L) {
      stop(errorCondition(
        paste0(
          "Edge-list cluster column assigns nodes to multiple groups: ",
          paste(utils::head(conflicts, 5L), collapse = ", "),
          ". Pass clusters as a named list or node-group data.frame instead."),
        class = "nestimate_bad_input", call = NULL))
    }

    full_map <- setNames(node_group$group, node_group$node)

    # Build cluster_list
    cluster_list <- split(names(full_map), unname(full_map))
    cluster_list <- lapply(cluster_list, sort)
    cluster_lookup <- full_map

    # Re-derive all_nodes from the lookup
    all_nodes <- sort(names(cluster_lookup))
  } else {
    # List or vector clusters
    if (is.null(clusters)) {
      stop(errorCondition("clusters argument is required for data.frame input",
        class = "nestimate_bad_input", call = NULL))
    }

    if (is.list(clusters) && !is.data.frame(clusters)) {
      cluster_list <- clusters
    } else {
      # Membership vector
      cluster_list <- .normalize_clusters(clusters, all_nodes)
    }

    cluster_lookup <- .build_cluster_lookup(cluster_list, all_nodes)
  }

  attr(df, "source") <- "edgelist"

  .build_from_transitions(from_vals, to_vals, weights,
                            cluster_lookup, cluster_list,
                            method, type, directed, compute_within,
                            data = df)
}

#' Build MCML from sequence data.frame
#' @noRd
.build_mcml_sequence <- function(df, clusters, method, type,
                                  directed, compute_within) {

  if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)

  stopifnot(is.data.frame(df))

  nc <- ncol(df)
  if (nc < 2) {
    stop(errorCondition(
      "Sequence data must have at least 2 columns (time steps)",
      class = "nestimate_bad_input", call = NULL))
  }

  # Extract consecutive pairs: (col[t], col[t+1]) for all rows
  pairs <- lapply(seq_len(nc - 1), function(t) {
    from_t <- as.character(df[[t]])
    to_t <- as.character(df[[t + 1]])
    data.frame(from = from_t, to = to_t, stringsAsFactors = FALSE)
  })
  pairs <- as.data.frame(data.table::rbindlist(pairs))

  # Remove NA pairs
  valid <- !is.na(pairs$from) & !is.na(pairs$to)
  from_vals <- pairs$from[valid]
  to_vals <- pairs$to[valid]
  weights <- rep(1, length(from_vals))

  observed_nodes <- as.character(unlist(df, use.names = FALSE))
  observed_nodes <- observed_nodes[!is.na(observed_nodes)]
  all_nodes <- sort(unique(observed_nodes))

  if (is.null(clusters)) {
    stop(errorCondition("clusters argument is required for sequence data",
      class = "nestimate_bad_input", call = NULL))
  }

  if (is.list(clusters) && !is.data.frame(clusters)) {
    cluster_list <- clusters
  } else {
    cluster_list <- .normalize_clusters(clusters, all_nodes)
  }

  cluster_lookup <- .build_cluster_lookup(cluster_list, all_nodes)

  .build_from_transitions(from_vals, to_vals, weights,
                            cluster_lookup, cluster_list,
                            method, type, directed, compute_within,
                            data = df)
}

#' Process weights based on type
#' @noRd
.process_weights <- function(raw_weights, type, directed = TRUE) {
  if (!is.matrix(raw_weights) || !is.numeric(raw_weights)) {
    stop("'raw_weights' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(raw_weights) != ncol(raw_weights)) {
    stop("'raw_weights' must be a square matrix.", call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single non-missing character value.", call. = FALSE)
  }
  valid_types <- c("raw", "frequency", "cooccurrence", "tna")
  if (!type %in% valid_types) {
    stop("'type' must be one of: ",
         paste(valid_types, collapse = ", "), call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!isTRUE(directed) || type == "cooccurrence") {
    raw_weights <- (raw_weights + t(raw_weights)) / 2
  }

  if (type == "raw" || type == "frequency") {
    return(raw_weights)
  }

  if (type == "cooccurrence") {
    return(raw_weights)
  }

  if (type == "tna") {
    # Row-normalize so rows sum to 1
    rs <- rowSums(raw_weights, na.rm = TRUE)
    processed <- raw_weights / ifelse(rs == 0 | is.na(rs), 1, rs)
    processed[is.na(processed)] <- 0
    return(processed)
  }

  raw_weights # nocov
}

#' Promote the Layers of an mcml to Networks
#'
#' @description
#' `as_tna()` turns an `mcml` object into a `netobject_group`: one network
#' for the cluster-level (macro) layer, followed by one per cluster for the
#' within-cluster layers. The stored weights are carried over as they are --
#' nothing is re-normalised -- so the networks hold the aggregation and
#' `type` chosen when the `mcml` was built. Despite its name it returns
#' Nestimate networks, not `tna` models; a `tna` model passed to it is
#' returned unchanged.
#'
#' This is the step that lets an MCML result reach the verbs that take a
#' group of networks: printing, the network metrics of [summary()],
#' [bootstrap_network()] (for an `mcml` built from sequences, whose layers
#' keep them; the layers of a [cluster_summary()] result carry none) and
#' drawing with cograph.
#'
#' @param x An `mcml` object from [cluster_summary()] or [build_mcml()]. The
#'   default method also takes a `tna` model.
#' @param expand For the `mcml` method: clusters whose member states replace
#'   the collapsed cluster node in the macro network, as in
#'   [macro_network()]. Default `NULL` (the macro network stays fully
#'   collapsed). The within-cluster networks are not affected. It needs an
#'   `mcml` that [build_mcml()] built from a wide sequence frame, a sequence
#'   matrix or an event log; one from [cluster_summary()], an edge list, a
#'   `netobject` or a `tna` model raises an error of class
#'   `nestimate_no_expand_source`.
#' @param ... Not used.
#'
#' @return A `netobject_group`: a named list whose first element, `macro`, is
#'   the k x k cluster-level network, followed by one network per cluster.
#'   Print it for one row per layer (nodes, edges, weight range), call
#'   [summary()] for the network metrics of each layer, and draw a layer
#'   with `cograph::splot()`. Each
#'   element is a `netobject` (class `c("netobject", "cograph_network")`)
#'   with the layer's weights, initial probabilities, node and edge tables
#'   and data, and records `method = "relative"` when the `mcml` holds
#'   row-normalised weights (`type = "tna"`), `method = "co_occurrence"`
#'   for symmetrised counts (`type = "cooccurrence"`), and
#'   `method = "frequency"` for counts (`type` `"raw"` or `"frequency"`) and
#'   for the aggregates of [cluster_summary()]. With
#'   `expand`, the `macro` element is the mixed-resolution network of
#'   [macro_network()] instead of the fully collapsed one. The default
#'   method returns a `tna` model unchanged and raises an error for any
#'   other class.
#'
#' @details
#' ## Workflow
#'
#' \preformatted{
#' # Full MCML workflow
#' net <- build_network(data, method = "relative")
#' cs   <- cluster_summary(net, clusters = group_assignments)
#' nets <- as_tna(cs)
#'
#' # Every layer is an ordinary netobject
#' print(nets)      # one line per layer
#' summary(nets)    # network metrics per layer
#' }
#'
#' ## Sink nodes
#'
#' Every cluster is returned, whatever its row sums. A node with no outgoing
#' weight is a legitimate sink (a terminal state), and its row stays all
#' zero in the network, for `method = "relative"` and `"frequency"` alike:
#' the stored weights are never re-normalised, so a sink row needs no
#' special handling.
#'
#' @export
#' @seealso
#'   \code{\link{cluster_summary}} and \code{\link{build_mcml}} to create the
#'   input object,
#'   \code{\link{macro_network}} for a macro layer with one cluster expanded,
#'   \code{\link{as_networks}} for the psychometric-network counterpart
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(runif(36), 6, 6)
#' rownames(mat) <- colnames(mat) <- LETTERS[1:6]
#' clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"), G3 = c("E", "F"))
#' cs <- cluster_summary(mat, clusters)
#' nets <- as_tna(cs)
#' nets
#' summary(nets)
as_tna <- function(x, ...) {
  UseMethod("as_tna")
}

#' @rdname as_tna
#' @export
as_tna.mcml <- function(x, expand = NULL, ...) {
  # `expand` replaces the collapsed macro layer with the mixed-resolution
  # network from macro_network(): the named clusters shown as their member
  # states, every other cluster still one node. The per-cluster layers are
  # untouched, so only the macro changes.
  if (!is.null(expand)) {
    macro_net <- macro_network(x, expand = expand,
                               method = if (identical(x$meta$type, "tna")) {
                                 "relative"
                               } else {
                                 "frequency"
                               })
  }
  # Determine the method to record on the wrapped netobjects:
  #   * "raw"/"frequency"/"aggregate" -> "frequency" (matrix-path or counts)
  #   * "tna" (already row-normalised) -> "relative"
  #   * "cooccurrence" (symmetrised counts, rows do not sum to 1, the
  #     "undirected co-occurrence" of build_mcml()) -> "co_occurrence";
  #     recording "relative" made downstream verbs read counts as
  #     probabilities.
  meta_type <- x$meta$type
  net_method <- if (is.null(meta_type) ||
                    meta_type %in% c("raw", "frequency", "aggregate")) {
    "frequency"
  } else if (identical(meta_type, "cooccurrence")) {
    "co_occurrence"
  } else {
    "relative"
  }
  directed <- if (!is.null(x$meta$directed)) x$meta$directed else TRUE

  # Macro
  if (is.null(expand)) {
    macro_net <- .wrap_netobject(x$macro$weights, data = x$macro$data,
                                 method = net_method, directed = directed,
                                 inits = x$macro$inits)
  }

  # Per-cluster. Every cluster is kept. A zero-sum row is a legitimate sink
  # (a terminal state with no out-transitions): the weights are stored as-is
  # by .wrap_netobject() -- which never re-normalises -- so a sink row simply
  # stays all-zero, exactly as the `relative` estimator and the frequency
  # path already keep it. Dropping the whole cluster over one sink node was
  # over-conservative and is no longer done.
  cluster_nets <- if (!is.null(x$clusters)) {
    nets <- lapply(names(x$clusters), function(cl) {
      .wrap_netobject(x$clusters[[cl]]$weights, data = x$clusters[[cl]]$data,
                      method = net_method, directed = directed,
                      inits = x$clusters[[cl]]$inits)
    })
    names(nets) <- names(x$clusters)
    nets
  } else {
    list()
  }

  result <- c(list(macro = macro_net), cluster_nets)
  class(result) <- "netobject_group"
  result
}


#' Wrap a weight matrix + optional data into a minimal netobject
#' @noRd
.wrap_netobject <- function(mat, data = NULL, method = "relative",
                            directed = TRUE, inits = NULL) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("'mat' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("'mat' must be a square matrix.", call. = FALSE)
  }
  .validate_mcml_matrix(mat)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)) {
    stop("'method' must be a single non-missing character value.",
         call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }
  states <- rownames(mat)
  if (is.null(states)) {
    states <- as.character(seq_len(nrow(mat)))
    dimnames(mat) <- list(states, states)
  }
  edges <- .extract_edges_from_matrix(mat, directed = directed)
  nodes_df <- data.frame(
    id = seq_along(states), label = states, name = states,
    x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
  )
  if (!is.null(inits)) {
    if (!is.numeric(inits) || length(inits) != length(states) ||
        any(is.na(inits) | !is.finite(inits))) {
      stop("'inits' must be a finite numeric vector with one value per state.",
           call. = FALSE)
    }
    if (!is.null(names(inits))) {
      missing_inits <- setdiff(states, names(inits))
      extra_inits <- setdiff(names(inits), states)
      if (length(missing_inits) > 0L || length(extra_inits) > 0L) {
        stop("'inits' names must match matrix state names.", call. = FALSE)
      }
      inits <- inits[states]
    } else {
      names(inits) <- states
    }
  }

  structure(
    list(
      data       = data,
      weights    = mat,
      inits      = inits,
      nodes      = nodes_df,
      edges      = edges,
      directed   = directed,
      method     = method,
      params     = list(),
      scaling    = NULL,
      threshold  = 0,
      n_nodes    = length(states),
      n_edges    = nrow(edges),
      level      = NULL,
      meta       = list(source = "nestimate", layout = NULL,
                        tna = list(method = method)),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}

#' @rdname as_tna
#' @export
as_tna.default <- function(x, ...) {
 if (inherits(x, "tna")) {
    return(x)
  }
  stop(errorCondition(
    paste0("Cannot convert object of class '", class(x)[1], "' to tna"),
    class = "nestimate_bad_input", call = NULL))
}


# ==============================================================================
# as_htna -- Full node-level network grouped by cluster (for cograph::plot_htna)
# ==============================================================================

#' Build a grouped node-level network (htna) from data and a clustering
#'
#' @description
#' `as_htna()` builds the full node-level network from the original data
#' with [build_network()] and attaches a cluster to every node, returning a
#' single `htna` network: every state is a node, and the clusters group the
#' nodes. It is the node-level counterpart of [build_mcml()]: where
#' `build_mcml()` collapses the network to a cluster-level (macro) summary,
#' `as_htna()` keeps every node and every transition, including the
#' between-cluster transitions an `mcml` keeps only in aggregate. It takes
#' data with a clustering, an `mcml` (which supplies the clustering, and
#' often the data), or a fitted `net_mmm`.
#'
#' The result is an ordinary `netobject`, so the inference verbs apply to it
#' ([bootstrap_network()], [net_centrality()], [permutation()]), and
#' cograph draws it as a grouped network:
#' `cograph::plot_htna(as_htna(data, clusters))`.
#'
#' @param x The data, in any form [build_network()] accepts (a sequence data
#'   frame, an event log, an edge list, a transition matrix, a `netobject`
#'   or a `tna` model); or an `mcml` object, which provides the cluster
#'   membership and, when it was built from a wide sequence frame or an
#'   event log, the data as well (see `data`); or a `net_mmm` fitted by
#'   [cluster_mmm()] on an HTNA network, which is materialised into one HTNA
#'   per sequence cluster using its preserved actor partition.
#' @param clusters The cluster of each node: a named list of node-name
#'   vectors, a membership vector with one value per node (in the order of
#'   the built network's nodes, alphabetical for sequence data, or named by
#'   node), or a two-column data frame of node and cluster names. Default `NULL`, under which the
#'   node groups of the network (`node_groups` with `node` and `group`
#'   columns) or an `mcml`'s clusters are used; with neither, an error is
#'   raised. Every node must belong to a cluster, and unassigned nodes are
#'   named in an error. A fitted `net_mmm` takes no `clusters`: passing them
#'   raises an error.
#' @param method The estimator passed to [build_network()]. Default
#'   `"relative"` (row-normalised transitions).
#' @param ... Further arguments for [build_network()], such as `actor`,
#'   `action` and `time` for an event log. A `group` argument raises an
#'   error, because a single network is expected.
#' @return For data and `mcml` input, a single `htna` network (class
#'   `c("htna", "netobject", "cograph_network")`) over all nodes. Print it for
#'   an overview, call [summary()] for its network metrics, and draw it with
#'   `cograph::plot_htna()`. The clusters are stored as a factor in the
#'   `groups` column of the node table and as character values in the
#'   `group` column of `node_groups`; `actor_levels` records their order and
#'   is also attached to `node_groups`, so the partition survives a round
#'   trip. For compatibility, the node table also keeps a `cluster` column,
#'   and the membership is kept in the `"cluster_members"` attribute. A
#'   fitted `net_mmm` returns an `htna_group`, one HTNA network per sequence
#'   cluster, with the MMM diagnostics kept; a `net_mmm` fitted on data
#'   without an HTNA actor partition raises an error.
#' @details
#' ## Why it rebuilds from data
#'
#' It rebuilds from data because an `mcml` holds no faithful node-level
#' network: its macro sequences are recoded to cluster labels, and its
#' per-cluster data keep only that cluster's states, so the node-level
#' transitions between clusters survive only in the original data. An
#' `mcml` therefore supplies the cluster membership, and either the source
#' it retained or the `data` passed supplies the transitions.
#' @seealso \code{\link{build_mcml}}, \code{\link{build_network}}; plot with
#'   \code{cograph::plot_htna()}.
#' @export
#' @examples
#' seqs <- data.frame(
#'   t1 = c("A", "C", "E", "B"), t2 = c("B", "D", "F", "A"),
#'   t3 = c("C", "A", "E", "D"), stringsAsFactors = FALSE
#' )
#' clusters <- list(C1 = c("A", "B"), C2 = c("C", "D"), C3 = c("E", "F"))
#' net <- as_htna(seqs, clusters)
#' net
#' if (requireNamespace("cograph", quietly = TRUE)) cograph::plot_htna(net)
as_htna <- function(x, clusters = NULL, method = "relative", ...) {
  UseMethod("as_htna")
}

#' @rdname as_htna
#' @param data For the `mcml` method: the original data the `mcml` was built
#'   from. Optional when it was built from a wide sequence frame, a
#'   character sequence matrix, or an event log (which
#'   `build_mcml()` turns into sequences first): that source is kept on the
#'   object, so `as_htna(mcml)` works on its own, except for an `mcml` built
#'   with `labels`, whose kept source holds the original state names (pass
#'   `clusters` with those names then). Required for
#'   an `mcml` built from a matrix, an edge list, a `netobject` or a `tna`
#'   model, which keep no usable node-level source; without it an error is
#'   raised. Default `NULL`.
#' @export
as_htna.mcml <- function(x, clusters = NULL, method = "relative",
                         data = NULL, ...) {
  # An mcml built from sequence/edgelist data stashes that source (see
  # build_mcml), so as_htna(mcml) works on its own; reuse it when `data`
  # is not supplied.
  src <- attr(x, "htna_source")
  if (is.null(data) && !is.null(src)) data <- src
  if (is.null(data)) {
    stop(errorCondition(
      paste0(
        "This mcml carries no expandable node-level source (it was built ",
        "from a matrix, aggregate, or edge list), so it cannot be expanded ",
        "to an htna on its own. Pass the original data:\n",
        "  as_htna(mcml, data = <original data>)   # reuses mcml$cluster_members\n",
        "or equivalently as_htna(<original data>, clusters = mcml$cluster_members)."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.null(clusters)) clusters <- x$cluster_members
  as_htna.default(data, clusters = clusters, method = method, ...)
}

#' @rdname as_htna
#' @export
as_htna.net_mmm <- function(x, clusters = NULL, method = "relative", ...) {
  if (!is.null(clusters)) {
    stop(errorCondition(
      paste0(
        "`clusters` is not used for a fitted `net_mmm`: its sequence ",
        "clusters are already fitted and its HTNA actor partition is ",
        "preserved separately."),
      class = "nestimate_bad_input", call = NULL))
  }

  nets <- build_network(x, method = method, ...)
  if (!inherits(nets, "htna_group")) {
    stop(errorCondition(
      paste0(
        "The fitted `net_mmm` does not carry an HTNA actor partition. ",
        "Fit `cluster_mmm()` on an HTNA object, or materialize ordinary ",
        "cluster networks with `build_network()`."),
      class = "nestimate_bad_input", call = NULL))
  }
  nets
}

#' @rdname as_htna
#' @export
as_htna.default <- function(x, clusters = NULL, method = "relative", ...) {
  net <- build_network(x, method = method, ...)
  if (!inherits(net, "netobject")) {
    stop(errorCondition(
      paste0(
        "as_htna() expects a single node-level network; build_network() ",
        "returned class '",
        class(net)[1],
        "'. Provide ungrouped data ",
        "(do not pass a 'group' argument)."),
      class = "nestimate_bad_input", call = NULL))
  }
  node_names <- net$nodes$label

  # Derive clusters from carried node groups when not supplied.
  if (is.null(clusters)) {
    ng <- net$node_groups
    if (is.data.frame(ng) && all(c("node", "group") %in% names(ng))) {
      clusters <- stats::setNames(ng$group, ng$node)
    } else {
      stop(errorCondition(
        paste0(
          "'clusters' is required: a named list of node vectors, a per-node ",
          "membership vector, or a two-column data frame."),
        class = "nestimate_bad_input", call = NULL))
    }
  }
  members <- .normalize_clusters(clusters, node_names)

  # node -> cluster label map; every node must belong to exactly one cluster.
  map <- stats::setNames(
    rep(names(members), lengths(members)),
    unlist(members, use.names = FALSE)
  )
  missing_nodes <- setdiff(node_names, names(map))
  if (length(missing_nodes) > 0L) {
    stop(errorCondition(paste0("These nodes are not assigned to any cluster: ",
         paste(missing_nodes, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  cluster_labels <- unname(map[node_names])
  actor_levels <- names(members)

  # Complete the same partition contract produced by htna::build_htna().
  # Keep the older `cluster` column and attribute below for callers that used
  # as_htna() before it carried the formal htna class.
  net$nodes$cluster <- cluster_labels
  net$nodes$groups <- factor(cluster_labels, levels = actor_levels)
  net$node_groups <- data.frame(
    node = node_names, group = cluster_labels,
    stringsAsFactors = FALSE
  )
  attr(net$node_groups, "actor_levels") <- actor_levels
  net$actor_levels <- actor_levels
  attr(net, "cluster_members") <- members
  class(net) <- unique(c("htna", "netobject", "cograph_network", class(net)))
  net
}

#' Normalize cluster specification to list format
#' @keywords internal
#' @noRd
.normalize_clusters <- function(clusters, node_names) {
  if (is.data.frame(clusters)) {
    # Data frame with node and group columns
    if (ncol(clusters) < 2L) {
      stop(errorCondition("clusters data.frame must have at least two columns.",
        class = "nestimate_bad_input", call = NULL))
    }
    nodes <- as.character(clusters[[1]])
    groups <- as.character(clusters[[2]])
    if (any(is.na(nodes) | !nzchar(nodes))) {
      stop(errorCondition(
        "clusters data.frame node column must not contain missing or empty values.",
        class = "nestimate_bad_input", call = NULL))
    }
    if (any(is.na(groups) | !nzchar(groups))) {
      stop(errorCondition(
        "clusters data.frame group column must not contain missing or empty values.",
        class = "nestimate_bad_input", call = NULL))
    }
    clusters <- split(nodes, groups)
  }

  if (is.list(clusters)) {
    # Already a list - validate node names
    .validate_cluster_partition(clusters, node_names)
    return(clusters)
  }

  if (is.vector(clusters) && (is.numeric(clusters) || is.integer(clusters))) {
    # Membership vector
    clusters <- .align_cluster_membership(clusters, node_names)
    if (any(is.na(clusters) | !is.finite(clusters))) {
      stop(errorCondition(
        "Membership vector must not contain missing or non-finite values.",
        class = "nestimate_bad_input", call = NULL))
    }
    # Convert to list
    unique_clusters <- sort(unique(clusters))
    cluster_list <- lapply(unique_clusters, function(k) {
      node_names[clusters == k]
    })
    names(cluster_list) <- as.character(unique_clusters)
    return(cluster_list)
  }

  if (is.factor(clusters) || is.character(clusters)) {
    # Named membership
    clusters <- .align_cluster_membership(clusters, node_names)
    clusters <- as.character(clusters)
    if (any(is.na(clusters) | !nzchar(clusters))) {
      stop(errorCondition(
        "Membership vector must not contain missing or empty values.",
        class = "nestimate_bad_input", call = NULL))
    }
    unique_clusters <- unique(clusters)
    cluster_list <- lapply(unique_clusters, function(k) {
      node_names[clusters == k]
    })
    names(cluster_list) <- unique_clusters
    return(cluster_list)
  }

  stop(errorCondition("clusters must be a list, numeric vector, or factor",
    class = "nestimate_bad_input", call = NULL))
}

.align_cluster_membership <- function(clusters, node_names) {
  if (length(clusters) != length(node_names)) {
    stop(errorCondition(paste0("Membership vector length (", length(clusters),
         ") must equal number of nodes (", length(node_names), ")"),
      class = "nestimate_bad_input", call = NULL))
  }
  nm <- names(clusters)
  if (is.null(nm)) {
    return(clusters)
  }
  if (any(is.na(nm) | !nzchar(nm))) {
    stop(errorCondition(
      "Named membership vector names must not contain missing or empty values.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (anyDuplicated(nm)) {
    stop(errorCondition("Named membership vector names must be unique.",
      class = "nestimate_bad_input", call = NULL))
  }
  missing_nodes <- setdiff(node_names, nm)
  extra_nodes <- setdiff(nm, node_names)
  if (length(missing_nodes) > 0L || length(extra_nodes) > 0L) {
    stop(errorCondition(
      "Named membership vector names must match node names exactly.",
      class = "nestimate_bad_input", call = NULL))
  }
  clusters[node_names]
}

.validate_cluster_partition <- function(clusters, node_names) {
  if (!is.list(clusters) || is.data.frame(clusters)) {
    stop(errorCondition("clusters must be a named list.",
      class = "nestimate_bad_input", call = NULL))
  }
  cluster_names <- names(clusters)
  if (is.null(cluster_names) || any(!nzchar(cluster_names)) ||
      any(is.na(cluster_names))) {
    stop(errorCondition("clusters list must have non-empty cluster names.",
      class = "nestimate_bad_input", call = NULL))
  }
  cluster_sizes <- vapply(clusters, length, integer(1L))
  if (any(cluster_sizes == 0L)) {
    stop(errorCondition(paste0("clusters list contains empty clusters: ",
         paste(names(clusters)[cluster_sizes == 0L], collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  nodes <- unlist(clusters, use.names = FALSE)
  nodes <- as.character(nodes)
  if (any(is.na(nodes)) || any(!nzchar(nodes))) {
    stop(errorCondition("clusters list contains missing or empty node names.",
      class = "nestimate_bad_input", call = NULL))
  }

  unknown <- setdiff(nodes, node_names)
  if (length(unknown) > 0L) {
    stop(errorCondition(paste0("Unknown nodes in clusters: ",
         paste(utils::head(unknown, 5L), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  duplicated_nodes <- unique(nodes[duplicated(nodes)])
  if (length(duplicated_nodes) > 0L) {
    stop(errorCondition(paste0("Nodes assigned to multiple clusters: ",
         paste(utils::head(duplicated_nodes, 5L), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  unmapped <- setdiff(node_names, nodes)
  if (length(unmapped) > 0L) {
    stop(errorCondition(paste0("Unmapped nodes: ",
         paste(utils::head(unmapped, 5L), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  invisible(TRUE)
}

# ==============================================================================
# S3 Methods
# ==============================================================================

#' @rdname build_mcml
#' @export
print.mcml <- function(x, ...) {
  .mcml_check_unused_dots("print.mcml", ...)
  n_clusters <- x$meta$n_clusters
  n_nodes <- x$meta$n_nodes
  cluster_sizes <- x$meta$cluster_sizes

  cat("MCML Network\n")
  cat("============\n")
  if (!is.null(x$meta$type)) {
    cat("Type:", x$meta$type, " | Method:", x$meta$method, "\n")
  } else {
    cat("Method:", x$meta$method, " (matrix path, aggregation-only)\n")
  }
  cat("Nodes:", n_nodes, " | Clusters:", n_clusters, "\n")
  if (!is.null(x$edges)) {
    cat("Transitions:", nrow(x$edges), "\n")
    n_between <- sum(x$edges$type == "between")
    n_within <- sum(x$edges$type == "within")
    cat("  Macro:", n_between, " | Per-cluster:", n_within, "\n\n")
  } else {
    cat("\n")
  }

  cat("Clusters:\n")
  for (cl in names(x$cluster_members)) {
    cat("  ", cl, " (", cluster_sizes[cl], "): ",
        paste(x$cluster_members[[cl]], collapse = ", "), "\n", sep = "")
  }

  cat("\nMacro (cluster-level) weights:\n")
  print(round(x$macro$weights, 4))

  invisible(x)
}

#' @rdname build_mcml
#' @export
summary.mcml <- function(object, ...) {
  .mcml_check_unused_dots("summary.mcml", ...)
  as_mat <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.matrix(x) && is.numeric(x)) return(x)
    if (is.list(x) && is.numeric(x$weights)) return(x$weights)
    if (is.list(x) && is.numeric(x$matrix))  return(x$matrix)
    NULL
  }

  macro    <- as_mat(object$macro)
  cluster_mats <- lapply(object$clusters, as_mat)
  # Clusters are named by the partition, which exists even when no
  # within-cluster layers were computed (compute_within = FALSE).
  cluster_names <- names(object$cluster_members) %||% rownames(macro) %||%
    names(cluster_mats)
  sizes <- object$meta$cluster_sizes
  if (is.null(sizes)) sizes <- vapply(object$cluster_members, length, integer(1L))

  # Total within-cluster weight, self-transitions included (they are part of
  # the macro diagonal); NA when the layers were not computed.
  within_total <- vapply(cluster_names, function(nm) {
    mat <- cluster_mats[[nm]]
    if (is.null(mat)) NA_real_ else sum(abs(mat))
  }, numeric(1L))

  directed <- isTRUE(object$meta$directed)
  if (!is.null(macro)) {
    diag0 <- macro
    diag(diag0) <- 0
    if (directed) {
      between_out <- rowSums(abs(diag0))
      between_in  <- colSums(abs(diag0))
    } else {
      # For undirected macro networks, in/out split is not meaningful;
      # report total incident weight under between_out and leave between_in NA.
      between_out <- rowSums(abs(diag0))
      between_in  <- rep(NA_real_, length(cluster_names))
    }
  } else {
    between_out <- rep(NA_real_, length(cluster_names))
    between_in  <- rep(NA_real_, length(cluster_names))
  }

  data.frame(
    cluster      = cluster_names,
    size         = as.integer(sizes),
    within_total = as.numeric(within_total),
    between_out  = as.numeric(between_out),
    between_in   = as.numeric(between_in),
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
}

.mcml_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(paste0(
    method, "() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", ")),
    class = "nestimate_bad_input", call = NULL))
}
