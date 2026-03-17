# ---- Multi-Cluster Multi-Layer Networks ----

#' Build Multi-Cluster Multi-Layer Network
#'
#' @description
#' Takes sequence data, a pre-computed matrix, or a network object, plus a
#' cluster mapping, and builds three levels of network:
#' \enumerate{
#'   \item Full node-level transition network
#'   \item Between-cluster network (states recoded to cluster labels)
#'   \item Within-cluster networks (one per cluster, non-member states filtered)
#' }
#'
#' When raw sequence data is provided, between and within networks are estimated
#' directly from the sequences using \code{\link{build_network}} — not by
#' post-hoc aggregation of the node-level matrix. This preserves actual
#' transition patterns including self-loops.
#'
#' Supports 6 data input formats and 5 cluster formats (see Parameters).
#'
#' @param data One of:
#'   \itemize{
#'     \item Data frame of sequence data (wide or long format)
#'     \item Pre-computed square matrix with matching row/column names
#'     \item Edge list data.frame with \code{from}, \code{to}, \code{weight}
#'     \item A \code{tna} object (from the tna package)
#'     \item A \code{cograph_network} object
#'     \item A \code{netobject} (from \code{\link{build_network}})
#'   }
#' @param clusters Cluster mapping in one of these formats:
#'   \itemize{
#'     \item Named list: \code{list(Group1 = c("A", "B"), Group2 = c("C", "D"))}
#'     \item Named vector: \code{c(A = "Group1", B = "Group1", C = "Group2")}
#'     \item Plain vector (length == n_nodes): \code{c("G1", "G1", "G2", "G2")}
#'     \item Data.frame with 2 columns (node, group)
#'     \item \code{NULL}: auto-detect from cograph_network node attributes
#'   }
#' @param method Character. Estimation method for the node-level network.
#'   One of \code{"relative"} (default), \code{"frequency"}, or
#'   \code{"co_occurrence"}. Aliases accepted. Only used when \code{data} is
#'   a raw data frame.
#' @param aggregation Character. How to aggregate edge weights when building
#'   between/within networks from a pre-computed matrix (not used for sequence
#'   data). One of \code{"sum"}, \code{"mean"} (default), \code{"median"},
#'   \code{"max"}, \code{"min"}.
#' @param params Named list. Additional parameters passed to the estimator.
#' @param scaling Character vector or NULL. Post-estimation scaling.
#' @param threshold Numeric. Edges below this absolute value are zeroed.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return An object of class \code{"mcml_network"} containing:
#' \describe{
#'   \item{matrix}{The full node-level transition matrix.}
#'   \item{between}{Between-cluster transition matrix.}
#'   \item{within}{Named list of within-cluster transition matrices.}
#'   \item{edges}{Data.frame with from, to, weight, cluster_from,
#'     cluster_to, type columns.}
#'   \item{clusters}{The cluster mapping (named list).}
#'   \item{cluster_summary}{The \code{cluster_summary} object from cograph
#'     (for \code{plot_mcml} compatibility).}
#'   \item{method}{The estimation method used.}
#'   \item{aggregation}{The aggregation method used.}
#'   \item{nodes}{Character vector of all node names.}
#'   \item{cluster_names}{Character vector of cluster names.}
#'   \item{n_nodes}{Total number of nodes.}
#'   \item{n_clusters}{Number of clusters.}
#' }
#'
#' @examples
#' \dontrun{
#' # From wide sequence data
#' mc <- build_mcml(
#'   group_regulation,
#'   clusters = list(
#'     Regulation = c("plan", "adapt", "monitor"),
#'     Social = c("discuss", "consensus", "cohesion"),
#'     Emotion = c("emotion", "coregulate", "synthesis")
#'   )
#' )
#' mc
#' mc$between
#' mc$within$Social
#' summary(mc)
#' plot(mc)
#'
#' # From pre-computed matrix
#' sim <- simulate_mtna(n_nodes = 4, n_types = 3, seed = 42)
#' mc <- build_mcml(sim$matrix, clusters = sim$node_types)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_mcml}}
#'
#' @export
build_mcml <- function(data,
                 clusters = NULL,
                 method = "relative",
                 aggregation = "mean",
                 params = list(),
                 scaling = NULL,
                 threshold = 0,
                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # --- Detect whether input is raw sequences or pre-computed ---
  is_sequences <- .is_sequence_data(data)

  # --- Resolve data to node-level matrix ---
  mat <- .resolve_mcml_data(data, method, params, scaling, threshold)

  # --- Resolve clusters ---
  clusters <- .resolve_clusters(clusters, rownames(mat), data)

  # --- Validate clusters ---
  stopifnot(is.list(clusters), length(clusters) >= 2)
  if (is.null(names(clusters)) || any(names(clusters) == "")) {
    stop("'clusters' must be a named list.", call. = FALSE)
  }

  all_nodes <- unlist(clusters, use.names = FALSE)
  if (anyDuplicated(all_nodes)) {
    dups <- all_nodes[duplicated(all_nodes)]
    stop("Nodes appear in multiple clusters: ",
         paste(unique(dups), collapse = ", "), call. = FALSE)
  }

  missing <- setdiff(all_nodes, rownames(mat))
  if (length(missing) > 0) {
    stop("Nodes not found in network: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # --- Build between and within networks ---
  if (is_sequences) {
    # Estimate directly from sequences
    bw <- .build_between_from_sequences(data, clusters, method, params,
                                         scaling, threshold)
    wn <- .build_within_from_sequences(data, clusters, method, params,
                                        scaling, threshold)
  } else {
    # Aggregate from pre-computed matrix
    bw <- .build_between_from_matrix(mat, clusters, aggregation)
    wn <- .build_within_from_matrix(mat, clusters, aggregation)
  }

  # --- cograph cluster_summary for plot_mcml compatibility ---
  cs <- tryCatch(
    cograph::cluster_summary(
      mat, clusters = clusters, method = aggregation,
      type = "tna", compute_within = TRUE
    ),
    error = function(e) NULL
  )

  # --- Build enriched edge list ---
  edges <- .build_mcml_edges(mat, clusters)

  # --- Build result ---
  result <- list(
    matrix = mat,
    between = bw,
    within = wn,
    edges = edges,
    clusters = clusters,
    cluster_summary = cs,
    method = method,
    aggregation = aggregation,
    nodes = rownames(mat),
    cluster_names = names(clusters),
    n_nodes = nrow(mat),
    n_clusters = length(clusters)
  )
  class(result) <- "mcml_network"
  result
}


# ---- Input type detection ----

#' Check if data is raw sequence data (vs pre-computed network)
#' @noRd
.is_sequence_data <- function(data) {
  if (inherits(data, "netobject")) return(FALSE)
  if (inherits(data, "tna")) return(FALSE)
  if (inherits(data, "cograph_network")) return(FALSE)
  if (is.matrix(data)) return(FALSE)
  # Edge list data.frame
  if (is.data.frame(data) &&
      all(c("from", "to", "weight") %in% names(data))) {
    return(FALSE)
  }
  # Remaining data.frames are sequence data
  if (is.data.frame(data)) return(TRUE)
  FALSE
}


# ---- Build between/within from sequences ----

#' Build between-cluster network by recoding states to cluster labels
#'
#' Recodes all states to their cluster labels and estimates transitions
#' via build_network(). This preserves all transitions including self-loops
#' (same-cluster to same-cluster), matching TNA behavior exactly.
#' @noRd
.build_between_from_sequences <- function(data, clusters, method, params,
                                           scaling, threshold) {
  recoded <- .recode_to_clusters(data, clusters, params)
  net <- build_network(recoded, method = method, params = params,
                       scaling = scaling, threshold = threshold)
  net$matrix
}

#' Build within-cluster networks by filtering to each cluster's states
#' @noRd
.build_within_from_sequences <- function(data, clusters, method, params,
                                          scaling, threshold) {
  result <- lapply(clusters, function(states) {
    if (length(states) < 2L) return(NULL)
    filtered <- .filter_to_cluster(data, states, params)
    tryCatch({
      net <- build_network(filtered, method = method, params = params,
                           scaling = scaling, threshold = threshold)
      net$matrix
    }, error = function(e) NULL)
  })
  # Remove NULLs (single-node or failed)
  result[!vapply(result, is.null, logical(1))]
}


# ---- Build between/within from pre-computed matrix ----

#' Build between-cluster matrix by aggregating node-level blocks
#'
#' Aggregates all blocks including diagonal (self-loops within same cluster).
#' Matches TNA behavior: all transitions are counted.
#' @noRd
.build_between_from_matrix <- function(mat, clusters, aggregation) {
  cl_names <- names(clusters)
  n_cl <- length(cl_names)
  bw <- matrix(0, n_cl, n_cl, dimnames = list(cl_names, cl_names))

  agg_fn <- switch(aggregation,
    sum = sum,
    mean = mean,
    median = stats::median,
    max = max,
    min = min,
    mean  # default
  )

  for (i in seq_len(n_cl)) {
    for (j in seq_len(n_cl)) {
      block <- mat[clusters[[i]], clusters[[j]], drop = FALSE]
      vals <- as.vector(block)
      vals <- vals[vals != 0]
      bw[i, j] <- if (length(vals) > 0) agg_fn(vals) else 0
    }
  }

  # Row-normalize
  rs <- rowSums(bw)
  nz <- rs > 0
  bw[nz, ] <- bw[nz, ] / rs[nz]
  bw
}

#' Build within-cluster matrices by extracting submatrices
#' @noRd
.build_within_from_matrix <- function(mat, clusters, aggregation) {
  result <- lapply(clusters, function(states) {
    if (length(states) < 2L) return(NULL)
    sub <- mat[states, states, drop = FALSE]
    # Row-normalize
    rs <- rowSums(sub)
    nz <- rs > 0
    sub[nz, ] <- sub[nz, ] / rs[nz]
    sub
  })
  result[!vapply(result, is.null, logical(1))]
}


# ---- Input resolution helpers ----

#' Resolve data input to a square matrix
#' @noRd
.resolve_mcml_data <- function(data, method, params, scaling, threshold) {
  # 1. netobject
  if (inherits(data, "netobject")) {
    return(data$matrix)
  }

  # 2. tna object
  if (inherits(data, "tna")) {
    return(data$weights)
  }

  # 3. cograph_network
  if (inherits(data, "cograph_network")) {
    return(cograph::to_matrix(data))
  }

  # 4. Pre-computed square matrix
  if (is.matrix(data) && nrow(data) == ncol(data) &&
      !is.null(rownames(data)) && !is.null(colnames(data)) &&
      setequal(rownames(data), colnames(data))) {
    return(data)
  }

  # 5. Edge list data.frame (has from, to, weight)
  if (is.data.frame(data) &&
      all(c("from", "to", "weight") %in% names(data))) {
    g <- cograph::as_cograph(data, directed = TRUE)
    return(cograph::to_matrix(g))
  }

  # 6. Wide/long sequence data.frame — build via build_network
  if (is.data.frame(data)) {
    net <- build_network(data, method = method, params = params,
                         scaling = scaling, threshold = threshold)
    return(net$matrix)
  }

  stop("Unsupported data input type for build_mcml().", call. = FALSE)
}


#' Resolve cluster input to a named list
#'
#' @param clusters Raw cluster input (list, vector, data.frame, or NULL).
#' @param node_names Character vector of node names from the matrix.
#' @param x Original data object (for auto-detection from cograph_network).
#' @return Named list mapping cluster names to node name vectors.
#' @noRd
.resolve_clusters <- function(clusters, node_names, x = NULL) {
  if (is.null(clusters)) {
    # Auto-detect from cograph_network
    if (inherits(x, "cograph_network") && !is.null(x$nodes)) {
      cl_col <- intersect(
        c("cluster", "clusters", "group", "groups"),
        names(x$nodes)
      )
      if (length(cl_col) > 0) {
        node_labels <- x$nodes$label %||% x$nodes$name
        grp_vals <- x$nodes[[cl_col[1]]]
        return(split(node_labels, grp_vals))
      }
    }
    stop("'clusters' is NULL and cannot be auto-detected from data.",
         call. = FALSE)
  }

  # Named list — use directly
  if (is.list(clusters) && !is.data.frame(clusters)) {
    return(clusters)
  }

  # Data.frame with 2 columns (node, group)
  if (is.data.frame(clusters)) {
    stopifnot(ncol(clusters) >= 2)
    return(split(as.character(clusters[[1]]), clusters[[2]]))
  }

  # Named vector (names = nodes, values = group labels)
  if ((is.character(clusters) || is.factor(clusters)) && !is.null(names(clusters))) {
    return(split(names(clusters), as.character(clusters)))
  }

  # Plain vector (length == n_nodes, values = group labels)
  if ((is.character(clusters) || is.factor(clusters)) && is.null(names(clusters))) {
    if (length(clusters) != length(node_names)) {
      stop("Plain cluster vector length (", length(clusters),
           ") must match number of nodes (", length(node_names), ").",
           call. = FALSE)
    }
    return(split(node_names, as.character(clusters)))
  }

  stop("Unsupported cluster format. Use a named list, named vector, ",
       "plain vector, or data.frame.", call. = FALSE)
}


#' Build enriched edge list from node-level matrix + cluster mapping
#' @noRd
.build_mcml_edges <- function(mat, clusters) {
  # Build node -> cluster lookup
  node_to_cluster <- character(0)
  for (cl_name in names(clusters)) {
    nodes <- clusters[[cl_name]]
    node_to_cluster[nodes] <- cl_name
  }

  # Extract all non-zero, non-diagonal edges
  idx <- which(mat != 0 & row(mat) != col(mat), arr.ind = TRUE)

  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0), weight = numeric(0),
      cluster_from = character(0), cluster_to = character(0),
      type = character(0), stringsAsFactors = FALSE
    ))
  }

  nms <- rownames(mat)
  from_nodes <- nms[idx[, 1]]
  to_nodes <- nms[idx[, 2]]
  cf <- node_to_cluster[from_nodes]
  ct <- node_to_cluster[to_nodes]

  edge_type <- ifelse(is.na(cf) | is.na(ct), NA_character_,
                       ifelse(cf == ct, "within", "between"))

  data.frame(
    from = from_nodes,
    to = to_nodes,
    weight = mat[idx],
    cluster_from = cf,
    cluster_to = ct,
    type = edge_type,
    stringsAsFactors = FALSE
  )
}


# ---- S3 methods for mcml_network ----

#' @export
print.mcml_network <- function(x, ...) {
  cat("Multi-Cluster Multi-Layer Network\n")
  cat(strrep("-", 38), "\n")
  cat("Nodes:      ", x$n_nodes, "\n")
  cat("Clusters:   ", x$n_clusters, " (",
      paste(x$cluster_names, collapse = ", "), ")\n")
  cat("Method:     ", x$method, "\n")
  cat("Aggregation:", x$aggregation, "\n\n")

  sizes <- vapply(x$clusters, length, integer(1))
  cat("Cluster sizes:\n")
  for (i in seq_along(x$clusters)) {
    cat("  ", x$cluster_names[i], ": ", sizes[i], " nodes (",
        paste(x$clusters[[i]], collapse = ", "), ")\n", sep = "")
  }

  cat("\nBetween-cluster weights:\n")
  print(round(x$between, 3))
  invisible(x)
}


#' @export
summary.mcml_network <- function(object, ...) {
  cat("Multi-Cluster Multi-Layer Network Summary\n")
  cat(strrep("=", 42), "\n\n")

  # Between-cluster summary
  cat("BETWEEN-CLUSTER NETWORK\n")
  cat(strrep("-", 24), "\n")
  bw <- object$between
  cat("Transition matrix (", nrow(bw), "x", ncol(bw), "):\n")
  print(round(bw, 3))

  # Strongest between-cluster transitions
  bw_nodiag <- bw
  diag(bw_nodiag) <- 0
  top_idx <- order(as.vector(bw_nodiag), decreasing = TRUE)
  n_show <- min(5, sum(bw_nodiag > 0))
  if (n_show > 0) {
    cat("\nStrongest transitions:\n")
    for (k in seq_len(n_show)) {
      i <- ((top_idx[k] - 1) %% nrow(bw)) + 1
      j <- ((top_idx[k] - 1) %/% nrow(bw)) + 1
      cat(sprintf("  %s -> %s: %.3f\n",
                  rownames(bw)[i], colnames(bw)[j], bw[i, j]))
    }
  }

  # Within-cluster summaries
  cat("\n\nWITHIN-CLUSTER NETWORKS\n")
  cat(strrep("-", 24), "\n")
  for (cl_name in names(object$within)) {
    w <- object$within[[cl_name]]
    n_edges <- sum(w > 0 & row(w) != col(w))
    cat(sprintf("\n%s (%d nodes, %d edges):\n",
                cl_name, nrow(w), n_edges))
    print(round(w, 3))
  }

  invisible(object)
}


#' @export
plot.mcml_network <- function(x, type = c("mcml", "between", "within"), ...) {
  type <- match.arg(type)

  if (type == "mcml") {
    if (!is.null(x$cluster_summary)) {
      cograph::plot_mcml(x$cluster_summary, ...)
    } else {
      message("plot_mcml requires cograph::cluster_summary (not available). ",
              "Falling back to between plot.")
      cograph::splot(x$between, ...)
    }
  } else if (type == "between") {
    cograph::splot(x$between, ...)
  } else if (type == "within") {
    n_within <- length(x$within)
    n_col <- min(n_within, 3)
    n_row <- ceiling(n_within / n_col)
    old_par <- graphics::par(mfrow = c(n_row, n_col))
    on.exit(graphics::par(old_par), add = TRUE)
    for (cl_name in names(x$within)) {
      cograph::splot(x$within[[cl_name]],
                     title = cl_name, ...)
    }
  }

  invisible(x)
}


# ============================================================================
# bootstrap_mcml()
# ============================================================================

#' Bootstrap Multi-Cluster Multi-Layer Network
#'
#' @description
#' Performs two independent bootstraps using \code{\link{bootstrap_network}}:
#' \enumerate{
#'   \item \strong{Between-cluster}: Recodes states to cluster labels, then
#'     bootstraps the cluster-level transition network.
#'   \item \strong{Within-cluster}: For each cluster with >= 2 nodes, filters
#'     data to keep only that cluster's states (others become NA), then
#'     bootstraps the within-cluster transition network.
#' }
#'
#' Requires raw wide-format sequence data (rows are sequences to resample).
#'
#' @param data Data frame of wide-format sequence data.
#' @param clusters Cluster mapping (any format accepted by
#'   \code{\link{build_mcml}}).
#' @param method Character. Estimation method (default: \code{"relative"}).
#' @param params Named list. Additional parameters for the estimator.
#' @param iter Integer. Number of bootstrap iterations (default: 1000).
#' @param ci_level Numeric. Significance level (default: 0.05).
#' @param inference Character. \code{"stability"} or \code{"threshold"}.
#' @param consistency_range Numeric vector of length 2 for stability inference.
#' @param seed Integer or NULL. Random seed.
#'
#' @return An object of class \code{"mcml_bootstrap"} containing:
#' \describe{
#'   \item{between}{\code{saqr_bootstrap} object for between-cluster network.}
#'   \item{within}{Named list of \code{saqr_bootstrap} objects (one per cluster
#'     with >= 2 nodes).}
#'   \item{original}{\code{mcml_network} object.}
#'   \item{clusters, method, iter, ci_level, inference}{Configuration fields.}
#' }
#'
#' @examples
#' \dontrun{
#' boot <- bootstrap_mcml(
#'   group_regulation,
#'   clusters = list(
#'     Social = c("discuss", "consensus", "cohesion"),
#'     Cognitive = c("plan", "adapt", "monitor"),
#'     Meta = c("emotion", "coregulate", "synthesis")
#'   ),
#'   iter = 500, seed = 42
#' )
#' print(boot)
#' summary(boot, level = "between")
#' }
#'
#' @seealso \code{\link{build_mcml}}, \code{\link{bootstrap_network}}
#'
#' @export
bootstrap_mcml <- function(data,
                           clusters,
                           method = "relative",
                           params = list(),
                           iter = 1000L,
                           ci_level = 0.05,
                           inference = "stability",
                           consistency_range = c(0.75, 1.25),
                           seed = NULL) {
  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("bootstrap_mcml() requires raw wide-format sequence data (data.frame).",
         call. = FALSE)
  }
  iter <- as.integer(iter)
  inference <- match.arg(inference, c("stability", "threshold"))

  # Build the original mcml to resolve clusters and get the full network
  original <- build_mcml(
    data, clusters = clusters, method = method,
    params = params, seed = seed
  )
  clusters <- original$clusters

  if (length(clusters) < 2) {
    stop("At least 2 clusters are required for bootstrap_mcml().", call. = FALSE)
  }

  # --- Between-cluster bootstrap ---
  between_data <- .recode_to_clusters(data, clusters, params)
  between_net <- build_network(between_data, method = method, params = params)
  boot_between <- bootstrap_network(
    between_net,
    iter = iter,
    ci_level = ci_level,
    inference = inference,
    consistency_range = consistency_range,
    seed = seed
  )

  # --- Within-cluster bootstraps ---
  multi_node_clusters <- clusters[vapply(clusters, length, integer(1)) >= 2L]

  boot_within <- lapply(multi_node_clusters, function(states) {
    filtered <- .filter_to_cluster(data, states, params)
    tryCatch(
      bootstrap_network(
        build_network(filtered, method = method, params = params),
        iter = iter,
        ci_level = ci_level,
        inference = inference,
        consistency_range = consistency_range
      ),
      error = function(e) NULL
    )
  })
  boot_within <- boot_within[!vapply(boot_within, is.null, logical(1))]

  # --- Assemble result ---
  result <- list(
    between = boot_between,
    within = boot_within,
    original = original,
    clusters = clusters,
    method = method,
    iter = iter,
    ci_level = ci_level,
    inference = inference,
    consistency_range = consistency_range
  )
  class(result) <- "mcml_bootstrap"
  result
}


# ---- Bootstrap helpers ----

#' Recode sequence data to cluster labels
#' @noRd
.recode_to_clusters <- function(data, clusters, params = list()) {
  lookup <- character(0)
  for (cl_name in names(clusters)) {
    for (state in clusters[[cl_name]]) {
      lookup[state] <- cl_name
    }
  }

  id_col <- params$id %||% params$id_col
  cols <- params$cols
  state_cols <- .select_state_cols(data, id_col, cols)

  result <- data
  result[, state_cols] <- lapply(data[, state_cols, drop = FALSE], function(col) {
    recoded <- lookup[as.character(col)]
    recoded[is.na(recoded) & !is.na(col)] <- NA_character_
    recoded
  })
  result
}


#' Filter data to keep only a cluster's states (others -> NA)
#' @noRd
.filter_to_cluster <- function(data, states, params = list()) {
  id_col <- params$id %||% params$id_col
  cols <- params$cols
  state_cols <- .select_state_cols(data, id_col, cols)

  result <- data
  result[, state_cols] <- lapply(data[, state_cols, drop = FALSE], function(col) {
    vals <- as.character(col)
    vals[!vals %in% states] <- NA_character_
    vals
  })
  result
}


# ---- S3 methods for mcml_bootstrap ----

#' @export
print.mcml_bootstrap <- function(x, ...) {
  cat("MCML Bootstrap\n")
  cat(strrep("-", 30), "\n")
  cat(sprintf("Iterations: %d  |  CI level: %.2f  |  Inference: %s\n",
              x$iter, x$ci_level, x$inference))
  cat(sprintf("Clusters: %d (%s)\n",
              length(x$clusters),
              paste(names(x$clusters), collapse = ", ")))

  cat("\nBetween-cluster:\n")
  cat(sprintf("  Original edges: %d  |  Significant: %d\n",
              x$between$original$n_edges, x$between$model$n_edges))

  if (length(x$within) > 0) {
    cat("\nWithin-cluster:\n")
    for (cl_name in names(x$within)) {
      w <- x$within[[cl_name]]
      cat(sprintf("  %s: %d original, %d significant\n",
                  cl_name, w$original$n_edges, w$model$n_edges))
    }
  }

  skipped <- setdiff(names(x$clusters), names(x$within))
  single <- skipped[vapply(x$clusters[skipped], length, integer(1)) < 2L]
  if (length(single) > 0) {
    cat(sprintf("\nSkipped (single-node): %s\n",
                paste(single, collapse = ", ")))
  }

  invisible(x)
}


#' @export
summary.mcml_bootstrap <- function(object,
                                    level = c("between", "within"),
                                    cluster = NULL, ...) {
  level <- match.arg(level)

  if (level == "between") {
    return(summary(object$between))
  }

  if (is.null(cluster)) {
    return(lapply(object$within, summary))
  }

  if (!cluster %in% names(object$within)) {
    stop("Cluster '", cluster, "' not found in within-cluster bootstraps.",
         call. = FALSE)
  }
  summary(object$within[[cluster]])
}


#' @export
plot.mcml_bootstrap <- function(x,
                                 type = c("between", "within"),
                                 cluster = NULL, ...) {
  type <- match.arg(type)

  if (type == "between") {
    plot(x$between, ...)
  } else {
    if (is.null(cluster)) {
      if (length(x$within) == 0) {
        stop("No within-cluster bootstraps available.", call. = FALSE)
      }
      cluster <- names(x$within)[1]
      message("Plotting within-cluster bootstrap for: ", cluster)
    }
    if (!cluster %in% names(x$within)) {
      stop("Cluster '", cluster, "' not found.", call. = FALSE)
    }
    plot(x$within[[cluster]], ...)
  }

  invisible(x)
}
