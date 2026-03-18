# Cluster Metrics for Network Analysis
# Summary measures for between/within clusters and multilayer networks

# ==============================================================================
# 1. Edge Weight Aggregation
# ==============================================================================

#' Aggregate Edge Weights
#'
#' Aggregates a vector of edge weights using various methods.
#' Compatible with igraph's edge.attr.comb parameter.
#'
#' @param w Numeric vector of edge weights
#' @param method Aggregation method: "sum", "mean", "median", "max", "min",
#'   "prod", "density", "geomean"
#' @param n_possible Number of possible edges (for density calculation)
#' @return Single aggregated value
#' @export
#' @examples
#' w <- c(0.5, 0.8, 0.3, 0.9)
#' aggregate_weights(w, "sum")   # 2.5
#' aggregate_weights(w, "mean")  # 0.625
#' aggregate_weights(w, "max")   # 0.9
aggregate_weights <- function(w, method = "sum", n_possible = NULL) {
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

#' @rdname aggregate_weights
#' @export
wagg <- aggregate_weights

# ==============================================================================
# 2. Cluster Summary (Between/Within Aggregates)
# ==============================================================================

#' Cluster Summary Statistics
#'
#' Aggregates node-level network weights to cluster-level summaries. Computes
#' both between-cluster transitions (how clusters connect to each other) and
#' within-cluster transitions (how nodes connect within each cluster).
#'
#' This is the core function for Multi-Cluster Multi-Level (MCML) analysis.
#' Use \code{as_tna()} to convert results to tna objects for further
#' analysis with the tna package.
#'
#' @param x Network input. Accepts multiple formats:
#'   \describe{
#'     \item{matrix}{Numeric adjacency/weight matrix. Row and column names are
#'       used as node labels. Values represent edge weights (e.g., transition
#'       counts, co-occurrence frequencies, or probabilities).}
#'     \item{netobject}{A cograph network object. The function extracts
#'       the weight matrix from \code{x$weights} or converts via
#'       \code{to_matrix()}. Clusters can be auto-detected from node attributes.}
#'     \item{tna}{A tna object from the tna package. Extracts \code{x$weights}.}
#'     \item{cluster_summary}{If already a cluster_summary, returns unchanged.}
#'   }
#'
#' @param clusters Cluster/group assignments for nodes. Accepts multiple formats:
#'   \describe{
#'     \item{NULL}{(default) Auto-detect from netobject. Looks for columns
#'       named 'clusters', 'cluster', 'groups', or 'group' in \code{x$nodes}.
#'       Throws an error if no cluster column is found.
#'       This option only works when \code{x} is a netobject.}
#'     \item{vector}{Cluster membership for each node, in the same order as the

#'       matrix rows/columns. Can be numeric (1, 2, 3) or character ("A", "B").
#'       Cluster names will be derived from unique values.
#'       Example: \code{c(1, 1, 2, 2, 3, 3)} assigns first two nodes to cluster 1.}
#'     \item{data.frame}{A data frame where the first column contains node names
#'       and the second column contains group/cluster names.
#'       Example: \code{data.frame(node = c("A", "B", "C"), group = c("G1", "G1", "G2"))}}
#'     \item{named list}{Explicit mapping of cluster names to node labels.
#'       List names become cluster names, values are character vectors of node
#'       labels that must match matrix row/column names.
#'       Example: \code{list(Alpha = c("A", "B"), Beta = c("C", "D"))}}
#'   }
#'
#' @param method Aggregation method for combining edge weights within/between
#'   clusters. Controls how multiple node-to-node edges are summarized:
#'   \describe{
#'     \item{"sum"}{(default) Sum of all edge weights. Best for count data
#'       (e.g., transition frequencies). Preserves total flow.}
#'     \item{"mean"}{Average edge weight. Best when cluster sizes differ and
#'       you want to control for size. Note: when input is already a transition
#'       matrix (rows sum to 1), "mean" avoids size bias.
#'       Example: cluster with 5 nodes won't have 5x the weight of cluster with 1 node.}
#'     \item{"median"}{Median edge weight. Robust to outliers.}
#'     \item{"max"}{Maximum edge weight. Captures strongest connection.}
#'     \item{"min"}{Minimum edge weight. Captures weakest connection.}
#'     \item{"density"}{Sum divided by number of possible edges. Normalizes
#'       by cluster size combinations.}
#'     \item{"geomean"}{Geometric mean of positive weights. Useful for
#'       multiplicative processes.}
#'   }
#'
#' @param type Post-processing applied to aggregated weights. Determines the
#'   interpretation of the resulting matrices:
#'   \describe{
#'     \item{"tna"}{(default) Row-normalize so each row sums to 1. Creates
#'       transition probabilities suitable for Markov chain analysis.
#'       Interpretation: "Given I'm in cluster A, what's the probability
#'       of transitioning to cluster B?"
#'       Required for use with tna package functions.
#'       Diagonal represents within-cluster transition probability.}
#'     \item{"raw"}{No normalization. Returns aggregated counts/weights as-is.
#'       Use for frequency analysis or when you need raw counts.
#'       Compatible with igraph's contract + simplify output.}
#'     \item{"cooccurrence"}{Symmetrize the matrix: (A + t(A)) / 2.
#'       For undirected co-occurrence analysis.}
#'     \item{"semi_markov"}{Row-normalize with duration weighting.
#'       For semi-Markov process analysis.}
#'   }
#'
#' @param directed Logical. If \code{TRUE} (default), treat network as directed.
#'   A->B and B->A are separate edges. If \code{FALSE}, edges are undirected
#'   and the matrix is symmetrized before processing.
#'
#' @param compute_within Logical. If \code{TRUE} (default), compute within-cluster
#'   transition matrices for each cluster. Each cluster gets its own n_i x n_i
#'   matrix showing internal node-to-node transitions.
#'   Set to \code{FALSE} to skip this computation for better performance when
#'   only between-cluster summary is needed.
#'
#' @return A \code{cluster_summary} object (S3 class) containing:
#'   \describe{
#'     \item{between}{List with two elements:
#'       \describe{
#'         \item{weights}{k x k matrix of cluster-to-cluster weights, where k is
#'           the number of clusters. Row i, column j contains the aggregated
#'           weight from cluster i to cluster j. Diagonal contains within-cluster
#'           totals. Processing depends on \code{type}.}
#'         \item{inits}{Numeric vector of length k. Initial state distribution
#'           across clusters, computed from column sums of the original matrix.
#'           Represents the proportion of incoming edges to each cluster.}
#'       }
#'     }
#'     \item{within}{Named list with one element per cluster. Each element contains:
#'       \describe{
#'         \item{weights}{n_i x n_i matrix for nodes within that cluster.
#'           Shows internal transitions between nodes in the same cluster.}
#'         \item{inits}{Initial distribution within the cluster.}
#'       }
#'       NULL if \code{compute_within = FALSE}.}
#'     \item{clusters}{Named list mapping cluster names to their member node labels.
#'       Example: \code{list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))}}
#'     \item{meta}{List of metadata:
#'       \describe{
#'         \item{type}{The \code{type} argument used ("tna", "raw", etc.)}
#'         \item{method}{The \code{method} argument used ("sum", "mean", etc.)}
#'         \item{directed}{Logical, whether network was treated as directed}
#'         \item{n_nodes}{Total number of nodes in original network}
#'         \item{n_clusters}{Number of clusters}
#'         \item{cluster_sizes}{Named vector of cluster sizes}
#'       }
#'     }
#'   }
#'
#' @details
#' ## Workflow
#'
#' Typical MCML analysis workflow:
#' \preformatted{
#' # 1. Create network
#' net <- build_network(data, method = "relative")
#' net$nodes$clusters <- group_assignments
#'
#' # 2. Compute cluster summary
#' cs <- cluster_summary(net, type = "tna")
#'
#' # 3. Convert to tna models
#' tna_models <- as_tna(cs)
#'
#' # 4. Analyze/visualize
#' plot(tna_models$macro)
#' tna::centralities(tna_models$macro)
#' }
#'
#' ## Between-Cluster Matrix Structure
#'
#' The \code{macro$weights} matrix has clusters as both rows and columns:
#' \itemize{
#'   \item Off-diagonal (row i, col j): Aggregated weight from cluster i to cluster j
#'   \item Diagonal (row i, col i): Within-cluster total (sum of internal edges in cluster i)
#' }
#'
#' When \code{type = "tna"}, rows sum to 1 and diagonal values represent
#' "retention rate" - the probability of staying within the same cluster.
#'
#' ## Choosing method and type
#'
#' \tabular{lll}{
#'   \strong{Input data} \tab \strong{Recommended} \tab \strong{Reason} \cr
#'   Edge counts \tab method="sum", type="tna" \tab Preserves total flow, normalizes to probabilities \cr
#'   Transition matrix \tab method="mean", type="tna" \tab Avoids cluster size bias \cr
#'   Frequencies \tab method="sum", type="raw" \tab Keep raw counts for analysis \cr
#'   Correlation matrix \tab method="mean", type="raw" \tab Average correlations \cr
#' }
#'
#' @export
#' @seealso
#'   \code{as_tna()} to convert results to tna objects,
#'   \code{plot()} for two-layer visualization,
#'   \code{plot()} for flat cluster visualization
#'
#' @examples
#' # -----------------------------------------------------
#' # Basic usage with matrix and cluster vector
#' # -----------------------------------------------------
#' mat <- matrix(runif(100), 10, 10)
#' rownames(mat) <- colnames(mat) <- LETTERS[1:10]
#'
#' clusters <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
#' cs <- cluster_summary(mat, clusters)
#'
#' # Access results
#' cs$macro$weights    # 3x3 cluster transition matrix
#' cs$macro$inits      # Initial distribution
#' cs$clusters$`1`$weights # Within-cluster 1 transitions
#' cs$meta               # Metadata
#'
#' # -----------------------------------------------------
#' # Named list clusters (more readable)
#' # -----------------------------------------------------
#' clusters <- list(
#'   Alpha = c("A", "B", "C"),
#'   Beta = c("D", "E", "F"),
#'   Gamma = c("G", "H", "I", "J")
#' )
#' cs <- cluster_summary(mat, clusters, type = "tna")
#' cs$macro$weights    # Rows/cols named Alpha, Beta, Gamma
#' cs$clusters$Alpha       # Within Alpha cluster
#'
#' # -----------------------------------------------------
#' # Auto-detect clusters from netobject
#' # -----------------------------------------------------
#' \dontrun{
#' net <- build_network(mat, method = "relative")
#' net$nodes$clusters <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
#' cs <- cluster_summary(net)  # No clusters argument needed
#' }
#'
#' # -----------------------------------------------------
#' # Different aggregation methods
#' # -----------------------------------------------------
#' cs_sum <- cluster_summary(mat, clusters, method = "sum")   # Total flow
#' cs_mean <- cluster_summary(mat, clusters, method = "mean") # Average
#' cs_max <- cluster_summary(mat, clusters, method = "max")   # Strongest
#'
#' # -----------------------------------------------------
#' # Raw counts vs TNA probabilities
#' # -----------------------------------------------------
#' cs_raw <- cluster_summary(mat, clusters, type = "raw")
#' cs_tna <- cluster_summary(mat, clusters, type = "tna")
#'
#' rowSums(cs_raw$macro$weights)  # Various sums
#' rowSums(cs_tna$macro$weights)  # All equal to 1
#'
#' # -----------------------------------------------------
#' # Skip within-cluster computation for speed
#' # -----------------------------------------------------
#' cs_fast <- cluster_summary(mat, clusters, compute_within = FALSE)
#' cs_fast$clusters  # NULL
#'
#' # -----------------------------------------------------
#' # Convert to tna objects for tna package
#' # -----------------------------------------------------
#' cs <- cluster_summary(mat, clusters, type = "tna")
#' tna_models <- as_tna(cs)
#' # tna_models$macro      # tna object
#' # tna_models$clusters$Alpha # tna object
cluster_summary <- function(x,
                            clusters = NULL,
                            method = c("sum", "mean", "median", "max",
                                       "min", "density", "geomean"),
                            type = c("tna", "cooccurrence", "semi_markov", "raw"),
                            directed = TRUE,
                            compute_within = TRUE) {

  # If already an mcml object, return as-is
  if (inherits(x, "mcml")) {
    return(x)
  }

  type <- match.arg(type)
  method <- match.arg(method)

  # Store original for cluster extraction
  x_orig <- x

  # Extract matrix from various input types
  if (inherits(x, "cograph_network")) {
    mat <- x$weights
  } else if (inherits(x, "netobject") || inherits(x, "netobject_ml")) {
    mat <- x$weights
  } else if (inherits(x, "tna")) {
    mat <- x$weights
  } else {
    mat <- x
  }

  if (is.null(clusters)) {
    stop("clusters argument is required", call. = FALSE)
  }

  # Validate input matrix
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("x must be a netobject, tna object, or numeric matrix", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("x must be a square matrix", call. = FALSE)
  }

  n <- nrow(mat)
  node_names <- rownames(mat)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))

  # Convert clusters to list format
  cluster_list <- .normalize_clusters(clusters, node_names)
  n_clusters <- length(cluster_list)
  cluster_names <- names(cluster_list)
  if (is.null(cluster_names)) cluster_names <- as.character(seq_len(n_clusters))
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
      between_raw[i, j] <- aggregate_weights(as.vector(w_ij), method, n_possible)
    }
  }

  # Process based on type
  between_weights <- .process_weights(between_raw, type, directed)

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
  between <- list(
    weights = between_weights,
    inits = between_inits,
    labels = cluster_names,
    data = NULL
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

        # Process based on type
        within_weights_i <- .process_weights(within_raw, type, directed)

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

      list(
        weights = within_weights_i,
        inits = within_inits_i,
        labels = cl_nodes,
        data = NULL
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
      meta = list(
        type = type,
        method = method,
        directed = directed,
        n_nodes = n,
        n_clusters = n_clusters,
        cluster_sizes = vapply(cluster_list, length, integer(1))
      )
    ),
    class = "mcml"
  )

  result
}

#' @rdname cluster_summary
#' @return See \code{\link{cluster_summary}}.
#' @export
#' @examples
#' \dontrun{
#' csum(matrix_data, clusters)
#' }
csum <- cluster_summary

# ==============================================================================
# 2b. Build MCML from Raw Transition Data
# ==============================================================================

#' Build MCML from Raw Transition Data
#'
#' Builds a Multi-Cluster Multi-Level (MCML) model from raw transition data
#' (edge lists or sequences) by recoding node labels to cluster labels and
#' counting actual transitions. Unlike \code{\link{cluster_summary}} which
#' aggregates a pre-computed weight matrix, this function works from the
#' original transition data to produce the TRUE Markov chain over cluster states.
#'
#' @param x Input data. Accepts multiple formats:
#'   \describe{
#'     \item{data.frame with from/to columns}{Edge list. Columns named
#'       from/source/src/v1/node1/i and to/target/tgt/v2/node2/j are
#'       auto-detected. Optional weight column (weight/w/value/strength).}
#'     \item{data.frame without from/to columns}{Sequence data. Each row is a
#'       sequence, columns are time steps. Consecutive pairs (t, t+1) become
#'       transitions.}
#'     \item{tna object}{If \code{x$data} is non-NULL, uses sequence path on
#'       the raw data. Otherwise falls back to \code{\link{cluster_summary}}.}
#'     \item{netobject}{If \code{x$data} is non-NULL, detects edge list
#'       vs sequence data. Otherwise falls back to \code{\link{cluster_summary}}.}
#'     \item{cluster_summary}{Returns as-is.}
#'     \item{square numeric matrix}{Falls back to \code{\link{cluster_summary}}.}
#'     \item{non-square or character matrix}{Treated as sequence data.}
#'   }
#'
#' @param clusters Cluster/group assignments. Accepts:
#'   \describe{
#'     \item{named list}{Direct mapping. List names = cluster names, values =
#'       character vectors of node labels.
#'       Example: \code{list(A = c("N1","N2"), B = c("N3","N4"))}}
#'     \item{data.frame}{A data frame where the first column contains node names
#'       and the second column contains group/cluster names.
#'       Example: \code{data.frame(node = c("N1","N2","N3"), group = c("A","A","B"))}}
#'     \item{membership vector}{Character or numeric vector. Node names are
#'       extracted from the data.
#'       Example: \code{c("A","A","B","B")}}
#'     \item{column name string}{For edge list data.frames, the name of a
#'       column containing cluster labels. The mapping is built from unique
#'       (node, group) pairs in both from and to columns.}
#'     \item{NULL}{Auto-detect from \code{netobject$nodes} or
#'       \code{$node_groups} (same logic as \code{\link{cluster_summary}}).}
#'   }
#'
#' @param method Aggregation method for combining edge weights: "sum", "mean",
#'   "median", "max", "min", "density", "geomean". Default "sum".
#' @param type Post-processing: "tna" (row-normalize), "cooccurrence"
#'   (symmetrize), "semi_markov", or "raw". Default "tna".
#' @param directed Logical. Treat as directed network? Default TRUE.
#' @param compute_within Logical. Compute within-cluster matrices? Default TRUE.
#'
#' @return A \code{cluster_summary} object with \code{meta$source = "transitions"},
#'   fully compatible with \code{plot()}, \code{as_tna()}, and
#'   \code{plot()}.
#'
#' @export
#' @seealso \code{\link{cluster_summary}} for matrix-based aggregation,
#'   \code{as_tna()} to convert to tna objects,
#'   \code{plot()} for visualization
#'
#' @examples
#' # Edge list with clusters
#' edges <- data.frame(
#'   from = c("A", "A", "B", "C", "C", "D"),
#'   to   = c("B", "C", "A", "D", "D", "A"),
#'   weight = c(1, 2, 1, 3, 1, 2)
#' )
#' clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
#' cs <- build_mcml(edges, clusters)
#' cs$macro$weights
#'
#' # Sequence data with clusters
#' seqs <- data.frame(
#'   T1 = c("A", "C", "B"),
#'   T2 = c("B", "D", "A"),
#'   T3 = c("C", "C", "D"),
#'   T4 = c("D", "A", "C")
#' )
#' cs <- build_mcml(seqs, clusters, type = "raw")
#' cs$macro$weights
build_mcml <- function(x,
                       clusters = NULL,
                       method = c("sum", "mean", "median", "max",
                                  "min", "density", "geomean"),
                       type = c("tna", "frequency", "cooccurrence",
                                "semi_markov", "raw"),
                       directed = TRUE,
                       compute_within = TRUE) {

  # If already an mcml object, return as-is
  if (inherits(x, "mcml")) {
    return(x)
  }

  type <- match.arg(type)
  method <- match.arg(method)

  # Coerce cograph_network so downstream branches see a netobject
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)

  input_type <- .detect_mcml_input(x)

  switch(input_type,
    "edgelist" = .build_mcml_edgelist(x, clusters, method, type,
                                       directed, compute_within),
    "sequence" = .build_mcml_sequence(x, clusters, method, type,
                                       directed, compute_within),
    "tna_data" = .build_mcml_sequence(x$data, clusters, method, type,
                                       directed, compute_within),
    "tna_matrix" = cluster_summary(x, clusters, method = method, type = type,
                                    directed = directed,
                                    compute_within = compute_within),
    "netobject_data" = {
      data <- x$data
      # Auto-detect clusters from network if not provided
      if (is.null(clusters)) {
        clusters <- .auto_detect_clusters(x)
      }
      sub_type <- .detect_mcml_input(data)
      if (sub_type == "edgelist") {
        .build_mcml_edgelist(data, clusters, method, type,
                              directed, compute_within)
      } else {
        .build_mcml_sequence(data, clusters, method, type,
                              directed, compute_within)
      }
    },
    "netobject_matrix" = {
      if (is.null(clusters)) {
        clusters <- .auto_detect_clusters(x)
      }
      cluster_summary(x, clusters, method = method, type = type,
                       directed = directed, compute_within = compute_within)
    },
    "matrix" = cluster_summary(x, clusters, method = method, type = type,
                                directed = directed,
                                compute_within = compute_within),
    stop("Cannot build MCML from input of class '", class(x)[1], "'",
         call. = FALSE)
  )
}

#' Detect input type for build_mcml
#' @keywords internal
.detect_mcml_input <- function(x) {
  if (inherits(x, "tna")) {
    if (!is.null(x$data)) return("tna_data")
    return("tna_matrix")
  }

  if (inherits(x, "netobject") || inherits(x, "netobject_ml")) {
    if (!is.null(x$data)) return("netobject_data")
    return("netobject_matrix")
  }

  if (is.data.frame(x)) {
    col_names <- tolower(names(x))
    from_cols <- c("from", "source", "src", "v1", "node1", "i")
    to_cols <- c("to", "target", "tgt", "v2", "node2", "j")
    has_from <- any(from_cols %in% col_names)
    has_to <- any(to_cols %in% col_names)
    if (has_from && has_to) return("edgelist")
    return("sequence")
  }

  if (is.matrix(x)) {
    if (is.numeric(x) && nrow(x) == ncol(x)) return("matrix")
    return("sequence")
  }

  "unknown"
}

#' Auto-detect clusters from netobject
#' @keywords internal
.auto_detect_clusters <- function(x) {
  clusters <- NULL
  if (!is.null(x$nodes)) {
    cluster_cols <- c("clusters", "cluster", "groups", "group")
    for (col in cluster_cols) {
      if (col %in% names(x$nodes)) {
        clusters <- x$nodes[[col]]
        break
      }
    }
  }
  if (is.null(clusters) && !is.null(x$node_groups)) {
    ng <- x$node_groups
    cluster_col <- intersect(c("cluster", "group", "layer"), names(ng))
    if (length(cluster_col) > 0) {
      clusters <- ng[[cluster_col[1]]]
    }
  }
  if (is.null(clusters)) {
    stop("No clusters found in netobject. ",
         "Add a 'clusters' column to nodes or provide clusters argument.",
         call. = FALSE)
  }
  clusters
}

#' Build node-to-cluster lookup from cluster specification
#' @keywords internal
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

#' Build cluster_summary from transition vectors
#' @keywords internal
.build_from_transitions <- function(from_nodes, to_nodes, weights,
                                     cluster_lookup, cluster_list,
                                     method, type, directed,
                                     compute_within, data = NULL) {

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

  # Include ALL transitions — node-level self-loops (A->A) are valid
  # cluster-level self-loops (e.g. discuss->discuss = Social->Social)
  b_from <- from_clusters
  b_to <- to_clusters
  b_w <- weights

  if (length(b_from) > 0) {
    # Build pair keys and aggregate
    pair_keys <- paste(b_from, b_to, sep = "\t")
    names(b_w) <- pair_keys
    agg_vals <- tapply(b_w, pair_keys, function(w) {
      n_possible <- NULL
      if (method == "density") {
        parts <- strsplit(names(w)[1], "\t")[[1]]
        n_i <- length(cluster_list[[parts[1]]])
        n_j <- length(cluster_list[[parts[2]]])
        n_possible <- n_i * n_j
      }
      aggregate_weights(w, method, n_possible)
    })

    for (key in names(agg_vals)) {
      parts <- strsplit(key, "\t")[[1]]
      between_raw[parts[1], parts[2]] <- agg_vals[[key]]
    }
  }

  # Process between weights
  between_weights <- .process_weights(between_raw, type, directed)

  # Compute inits from column sums
  col_sums <- colSums(between_raw)
  total <- sum(col_sums)
  if (total > 0) {
    between_inits <- col_sums / total
  } else {
    between_inits <- rep(1 / n_clusters, n_clusters)
  }
  names(between_inits) <- cluster_names

  # ---- Build recoded sequence data if input is sequences ----
  between_seq_data <- NULL
  within_seq_data_list <- NULL
  is_seq <- is.data.frame(data) && !any(tolower(names(data)) %in%
    c("from", "source", "src", "v1", "node1", "i",
      "to", "target", "tgt", "v2", "node2", "j"))

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
  }

  between <- list(
    weights = between_weights,
    inits = between_inits,
    labels = cluster_names,
    data = between_seq_data
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
        # Single node: count self-loops
        keep_self <- w_from %in% cl_nodes & w_to %in% cl_nodes
        self_weight <- if (any(keep_self)) {
          aggregate_weights(w_w[keep_self], method)
        } else { 0 }
        within_weights_i <- matrix(self_weight, 1, 1,
                                    dimnames = list(cl_nodes, cl_nodes))
        within_inits_i <- setNames(1, cl_nodes)
      } else {
        # Filter transitions for this cluster (keep self-loops)
        keep <- w_from %in% cl_nodes & w_to %in% cl_nodes
        cf <- w_from[keep]
        ct <- w_to[keep]
        cw <- w_w[keep]

        within_raw <- matrix(0, n_i, n_i,
                              dimnames = list(cl_nodes, cl_nodes))

        if (length(cf) > 0) {
          pair_keys <- paste(cf, ct, sep = "\t")
          agg_vals <- tapply(cw, pair_keys, function(w) {
            aggregate_weights(w, method)
          })
          for (key in names(agg_vals)) {
            parts <- strsplit(key, "\t")[[1]]
            within_raw[parts[1], parts[2]] <- agg_vals[[key]]
          }
        }

        within_weights_i <- .process_weights(within_raw, type, directed)

        col_sums_w <- colSums(within_raw, na.rm = TRUE)
        total_w <- sum(col_sums_w, na.rm = TRUE)
        within_inits_i <- if (!is.na(total_w) && total_w > 0) {
          col_sums_w / total_w
        } else {
          rep(1 / n_i, n_i)
        }
        names(within_inits_i) <- cl_nodes
      }

      # Attach filtered sequence data for this cluster
      cl_seq_data <- if (!is.null(within_seq_data_list)) {
        within_seq_data_list[[cl_name]]
      } else {
        NULL
      }

      list(
        weights = within_weights_i,
        inits = within_inits_i,
        labels = cl_nodes,
        data = cl_seq_data
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
  all_nodes <- sort(unique(c(from_nodes, to_nodes)))
  n_nodes <- length(all_nodes)

  structure(
    list(
      macro = between,
      clusters = within_data,
      cluster_members = cluster_list,
      edges = edges,
      data = data,
      meta = list(
        type = type,
        method = method,
        directed = directed,
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
#' @keywords internal
.build_mcml_edgelist <- function(df, clusters, method, type,
                                  directed, compute_within) {

  col_names <- tolower(names(df))

  # Detect from/to columns
  from_col <- which(col_names %in% c("from", "source", "src",
                                       "v1", "node1", "i"))[1]
  if (is.na(from_col)) from_col <- 1L

  to_col <- which(col_names %in% c("to", "target", "tgt",
                                     "v2", "node2", "j"))[1]
  if (is.na(to_col)) to_col <- 2L

  # Detect weight column
  weight_col <- which(col_names %in% c("weight", "w", "value", "strength"))[1]
  has_weight <- !is.na(weight_col)

  from_vals <- as.character(df[[from_col]])
  to_vals <- as.character(df[[to_col]])
  weights <- if (has_weight) as.numeric(df[[weight_col]]) else rep(1, nrow(df))

  # Remove rows with NA in from/to
  valid <- !is.na(from_vals) & !is.na(to_vals)
  from_vals <- from_vals[valid]
  to_vals <- to_vals[valid]
  weights <- weights[valid]

  all_nodes <- sort(unique(c(from_vals, to_vals)))

  # Handle clusters parameter
  if (is.character(clusters) && length(clusters) == 1 &&
      clusters %in% names(df)) {
    # Column name: build lookup from both from+group and to+group
    group_col <- df[[clusters]]
    group_col <- as.character(group_col[valid])

    # Build mapping from from-side
    from_map <- setNames(group_col, from_vals)
    # Build mapping from to-side
    to_map <- setNames(group_col, to_vals)
    # Merge (from takes priority if conflicting, but shouldn't)
    full_map <- c(to_map, from_map)
    # Keep unique node -> cluster mapping
    full_map <- full_map[!duplicated(names(full_map))]

    # Build cluster_list
    cluster_list <- split(names(full_map), unname(full_map))
    cluster_list <- lapply(cluster_list, sort)
    cluster_lookup <- full_map

    # Re-derive all_nodes from the lookup
    all_nodes <- sort(names(cluster_lookup))
  } else {
    # List or vector clusters
    if (is.null(clusters)) {
      stop("clusters argument is required for data.frame input", call. = FALSE)
    }

    if (is.list(clusters) && !is.data.frame(clusters)) {
      cluster_list <- clusters
    } else {
      # Membership vector
      cluster_list <- .normalize_clusters(clusters, all_nodes)
    }

    cluster_lookup <- .build_cluster_lookup(cluster_list, all_nodes)
  }

  .build_from_transitions(from_vals, to_vals, weights,
                            cluster_lookup, cluster_list,
                            method, type, directed, compute_within,
                            data = df)
}

#' Build MCML from sequence data.frame
#' @keywords internal
.build_mcml_sequence <- function(df, clusters, method, type,
                                  directed, compute_within) {

  if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)

  stopifnot(is.data.frame(df))

  nc <- ncol(df)
  if (nc < 2) {
    stop("Sequence data must have at least 2 columns (time steps)",
         call. = FALSE)
  }

  # Extract consecutive pairs: (col[t], col[t+1]) for all rows
  pairs <- lapply(seq_len(nc - 1), function(t) {
    from_t <- as.character(df[[t]])
    to_t <- as.character(df[[t + 1]])
    data.frame(from = from_t, to = to_t, stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, pairs)

  # Remove NA pairs
  valid <- !is.na(pairs$from) & !is.na(pairs$to)
  from_vals <- pairs$from[valid]
  to_vals <- pairs$to[valid]
  weights <- rep(1, length(from_vals))

  all_nodes <- sort(unique(c(from_vals, to_vals)))

  if (is.null(clusters)) {
    stop("clusters argument is required for sequence data", call. = FALSE)
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
#' @keywords internal
.process_weights <- function(raw_weights, type, directed = TRUE) {
  if (type == "raw" || type == "frequency") {
    return(raw_weights)
  }

  if (type == "cooccurrence") {
    # Symmetrize
    return((raw_weights + t(raw_weights)) / 2)
  }

  if (type == "tna" || type == "semi_markov") {
    # Row-normalize so rows sum to 1
    rs <- rowSums(raw_weights, na.rm = TRUE)
    processed <- raw_weights / ifelse(rs == 0 | is.na(rs), 1, rs)
    processed[is.na(processed)] <- 0
    return(processed)
  }

  # Default: return as-is
  raw_weights # nocov
}

#' Convert cluster_summary to tna Objects
#'
#' Converts a \code{cluster_summary} object to proper tna objects that can be
#' used with all functions from the tna package. Creates both a between-cluster
#' tna model (cluster-level transitions) and within-cluster tna models (internal
#' transitions within each cluster).
#'
#' This is the final step in the MCML workflow, enabling full integration with
#' the tna package for centrality analysis, bootstrap validation, permutation
#' tests, and visualization.
#'
#' @param x A \code{cluster_summary} object created by \code{\link{cluster_summary}}.
#'   The cluster_summary should typically be created with \code{type = "tna"} to
#'   ensure row-normalized transition probabilities. If created with
#'   \code{type = "raw"}, the raw counts will be passed to \code{tna::tna()}
#'   which will normalize them.
#'
#' @return A \code{cluster_tna} object (S3 class) containing:
#'   \describe{
#'     \item{between}{A tna object representing cluster-level transitions.
#'       Contains \code{$weights} (k x k transition matrix), \code{$inits}
#'       (initial distribution), and \code{$labels} (cluster names).
#'       Use this for analyzing how learners/entities move between high-level
#'       groups or phases.}
#'     \item{within}{Named list of tna objects, one per cluster. Each tna object
#'       represents internal transitions within that cluster. Contains
#'       \code{$weights} (n_i x n_i matrix), \code{$inits} (initial distribution),
#'       and \code{$labels} (node labels). Clusters with single nodes or zero-row
#'       nodes are excluded (tna requires positive row sums).}
#'   }
#'
#' @details
#' ## Requirements
#'
#' The tna package must be installed. If not available, the function throws
#' an error with installation instructions.
#'
#' ## Workflow
#'
#' \preformatted{
#' # Full MCML workflow
#' net <- build_network(data, method = "relative")
#' net$nodes$clusters <- group_assignments
#' cs <- cluster_summary(net, type = "tna")
#' tna_models <- as_tna(cs)
#'
#' # Now use tna package functions
#' plot(tna_models$macro)
#' tna::centralities(tna_models$macro)
#' tna::bootstrap(tna_models$macro, iter = 1000)
#'
#' # Analyze within-cluster patterns
#' plot(tna_models$clusters$ClusterA)
#' tna::centralities(tna_models$clusters$ClusterA)
#' }
#'
#' ## Excluded Clusters
#'
#' A within-cluster tna cannot be created when:
#' \itemize{
#'   \item The cluster has only 1 node (no internal transitions possible)
#'   \item Some nodes in the cluster have no outgoing edges (row sums to 0)
#' }
#'
#' These clusters are silently excluded from \code{$clusters}. The between-cluster
#' model still includes all clusters.
#'
#' @export
#' @seealso
#'   \code{\link{cluster_summary}} to create the input object,
#'   \code{plot()} for visualization without conversion,
#'   \code{tna::tna} for the underlying tna constructor
#'
#' @examples
#' # -----------------------------------------------------
#' # Basic usage
#' # -----------------------------------------------------
#' mat <- matrix(runif(36), 6, 6)
#' rownames(mat) <- colnames(mat) <- LETTERS[1:6]
#'
#' clusters <- list(
#'   G1 = c("A", "B"),
#'   G2 = c("C", "D"),
#'   G3 = c("E", "F")
#' )
#'
#' cs <- cluster_summary(mat, clusters, type = "tna")
#' tna_models <- as_tna(cs)
#'
#' # Print summary
#' tna_models
#'
#' # -----------------------------------------------------
#' # Access components
#' # -----------------------------------------------------
#' # Between-cluster tna
#' tna_models$macro
#' tna_models$macro$weights  # 3x3 transition matrix
#' tna_models$macro$inits    # Initial distribution
#' tna_models$macro$labels   # c("G1", "G2", "G3")
#'
#' # Within-cluster tnas
#' names(tna_models$clusters)    # Which clusters have within models
#' tna_models$clusters$G1        # tna for cluster G1
#' tna_models$clusters$G1$weights  # 2x2 matrix (A, B)
#'
#' # -----------------------------------------------------
#' # Use with tna package (requires tna)
#' # -----------------------------------------------------
#' \dontrun{
#' # Plot
#' plot(tna_models$macro)
#' plot(tna_models$clusters$G1)
#'
#' # Centrality analysis
#' tna::centralities(tna_models$macro)
#' tna::centralities(tna_models$clusters$G1)
#' tna::centralities(tna_models$clusters$G2)
#' }
#'
#' \dontrun{
#' # Bootstrap validation (requires tna built from sequence data)
#' boot <- tna::bootstrap(tna_models$macro, iter = 1000)
#' summary(boot)
#' }
#'
#' # -----------------------------------------------------
#' # Check which within-cluster models were created
#' # -----------------------------------------------------
#' cs <- cluster_summary(mat, clusters, type = "tna")
#' tna_models <- as_tna(cs)
#'
#' # All cluster names
#' names(cs$clusters)
#'
#' # Clusters with valid within-models
#' names(tna_models$clusters)
#'
#' # Clusters excluded (single node or zero rows)
#' setdiff(names(cs$clusters), names(tna_models$clusters))
as_tna <- function(x) {
  UseMethod("as_tna")
}

#' @rdname as_tna
#' @return A \code{netobject_group} with data preserved from each sub-network.
#' @export
as_tna.mcml <- function(x) {
  # Determine method from meta
  meta_type <- x$meta$type
  net_method <- if (!is.null(meta_type) && meta_type %in% c("raw", "frequency")) {
    "frequency"
  } else {
    "relative"
  }
  directed <- if (!is.null(x$meta$directed)) x$meta$directed else TRUE

  # Macro
  macro_net <- .wrap_netobject(x$macro$weights, data = x$macro$data,
                               method = net_method, directed = directed)

  # Per-cluster
  cluster_nets <- if (!is.null(x$clusters)) {
    nets <- lapply(names(x$clusters), function(cl) {
      w <- x$clusters[[cl]]$weights
      d <- x$clusters[[cl]]$data
      if (any(rowSums(w) == 0)) return(NULL)
      .wrap_netobject(w, data = d, method = net_method, directed = directed)
    })
    names(nets) <- names(x$clusters)
    nets[!vapply(nets, is.null, logical(1))]
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
                            directed = TRUE) {
  states <- rownames(mat)
  edges <- .extract_edges_from_matrix(mat, directed = directed)
  nodes_df <- data.frame(
    id = seq_along(states), label = states, name = states,
    x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
  )

  structure(
    list(
      data       = data,
      weights    = mat,
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
#' @return A \code{tna} object constructed from the input.
#' @export
as_tna.default <- function(x) {
 if (inherits(x, "tna")) {
    return(x)
  }
  stop("Cannot convert object of class '", class(x)[1], "' to tna", call. = FALSE)
}

#' Normalize cluster specification to list format
#' @keywords internal
#' @noRd
.normalize_clusters <- function(clusters, node_names) {
  if (is.data.frame(clusters)) {
    # Data frame with node and group columns
    stopifnot(ncol(clusters) >= 2)
    nodes <- as.character(clusters[[1]])
    groups <- as.character(clusters[[2]])
    clusters <- split(nodes, groups)
  }

  if (is.list(clusters)) {
    # Already a list - validate node names
    all_nodes <- unlist(clusters)
    if (!all(all_nodes %in% node_names)) {
      missing <- setdiff(all_nodes, node_names)
      stop("Unknown nodes in clusters: ",
           paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
    }
    return(clusters)
  }

  if (is.vector(clusters) && (is.numeric(clusters) || is.integer(clusters))) {
    # Membership vector
    if (length(clusters) != length(node_names)) {
      stop("Membership vector length (", length(clusters),
           ") must equal number of nodes (", length(node_names), ")",
           call. = FALSE)
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
    if (length(clusters) != length(node_names)) {
      stop("Membership vector length must equal number of nodes", call. = FALSE)
    }
    clusters <- as.character(clusters)
    unique_clusters <- unique(clusters)
    cluster_list <- lapply(unique_clusters, function(k) {
      node_names[clusters == k]
    })
    names(cluster_list) <- unique_clusters
    return(cluster_list)
  }

  stop("clusters must be a list, numeric vector, or factor", call. = FALSE)
}

# ==============================================================================
# S3 Methods
# ==============================================================================

#' Print Method for mcml
#'
#' @param x An \code{mcml} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.mcml <- function(x, ...) {
  n_clusters <- x$meta$n_clusters
  n_nodes <- x$meta$n_nodes
  cluster_sizes <- x$meta$cluster_sizes

  cat("MCML Network\n")
  cat("============\n")
  cat("Type:", x$meta$type, " | Method:", x$meta$method, "\n")
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

#' Summary Method for mcml
#'
#' @param object An \code{mcml} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
summary.mcml <- function(object, ...) {
  print(object, ...)
}
