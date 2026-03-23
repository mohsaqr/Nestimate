#' Temporal Network Analysis
#'
#' @description
#' Constructs and analyzes temporal networks from edge-timestamped data.
#' Implements time-respecting path algorithms, temporal centrality measures,
#' edge dynamics, burstiness analysis, and time-varying degree computation.
#' All algorithms are implemented from scratch using only base R and igraph.
#'
#' @param edges Data frame with edge events. Must contain columns for source,
#'   target, onset time, and terminus time. Each row represents one edge event
#'   active during \code{[onset, terminus)}.
#' @param from_col Character. Column name for source vertex. Default: \code{"from"}.
#' @param to_col Character. Column name for target vertex. Default: \code{"to"}.
#' @param onset_col Character. Column name for onset time. Default: \code{"onset"}.
#' @param terminus_col Character. Column name for terminus time.
#'   Default: \code{"terminus"}.
#' @param directed Logical. Whether edges are directed. Default: \code{TRUE}.
#' @param time_interval Numeric. Bin width for time series aggregation.
#'   Default: \code{1}.
#'
#' @return An object of class \code{"temporal_network"} containing:
#' \describe{
#'   \item{edges}{Canonicalized edge list (from, to, onset, terminus).}
#'   \item{n_vertices}{Number of unique vertices.}
#'   \item{n_edges}{Number of edge events.}
#'   \item{directed}{Whether the network is directed.}
#'   \item{time_range}{Numeric vector of length 2: (min onset, max terminus).}
#'   \item{time_bins}{Numeric vector of time bin start points.}
#'   \item{time_interval}{Bin width used.}
#'   \item{vertex_names}{Character vector of vertex names.}
#'   \item{degree}{Matrix (vertices x bins): total degree over time.}
#'   \item{indegree}{Matrix (vertices x bins): in-degree (directed only).}
#'   \item{outdegree}{Matrix (vertices x bins): out-degree (directed only).}
#'   \item{reachability_fwd}{Integer vector: forward reachable set size per vertex.}
#'   \item{reachability_bkwd}{Integer vector: backward reachable set size per vertex.}
#'   \item{temporal_closeness}{Numeric vector: temporal closeness centrality.}
#'   \item{temporal_betweenness}{Numeric vector: temporal betweenness centrality.}
#'   \item{formation}{Integer vector: edges forming per time bin.}
#'   \item{dissolution}{Integer vector: edges dissolving per time bin.}
#'   \item{edge_durations}{Numeric vector: total duration per unique edge pair.}
#'   \item{iet_vertex}{List of inter-event time vectors per vertex.}
#'   \item{burstiness}{Numeric vector: per-vertex burstiness coefficient.}
#'   \item{temporal_density}{Numeric scalar: fraction of possible temporal edges.}
#'   \item{snapshots}{List of igraph objects per time bin.}
#'   \item{density_bins}{Numeric vector (n_bins): static edge density per bin.}
#'   \item{reciprocity}{Numeric vector (n_bins): reciprocity per bin.
#'     Directed only; NULL for undirected. NaN on empty graphs mapped to 0.}
#'   \item{mutuality}{Integer vector (n_bins): mutual dyad count per bin.
#'     Directed only; NULL for undirected.}
#'   \item{dyad_census}{Integer matrix (3 x n_bins): Mutual/Asymmetric/Null
#'     dyads per bin. Directed only; NULL for undirected.}
#'   \item{transitivity}{Numeric vector (n_bins): transitivity per bin.
#'     Directed: weak transitivity (fraction of directed two-paths that are
#'     closed). Undirected: global clustering coefficient. Vacuous case
#'     (no two-paths/triples) returns 1, matching \code{sna::gtrans()}.}
#'   \item{centralization_degree}{Numeric vector (n_bins): degree centralization.}
#'   \item{centralization_betweenness}{Numeric vector (n_bins): betweenness
#'     centralization.}
#'   \item{centralization_closeness}{Numeric vector (n_bins): closeness
#'     centralization.}
#'   \item{n_components}{Integer vector (n_bins): number of weakly connected
#'     components per bin.}
#'   \item{triad_census}{Integer matrix (16 x n_bins): triad census per bin
#'     using MAN labeling. Directed only; NULL for undirected.}
#'   \item{mean_distance}{Numeric vector (n_bins): mean shortest path length.}
#'   \item{diameter}{Numeric vector (n_bins): longest shortest path.}
#'   \item{assortativity}{Numeric vector (n_bins): degree assortativity.
#'     NA for bins with no edges.}
#'   \item{closeness_snapshot}{Numeric matrix (n_vertices x n_bins): static
#'     closeness centrality per snapshot. NaN mapped to 0.}
#'   \item{betweenness_snapshot}{Numeric matrix (n_vertices x n_bins): static
#'     betweenness centrality per snapshot.}
#'   \item{eigenvector}{Numeric matrix (n_vertices x n_bins): eigenvector
#'     centrality per snapshot.}
#'   \item{page_rank}{Numeric matrix (n_vertices x n_bins): PageRank per
#'     snapshot.}
#'   \item{hub_score}{Numeric matrix (n_vertices x n_bins): HITS hub score
#'     per snapshot.}
#'   \item{authority_score}{Numeric matrix (n_vertices x n_bins): HITS authority
#'     score per snapshot.}
#'   \item{constraint}{Numeric matrix (n_vertices x n_bins): Burt's structural
#'     constraint. NA for isolated vertices.}
#'   \item{coreness}{Integer matrix (n_vertices x n_bins): k-core decomposition.}
#' }
#'
#' @examples
#' \dontrun{
#' edges <- data.frame(
#'   from = c("A", "B", "A"),
#'   to   = c("B", "C", "C"),
#'   onset = c(1, 3, 5),
#'   terminus = c(4, 7, 6)
#' )
#' tn <- temporal_network(edges)
#' print(tn)
#' summary(tn)
#' plot(tn, type = "degree")
#' }
#'
#' @export
temporal_network <- function(edges,
                             from_col = "from",
                             to_col = "to",
                             onset_col = "onset",
                             terminus_col = "terminus",
                             directed = TRUE,
                             time_interval = 1) {

  # --- Validate inputs ---
  stopifnot(
    is.data.frame(edges),
    is.character(from_col), length(from_col) == 1L,
    is.character(to_col), length(to_col) == 1L,
    is.character(onset_col), length(onset_col) == 1L,
    is.character(terminus_col), length(terminus_col) == 1L,
    is.logical(directed), length(directed) == 1L,
    is.numeric(time_interval), length(time_interval) == 1L,
    time_interval > 0
  )

  required_cols <- c(from_col, to_col, onset_col, terminus_col)
  missing_cols <- setdiff(required_cols, names(edges))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  from_raw <- edges[[from_col]]
  to_raw <- edges[[to_col]]
  onset_raw <- edges[[onset_col]]
  terminus_raw <- edges[[terminus_col]]

  if (!is.numeric(onset_raw)) stop("onset column must be numeric")
  if (!is.numeric(terminus_raw)) stop("terminus column must be numeric")

  if (anyNA(from_raw) || anyNA(to_raw)) stop("from/to columns must not contain NA")
  if (anyNA(onset_raw) || anyNA(terminus_raw)) {
    stop("onset/terminus columns must not contain NA")
  }
  if (any(terminus_raw <= onset_raw)) {
    stop("terminus must be strictly greater than onset for every edge")
  }

  # --- Canonicalize ---
  # Sort numerically if all names are numeric, otherwise alphabetically
  all_names <- unique(c(as.character(from_raw), as.character(to_raw)))
  numeric_names <- suppressWarnings(as.numeric(all_names))
  if (!anyNA(numeric_names)) {
    vertex_names <- all_names[order(numeric_names)]
  } else {
    vertex_names <- sort(all_names)
  }
  n_vertices <- length(vertex_names)
  vertex_lookup <- setNames(seq_len(n_vertices), vertex_names)

  from_id <- vertex_lookup[as.character(from_raw)]
  to_id <- vertex_lookup[as.character(to_raw)]

  canon_edges <- data.frame(
    from = from_id,
    to = to_id,
    onset = onset_raw,
    terminus = terminus_raw,
    stringsAsFactors = FALSE
  )
  # Sort by onset
  ord <- order(canon_edges$onset)
  canon_edges <- canon_edges[ord, , drop = FALSE]
  rownames(canon_edges) <- NULL

  # --- Time bins ---
  t_min <- min(canon_edges$onset)
  t_max <- max(canon_edges$terminus)
  # Extend one bin past t_max so dissolution at t_max is captured (match tsna)
  time_bins <- seq(t_min, t_max + time_interval, by = time_interval)
  # Trim if overshoot is excessive (keep at most 1 bin past t_max)
  while (length(time_bins) > 2L && time_bins[length(time_bins) - 1L] > t_max) {
    time_bins <- time_bins[-length(time_bins)]
  }
  n_bins <- length(time_bins) - 1L

  # --- Time-varying degree ---
  degree_mat <- .compute_time_varying_degree(
    canon_edges, n_vertices, time_bins, n_bins, directed
  )

  # --- Temporal paths, reachability, closeness, betweenness ---
  path_results <- .compute_all_temporal_paths(
    canon_edges, n_vertices, directed
  )

  # --- Edge formation / dissolution ---
  form_diss <- .compute_formation_dissolution(
    canon_edges, time_bins, n_bins
  )

  # --- Edge durations ---
  edge_durations <- .compute_edge_durations(canon_edges, n_vertices)

  # --- Inter-event times and burstiness ---
  iet_burst <- .compute_iet_burstiness(
    canon_edges, n_vertices, directed
  )

  # --- Temporal density ---
  temporal_density <- .compute_temporal_density(
    canon_edges, n_vertices, directed, t_min, t_max
  )

  # --- Static snapshots ---
  snapshots <- .build_static_snapshots(
    canon_edges, n_vertices, vertex_names, directed, time_bins, n_bins
  )

  # --- Snapshot-based metrics ---
  snap_metrics <- .compute_snapshot_metrics(
    snapshots, n_vertices, vertex_names, directed, n_bins
  )

  # --- Assemble result ---
  result <- list(
    edges              = canon_edges,
    n_vertices         = n_vertices,
    n_edges            = nrow(canon_edges),
    directed           = directed,
    time_range         = c(t_min, t_max),
    time_bins          = time_bins,
    time_interval      = time_interval,
    vertex_names       = vertex_names,
    degree             = degree_mat$total,
    indegree           = degree_mat$indeg,
    outdegree          = degree_mat$outdeg,
    reachability_fwd   = path_results$reachability_fwd,
    reachability_bkwd  = path_results$reachability_bkwd,
    temporal_closeness = path_results$temporal_closeness,
    temporal_betweenness = path_results$temporal_betweenness,
    formation          = form_diss$formation,
    dissolution        = form_diss$dissolution,
    edge_durations     = edge_durations,
    iet_vertex         = iet_burst$iet_vertex,
    burstiness         = iet_burst$burstiness,
    temporal_density   = temporal_density,
    snapshots          = snapshots
  )
  result <- c(result, snap_metrics)

  class(result) <- "temporal_network"
  result
}


# ===========================================================================
# Internal: Time-varying degree
# ===========================================================================

#' @keywords internal
.compute_time_varying_degree <- function(edges, n_vertices, time_bins, n_bins,
                                         directed) {
  total_mat <- matrix(0L, nrow = n_vertices, ncol = n_bins)
  indeg_mat <- if (directed) matrix(0L, nrow = n_vertices, ncol = n_bins) else NULL
  outdeg_mat <- if (directed) matrix(0L, nrow = n_vertices, ncol = n_bins) else NULL

  # For each bin, compute point-in-time degree at bin_start (match tsna tDegree)
  vapply(seq_len(n_bins), function(k) {
    bin_start <- time_bins[k]
    # Edge active at time t: onset <= t AND terminus > t
    active <- which(edges$onset <= bin_start & edges$terminus > bin_start)

    if (length(active) > 0L) {
      from_active <- edges$from[active]
      to_active <- edges$to[active]

      if (directed) {
        outdeg_mat[, k] <<- tabulate(from_active, nbins = n_vertices)
        indeg_mat[, k] <<- tabulate(to_active, nbins = n_vertices)
        total_mat[, k] <<- outdeg_mat[, k] + indeg_mat[, k]
      } else {
        # Undirected: each edge contributes 1 to each endpoint
        both <- c(from_active, to_active)
        total_mat[, k] <<- tabulate(both, nbins = n_vertices)
      }
    }
    0L
  }, integer(1L))

  list(total = total_mat, indeg = indeg_mat, outdeg = outdeg_mat)
}


# ===========================================================================
# Internal: Temporal BFS (earliest-arrival)
# ===========================================================================

#' @keywords internal
.temporal_bfs <- function(edges, n_vertices, source, start_time = -Inf) {
  # edges must be sorted by onset
  earliest <- rep(Inf, n_vertices)
  previous <- rep(NA_integer_, n_vertices)
  earliest[source] <- start_time

  n_edges <- nrow(edges)
  e_from <- edges$from
  e_to <- edges$to
  e_onset <- edges$onset
  e_terminus <- edges$terminus

  # Repeat until convergence (for non-DAG temporal structures)
  changed <- TRUE
  max_iter <- n_vertices  # worst case: one vertex per iteration
  iter <- 0L


  while (changed && iter < max_iter) {
    changed <- FALSE
    iter <- iter + 1L

    for (i in seq_len(n_edges)) {
      arr_u <- earliest[e_from[i]]
      if (is.infinite(arr_u) && arr_u > 0) next  # unreachable

      # Can traverse if max(arr_u, onset) < terminus
      depart <- max(arr_u, e_onset[i])
      if (depart < e_terminus[i]) {
        if (depart < earliest[e_to[i]]) {
          earliest[e_to[i]] <- depart
          previous[e_to[i]] <- e_from[i]
          changed <- TRUE
        }
      }
    }
  }

  list(earliest = earliest, previous = previous)
}

#' @keywords internal
.temporal_bfs_backward <- function(edges, n_vertices, target) {
  # Backward: find latest-departure times to reach target

  # Reverse edges and negate times
  rev_edges <- data.frame(
    from = edges$to,
    to = edges$from,
    onset = -edges$terminus,
    terminus = -edges$onset,
    stringsAsFactors = FALSE
  )
  rev_edges <- rev_edges[order(rev_edges$onset), , drop = FALSE]

  result <- .temporal_bfs(rev_edges, n_vertices, target, start_time = -Inf)
  # Vertices reachable = those with finite earliest (in negated time)
  list(
    reachable = which(is.finite(result$earliest)),
    latest_departure = -result$earliest
  )
}


# ===========================================================================
# Internal: All temporal paths + centralities
# ===========================================================================

#' @keywords internal
.compute_all_temporal_paths <- function(edges, n_vertices, directed) {
  # For undirected networks, add reverse edges so BFS traverses both directions
  if (!directed) {
    rev_edges <- data.frame(
      from = edges$to, to = edges$from,
      onset = edges$onset, terminus = edges$terminus,
      stringsAsFactors = FALSE
    )
    edges <- rbind(edges, rev_edges)
    edges <- edges[order(edges$onset), , drop = FALSE]
    rownames(edges) <- NULL
  }

  reachability_fwd <- integer(n_vertices)
  reachability_bkwd <- integer(n_vertices)
  temporal_closeness <- numeric(n_vertices)
  temporal_betweenness <- numeric(n_vertices)
  t_min <- min(edges$onset)

  # Forward reachability + closeness + betweenness
  for (s in seq_len(n_vertices)) {
    bfs <- .temporal_bfs(edges, n_vertices, s, start_time = t_min)
    reachable <- which(is.finite(bfs$earliest) & seq_len(n_vertices) != s)
    reachability_fwd[s] <- length(reachable) + 1L  # +1 includes self (match tsna tReach)

    # Temporal closeness
    if (length(reachable) > 0L) {
      distances <- bfs$earliest[reachable] - t_min
      total_dist <- sum(distances)
      if (total_dist > 0) {
        temporal_closeness[s] <- length(reachable) / total_dist
      }
    }

    # Temporal betweenness: trace back paths
    for (t in reachable) {
      path <- .trace_path(bfs$previous, s, t)
      if (length(path) > 2L) {
        # Intermediate vertices (exclude source and target)
        intermediates <- path[2:(length(path) - 1L)]
        temporal_betweenness[intermediates] <-
          temporal_betweenness[intermediates] + 1
      }
    }
  }

  # Backward reachability
  for (t in seq_len(n_vertices)) {
    bkwd <- .temporal_bfs_backward(edges, n_vertices, t)
    # Target has -Inf (not finite), so is.finite already excludes it
    reachability_bkwd[t] <- length(bkwd$reachable) + 1L  # +1 includes self (match tsna tReach)
  }

  # Normalize betweenness
  if (directed && n_vertices > 2L) {
    normalization <- (n_vertices - 1L) * (n_vertices - 2L)
    temporal_betweenness <- temporal_betweenness / normalization
  } else if (!directed && n_vertices > 2L) {
    normalization <- (n_vertices - 1L) * (n_vertices - 2L) / 2
    temporal_betweenness <- temporal_betweenness / normalization
  }

  list(
    reachability_fwd = reachability_fwd,
    reachability_bkwd = reachability_bkwd,
    temporal_closeness = temporal_closeness,
    temporal_betweenness = temporal_betweenness
  )
}

#' @keywords internal
.trace_path <- function(previous, source, target) {
  path <- target
  current <- target
  max_steps <- length(previous)
  step <- 0L

  while (!is.na(previous[current]) && current != source && step < max_steps) {
    current <- previous[current]
    path <- c(current, path)
    step <- step + 1L
  }

  if (current == source) path else integer(0L)
}


# ===========================================================================
# Internal: Edge formation / dissolution
# ===========================================================================

#' @keywords internal
.compute_formation_dissolution <- function(edges, time_bins, n_bins) {
  formation <- integer(n_bins)
  dissolution <- integer(n_bins)

  vapply(seq_len(n_bins), function(k) {
    bin_start <- time_bins[k]
    bin_end <- time_bins[k + 1L]
    # Formation: edges whose onset falls in [bin_start, bin_end)
    formation[k] <<- sum(edges$onset >= bin_start & edges$onset < bin_end)
    # Dissolution: edges whose terminus falls in [bin_start, bin_end) (match tsna)
    dissolution[k] <<- sum(edges$terminus >= bin_start & edges$terminus < bin_end)
    0L
  }, integer(1L))

  list(formation = formation, dissolution = dissolution)
}


# ===========================================================================
# Internal: Edge durations
# ===========================================================================

#' @keywords internal
.compute_edge_durations <- function(edges, n_vertices) {
  # Per unique (from, to) pair: sum of durations across all spells
  durations <- edges$terminus - edges$onset
  pair_key <- paste(edges$from, edges$to, sep = "_")
  tapply(durations, pair_key, sum)
}


# ===========================================================================
# Internal: Inter-event times and burstiness
# ===========================================================================

#' @keywords internal
.compute_iet_burstiness <- function(edges, n_vertices, directed) {
  iet_vertex <- vector("list", n_vertices)
  burstiness <- rep(NA_real_, n_vertices)

  for (v in seq_len(n_vertices)) {
    # Event times for vertex v: onset of edges incident to v
    if (directed) {
      event_times <- sort(c(
        edges$onset[edges$from == v],
        edges$onset[edges$to == v]
      ))
    } else {
      event_times <- sort(c(
        edges$onset[edges$from == v],
        edges$onset[edges$to == v]
      ))
    }

    event_times <- unique(event_times)

    if (length(event_times) >= 2L) {
      iet <- diff(event_times)
      iet_vertex[[v]] <- iet
      mu <- mean(iet)
      sigma <- if (length(iet) >= 2L) stats::sd(iet) else 0
      denom <- sigma + mu
      if (!is.na(denom) && denom > 0) {
        burstiness[v] <- (sigma - mu) / denom
      } else {
        burstiness[v] <- 0
      }
    } else {
      iet_vertex[[v]] <- numeric(0L)
      burstiness[v] <- NA_real_
    }
  }

  list(iet_vertex = iet_vertex, burstiness = burstiness)
}


# ===========================================================================
# Internal: Temporal density
# ===========================================================================

#' @keywords internal
.compute_temporal_density <- function(edges, n_vertices, directed, t_min, t_max) {
  # Total possible directed edge slots: n*(n-1) * (t_max - t_min)
  # (or n*(n-1)/2 for undirected)
  time_span <- t_max - t_min
  if (time_span <= 0) return(0)

  if (directed) {
    max_possible <- n_vertices * (n_vertices - 1L) * time_span
  } else {
    max_possible <- n_vertices * (n_vertices - 1L) / 2 * time_span
  }

  if (max_possible <= 0) return(0)

  # Actual temporal volume: sum of all edge durations
  total_duration <- sum(edges$terminus - edges$onset)
  total_duration / max_possible
}


# ===========================================================================
# Internal: Static snapshots
# ===========================================================================

#' @keywords internal
.build_static_snapshots <- function(edges, n_vertices, vertex_names,
                                    directed, time_bins, n_bins) {
  lapply(seq_len(n_bins), function(k) {
    bin_start <- time_bins[k]
    # Point-in-time at bin_start (matches tsna network.collapse / tSnaStats)
    active <- which(edges$onset <= bin_start & edges$terminus > bin_start)

    if (length(active) > 0L) {
      el <- cbind(edges$from[active], edges$to[active])
      g <- igraph::graph_from_edgelist(el, directed = directed)
      # Ensure all vertices are present
      missing_v <- setdiff(seq_len(n_vertices), igraph::V(g))
      if (length(missing_v) > 0L) {
        g <- igraph::add_vertices(g, length(missing_v))
      }
      # Simplify: collapse multi-edges (matches network.collapse)
      g <- igraph::simplify(g)
      igraph::V(g)$name <- vertex_names[seq_len(igraph::vcount(g))]
      g
    } else {
      g <- igraph::make_empty_graph(n = n_vertices, directed = directed)
      igraph::V(g)$name <- vertex_names
      g
    }
  })
}


# ===========================================================================
# Internal: Snapshot-based graph & node metrics
# ===========================================================================

#' @keywords internal
.compute_snapshot_metrics <- function(snapshots, n_vertices, vertex_names,
                                      directed, n_bins) {
  # --- Graph-level vectors ---
  density_bins <- numeric(n_bins)
  transitivity_bins <- numeric(n_bins)
  centr_deg <- numeric(n_bins)
  centr_betw_vec <- numeric(n_bins)
  centr_clo_vec <- numeric(n_bins)
  n_comp <- integer(n_bins)
  mean_dist <- numeric(n_bins)
  diam_vec <- numeric(n_bins)
  assort_vec <- rep(NA_real_, n_bins)

  # Directed-only graph-level
  recip_vec <- if (directed) numeric(n_bins) else NULL
  mutual_vec <- if (directed) integer(n_bins) else NULL
  dyad_mat <- if (directed) {
    m <- matrix(0L, nrow = 3L, ncol = n_bins)
    rownames(m) <- c("Mutual", "Asymmetric", "Null")
    m
  } else NULL
  triad_mat <- if (directed) {
    m <- matrix(0L, nrow = 16L, ncol = n_bins)
    rownames(m) <- c(
      "003", "012", "102", "021D", "021U", "021C", "111D", "111U",
      "030T", "030C", "201", "120D", "120U", "120C", "210", "300"
    )
    m
  } else NULL

  # --- Node-level matrices ---
  close_snap <- matrix(0, nrow = n_vertices, ncol = n_bins)
  betw_snap <- matrix(0, nrow = n_vertices, ncol = n_bins)
  eig_mat <- matrix(0, nrow = n_vertices, ncol = n_bins)
  pr_mat <- matrix(0, nrow = n_vertices, ncol = n_bins)
  hub_mat <- matrix(0, nrow = n_vertices, ncol = n_bins)
  auth_mat <- matrix(0, nrow = n_vertices, ncol = n_bins)
  constr_mat <- matrix(NA_real_, nrow = n_vertices, ncol = n_bins)
  core_mat <- matrix(0L, nrow = n_vertices, ncol = n_bins)

  rownames(close_snap) <- vertex_names
  rownames(betw_snap) <- vertex_names
  rownames(eig_mat) <- vertex_names
  rownames(pr_mat) <- vertex_names
  rownames(hub_mat) <- vertex_names
  rownames(auth_mat) <- vertex_names
  rownames(constr_mat) <- vertex_names
  rownames(core_mat) <- vertex_names

  # --- Single pass over all snapshots ---
  invisible(vapply(seq_len(n_bins), function(k) {
    g <- snapshots[[k]]
    ne <- igraph::ecount(g)

    # Graph-level
    density_bins[k] <<- igraph::edge_density(g)

    # Transitivity matching sna::gtrans convention:
    # - Directed: weak transitivity = closed_two_paths / total_two_paths
    # - Undirected: global clustering = 3*triangles / connected_triples
    # - Vacuous case (no two-paths/triples): 1 (matches sna)
    if (directed) {
      A <- as.matrix(igraph::as_adjacency_matrix(g))
      diag(A) <- 0
      A2 <- A %*% A
      diag(A2) <- 0
      total_two_paths <- sum(A2)
      if (total_two_paths > 0) {
        transitivity_bins[k] <<- sum(A2 * A) / total_two_paths
      } else {
        transitivity_bins[k] <<- 1
      }
    } else {
      tv <- igraph::transitivity(g, type = "global")
      transitivity_bins[k] <<- if (is.nan(tv)) 1 else tv
    }

    if (ne > 0L && n_vertices >= 2L) {
      centr_deg[k] <<- igraph::centr_degree(g)$centralization
      centr_betw_vec[k] <<- igraph::centr_betw(g)$centralization
      centr_clo_vec[k] <<- tryCatch({
        val <- igraph::centr_clo(g)$centralization
        if (is.finite(val)) val else 0
      }, error = function(e) 0)
    }

    n_comp[k] <<- igraph::components(g)$no

    md <- igraph::mean_distance(g, unconnected = TRUE)
    mean_dist[k] <<- if (is.nan(md)) NA_real_ else md

    diam_vec[k] <<- igraph::diameter(g)

    if (ne > 0L) {
      av <- igraph::assortativity_degree(g)
      assort_vec[k] <<- if (is.nan(av)) NA_real_ else av
    }

    if (directed) {
      rv <- igraph::reciprocity(g)
      recip_vec[k] <<- if (is.nan(rv)) 0 else rv
      dc <- igraph::dyad_census(g)
      mutual_vec[k] <<- as.integer(dc$mut)
      dyad_mat[, k] <<- c(
        as.integer(dc$mut), as.integer(dc$asym), as.integer(dc$null)
      )
      triad_mat[, k] <<- igraph::triad_census(g)
    }

    # Node-level
    cl <- igraph::closeness(g, normalized = TRUE)
    cl[is.nan(cl)] <- 0
    close_snap[, k] <<- cl

    betw_snap[, k] <<- igraph::betweenness(g, normalized = TRUE)

    if (ne > 0L) {
      eig_mat[, k] <<- tryCatch(
        igraph::eigen_centrality(g)$vector,
        error = function(e) rep(0, n_vertices)
      )
      pr_mat[, k] <<- igraph::page_rank(g)$vector
      hits <- igraph::hits_scores(g)
      hub_mat[, k] <<- hits$hub
      auth_mat[, k] <<- hits$authority
      con <- igraph::constraint(g)
      con[is.nan(con)] <- NA_real_
      constr_mat[, k] <<- con
    } else {
      pr_mat[, k] <<- rep(1 / n_vertices, n_vertices)
    }

    core_mat[, k] <<- igraph::coreness(g)

    0L
  }, integer(1L)))

  list(
    density_bins              = density_bins,
    reciprocity               = recip_vec,
    mutuality                 = mutual_vec,
    dyad_census               = dyad_mat,
    transitivity              = transitivity_bins,
    centralization_degree     = centr_deg,
    centralization_betweenness = centr_betw_vec,
    centralization_closeness  = centr_clo_vec,
    n_components              = n_comp,
    triad_census              = triad_mat,
    mean_distance             = mean_dist,
    diameter                  = diam_vec,
    assortativity             = assort_vec,
    closeness_snapshot        = close_snap,
    betweenness_snapshot      = betw_snap,
    eigenvector               = eig_mat,
    page_rank                 = pr_mat,
    hub_score                 = hub_mat,
    authority_score           = auth_mat,
    constraint                = constr_mat,
    coreness                  = core_mat
  )
}


# ===========================================================================
# S3: print.temporal_network
# ===========================================================================

#' @export
print.temporal_network <- function(x, ...) {
  dir_label <- if (x$directed) "directed" else "undirected"
  n_bins <- length(x$time_bins) - 1L
  mean_deg <- mean(x$degree)
  range_deg <- range(x$degree)
  mean_fwd <- mean(x$reachability_fwd)
  pct_fwd <- round(100 * mean_fwd / x$n_vertices, 1)
  mean_close <- mean(x$temporal_closeness, na.rm = TRUE)
  mean_form <- mean(x$formation)
  sd_form <- stats::sd(x$formation)
  mean_dur <- mean(x$edge_durations, na.rm = TRUE)
  valid_burst <- x$burstiness[!is.na(x$burstiness)]
  mean_burst <- if (length(valid_burst) > 0L) mean(valid_burst) else NA_real_

  cat(sprintf("Temporal Network [%s]\n", dir_label))
  cat(sprintf(
    "  Vertices: %d  |  Edge events: %d\n",
    x$n_vertices, x$n_edges
  ))
  cat(sprintf(
    "  Time range: [%.1f, %.1f]  |  Interval: %g  |  Bins: %d\n",
    x$time_range[1], x$time_range[2], x$time_interval, n_bins
  ))
  cat("\n")
  cat(sprintf(
    "  Mean degree: %.1f (range: %.0f-%.0f)\n",
    mean_deg, range_deg[1], range_deg[2]
  ))
  cat(sprintf(
    "  Mean forward reachability: %.1f / %d (%.1f%%)\n",
    mean_fwd, x$n_vertices, pct_fwd
  ))
  cat(sprintf("  Mean temporal closeness: %.4f\n", mean_close))
  cat(sprintf(
    "  Edge formation rate: %.1f/bin (sd: %.1f)\n",
    mean_form, sd_form
  ))
  cat(sprintf("  Mean edge duration: %.1f bins\n", mean_dur))
  if (!is.na(mean_burst)) {
    cat(sprintf("  Mean burstiness: %.3f\n", mean_burst))
  }
  mean_trans <- mean(x$transitivity, na.rm = TRUE)
  cat(sprintf("  Mean transitivity: %.3f\n", mean_trans))
  if (x$directed && !is.null(x$reciprocity)) {
    mean_recip <- mean(x$reciprocity, na.rm = TRUE)
    cat(sprintf("  Mean reciprocity: %.3f\n", mean_recip))
  }
  cat("\n")

  # Top temporal betweenness
  tw_between <- x$temporal_betweenness
  if (any(tw_between > 0)) {
    top_idx <- utils::head(order(tw_between, decreasing = TRUE), 3L)
    top_str <- vapply(top_idx, function(i) {
      sprintf("%s: %.4f", x$vertex_names[i], tw_between[i])
    }, character(1L))
    cat("  Top temporal betweenness:\n")
    cat(sprintf("    %s\n", paste(top_str, collapse = "  |  ")))
  }

  invisible(x)
}


# ===========================================================================
# S3: summary.temporal_network
# ===========================================================================

#' @export
summary.temporal_network <- function(object, ...) {
  cat("=== Temporal Network Summary ===\n\n")
  print.temporal_network(object)
  cat("\n--- Per-vertex centralities ---\n")

  centrality_df <- data.frame(
    vertex = object$vertex_names,
    reach_fwd = object$reachability_fwd,
    reach_bkwd = object$reachability_bkwd,
    closeness = round(object$temporal_closeness, 4),
    betweenness = round(object$temporal_betweenness, 4),
    burstiness = round(object$burstiness, 3),
    stringsAsFactors = FALSE
  )
  print(centrality_df, row.names = FALSE)

  cat("\n--- Degree time series (mean per bin) ---\n")
  mean_per_bin <- colMeans(object$degree)
  n_show <- min(length(mean_per_bin), 20L)
  cat(sprintf(
    "  First %d bins: %s\n",
    n_show,
    paste(round(mean_per_bin[seq_len(n_show)], 2), collapse = ", ")
  ))

  cat("\n--- Edge dynamics ---\n")
  cat(sprintf("  Total formation events: %d\n", sum(object$formation)))
  cat(sprintf("  Total dissolution events: %d\n", sum(object$dissolution)))
  cat(sprintf("  Unique edge pairs: %d\n", length(object$edge_durations)))
  cat(sprintf("  Mean edge duration: %.2f\n", mean(object$edge_durations)))
  cat(sprintf("  Median edge duration: %.2f\n", stats::median(object$edge_durations)))
  cat(sprintf("  SD edge duration: %.2f\n", stats::sd(object$edge_durations)))

  cat(sprintf("\n  Temporal density: %.4f\n", object$temporal_density))

  cat("\n--- Graph-level time series ---\n")
  cat(sprintf("  Mean density (per bin): %.4f\n", mean(object$density_bins)))
  cat(sprintf("  Mean transitivity: %.4f\n", mean(object$transitivity)))
  if (object$directed && !is.null(object$reciprocity)) {
    cat(sprintf("  Mean reciprocity: %.4f\n", mean(object$reciprocity)))
    cat(sprintf("  Mean mutuality: %.1f\n", mean(object$mutuality)))
  }
  cat(sprintf(
    "  Mean centralization (degree): %.4f\n",
    mean(object$centralization_degree)
  ))
  cat(sprintf("  Mean n_components: %.1f\n", mean(object$n_components)))

  invisible(object)
}


# ===========================================================================
# S3: plot.temporal_network
# ===========================================================================

#' @export
plot.temporal_network <- function(x, type = "degree", ...) {
  type <- match.arg(type, c(
    "degree", "formation", "reachability", "centrality",
    "burstiness", "duration", "iet", "snapshot",
    "reciprocity", "centralization", "eigenvector", "dyad_census",
    "proximity"
  ))

  switch(type,
    degree = .plot_degree(x, ...),
    formation = .plot_formation(x, ...),
    reachability = .plot_reachability(x, ...),
    centrality = .plot_centrality(x, ...),
    burstiness = .plot_burstiness(x, ...),
    duration = .plot_duration(x, ...),
    iet = .plot_iet(x, ...),
    snapshot = .plot_snapshot(x, ...),
    reciprocity = .plot_reciprocity(x, ...),
    centralization = .plot_centralization_ts(x, ...),
    eigenvector = .plot_eigenvector_ts(x, ...),
    dyad_census = .plot_dyad_census(x, ...),
    proximity = .plot_proximity_timeline(x, ...)
  )
}

#' @keywords internal
.plot_degree <- function(x, ...) {
  n_bins <- ncol(x$degree)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2
  mean_deg <- colMeans(x$degree)

  df <- data.frame(time = bin_mids, mean_degree = mean_deg)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$mean_degree)) +
    ggplot2::geom_line(linewidth = 0.8, color = "#2c7fb8") +
    ggplot2::geom_point(size = 1.5, color = "#2c7fb8") +
    ggplot2::labs(
      title = "Mean Degree Over Time",
      x = "Time", y = "Mean Degree"
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_formation <- function(x, ...) {
  n_bins <- length(x$formation)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2

  df <- data.frame(
    time = rep(bin_mids, 2L),
    count = c(x$formation, x$dissolution),
    type = rep(c("Formation", "Dissolution"), each = n_bins)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$time, y = .data$count, color = .data$type
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = c(
      "Formation" = "#31a354", "Dissolution" = "#e6550d"
    )) +
    ggplot2::labs(
      title = "Edge Formation & Dissolution",
      x = "Time", y = "Count", color = ""
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_reachability <- function(x, ...) {
  df <- data.frame(
    reachability = x$reachability_fwd,
    vertex = x$vertex_names
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$reachability)) +
    ggplot2::geom_histogram(
      bins = max(5L, min(30L, x$n_vertices)),
      fill = "#756bb1", color = "white"
    ) +
    ggplot2::labs(
      title = "Forward Reachability Distribution",
      x = "Reachable Vertices", y = "Count"
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_centrality <- function(x, ...) {
  # Bar chart ordered by betweenness
  ord <- order(x$temporal_betweenness, decreasing = TRUE)
  n_show <- min(x$n_vertices, 20L)
  top_idx <- ord[seq_len(n_show)]

  df <- data.frame(
    vertex = factor(
      x$vertex_names[top_idx],
      levels = x$vertex_names[top_idx]
    ),
    betweenness = x$temporal_betweenness[top_idx],
    closeness = x$temporal_closeness[top_idx]
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$vertex, y = .data$betweenness
  )) +
    ggplot2::geom_col(fill = "#de2d26") +
    ggplot2::labs(
      title = "Temporal Betweenness Centrality",
      x = "Vertex", y = "Betweenness"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  p
}

#' @keywords internal
.plot_burstiness <- function(x, ...) {
  valid <- x$burstiness[!is.na(x$burstiness)]
  if (length(valid) == 0L) {
    message("No burstiness values to plot (need vertices with >= 2 events)")
    return(invisible(NULL))
  }
  df <- data.frame(burstiness = valid)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$burstiness)) +
    ggplot2::geom_histogram(
      bins = max(5L, min(30L, length(valid))),
      fill = "#fd8d3c", color = "white"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      title = "Burstiness Distribution",
      subtitle = "B < 0: periodic | B = 0: Poisson | B > 0: bursty",
      x = "Burstiness Coefficient", y = "Count"
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_duration <- function(x, ...) {
  df <- data.frame(duration = as.numeric(x$edge_durations))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$duration)) +
    ggplot2::geom_histogram(
      bins = max(5L, min(30L, length(x$edge_durations))),
      fill = "#3182bd", color = "white"
    ) +
    ggplot2::labs(
      title = "Edge Duration Distribution",
      x = "Duration", y = "Count"
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_iet <- function(x, ...) {
  all_iet <- unlist(x$iet_vertex)
  if (length(all_iet) == 0L) {
    message("No inter-event times to plot")
    return(invisible(NULL))
  }
  df <- data.frame(iet = all_iet)

  # Log-log if heavy-tailed (range > 2 orders of magnitude)
  use_log <- (max(all_iet) / min(all_iet[all_iet > 0])) > 100

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iet)) +
    ggplot2::geom_histogram(
      bins = max(5L, min(50L, length(all_iet))),
      fill = "#74c476", color = "white"
    ) +
    ggplot2::labs(
      title = "Inter-Event Time Distribution",
      x = "Inter-Event Time", y = "Count"
    ) +
    ggplot2::theme_minimal()

  if (use_log) {
    p <- p +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "Inter-Event Time (log scale)")
  }
  p
}

#' @keywords internal
.plot_snapshot <- function(x, n_snapshots = 4L, ...) {
  n_bins <- length(x$snapshots)
  if (n_bins == 0L) {
    message("No snapshots available")
    return(invisible(NULL))
  }

  # Select evenly spaced snapshots
  idx <- unique(round(seq(1, n_bins, length.out = min(n_snapshots, n_bins))))

  old_par <- graphics::par(
    mfrow = c(ceiling(length(idx) / 2), 2L),
    mar = c(1, 1, 2, 1)
  )
  on.exit(graphics::par(old_par), add = TRUE)

  for (i in idx) {
    g <- x$snapshots[[i]]
    bin_label <- sprintf("t = [%.1f, %.1f)", x$time_bins[i], x$time_bins[i + 1L])
    igraph::plot.igraph(
      g,
      main = bin_label,
      vertex.size = 15,
      vertex.label.cex = 0.7,
      edge.arrow.size = 0.4,
      ...
    )
  }
  invisible(NULL)
}

#' @keywords internal
.plot_reciprocity <- function(x, ...) {
  n_bins <- length(x$transitivity)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2

  df_trans <- data.frame(
    time = bin_mids,
    value = x$transitivity,
    metric = rep("Transitivity", n_bins)
  )

  if (x$directed && !is.null(x$reciprocity)) {
    df_recip <- data.frame(
      time = bin_mids,
      value = x$reciprocity,
      metric = rep("Reciprocity", n_bins)
    )
    df <- rbind(df_trans, df_recip)
  } else {
    df <- df_trans
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$time, y = .data$value, color = .data$metric
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = c(
      "Reciprocity" = "#e41a1c", "Transitivity" = "#377eb8"
    )) +
    ggplot2::labs(
      title = "Reciprocity & Transitivity Over Time",
      x = "Time", y = "Value", color = ""
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_centralization_ts <- function(x, ...) {
  n_bins <- length(x$centralization_degree)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2

  df <- data.frame(
    time = rep(bin_mids, 3L),
    value = c(
      x$centralization_degree,
      x$centralization_betweenness,
      x$centralization_closeness
    ),
    measure = rep(
      c("Degree", "Betweenness", "Closeness"), each = n_bins
    )
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$time, y = .data$value, color = .data$measure
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = c(
      "Degree" = "#e41a1c", "Betweenness" = "#377eb8",
      "Closeness" = "#4daf4a"
    )) +
    ggplot2::labs(
      title = "Centralization Over Time",
      x = "Time", y = "Centralization", color = ""
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_eigenvector_ts <- function(x, ...) {
  n_bins <- ncol(x$eigenvector)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2
  mean_eig <- colMeans(x$eigenvector)

  df <- data.frame(time = bin_mids, mean_eigenvector = mean_eig)
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$time, y = .data$mean_eigenvector
  )) +
    ggplot2::geom_line(linewidth = 0.8, color = "#984ea3") +
    ggplot2::geom_point(size = 1.5, color = "#984ea3") +
    ggplot2::labs(
      title = "Mean Eigenvector Centrality Over Time",
      x = "Time", y = "Mean Eigenvector Centrality"
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_dyad_census <- function(x, ...) {
  if (!x$directed || is.null(x$dyad_census)) {
    message("Dyad census is only available for directed networks")
    return(invisible(NULL))
  }

  n_bins <- ncol(x$dyad_census)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2

  df <- data.frame(
    time = rep(bin_mids, 3L),
    count = c(x$dyad_census["Mutual", ],
              x$dyad_census["Asymmetric", ],
              x$dyad_census["Null", ]),
    type = rep(c("Mutual", "Asymmetric", "Null"), each = n_bins)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$time, y = .data$count, fill = .data$type
  )) +
    ggplot2::geom_area(alpha = 0.7) +
    ggplot2::scale_fill_manual(values = c(
      "Mutual" = "#2ca02c", "Asymmetric" = "#ff7f0e", "Null" = "#d62728"
    )) +
    ggplot2::labs(
      title = "Dyad Census Over Time",
      x = "Time", y = "Count", fill = ""
    ) +
    ggplot2::theme_minimal()
  p
}


# ===========================================================================
# Internal: proximity timeline plot
# ===========================================================================

#' @keywords internal
.resolve_proximity_metric <- function(x, metric) {
  metric <- match.arg(metric, c(
    "eigenvector", "degree", "closeness", "betweenness",
    "page_rank", "hub_score", "authority_score", "proximity"
  ))

  metric_map <- list(
    eigenvector      = "eigenvector",
    degree           = "degree",
    closeness        = "closeness_snapshot",
    betweenness      = "betweenness_snapshot",
    page_rank        = "page_rank",
    hub_score        = "hub_score",
    authority_score  = "authority_score"
  )

  if (metric == "proximity") {
    return(.compute_proximity_mds(x))
  }

  field <- metric_map[[metric]]
  mat <- x[[field]]
  if (is.null(mat)) {
    stop(sprintf("Metric '%s' (field '%s') not found in temporal_network object",
                 metric, field))
  }
  mat
}


#' @keywords internal
.compute_proximity_mds <- function(x) {
  n_v <- x$n_vertices
  n_bins <- length(x$snapshots)
  pos_mat <- matrix(0, nrow = n_v, ncol = n_bins)
  rownames(pos_mat) <- x$vertex_names

  # Compute strength (weighted degree) at each bin
  str_mat <- vapply(seq_len(n_bins), function(k) {
    g <- x$snapshots[[k]]
    if (igraph::ecount(g) == 0L) return(rep(0, n_v))
    igraph::strength(g, mode = "all")
  }, numeric(n_v))
  if (is.null(dim(str_mat))) str_mat <- matrix(str_mat, nrow = n_v)
  rownames(str_mat) <- x$vertex_names

  # Interleaved slot positions: 0, +s, -s, +2s, -2s, ...
  n_half <- ceiling((n_v - 1L) / 2)
  step <- if (n_half > 0) 1 / n_half else 1
  slots <- numeric(n_v)
  slots[1L] <- 0
  idx <- 2L
  for (i in seq_len(n_half)) {
    slots[idx] <- i * step
    idx <- idx + 1L
    if (idx <= n_v) {
      slots[idx] <- -i * step
      idx <- idx + 1L
    }
  }

  # At each bin: sort by strength, assign to interleaved slots
  vapply(seq_len(n_bins), function(k) {
    str_k <- str_mat[, k]
    sorted_idx <- order(-str_k, method = "radix")
    pos_mat[sorted_idx, k] <<- slots
    0L
  }, integer(1L))

  pos_mat
}


#' Spline-interpolate a metric matrix to higher temporal resolution.
#' @param metric_mat n_vertices x n_bins matrix
#' @param bin_mids  numeric vector of original time points (length n_bins)
#' @param fine_times numeric vector of target time points
#' @return n_vertices x length(fine_times) matrix
#' @keywords internal
.interpolate_metric <- function(metric_mat, bin_mids, fine_times) {
  n_v <- nrow(metric_mat)
  n_bins <- ncol(metric_mat)

  # Need >= 4 bins for cubic spline; fall back to linear for fewer
  method <- if (n_bins >= 4L) "natural" else "linear"

  fine_mat <- vapply(seq_len(n_v), function(i) {
    y <- metric_mat[i, ]
    # If all identical (zero variance), return constant
    if (max(y) == min(y)) return(rep(y[1L], length(fine_times)))
    if (method == "linear") {
      stats::approx(bin_mids, y, xout = fine_times, rule = 2L)$y
    } else {
      stats::spline(bin_mids, y, xout = fine_times, method = "natural")$y
    }
  }, numeric(length(fine_times)))

  # vapply returns n_fine x n_v; transpose to n_v x n_fine
  t(fine_mat)
}


#' @keywords internal
.plot_proximity_timeline <- function(x,
                                     metric = "proximity",
                                     size_metric = "degree",
                                     vertex_group = NULL,
                                     vertex_color = NULL,
                                     vertex_label = NULL,
                                     labels_at = NULL,
                                     label_size = 3,
                                     line_width = c(0.5, 3),
                                     alpha = 0.4,
                                     highlight = NULL,
                                     n_slices = NULL,
                                     spline = TRUE,
                                     ...) {
  # --- Resolve metric matrices ---
  metric_mat <- .resolve_proximity_metric(x, metric)
  n_v <- nrow(metric_mat)
  n_bins <- ncol(metric_mat)
  bin_mids <- x$time_bins[seq_len(n_bins)] + x$time_interval / 2
  vnames <- x$vertex_names

  has_size <- !(is.null(size_metric) || identical(size_metric, FALSE))
  if (has_size) {
    size_mat <- .resolve_proximity_metric(x, size_metric)
  }

  # --- Interpolate to high-resolution time series ---
  if (is.null(n_slices)) n_slices <- max(200L, n_bins * 10L)

  do_spline <- spline && n_bins >= 3L && n_slices > n_bins
  if (do_spline) {
    fine_times <- seq(min(bin_mids), max(bin_mids), length.out = n_slices)
    plot_mat <- .interpolate_metric(metric_mat, bin_mids, fine_times)
    rownames(plot_mat) <- vnames
    # Clamp spline overshoot to observed range per vertex
    metric_range <- range(metric_mat, na.rm = TRUE)
    plot_mat[plot_mat < metric_range[1L]] <- metric_range[1L]
    plot_mat[plot_mat > metric_range[2L]] <- metric_range[2L]
    if (has_size) {
      size_fine <- .interpolate_metric(size_mat, bin_mids, fine_times)
      size_fine[size_fine < 0] <- 0
    }
    plot_times <- fine_times
    n_plot <- n_slices
  } else {
    plot_times <- bin_mids
    plot_mat <- metric_mat
    if (has_size) size_fine <- size_mat
    n_plot <- n_bins
  }

  # --- Default labels ---
  if (is.null(vertex_label)) vertex_label <- vnames

  # --- Build segment data (consecutive point pairs per vertex) ---
  seg_list <- lapply(seq_len(n_v), function(i) {
    n_seg <- n_plot - 1L
    if (n_seg < 1L) return(NULL)
    seg <- data.frame(
      x    = plot_times[-n_plot],
      xend = plot_times[-1L],
      y    = plot_mat[i, -n_plot],
      yend = plot_mat[i, -1L],
      vertex = vnames[i],
      label  = vertex_label[i],
      stringsAsFactors = FALSE
    )
    if (has_size) {
      seg$thickness <- (size_fine[i, -n_plot] + size_fine[i, -1L]) / 2
    }
    seg
  })
  seg_df <- do.call(rbind, seg_list)

  # Handle edge case: no segments (single time bin)
  if (is.null(seg_df) || nrow(seg_df) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Proximity Timeline", x = "Time", y = metric)
    return(p)
  }

  # --- Pre-scale linewidth (base R, no scales dependency) ---
  if (has_size) {
    rng <- range(seg_df$thickness, na.rm = TRUE)
    if (rng[1L] == rng[2L]) {
      seg_df$lw <- rep(mean(line_width), nrow(seg_df))
    } else {
      seg_df$lw <- (seg_df$thickness - rng[1L]) / (rng[2L] - rng[1L]) *
        (line_width[2L] - line_width[1L]) + line_width[1L]
    }
  } else {
    seg_df$lw <- rep(
      if (length(line_width) == 2L) mean(line_width) else line_width[1L],
      nrow(seg_df)
    )
  }

  # --- Color palette: Okabe-Ito default ---
  okabe_ito <- c(
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7",
    "#56B4E9", "#882255", "#332288", "#117733", "#AA4499"
  )

  if (!is.null(vertex_group)) {
    stopifnot(length(vertex_group) == n_v)
    grp_map <- stats::setNames(as.character(vertex_group), vnames)
    seg_df$color_group <- grp_map[seg_df$vertex]
    unique_groups <- unique(seg_df$color_group)
    n_groups <- length(unique_groups)
    if (!is.null(vertex_color)) {
      stopifnot(length(vertex_color) == n_groups)
      color_map <- stats::setNames(vertex_color, unique_groups)
    } else {
      color_map <- stats::setNames(rep_len(okabe_ito, n_groups), unique_groups)
    }
  } else {
    seg_df$color_group <- seg_df$vertex
    if (!is.null(vertex_color)) {
      stopifnot(length(vertex_color) == n_v)
      color_map <- stats::setNames(vertex_color, vnames)
    } else {
      color_map <- stats::setNames(rep_len(okabe_ito, n_v), vnames)
    }
  }

  # Pre-map colors into data.frame
  seg_df$color <- color_map[seg_df$color_group]

  # --- Highlight: set per-segment alpha ---
  if (!is.null(highlight)) {
    seg_df$seg_alpha <- ifelse(seg_df$vertex %in% highlight, alpha, 0.15)
  } else {
    seg_df$seg_alpha <- alpha
  }

  # --- Build ggplot ---
  p <- ggplot2::ggplot(seg_df) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$x, xend = .data$xend,
        y = .data$y, yend = .data$yend
      ),
      color = seg_df$color,
      linewidth = seg_df$lw,
      alpha = seg_df$seg_alpha,
      lineend = "round"
    )

  # --- Labels ---
  # Resolve color keys for label data
  .label_colors <- function(vnames_vec) {
    if (!is.null(vertex_group)) {
      grp <- stats::setNames(as.character(vertex_group), vnames)[vnames_vec]
      color_map[grp]
    } else {
      color_map[vnames_vec]
    }
  }

  if (is.null(labels_at)) {
    # Direct endpoint labels at last time point
    last_idx <- n_plot
    label_df <- data.frame(
      x = plot_times[last_idx],
      y = plot_mat[, last_idx],
      label = vertex_label,
      stringsAsFactors = FALSE
    )
    label_df$color <- .label_colors(vnames)
    label_df$lbl_alpha <- if (!is.null(highlight)) {
      ifelse(vnames %in% highlight, 1, 0.3)
    } else {
      rep(1, n_v)
    }
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      color = label_df$color,
      alpha = label_df$lbl_alpha,
      size = label_size, hjust = -0.15, vjust = 0.5,
      show.legend = FALSE
    )
  } else {
    # Labels at specified bins (suppresses endpoint labels)
    labels_at <- labels_at[labels_at >= 1L & labels_at <= n_bins]
    if (length(labels_at) > 0L) {
      lbl_list <- lapply(labels_at, function(k) {
        data.frame(
          x     = bin_mids[k],
          y     = metric_mat[, k],
          label = vertex_label,
          color = .label_colors(vnames),
          stringsAsFactors = FALSE
        )
      })
      label_df <- do.call(rbind, lbl_list)
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        color = label_df$color,
        size = label_size, hjust = -0.15, vjust = 0.5,
        show.legend = FALSE
      )
    }
  }

  # --- Axis labels and theme ---
  metric_labels <- c(
    eigenvector     = "Eigenvector Centrality",
    degree          = "Degree",
    closeness       = "Closeness Centrality",
    betweenness     = "Betweenness Centrality",
    page_rank       = "PageRank",
    hub_score       = "Hub Score",
    authority_score = "Authority Score",
    proximity       = "Network Proximity"
  )
  y_label <- metric_labels[[metric]]

  p <- p +
    ggplot2::labs(
      title = sprintf("Proximity Timeline (%s)", y_label),
      x = "Time", y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p
}


#' @keywords internal
.proximity_size_label <- function(size_metric) {
  labels <- c(
    eigenvector     = "Eigenvector",
    degree          = "Degree",
    closeness       = "Closeness",
    betweenness     = "Betweenness",
    page_rank       = "PageRank",
    hub_score       = "Hub Score",
    authority_score = "Authority",
    proximity       = "Proximity"
  )
  labels[[size_metric]]
}


# ===========================================================================
# Exported: temporal_paths
# ===========================================================================

#' Compute Temporal Paths from a Source Vertex
#'
#' @description
#' Runs temporal BFS from a specified vertex and returns earliest-arrival
#' information for all reachable vertices.
#'
#' @param x A \code{temporal_network} object.
#' @param from Character or integer. Source vertex name or ID.
#' @param start Numeric or NULL. Start time for BFS. Default: minimum onset.
#'
#' @return A data frame of class \code{"temporal_paths"} with columns:
#'   \code{vertex}, \code{arrival_time}, \code{previous}, \code{n_hops}.
#'   Attributes \code{source} and \code{start_time} record the BFS origin.
#'
#' @export
temporal_paths <- function(x, from, start = NULL) {
  stopifnot(inherits(x, "temporal_network"))

  # Resolve vertex
  if (is.character(from)) {
    from_id <- match(from, x$vertex_names)
    if (is.na(from_id)) stop(sprintf("Vertex '%s' not found", from))
  } else {
    from_id <- as.integer(from)
    stopifnot(from_id >= 1L, from_id <= x$n_vertices)
  }

  if (is.null(start)) start <- x$time_range[1]

  # For undirected networks, add reverse edges so BFS traverses both directions
  bfs_edges <- x$edges
  if (!x$directed) {
    rev_edges <- data.frame(
      from = bfs_edges$to, to = bfs_edges$from,
      onset = bfs_edges$onset, terminus = bfs_edges$terminus,
      stringsAsFactors = FALSE
    )
    bfs_edges <- rbind(bfs_edges, rev_edges)
    bfs_edges <- bfs_edges[order(bfs_edges$onset), , drop = FALSE]
    rownames(bfs_edges) <- NULL
  }

  bfs <- .temporal_bfs(bfs_edges, x$n_vertices, from_id, start_time = start)

  # Compute n_hops by tracing paths
  n_hops <- integer(x$n_vertices)
  for (v in seq_len(x$n_vertices)) {
    if (v == from_id) {
      n_hops[v] <- 0L
      next
    }
    if (is.infinite(bfs$earliest[v])) {
      n_hops[v] <- NA_integer_
      next
    }
    path <- .trace_path(bfs$previous, from_id, v)
    n_hops[v] <- if (length(path) > 0L) length(path) - 1L else NA_integer_
  }

  result <- data.frame(
    vertex = x$vertex_names,
    arrival_time = bfs$earliest,
    previous = ifelse(
      is.na(bfs$previous), NA_character_,
      x$vertex_names[bfs$previous]
    ),
    n_hops = n_hops,
    stringsAsFactors = FALSE
  )
  attr(result, "source") <- x$vertex_names[from_id]
  attr(result, "start_time") <- start
  class(result) <- c("temporal_paths", "data.frame")
  result
}


# ===========================================================================
# S3: plot.temporal_paths
# ===========================================================================

#' Plot Temporal Path Tree
#'
#' @description
#' Visualizes the earliest-arrival path tree from a temporal BFS as a
#' network graph. Uses a tree layout with the source vertex as root.
#' Edges are drawn as thick directed arrows with arrival-time labels.
#' The source vertex is highlighted with a distinct border.
#'
#' @param x A \code{temporal_paths} object returned by \code{temporal_paths()}.
#' @param ... Additional arguments passed to \code{igraph::plot.igraph()}.
#'
#' @return Invisibly returns the plot coordinates (matrix with x, y columns).
#'
#' @export
plot.temporal_paths <- function(x, ...) {
  source_name <- attr(x, "source")

  # Separate reachable (finite arrival) from unreachable
  reachable <- x[is.finite(x$arrival_time), , drop = FALSE]

  if (nrow(reachable) == 0L) {
    message("No reachable vertices to plot")
    return(invisible(NULL))
  }

  # Build edge list from predecessor links (predecessor -> vertex)
  has_pred <- !is.na(reachable$previous)
  if (any(has_pred)) {
    el <- cbind(
      reachable$previous[has_pred],
      reachable$vertex[has_pred]
    )
    edge_labels <- round(reachable$arrival_time[has_pred], 1)
  } else {
    el <- matrix(character(0), ncol = 2)
    edge_labels <- character(0)
  }

  # Create igraph tree
  g <- igraph::make_empty_graph(directed = TRUE)
  g <- igraph::add_vertices(g, nrow(reachable), name = reachable$vertex)

  if (nrow(el) > 0L) {
    g <- igraph::add_edges(g, as.vector(t(el)))
  }

  # Tree layout with source as root
  root_idx <- which(igraph::V(g)$name == source_name)
  coords <- igraph::layout_as_tree(g, root = root_idx, mode = "out")

  # Vertex styling
  is_source <- igraph::V(g)$name == source_name
  v_color <- ifelse(is_source, "#e41a1c", "#2c7fb8")
  v_frame <- ifelse(is_source, "#e41a1c", "#2c7fb8")
  v_size <- ifelse(is_source, 20, 15)
  v_label_color <- "white"

  # Default arguments, user can override via ...
  defaults <- list(
    x = g,
    layout = coords,
    vertex.color = v_color,
    vertex.frame.color = v_frame,
    vertex.size = v_size,
    vertex.label.color = v_label_color,
    vertex.label.cex = 0.9,
    vertex.label.font = 2,
    edge.color = "#e41a1c",
    edge.width = 3,
    edge.arrow.size = 0.6,
    edge.label = edge_labels,
    edge.label.color = "#e41a1c",
    edge.label.cex = 0.8,
    main = sprintf("Temporal Paths from %s", source_name),
    margin = c(0, 0, 0, 0)
  )

  # Merge user args (override defaults)
  user_args <- list(...)
  final_args <- modifyList(defaults, user_args)

  do.call(igraph::plot.igraph, final_args)
  invisible(coords)
}


# ===========================================================================
# Exported: extract_snapshot
# ===========================================================================

#' Extract Static Network Snapshot at a Time Point
#'
#' @description
#' Extracts a static igraph network containing all edges active at a given
#' time point.
#'
#' @param x A \code{temporal_network} object.
#' @param at Numeric. Time point at which to extract the snapshot.
#'
#' @return An \code{igraph} object with edges active at time \code{at}.
#'
#' @export
extract_snapshot <- function(x, at) {
  stopifnot(inherits(x, "temporal_network"))
  stopifnot(is.numeric(at), length(at) == 1L)

  # Find active edges at time `at`: onset <= at < terminus
  active <- which(x$edges$onset <= at & x$edges$terminus > at)

  if (length(active) > 0L) {
    el <- cbind(x$edges$from[active], x$edges$to[active])
    g <- igraph::graph_from_edgelist(el, directed = x$directed)
    # Ensure all vertices present
    missing_v <- setdiff(seq_len(x$n_vertices), igraph::V(g))
    if (length(missing_v) > 0L) {
      g <- igraph::add_vertices(g, length(missing_v))
    }
    igraph::V(g)$name <- x$vertex_names[seq_len(igraph::vcount(g))]
  } else {
    g <- igraph::make_empty_graph(n = x$n_vertices, directed = x$directed)
    igraph::V(g)$name <- x$vertex_names
  }

  g
}
