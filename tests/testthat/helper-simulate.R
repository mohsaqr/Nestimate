# Minimal simulation helpers for Nestimate tests
# These replace Saqrlab simulation functions with standalone versions

simulate_data <- function(type = "mlvar", seed = NULL, n_subjects = 20,
                          d = 6, n_obs = 40, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (type == "mlvar") {
    vars <- paste0("V", seq_len(d))
    rows <- n_subjects * n_obs

    # Generate known temporal and contemporaneous matrices
    true_B <- matrix(0, d, d)
    diag(true_B) <- runif(d, 0.2, 0.5)
    # Add a few off-diagonal temporal effects
    for (i in seq_len(min(d - 1, 3))) {
      true_B[i, i + 1] <- runif(1, 0.1, 0.3)
    }

    true_Theta <- diag(d)
    for (i in seq_len(min(d - 1, 2))) {
      true_Theta[i, i + 1] <- true_Theta[i + 1, i] <- runif(1, 0.1, 0.3)
    }

    # Generate VAR data
    between_means <- matrix(rnorm(n_subjects * d, 0, 1), n_subjects, d)
    df <- data.frame(
      id = rep(seq_len(n_subjects), each = n_obs),
      day = rep(rep(seq_len(n_obs %/% 10), each = 10), n_subjects),
      beep = rep(rep(seq_len(10), n_obs %/% 10), n_subjects)
    )
    y_mat <- matrix(0, rows, d)
    for (s in seq_len(n_subjects)) {
      idx <- ((s - 1) * n_obs + 1):(s * n_obs)
      y_mat[idx[1], ] <- between_means[s, ] + rnorm(d)
      for (t in 2:n_obs) {
        innovation <- rnorm(d)
        y_mat[idx[t], ] <- between_means[s, ] +
          true_B %*% (y_mat[idx[t - 1], ] - between_means[s, ]) + innovation
      }
    }
    for (j in seq_len(d)) df[[vars[j]]] <- y_mat[, j]

    rownames(true_B) <- colnames(true_B) <- vars
    rownames(true_Theta) <- colnames(true_Theta) <- vars

    attr(df, "vars") <- vars
    attr(df, "type") <- "mlvar"
    attr(df, "true_temporal") <- true_B
    attr(df, "true_contemporaneous") <- true_Theta
    attr(df, "info") <- list()
    return(df)
  }
  stop(sprintf("simulate_data type '%s' not supported in Nestimate test helper", type))
}

#' Generate random wide-format sequence data
#' @noRd
simulate_sequences <- function(n_actors = 10, n_states = 5,
                               seq_length = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- lapply(seq_len(n_actors), function(i) {
    sample(states, seq_length, replace = TRUE)
  })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_length))
  df
}

#' Generate random continuous data for association-network tests
#' @noRd
simulate_continuous <- function(n = 100, p = 5, rho = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- diag(p)
  for (i in seq_len(p - 1L)) {
    Sigma[i, i + 1L] <- rho
    Sigma[i + 1L, i] <- rho
  }
  L <- chol(Sigma)
  mat <- matrix(rnorm(n * p), n, p) %*% L
  df <- as.data.frame(mat)
  colnames(df) <- paste0("V", seq_len(p))
  df
}

#' Generate a random binary incidence matrix for hypergraph equivalence tests.
#'
#' Columns are hyperedges, rows are nodes. Each hyperedge's size is drawn
#' uniformly from `edge_size_range`, members are sampled without replacement
#' from the node set. `density` adds per-edge jitter: with probability
#' (1 - density) an edge is skipped (returned as a degenerate size-1 edge, or
#' dropped if too small), ensuring size variability across runs.
#'
#' @param n_nodes Integer >= 3. Number of nodes.
#' @param n_edges Integer >= 1. Number of hyperedges.
#' @param edge_size_range Integer vector length 2. Min/max hyperedge size.
#'   Must satisfy 2 <= min <= max <= n_nodes.
#' @param density Numeric in (0, 1]. Probability that a given edge is kept.
#' @param seed Integer or NULL. Random seed.
#' @return list(incidence = matrix(n_nodes, n_edges_kept),
#'              hyperedges = list of integer node indices,
#'              nodes = character vector of node names)
#' @noRd
simulate_hypergraph_incidence <- function(n_nodes = 15, n_edges = 20,
                                          edge_size_range = c(2L, 4L),
                                          density = 0.8,
                                          seed = NULL) {
  stopifnot(n_nodes >= 3L, n_edges >= 1L,
            length(edge_size_range) == 2L,
            edge_size_range[1] >= 2L,
            edge_size_range[2] >= edge_size_range[1],
            edge_size_range[2] <= n_nodes,
            density > 0, density <= 1)
  if (!is.null(seed)) set.seed(seed)
  nodes <- paste0("n", seq_len(n_nodes))
  sizes <- sample(edge_size_range[1]:edge_size_range[2], n_edges, replace = TRUE)
  keep <- stats::runif(n_edges) < density
  edges <- mapply(function(sz, k) {
    if (!k) return(integer(0))
    sort(sample.int(n_nodes, sz))
  }, sizes, keep, SIMPLIFY = FALSE)
  edges <- edges[vapply(edges, length, integer(1)) >= 2L]
  if (length(edges) == 0L) {
    edges <- list(sort(sample.int(n_nodes, edge_size_range[1])))
  }
  inc <- matrix(0L, nrow = n_nodes, ncol = length(edges),
                dimnames = list(nodes, paste0("e", seq_along(edges))))
  for (j in seq_along(edges)) inc[edges[[j]], j] <- 1L
  list(incidence = inc, hyperedges = edges, nodes = nodes)
}
