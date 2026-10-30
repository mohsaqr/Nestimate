# ---- Hypergraph eigenvector centralities (HON-5) -------------------------
# Three variants from Benson (2019, arXiv:1807.09644):
#   * clique-motif ("CEC") - standard eigenvector centrality on the
#                            clique-expanded pairwise graph
#   * Z-eigenvector ("Z")  - linear tensor eigenvector
#   * H-eigenvector ("H")  - H-eigenvector (power-k-1 recurrence)
#
# CEC is the Perron vector of each connected component of the clique
# adjacency, solved directly (power iteration oscillates on bipartite
# components). ZEC and HEC are shifted tensor power iterations, defined
# for uniform hypergraphs only.

#' Hypergraph eigenvector centralities
#'
#' Computes one or more eigenvector-style centralities on a
#' [net_hypergraph][build_hypergraph]: *clique-motif* (CEC),
#' *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC). Each variant
#' captures influence differently - CEC flattens group structure via
#' clique expansion, while ZEC and HEC propagate through the
#' higher-order groups directly.
#'
#' @param hg A `net_hypergraph` (from [build_hypergraph()] or
#'   [bipartite_groups()]).
#' @param type Character vector, any subset of `c("clique", "Z", "H")`.
#'   Default `"clique"`. `"Z"` and `"H"` need a uniform hypergraph (every
#'   hyperedge with two or more members has the same size).
#' @param max_iter Maximum number of power-iteration steps of the `"Z"`
#'   and `"H"` iterations. Default `1000`. `"clique"` is solved directly.
#' @param tol Convergence tolerance on the L1 change between successive
#'   iterates of `"Z"` and `"H"`. Default `1e-8`. An iteration that does not
#'   reach it within `max_iter` steps warns with class
#'   `nestimate_no_converge`.
#' @param normalize Logical. If `TRUE` (default), each returned
#'   centrality vector is L2-normalized to unit norm (compatible with
#'   `igraph::eigen_centrality()`'s scale for type `"clique"`). If
#'   `FALSE`, the vector is rescaled so its largest absolute entry is 1.
#'
#' @return A named list, one component per requested `type`, in the order
#'   given by `type`. Each component is a numeric vector of length
#'   `hg$n_nodes` named by `hg$nodes`. A hypergraph with no hyperedges
#'   yields all-zero vectors. Requesting `"Z"` or `"H"` on a non-uniform
#'   hypergraph raises an error of class `nestimate_nonuniform_hypergraph`.
#'
#' @details
#' **Clique-motif eigenvector centrality (CEC)**: forms the
#' clique-expanded pairwise graph \eqn{W} where
#' \eqn{W_{ij} = |\{e : i, j \in e\}|} and returns the leading
#' eigenvector of \eqn{W}. Equivalent to running
#' `igraph::eigen_centrality()` on [clique_expansion()] output. The vector
#' is solved directly from a symmetric eigendecomposition of each connected
#' component (Perron-Frobenius makes the leading eigenvalue of a connected
#' component simple), not by power iteration, which oscillates when the
#' clique expansion is bipartite (for example a star of pairs). When several
#' components share the largest eigenvalue, the result is the projection of
#' the uniform vector onto their joint leading eigenspace; components with a
#' smaller leading eigenvalue score zero.
#'
#' **Z-eigenvector centrality (ZEC)**: solves the eigen-equation on the
#' hyperedge tensor,
#' \deqn{\lambda\, x_i \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j,}
#' by shifted power iteration. Tensor products and normalization are
#' evaluated in log space to retain finite, normalized centralities for
#' large hyperedges.
#'
#' **H-eigenvector centrality (HEC)**: solves
#' \deqn{\lambda\, x_i^{k-1} \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j.}
#'
#' **Uniform hypergraphs only.** Benson (2019) defines ZEC and HEC for
#' \eqn{k}-uniform hypergraphs. With mixed sizes the right-hand side mixes
#' polynomial degrees, so the equation has no scale-free solution: the
#' vector satisfies it only at the normalisation chosen inside the
#' iteration. Hyperedges with one member contribute nothing and are ignored
#' by the uniformity check. For mixed sizes use `"clique"`, or keep one
#' size with [build_hypergraph()] (`include_pairwise = FALSE`).
#'
#' @seealso [build_hypergraph()], [clique_expansion()],
#'   [hypergraph_measures()].
#'
#' @examples
#' df <- data.frame(
#'   player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
#' )
#' hg <- bipartite_groups(df, "player", "session")
#' hypergraph_centrality(hg)
#'
#' # Z and H need a uniform hypergraph: here every session has three players
#' df3 <- data.frame(
#'   player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S2", "S3", "S3", "S3")
#' )
#' hypergraph_centrality(bipartite_groups(df3, "player", "session"),
#'                       type = c("clique", "Z", "H"))
#'
#' @references
#' Benson, A. R. (2019). Three hypergraph eigenvector centralities.
#' \emph{SIAM Journal on Mathematics of Data Science} 1(2), 293-312.
#' arXiv:1807.09644.
#'
#' @note The `"clique"` (CEC) variant is validated against
#'   `igraph::eigen_centrality` (cosine ~ 1). The `"Z"` and `"H"` variants are
#'   **(experimental)** - validated only against a clean-room list-based
#'   tensor power iteration (same operator, different loop structure); no
#'   R package exposes tensor eigenvectors as a primitive for independent
#'   comparison.
#'
#' @export
hypergraph_centrality <- function(hg,
                                   type     = "clique",
                                   max_iter = 1000L,
                                   tol      = 1e-8,
                                   normalize = TRUE) {
  stopifnot(
    inherits(hg, "net_hypergraph"),
    is.numeric(max_iter), length(max_iter) == 1L, max_iter > 0,
    is.numeric(tol), length(tol) == 1L, tol > 0,
    is.logical(normalize), length(normalize) == 1L
  )
  type <- match.arg(type, choices = c("clique", "Z", "H"), several.ok = TRUE)

  n     <- hg$n_nodes
  nodes <- hg$nodes

  # Degenerate: no nodes
  if (n == 0L) {
    out <- stats::setNames(rep(list(numeric(0L)), length(type)), type)
    return(out)
  }

  # Degenerate: no hyperedges -> zero centralities everywhere
  if (hg$n_hyperedges == 0L) {
    zero_vec <- stats::setNames(rep(0, n), nodes)
    out <- stats::setNames(rep(list(zero_vec), length(type)), type)
    return(out)
  }

  # Shared initialization: uniform positive vector
  x0 <- rep(1 / sqrt(n), n)

  # Hyperedges as integer vector lists; pre-compute for Z/H
  hyperedges <- hg$hyperedges
  edge_sizes <- vapply(hyperedges, length, integer(1L))
  k_max      <- max(edge_sizes)

  # Z and H are defined for uniform hypergraphs; singletons contribute nothing
  tensor_sizes <- unique(edge_sizes[edge_sizes >= 2L])
  if (any(c("Z", "H") %in% type) && length(tensor_sizes) > 1L) {
    stop(errorCondition(
      sprintf(paste0(
        "`type = \"Z\"` and `\"H\"` need a uniform hypergraph; ",
        "hyperedge sizes here are %s. Use `type = \"clique\"`, or build ",
        "the hypergraph with one hyperedge size."),
        paste(sort(tensor_sizes), collapse = ", ")),
      class = "nestimate_nonuniform_hypergraph", call = NULL
    ))
  }

  out <- list()

  # ---- CEC: Perron vector of the clique expansion W, per component ----
  if ("clique" %in% type) {
    W <- tcrossprod((hg$incidence != 0) * 1.0)
    diag(W) <- 0
    x <- .hg_perron_vector(.hg_component_spectra(W), x0)
    if (!normalize && any(x != 0)) x <- x / max(abs(x))
    out$clique <- stats::setNames(x, nodes)
  }

  # ---- ZEC: lambda x = sum_{e contains i} prod_{j in e,j!=i} x_j ----
  if ("Z" %in% type) {
    out$Z <- stats::setNames(
      if (k_max >= 2L) {
        .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                              exponent = 1L, max_iter = max_iter, tol = tol,
                              x0 = x0, normalize = normalize, label = "Z")
      } else {
        rep(0, n)  # singletons only: no tensor entry, zero by convention
      },
      nodes
    )
  }

  # ---- HEC: lambda x^{k-1} = sum_{e contains i} prod_{j in e,j!=i} x_j ----
  if ("H" %in% type) {
    out$H <- stats::setNames(
      if (k_max >= 2L) {
        .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                              exponent = k_max - 1L,
                              max_iter = max_iter, tol = tol,
                              x0 = x0, normalize = normalize, label = "H")
      } else {
        rep(0, n)
      },
      nodes
    )
  }

  # Preserve user-requested order
  out[type]
}

# The tensor map f(x)_i = sum_{e contains i} prod_{j in e, j != i} x_j, with
# the (exponent)-th root taken for H. O(k) per hyperedge via total-product /
# x_j, with zeros handled explicitly.
#' @noRd
.hg_tensor_apply <- function(hyperedges, edge_sizes, x, n, exponent) {
  y <- numeric(n)
  # A loop over hyperedges: each adds into the entries of its own members.
  for (e_idx in seq_along(hyperedges)) {
    e  <- hyperedges[[e_idx]]
    ke <- edge_sizes[e_idx]
    if (ke < 2L) next
    x_e <- x[e]
    if (any(x_e == 0)) {
      nz_idx <- which(x_e != 0)
      if (length(nz_idx) == ke - 1L) {
        zero_pos <- setdiff(seq_len(ke), nz_idx)
        y[e[zero_pos]] <- y[e[zero_pos]] + prod(x_e[nz_idx])
      }
      # With >=2 zeros in x_e, every leave-one-out product is zero.
    } else {
      y[e] <- y[e] + prod(x_e) / x_e
    }
  }
  if (exponent > 1L) {
    # k-1-root, preserving sign (for non-negative hypergraphs all >= 0)
    y <- sign(y) * abs(y)^(1 / exponent)
  }
  y
}

# Unit-L2 direction of the tensor map, evaluated in log space. A common
# scale cancels in normalization (also after the H root), so neither tiny
# products nor their squared norm need be represented in floating point.
#' @noRd
.hg_tensor_direction <- function(hyperedges, edge_sizes, x, n, exponent) {
  contributions <- Map(function(e, size) {
    if (size < 2L) return(rep(-Inf, size))
    values <- x[e]
    positive <- values > 0
    if (all(positive)) {
      logs <- log(values)
      return(sum(logs) - logs)
    }
    out <- rep(-Inf, size)
    if (sum(positive) == size - 1L) out[!positive] <- sum(log(values[positive]))
    out
  }, hyperedges, edge_sizes)
  log_values <- unlist(contributions, use.names = FALSE)
  finite <- is.finite(log_values)
  if (!any(finite)) return(numeric(n))
  node_indices <- unlist(hyperedges, use.names = FALSE)[finite]
  groups <- split(log_values[finite], node_indices)
  log_totals <- vapply(groups, function(v) {
    largest <- max(v)
    largest + log(sum(exp(v - largest)))
  }, numeric(1L))
  y <- numeric(n)
  y[as.integer(names(log_totals))] <-
    exp((log_totals - max(log_totals)) / exponent)
  y / sqrt(sum(y^2))
}

# Shared tensor-power-iteration kernel (Kolda & Mayo 2011 SS-HOPM shift):
# x_{t+1} ~ f(x_t) / ||f(x_t)|| + shift * x_t. Its fixed points are exactly
# the eigenvectors f(x) = c x. The update is normalised before the shift is
# added: for large hyperedges f(x) of a unit vector is a product of many
# small numbers, and an unnormalised shift swamped it, so the iterate stopped
# moving at once and passed the step test without solving the equation (a
# 3-petal sunflower of 20-member hyperedges returned a hub/petal ratio of 1
# instead of sqrt(3)). Convergence therefore also requires the eigen-equation
# residual to be small; otherwise the result warns with class
# `nestimate_no_converge` (naming `label`).
#
# exponent = 1    => Z-eigenvector (no post-root)
# exponent = k-1  => H-eigenvector (k-1-root)
#' @noRd
.hg_tensor_power_iter <- function(hyperedges, edge_sizes, n, exponent,
                                   max_iter, tol, x0, normalize,
                                   shift = 1, label = "tensor") {
  x <- x0
  step_converged <- FALSE
  # A loop, not a vectorised map: every step needs the previous iterate.
  for (iter in seq_len(max_iter)) {
    fx <- .hg_tensor_direction(hyperedges, edge_sizes, x, n, exponent)
    nf <- sqrt(sum(fx^2))
    if (nf == 0) {
      # no tensor entry touches the support of x: nothing to iterate
      x <- fx
      step_converged <- TRUE
      break
    }
    y <- fx / nf + shift * x
    y <- y / sqrt(sum(y^2))
    if (sum(abs(y - x)) < tol) {
      x <- y
      step_converged <- TRUE
      break
    }
    x <- y
  }

  # Eigen-equation residual ||f(x) - (x'f(x)) x|| / ||f(x)|| at the result
  residual_ok <- TRUE
  if (any(x != 0)) {
    fx <- .hg_tensor_direction(hyperedges, edge_sizes, x, n, exponent)
    nf <- sqrt(sum(fx^2))
    if (nf > 0) {
      residual <- sqrt(sum((fx - sum(x * fx) * x)^2)) / nf
      residual_ok <- residual <= max(sqrt(tol), 1e-6)
    }
  }
  if (!step_converged || !residual_ok) {
    warning(warningCondition(
      sprintf(paste0("%s eigenvector centrality did not converge in %d ",
                     "iterations (tol = %g): the eigen-equation is not solved"),
              label, as.integer(max_iter), tol),
      class = "nestimate_no_converge", call = NULL
    ))
  }

  if (any(x != 0) && sum(x) < 0) x <- -x
  if (!normalize && any(x != 0)) x <- x / max(abs(x))
  x
}

# Connected components of a symmetric adjacency by min-label propagation;
# returns one integer label per node, numbered in order of first appearance.
#' @noRd
.hg_component_labels <- function(adjacency) {
  n <- nrow(adjacency)
  label <- seq_len(n)
  cells <- which(adjacency != 0, arr.ind = TRUE)
  if (nrow(cells) == 0L) return(label)
  row_of <- factor(cells[, 1L], levels = seq_len(n))
  repeat {
    neighbour <- vapply(
      split(label[cells[, 2L]], row_of),
      \(v) if (length(v)) min(v) else .Machine$integer.max, integer(1L)
    )
    nxt <- pmin(label, neighbour)
    if (identical(nxt, label)) break
    label <- nxt
  }
  match(label, unique(label))
}

# Symmetric eigendecomposition of each connected component of a clique
# adjacency: one element per component with `index` (its nodes), `values`
# (decreasing) and `vectors` (columns). By Perron-Frobenius the leading
# eigenvalue of a connected non-negative component is simple and its vector
# has one sign, so the leading vector is well defined component by component
# even where the whole matrix has a repeated leading eigenvalue (equal
# components) or a +/- pair (bipartite components).
#' @noRd
.hg_component_spectra <- function(adjacency) {
  labels <- .hg_component_labels(adjacency)
  lapply(split(seq_len(nrow(adjacency)), labels), \(index) {
    eig <- eigen(adjacency[index, index, drop = FALSE], symmetric = TRUE)
    list(index = index, values = eig$values, vectors = eig$vectors)
  })
}

# Unit-L2 leading eigenvector of a clique adjacency from its component
# spectra. Components whose leading eigenvalue equals the global one (to
# rounding) share the leading eigenspace; the result is the projection of
# the start vector `x0` onto that space -- the limit of power iteration from
# `x0` -- and every other node is 0. An all-zero adjacency returns zeros.
#' @noRd
.hg_perron_vector <- function(spectra, x0) {
  lead <- vapply(spectra, \(piece) piece$values[1L], numeric(1L))
  top_value <- max(lead)
  if (!(top_value > 0)) return(rep(0, length(x0)))
  near <- abs(lead - top_value) <= sqrt(.Machine$double.eps) * top_value
  x <- numeric(length(x0))
  parts <- lapply(spectra[near], \(piece) {
    v <- piece$vectors[, 1L]
    v <- if (sum(v) < 0) -v else v
    # a Perron vector is non-negative; clamp rounding noise below zero
    v <- pmax(v, 0)
    sum(v * x0[piece$index]) * v
  })
  x[unlist(lapply(spectra[near], `[[`, "index"))] <- unlist(parts)
  x / sqrt(sum(x^2))
}
