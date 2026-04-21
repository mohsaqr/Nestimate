# ---- Hypergraph eigenvector centralities (HON-5) -------------------------
# Three variants from Benson (2019, arXiv:1807.09644):
#   * clique-motif ("CEC") — standard eigenvector centrality on the
#                            clique-expanded pairwise graph
#   * Z-eigenvector ("Z")  — linear tensor eigenvector
#   * H-eigenvector ("H")  — H-eigenvector (power-k-1 recurrence)
#
# All three solved by power iteration on the hyperedge list.

#' Hypergraph eigenvector centralities
#'
#' Computes one or more eigenvector-style centralities on a
#' [net_hypergraph][build_hypergraph]: *clique-motif* (CEC),
#' *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC). Each variant
#' captures influence differently — CEC flattens group structure via
#' clique expansion, while ZEC and HEC propagate through the
#' higher-order groups directly.
#'
#' @param hg A `net_hypergraph` (from [build_hypergraph()] or
#'   [bipartite_groups()]).
#' @param type Character vector, any subset of `c("clique", "Z", "H")`.
#'   Default computes all three.
#' @param max_iter Maximum number of power-iteration steps. Default
#'   `1000`.
#' @param tol Convergence tolerance on the L1 change between successive
#'   iterates. Default `1e-8`.
#' @param normalize Logical. If `TRUE` (default), each returned
#'   centrality vector is L2-normalized to unit norm (compatible with
#'   `igraph::eigen_centrality()`'s scale for type `"clique"`).
#'
#' @return A named list; one component per requested `type`. Each
#'   component is a named numeric vector of length `hg$n_nodes`.
#'
#' @details
#' **Clique-motif eigenvector centrality (CEC)**: forms the
#' clique-expanded pairwise graph \eqn{W} where
#' \eqn{W_{ij} = |\{e : i, j \in e\}|} and returns the leading
#' eigenvector of \eqn{W}. Equivalent to running
#' `igraph::eigen_centrality()` on [clique_expansion()] output.
#'
#' **Z-eigenvector centrality (ZEC)**: solves the linear
#' eigen-equation on the hyperedge tensor,
#' \deqn{\lambda\, x_i \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j,}
#' via power iteration. Works for hypergraphs with mixed edge sizes.
#'
#' **H-eigenvector centrality (HEC)**: solves the power-k-1
#' eigen-equation,
#' \deqn{\lambda\, x_i^{k-1} \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j.}
#' For uniform hypergraphs (all hyperedges of size \eqn{k}), this is
#' equivalent to normalizing the ZEC update by the geometric-mean
#' exponent \eqn{1/(k-1)}. For mixed sizes, the effective exponent is
#' taken from the largest hyperedge; expect slightly different rankings
#' from ZEC in the mixed case.
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
#' cent <- hypergraph_centrality(hg)
#' # Compare rankings across the three variants
#' do.call(cbind, cent)
#'
#' @references
#' Benson, A. R. (2019). Three hypergraph eigenvector centralities.
#' \emph{SIAM Journal on Mathematics of Data Science} 1(2), 293-312.
#' arXiv:1807.09644.
#'
#' @note The `"clique"` (CEC) variant is validated against
#'   `igraph::eigen_centrality` (cosine ~ 1). The `"Z"` and `"H"` variants are
#'   **(experimental)** — validated only against a clean-room list-based
#'   tensor power iteration (same operator, different loop structure); no
#'   R package exposes tensor eigenvectors as a primitive for independent
#'   comparison.
#'
#' @export
hypergraph_centrality <- function(hg,
                                   type     = c("clique", "Z", "H"),
                                   max_iter = 1000L,
                                   tol      = 1e-8,
                                   normalize = TRUE) {
  stopifnot(
    inherits(hg, "net_hypergraph"),
    is.numeric(max_iter), length(max_iter) == 1L, max_iter > 0,
    is.numeric(tol), length(tol) == 1L, tol > 0,
    is.logical(normalize), length(normalize) == 1L
  )
  type <- match.arg(type, several.ok = TRUE)

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

  out <- list()

  # ---- CEC: power iteration on clique-expansion W ----
  if ("clique" %in% type) {
    B_bin <- (hg$incidence > 0) * 1.0
    W <- tcrossprod(B_bin)
    diag(W) <- 0
    x <- x0
    for (iter in seq_len(max_iter)) {
      y <- as.numeric(W %*% x)
      nrm <- sqrt(sum(y^2))
      if (nrm == 0) break
      y <- y / nrm
      if (sum(abs(y - x)) < tol) break
      x <- y
    }
    # Sign convention: positive entries (W is non-negative so Perron vec is positive)
    if (any(x != 0) && sum(x) < 0) x <- -x
    if (!normalize && any(x != 0)) x <- x / max(abs(x))
    out$clique <- stats::setNames(x, nodes)
  }

  # ---- ZEC: λ x = Σ_{e∋i} Π_{j∈e,j≠i} x_j ----
  if ("Z" %in% type) {
    out$Z <- stats::setNames(
      .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                            exponent = 1L, max_iter = max_iter, tol = tol,
                            x0 = x0, normalize = normalize),
      nodes
    )
  }

  # ---- HEC: λ x^{k-1} = Σ_{e∋i} Π_{j∈e,j≠i} x_j ----
  if ("H" %in% type) {
    out$H <- stats::setNames(
      .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                            exponent = k_max - 1L,
                            max_iter = max_iter, tol = tol,
                            x0 = x0, normalize = normalize),
      nodes
    )
  }

  # Preserve user-requested order
  out[type]
}

# Shared tensor-power-iteration kernel.
# Uses Kolda-Mayo SSHOPM shift: x_{k+1} ~ f(x_k) + shift * x_k
# which guarantees monotone convergence for non-negative tensors
# (Chang, Pearson & Zhang 2009 / Kolda & Mayo 2011).
#
# exponent = 1    ⇒ Z-eigenvector (no post-root)
# exponent = k-1  ⇒ H-eigenvector (k-1-root)
#' @noRd
.hg_tensor_power_iter <- function(hyperedges, edge_sizes, n, exponent,
                                   max_iter, tol, x0, normalize,
                                   shift = 1) {
  x <- x0
  for (iter in seq_len(max_iter)) {
    y <- numeric(n)
    for (e_idx in seq_along(hyperedges)) {
      e  <- hyperedges[[e_idx]]
      ke <- edge_sizes[e_idx]
      if (ke < 2L) next
      x_e <- x[e]
      # Product over all members excluding each i, done in O(k) per edge
      # via total-product / x_j (handling zeros explicitly)
      total <- prod(x_e)
      if (any(x_e == 0)) {
        nz_idx <- which(x_e != 0)
        if (length(nz_idx) == ke - 1L) {
          zero_pos <- setdiff(seq_len(ke), nz_idx)
          y[e[zero_pos]] <- y[e[zero_pos]] + prod(x_e[nz_idx])
        }
        # With >=2 zeros in x_e, every leave-one-out product is zero.
      } else {
        for (idx in seq_len(ke)) {
          y[e[idx]] <- y[e[idx]] + total / x_e[idx]
        }
      }
    }

    if (exponent > 1L) {
      # k-1-root, preserving sign (for non-negative hypergraphs all >= 0)
      y <- sign(y) * abs(y)^(1 / exponent)
    }

    # SSHOPM shift: bias towards x to stabilize oscillating tensor iterates
    y <- y + shift * x

    nrm <- sqrt(sum(y^2))
    if (nrm == 0) {
      x <- y
      break
    }
    y <- y / nrm

    if (sum(abs(y - x)) < tol) {
      x <- y
      break
    }
    x <- y
  }

  if (any(x != 0) && sum(x) < 0) x <- -x
  if (!normalize && any(x != 0)) x <- x / max(abs(x))
  x
}
