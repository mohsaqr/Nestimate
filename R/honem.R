# ---- HONEM: Higher-Order Network Embedding ----
#
# A HONEM-style embedding after Saebi et al. (2020): exponentially decaying
# matrix powers of the HON transition matrix, then truncated SVD. It is not
# the paper's construction, whose order-v matrices D^v over first-order nodes
# are not fully specified there (see the build_honem() documentation).
#
# NO random walks or skip-gram - direct matrix factorization.

# ---------------------------------------------------------------------------
# Internal: Build transition matrix from HON adjacency
# ---------------------------------------------------------------------------

#' Build row-normalized transition matrix from HON output
#'
#' @param mat Square weighted adjacency matrix from build_hon().
#' @return Row-stochastic transition matrix.
#' @noRd
.honem_transition_matrix <- function(mat) {
  row_sums <- rowSums(mat)
  nonzero <- row_sums > 0
  result <- mat
  result[nonzero, ] <- result[nonzero, ] / row_sums[nonzero]
  result
}

# ---------------------------------------------------------------------------
# Internal: Compute neighborhood matrix
# ---------------------------------------------------------------------------

#' Compute HONEM neighborhood matrix
#'
#' S = (1/Z) * sum_{k=0}^{L} exp(-k) * D^{k+1}
#' where D is the row-normalized transition matrix and Z is the normalization.
#'
#' @param D Row-stochastic transition matrix.
#' @param max_power Integer. Maximum power L (default: 10 or diameter).
#' @return Dense neighborhood matrix S.
#' @noRd
.honem_neighborhood_matrix <- function(D, max_power = 10L) {
  n <- nrow(D)

  # Normalization constant: Z = sum_{k=0}^{L} exp(-k)
  weights <- exp(-(0L:max_power))
  Z <- sum(weights)

  # Accumulate weighted matrix powers
  S <- matrix(0, nrow = n, ncol = n)
  D_power <- D  # D^1

  vapply(seq_along(weights), function(k_idx) {
    S <<- S + weights[k_idx] * D_power
    D_power <<- D_power %*% D  # D^{k+2}
    0
  }, numeric(1L))

  S / Z
}

# ---------------------------------------------------------------------------
# Internal: Truncated SVD for embeddings
# ---------------------------------------------------------------------------

#' Compute embeddings via truncated SVD
#'
#' @param S Neighborhood matrix.
#' @param dim Integer. Embedding dimension.
#' @return List with embeddings (n x dim matrix), singular_values, explained.
#' @noRd
.honem_svd <- function(S, dim) {
  n <- nrow(S)
  dim <- min(dim, n - 1L)

  sv <- svd(S, nu = dim, nv = dim)

  # Embedding: U * sqrt(Sigma)
  sigma_sqrt <- sqrt(sv$d[seq_len(dim)])
  embeddings <- sweep(sv$u[, seq_len(dim), drop = FALSE], 2, sigma_sqrt, `*`)
  rownames(embeddings) <- rownames(S)
  colnames(embeddings) <- paste0("dim_", seq_len(dim))

  total_var <- sum(sv$d^2)
  explained <- sum(sv$d[seq_len(dim)]^2) / total_var

  list(
    embeddings = embeddings,
    singular_values = sv$d[seq_len(dim)],
    explained_variance = explained
  )
}

# ---------------------------------------------------------------------------
# Main function: build_honem
# ---------------------------------------------------------------------------

#' Build HONEM Embeddings for Higher-Order Networks
#'
#' @description
#' `build_honem()` places the nodes of a higher-order network (HON) in a
#' low-dimensional space, after the HONEM approach of Saebi et al. (2020):
#' nodes whose multi-step transition profiles are alike (where they lead in
#' one, two or more steps) get nearby coordinates. It reads a `net_hon` from
#' [build_hon()], whose higher-order nodes (such as `"A -> B"`) are embedded
#' beside the first-order ones, or any square weighted adjacency matrix. The
#' coordinates come from a decayed multi-step proximity matrix (see
#' *Details*).
#'
#' The coordinates have no units and the sign of each dimension is
#' arbitrary; what carries meaning is how close nodes lie to each other.
#' The variance explained says how much of the proximity matrix the `dim`
#' kept dimensions reproduce; a plot of two of them shows less than that.
#'
#' @details
#' No random walks or skip-gram are used: the coordinates come from the
#' truncated singular value decomposition of the proximity matrix. With
#' \eqn{D} the row-normalised transition matrix over all HON nodes (rows
#' without outgoing edges stay zero), the proximity matrix is
#' \eqn{S = Z^{-1} \sum_{k=0}^{L} e^{-k} D^{k+1}} with
#' \eqn{Z = \sum_{k=0}^{L} e^{-k}} and \eqn{L} = `max_power`, where
#' \eqn{D^{k+1}} is a matrix power. Each node's embedding is its row of
#' \eqn{U_d \sqrt{\Sigma_d}}, from the first \eqn{d} = `dim` singular
#' vectors and values of \eqn{S}. `explained_variance` is the share of the
#' sum of squared singular values held by the first \eqn{d}; \eqn{S} is not
#' centred, so this is not a principal-component variance.
#'
#' \strong{How this differs from the paper.} Saebi et al. write the same
#' sum, but there \eqn{D^v} is a \emph{v-th-order neighbourhood matrix
#' over first-order nodes} (the superscript is an order, not a power),
#' filled from the HON edge weights of order-v dependencies, with \eqn{L}
#' the highest order present, rule extraction changed to keep every lower
#' order, and one embedding per first-order node. The paper gives no
#' formula for \eqn{D^v}, and its worked example (Supplementary Table S.1)
#' is internally inconsistent; no reference code is public. This function
#' therefore does not claim to reproduce the paper's embeddings.
#'
#' @param hon The network to embed: a `net_hon` object from [build_hon()],
#'   or a square weighted adjacency matrix with at least 2 rows. Anything
#'   else raises an error.
#' @param dim The number of embedding dimensions. Default `32`. It is capped
#'   at the number of nodes minus 1 without a warning; `print()` reports the
#'   dimension used.
#' @param max_power \eqn{L}, the highest matrix power in the decayed sum, at
#'   least 1. Default `10`. Higher values let longer paths shape the
#'   proximities.
#' @return A `net_honem` object. `print()` shows the number of nodes, the
#'   dimensions, `max_power` and the variance explained; `summary()` prints
#'   the leading singular values and the range of the coordinates, and
#'   returns one row per node; `plot()` draws the nodes on two dimensions.
#'   The object's fields:
#'   \describe{
#'     \item{embeddings}{Numeric matrix (nodes x `dim`) of coordinates, with
#'       the node names as row names and columns `dim_1`, `dim_2`, ...}
#'     \item{nodes}{Character vector of node names.}
#'     \item{singular_values}{The `dim` leading singular values.}
#'     \item{explained_variance}{Share of the squared singular values held
#'       by the kept dimensions.}
#'     \item{dim}{Embedding dimension used.}
#'     \item{max_power}{Highest matrix power used.}
#'     \item{n_nodes}{Number of nodes embedded.}
#'   }
#'
#' @references
#' Saebi, M., Ciampaglia, G. L., Kaplan, L. M., & Chawla, N. V. (2020).
#' HONEM: Learning Embedding for Higher Order Networks. \emph{Big Data},
#' 8(4), 255-269.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hem <- build_honem(build_hon(seqs, max_order = 2), dim = 2)
#'
#' \donttest{
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
#'               c("B","C","D","A"), c("C","D","A","B"))
#' hon <- build_hon(trajs, max_order = 2)
#' emb <- build_honem(hon, dim = 4)
#' print(emb)
#' plot(emb)
#' }
#'
#' @export
build_honem <- function(hon, dim = 32L, max_power = 10L) {
  dim <- as.integer(dim)
  max_power <- as.integer(max_power)
  .assert(
    "'dim' must be >= 1" = dim >= 1L,
    "'max_power' must be >= 1" = max_power >= 1L
  )

  # Extract adjacency matrix
  if (inherits(hon, "net_hon")) {
    mat <- hon$matrix
  } else if (is.matrix(hon) && nrow(hon) == ncol(hon)) {
    mat <- hon
  } else {
    stop(errorCondition("'hon' must be a net_hon object or a square matrix",
                        class = "nestimate_bad_input", call = NULL))
  }

  n <- nrow(mat)
  if (n < 2L) stop(errorCondition("Need at least 2 nodes for embedding",
                                  class = "nestimate_bad_input", call = NULL))

  dim <- min(dim, n - 1L)

  # Build transition matrix
  D <- .honem_transition_matrix(mat)

  # Compute neighborhood matrix
  S <- .honem_neighborhood_matrix(D, max_power)
  rownames(S) <- rownames(mat)
  colnames(S) <- colnames(mat)

  # Truncated SVD
  svd_result <- .honem_svd(S, dim)

  result <- list(
    embeddings = svd_result$embeddings,
    nodes = rownames(mat),
    singular_values = svd_result$singular_values,
    explained_variance = svd_result$explained_variance,
    dim = dim,
    max_power = max_power,
    n_nodes = n
  )

  class(result) <- "net_honem"
  result
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname build_honem
#' @param x For the `print()` and `plot()` methods: a `net_honem` object.
#' @param ... For the `plot()` method: further arguments passed to
#'   [graphics::plot()]. For the `print()` and `summary()` methods: ignored.
#' @param object For the `summary()` method: a `net_honem` object.
#' @param dims For the `plot()` method: the two dimensions to draw. Default
#'   `c(1, 2)`; a dimension the embedding does not have raises an error.
#'   Nodes are labelled when there are at most 50; an embedding with a single
#'   dimension draws nothing and says so in a message.
#' @return The `print()` and `plot()` methods return the input invisibly.
#'   `summary()` returns a data frame with one row per node: `node` (the node
#'   name), then `dim1`, `dim2`, ... with the coordinates; the summary text is
#'   printed as a side effect.
#' @export
print.net_honem <- function(x, ...) {
  cat("HONEM: Higher-Order Network Embedding\n")
  cat(sprintf("  Nodes:      %d\n", x$n_nodes))
  cat(sprintf("  Dimensions: %d\n", x$dim))
  cat(sprintf("  Max power:  %d\n", x$max_power))
  cat(sprintf("  Variance explained: %.1f%%\n", x$explained_variance * 100))
  invisible(x)
}

#' @rdname build_honem
#' @export
summary.net_honem <- function(object, ...) {
  cat("HONEM Summary\n\n")
  cat(sprintf("  Nodes: %d | Dimensions: %d\n", object$n_nodes, object$dim))
  cat(sprintf("  Variance explained: %.1f%%\n",
              object$explained_variance * 100))
  cat(sprintf("  Top singular values: %s\n",
              paste(round(object$singular_values[seq_len(min(5, object$dim))],
                          3), collapse = ", ")))
  cat(sprintf("  Embedding range: [%.3f, %.3f]\n",
              min(object$embeddings), max(object$embeddings)))

  emb <- object$embeddings
  node_labels <- rownames(emb)
  if (is.null(node_labels)) node_labels <- paste0("n", seq_len(nrow(emb)))
  dim_cols <- as.data.frame(emb)
  colnames(dim_cols) <- paste0("dim", seq_len(ncol(emb)))
  data.frame(node = node_labels, dim_cols,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @rdname build_honem
#' @export
plot.net_honem <- function(x, dims = c(1L, 2L), ...) {
  if (x$dim < 2L) {
    message("Need at least 2 dimensions to plot")
    return(invisible(x))
  }

  emb <- x$embeddings
  if (!is.numeric(dims) || length(dims) != 2L || anyNA(dims) ||
      any(dims != round(dims)) || any(dims < 1L) || any(dims > ncol(emb))) {
    stop(errorCondition(
      sprintf(paste0("`dims` must be two dimension numbers between 1 and %d ",
                     "(the embedding has %d dimensions)."),
              ncol(emb), ncol(emb)),
      class = "nestimate_bad_input", call = NULL))
  }
  d1 <- dims[1L]
  d2 <- dims[2L]

  plot(emb[, d1], emb[, d2], pch = 19, col = "steelblue",
       xlab = sprintf("Dimension %d", d1),
       ylab = sprintf("Dimension %d", d2),
       main = sprintf("HONEM Embedding (%d nodes, %.0f%% var)",
                       x$n_nodes, x$explained_variance * 100),
       ...)

  if (x$n_nodes <= 50) {
    graphics::text(emb[, d1], emb[, d2], labels = x$nodes,
                   pos = 3, cex = 0.7)
  }

  invisible(x)
}
