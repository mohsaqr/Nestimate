# ---- Diagram-level operations on persistent_homology objects ----
#
# bottleneck_distance() and persistence_landscape() take persistent_homology
# objects produced by persistent_homology() and operate on the persistence
# diagram only -- the boundary-matrix reduction lives in simplicial.R.

# =========================================================================
# Bottleneck distance
# =========================================================================

#' Bottleneck Distance Between Persistence Diagrams
#'
#' @description
#' `bottleneck_distance()` returns the bottleneck distance between two
#' persistence diagrams, one value per homology dimension. It takes two
#' [persistent_homology()] results, or two data frames of persistence
#' pairs. For the finite pairs, the distance is
#' \deqn{W_\infty(D_1, D_2) = \inf_{\gamma} \sup_{p \in D_1} \|p - \gamma(p)\|_\infty,}
#' where \eqn{\gamma} ranges over bijections \eqn{D_1 \cup \Delta \to
#' D_2 \cup \Delta} and \eqn{\Delta = \{(x,x)\}} is the diagonal. Each
#' point is matched to a point of the other diagram or to its projection
#' onto the diagonal, at cost \eqn{|d - b|/2}.
#'
#' The distance is the largest shift needed to match the two diagrams'
#' features: a small value means every feature of either diagram has a close
#' counterpart in the other or is short-lived enough to sit near the
#' diagonal. It is a distance between two observed diagrams and carries no
#' sampling uncertainty of its own; to compare diagrams on a common grid,
#' or to average several, use [persistence_landscape()].
#'
#' @details
#' The smallest feasible \eqn{\varepsilon} is found by binary search over
#' the candidate costs (every pairwise \eqn{L_\infty} distance and every
#' diagonal projection), checking each with a Kuhn bipartite matching.
#'
#' Essential classes (`death = Inf` in VR mode, `death = 0` with
#' `birth > 0` in clique mode, as [persistence_landscape()] reads them; a
#' `(0, 0)` pair is a zero-length finite point) are matched one-to-one within
#' each
#' dimension by their sorted births. When the two diagrams have different
#' numbers of essential classes in a dimension, the distance for that
#' dimension is `Inf`.
#'
#' @param d1,d2 The two diagrams: [persistent_homology()] results, or data
#'   frames with the columns `dimension`, `birth` and `death` (other
#'   columns are ignored). Any other input, or a data frame missing one of
#'   those columns, raises an error.
#' @param dimension The homology dimensions to compare, as whole numbers.
#'   Default `NULL`: every dimension that appears in either diagram. A
#'   dimension absent from both gives 0.
#' @param tol The tolerance added to each candidate cost when checking
#'   whether a matching exists, so that floating-point error does not miss
#'   an exact match. Default `.Machine$double.eps ^ 0.5`.
#'
#' @return A named numeric vector with one distance per dimension, named
#'   `"dim_<k>"` (empty when there is no dimension to compare). `Inf` marks
#'   a structural mismatch: different numbers of essential classes in that
#'   dimension. The distance of a diagram to itself is 0.
#'
#' @references
#' Edelsbrunner, H. & Harer, J. (2010). \emph{Computational Topology: An
#' Introduction}. AMS. Section VIII.
#'
#' @examples
#' mat1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
#' rownames(mat1) <- colnames(mat1) <- c("A","B","C")
#' ph1 <- persistent_homology(mat1, n_steps = 5)
#' bottleneck_distance(ph1, ph1)  # self-distance is 0
#'
#' @export
bottleneck_distance <- function(d1, d2, dimension = NULL,
                                 tol = .Machine$double.eps ^ 0.5) {
  df1 <- .ph_as_diagram(d1)
  df2 <- .ph_as_diagram(d2)
  dims <- if (is.null(dimension)) {
    sort(unique(c(df1$dimension, df2$dimension)))
  } else {
    .assert("`dimension` must be whole numbers" =
              is.numeric(dimension) && !anyNA(dimension) &&
                all(dimension == as.integer(dimension)))
    as.integer(dimension)
  }
  if (length(dims) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }
  out <- vapply(dims, function(k) {
    .bottleneck_one_dim(df1[df1$dimension == k, , drop = FALSE],
                       df2[df2$dimension == k, , drop = FALSE],
                       tol = tol)
  }, numeric(1))
  names(out) <- paste0("dim_", dims)
  out
}

#' @noRd
.ph_as_diagram <- function(x) {
  if (inherits(x, "persistent_homology")) return(x$persistence)
  if (is.data.frame(x)) {
    need <- c("dimension", "birth", "death")
    if (!all(need %in% names(x))) {
      stop(errorCondition(paste0("data.frame must have columns: ",
           paste(need, collapse = ", ")),
           class = "nestimate_bad_input", call = NULL))
    }
    return(x[, c("dimension", "birth", "death")])
  }
  stop(errorCondition("d1/d2 must be persistent_homology or data.frame.",
                      class = "nestimate_bad_input", call = NULL))
}

#' @noRd
.bottleneck_one_dim <- function(p1, p2, tol = .Machine$double.eps ^ 0.5) {
  # Classify rows as finite or essential. In clique mode essentials have
  # death == 0; in VR mode essentials have death == Inf. We treat both as
  # essential.
  ess1 <- .is_essential(p1)
  ess2 <- .is_essential(p2)

  fin1 <- as.matrix(p1[!ess1, c("birth", "death"), drop = FALSE])
  fin2 <- as.matrix(p2[!ess2, c("birth", "death"), drop = FALSE])
  ess_b1 <- p1$birth[ess1]
  ess_b2 <- p2$birth[ess2]

  # Essential half: must match one-to-one in birth values.
  if (length(ess_b1) != length(ess_b2)) {
    return(Inf)
  }
  ess_cost <- if (length(ess_b1) == 0L) 0 else {
    # Optimal bottleneck on a 1-D matching: pair sorted.
    max(abs(sort(ess_b1) - sort(ess_b2)))
  }

  if (nrow(fin1) == 0L && nrow(fin2) == 0L) {
    return(ess_cost)
  }

  # Candidate eps values: pairwise distances + diagonal projections.
  fin_cost <- .bottleneck_finite(fin1, fin2, tol = tol)
  max(ess_cost, fin_cost)
}

#' @noRd
.is_essential <- function(p) {
  # death 0 marks an essential class in the clique encoding only when the
  # class was born (birth > 0); a (0, 0) pair is a zero-length finite point,
  # as persistence_landscape() reads it
  is.infinite(p$death) | (p$death == 0 & p$birth > 0)
}

#' @noRd
.bottleneck_finite <- function(p1, p2, tol = .Machine$double.eps ^ 0.5) {
  n1 <- nrow(p1); n2 <- nrow(p2)
  if (n1 == 0L && n2 == 0L) return(0)
  if (n1 == 0L) return(max(abs(p2[, "death"] - p2[, "birth"]) / 2))
  if (n2 == 0L) return(max(abs(p1[, "death"] - p1[, "birth"]) / 2))

  # Pairwise L_inf distances p1 x p2
  d_p1_p2 <- pmax(
    abs(outer(p1[, "birth"], p2[, "birth"], "-")),
    abs(outer(p1[, "death"], p2[, "death"], "-"))
  )
  d_p1_diag <- abs(p1[, "death"] - p1[, "birth"]) / 2
  d_p2_diag <- abs(p2[, "death"] - p2[, "birth"]) / 2
  candidates <- sort(unique(c(0, as.vector(d_p1_p2), d_p1_diag, d_p2_diag)))

  if (.bottleneck_feasible(p1, p2, d_p1_p2, d_p1_diag, d_p2_diag,
                            candidates[length(candidates)] + tol)) {
    # find smallest eps via binary search (a loop: each step halves the
    # interval the previous one left)
    lo <- 1L; hi <- length(candidates)
    while (lo < hi) {
      mid <- (lo + hi) %/% 2L
      if (.bottleneck_feasible(p1, p2, d_p1_p2, d_p1_diag, d_p2_diag,
                                candidates[mid] + tol)) {
        hi <- mid
      } else {
        lo <- mid + 1L
      }
    }
    candidates[lo]
  } else {
    Inf
  }
}

#' @noRd
.bottleneck_feasible <- function(p1, p2, d_p1_p2, d_p1_diag, d_p2_diag, eps) {
  # Build bipartite adjacency on (n1 + n2) left, (n1 + n2) right.
  # Left: rows 1..n1 = p1 points, rows n1+1..n1+n2 = diagonal copies of p2
  # Right: cols 1..n2 = p2 points, cols n2+1..n2+n1 = diagonal copies of p1
  # Edges:
  #   p1[i] -> p2[j]            if d_p1_p2[i,j] <= eps
  #   p1[i] -> diag(p1[i])      if d_p1_diag[i] <= eps  (col n2 + i)
  #   diag(p2[j]) -> p2[j]      if d_p2_diag[j] <= eps  (row n1 + j -> col j)
  #   diag(p2[j]) -> diag(p1[i]) always (cost 0)        (row n1+j -> col n2+i)
  n1 <- nrow(p1); n2 <- nrow(p2)
  adj <- matrix(FALSE, n1 + n2, n1 + n2)
  if (n1 > 0L && n2 > 0L) adj[seq_len(n1), seq_len(n2)] <- d_p1_p2 <= eps
  if (n1 > 0L) {
    rr <- seq_len(n1); cc <- n2 + seq_len(n1)
    diag_block <- matrix(FALSE, n1, n1)
    diag(diag_block) <- d_p1_diag <= eps
    adj[rr, cc] <- diag_block
  }
  if (n2 > 0L) {
    rr <- n1 + seq_len(n2); cc <- seq_len(n2)
    diag_block <- matrix(FALSE, n2, n2)
    diag(diag_block) <- d_p2_diag <= eps
    adj[rr, cc] <- diag_block
  }
  if (n1 > 0L && n2 > 0L) {
    adj[(n1 + 1L):(n1 + n2), (n2 + 1L):(n2 + n1)] <- TRUE
  }
  .kuhn_match_size(adj) == nrow(adj)
}

#' @noRd
.kuhn_match_size <- function(adj) {
  # Kuhn bipartite cardinality matching via DFS augmenting paths.
  # Sequential by construction (augmenting paths mutate match state) --
  # same for-loop exception class as the boundary-matrix reduction.
  n_l <- nrow(adj); n_r <- ncol(adj)
  if (n_l == 0L || n_r == 0L) return(0L)
  match_r <- integer(n_r)
  count <- 0L

  try_augment <- function(u, visited_env) {
    candidates <- which(adj[u, ])
    for (v in candidates) {
      if (!visited_env$v[v]) {
        visited_env$v[v] <- TRUE
        if (match_r[v] == 0L || try_augment(match_r[v], visited_env)) {
          match_r[v] <<- u
          return(TRUE)
        }
      }
    }
    FALSE
  }

  for (u in seq_len(n_l)) {
    visited_env <- new.env(parent = emptyenv())
    visited_env$v <- logical(n_r)
    if (try_augment(u, visited_env)) count <- count + 1L
  }
  count
}

# =========================================================================
# Persistence landscapes (Bubenik 2015)
# =========================================================================

#' Persistence Landscape
#'
#' @description
#' `persistence_landscape()` computes the persistence landscape (Bubenik
#' 2015) of one homology dimension of a persistence diagram, evaluated on a
#' grid. Each (birth, death) pair contributes a tent function
#' \deqn{\Lambda_{(b,d)}(t) = \max(0, \min(t - b, d - t)),}
#' and the \eqn{k}-th landscape function \eqn{\lambda^{(k)}(t)} is the
#' \eqn{k}-th largest of \eqn{\{\Lambda_{(b_i,d_i)}(t)\}_i} at each
#' \eqn{t}. It takes a [persistent_homology()] result or a data frame of
#' persistence pairs.
#'
#' At each filtration value \eqn{t}, a tent is the distance from \eqn{t}
#' to the nearer end of a feature's interval (0 outside it), so the first
#' landscape function follows the feature that \eqn{t} lies deepest
#' inside, the second the feature with the second largest tent, and so
#' on. This need not be the most persistent feature alive: a long-lived
#' feature near its birth or death can sit below a shorter one near its
#' midpoint. Landscapes
#' are stable: the sup-norm distance between two landscapes is at most the
#' [bottleneck_distance()] between their diagrams. They also live in a
#' Banach space, where they can be averaged (Bubenik 2015).
#'
#' Essential classes are left out, because the tent of an infinitely lived
#' class is undefined on a finite grid: pairs with `death = Inf` (VR mode)
#' and pairs with `death = 0` but `birth > 0` (the clique-mode encoding of
#' an essential class) are dropped before the landscape is built. When no
#' finite pair remains in the requested dimension, every landscape
#' function is zero on the grid.
#'
#' @param ph The diagram: a [persistent_homology()] result, or a data frame
#'   with the columns `dimension`, `birth` and `death`.
#' @param k_max The number of landscape functions to compute, a positive
#'   whole number. Functions beyond the number of pairs are zero. Default
#'   `5`.
#' @param dimension The homology dimension, a single non-negative whole
#'   number. Default `1`.
#' @param t_grid The filtration values to evaluate the landscape at.
#'   Default `NULL`: 200 evenly spaced points from the smallest birth to the
#'   largest death of the remaining pairs (each pair is ordered so that
#'   birth comes first, since clique filtrations descend), or from 0 to 1
#'   when no pair remains. Pass the same grid to landscapes you want to
#'   compare.
#'
#' @return A `persistence_landscape` object. Printing it shows the
#'   dimension, `k_max`, the grid range and the sup norm (largest value) of
#'   each landscape function; `plot()` draws the functions as lines over
#'   the grid. The object holds:
#' \describe{
#'   \item{landscape}{A data frame with one row per landscape function and
#'     grid point, and the columns `k`, `t` and `value`.}
#'   \item{dimension}{The homology dimension.}
#'   \item{k_max}{The number of landscape functions.}
#'   \item{t_grid}{The evaluation grid.}
#' }
#'
#' @references
#' Bubenik, P. (2015). Statistical topological data analysis using
#' persistence landscapes. \emph{Journal of Machine Learning Research}
#' \strong{16}, 77-102.
#'
#' @examples
#' mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
#' rownames(mat) <- colnames(mat) <- c("A","B","C")
#' ph <- persistent_homology(mat, n_steps = 5)
#' pl <- persistence_landscape(ph, k_max = 3, dimension = 0)
#'
#' @export
persistence_landscape <- function(ph, k_max = 5L, dimension = 1L,
                                   t_grid = NULL) {
  .assert(
    "`k_max` must be a single positive whole number" =
      is.numeric(k_max) && length(k_max) == 1L &&
        !is.na(k_max) && k_max >= 1L && k_max == as.integer(k_max),
    "`dimension` must be a single non-negative whole number" =
      is.numeric(dimension) && length(dimension) == 1L &&
        !is.na(dimension) && dimension >= 0L &&
        dimension == as.integer(dimension)
  )
  df <- .ph_as_diagram(ph)
  sub <- df[df$dimension == as.integer(dimension), , drop = FALSE]
  # Drop essential pairs (death = Inf) for landscape computation; tent for an
  # infinite-lived class is undefined on a finite grid.
  sub <- sub[is.finite(sub$death), , drop = FALSE]
  # In clique mode, essential pairs have death = 0; treat as Inf-equivalent
  # only if birth > 0 (otherwise it's a genuine point at the origin).
  sub <- sub[!(sub$death == 0 & sub$birth > 0), , drop = FALSE]

  if (nrow(sub) == 0L) {
    if (is.null(t_grid)) t_grid <- seq(0, 1, length.out = 200L)
    out <- data.frame(
      k = rep(seq_len(k_max), each = length(t_grid)),
      t = rep(t_grid, k_max),
      value = 0, stringsAsFactors = FALSE
    )
    return(structure(list(
      landscape = out, dimension = as.integer(dimension),
      k_max = as.integer(k_max), t_grid = t_grid
    ), class = "persistence_landscape"))
  }

  bd <- as.matrix(sub[, c("birth", "death")])
  # Normalize orientation: birth <= death (some conventions store the other way)
  bd <- t(apply(bd, 1L, function(row) c(min(row), max(row))))

  if (is.null(t_grid)) {
    t_grid <- seq(min(bd[, 1L]), max(bd[, 2L]), length.out = 200L)
  }
  n_pairs <- nrow(bd); n_t <- length(t_grid)

  # tents[i, j] = tent value of pair i at t_grid[j]. Build via outer to keep
  # shape stable when n_pairs == 1 (vapply collapses to a vector otherwise).
  left  <- outer(bd[, 1L], t_grid, function(b, t) t - b)
  right <- outer(bd[, 2L], t_grid, function(d, t) d - t)
  tents <- pmax(0, pmin(left, right))
  dim(tents) <- c(n_pairs, n_t)

  # For each t-column, take the k_max largest values across pairs (rows).
  lambda <- vapply(seq_len(n_t), function(j) {
    v <- sort(tents[, j], decreasing = TRUE)
    k_eff <- min(k_max, length(v))
    c(v[seq_len(k_eff)], rep(0, k_max - k_eff))
  }, numeric(k_max))
  # lambda is k_max x n_t -- flatten ROW-MAJOR so values align with the
  # (k, t) iteration order below.
  dim(lambda) <- c(k_max, n_t)

  out <- data.frame(
    k = rep(seq_len(k_max), each = n_t),
    t = rep(t_grid, k_max),
    value = as.vector(t(lambda)),
    stringsAsFactors = FALSE
  )

  structure(list(
    landscape = out, dimension = as.integer(dimension),
    k_max = as.integer(k_max), t_grid = t_grid
  ), class = "persistence_landscape")
}

#' @rdname persistence_landscape
#' @param x For the `print()` and `plot()` methods: a
#'   `persistence_landscape` object.
#' @param ... Not used.
#' @return The `print()` method returns its input invisibly, and `plot()` a
#'   `ggplot` object with one line per landscape function.
#' @export
print.persistence_landscape <- function(x, ...) {
  cat(sprintf("Persistence Landscape (dimension %d, k_max = %d)\n",
              x$dimension, x$k_max))
  cat(sprintf("  Evaluated at %d grid points in [%.4f, %.4f]\n",
              length(x$t_grid), min(x$t_grid), max(x$t_grid)))
  k_norms <- vapply(seq_len(x$k_max), function(k) {
    sub <- x$landscape[x$landscape$k == k, ]
    max(sub$value)
  }, numeric(1))
  cat("  Sup norms by k:",
      paste(sprintf("k%d=%.4f", seq_len(x$k_max), k_norms),
            collapse = "  "), "\n")
  invisible(x)
}

#' @rdname persistence_landscape
#' @export
plot.persistence_landscape <- function(x, ...) {
  df <- x$landscape
  df$k <- factor(df$k, levels = seq_len(x$k_max))
  # Render legend labels via plotmath (lambda[k]) rather than a literal
  # Unicode "lambda" string: the latter is drawn by the graphics device font
  # and fails to transcode in non-UTF-8 check locales ("conversion failure
  # ... in 'mbcsToSbcs'"). plotmath uses R's symbol-font tables instead.
  lambda_labels <- parse(text = paste0("lambda[", seq_len(x$k_max), "]"))
  ggplot2::ggplot(df, ggplot2::aes(x = t, y = value, color = k)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_discrete(labels = lambda_labels) +
    ggplot2::labs(
      title = sprintf("Persistence Landscape (dimension %d)", x$dimension),
      subtitle = sprintf("Bubenik 2015; top %d landscape functions", x$k_max),
      x = "Filtration parameter t", y = expression(lambda^(k)),
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}
