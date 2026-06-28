# ---- Edge pruning workflow (prune / deprune / reprune / details) ----
#
# Mirrors tna::prune()/deprune()/reprune()/pruning_details() in pure base R
# (no new dependencies). Pruning is stored, non-destructively, as a "pruning"
# attribute on the netobject: the original weights, the pruned weights, an
# `active` flag, the method/parameters used, and a tidy table of removed
# edges. `net_deprune()`/`net_reprune()` toggle `active` and swap `$weights`
# back and forth without recomputing anything.

#' Weak connectivity of a weighted adjacency matrix
#'
#' TRUE if every node is reachable from node 1 ignoring edge direction.
#' Uses a boolean transitive closure by repeated squaring (no igraph).
#' @param W Square numeric weight matrix (zeros = no edge).
#' @return Logical scalar.
#' @noRd
.is_weakly_connected <- function(W) {
  n <- nrow(W)
  if (n <= 1L) return(TRUE)
  reach <- (W > 0) | (t(W) > 0)
  diag(reach) <- TRUE
  steps <- ceiling(log2(n)) + 1L
  closure <- Reduce(function(R, .) (R %*% R) > 0, seq_len(steps), reach)
  all(closure[1L, ])
}

#' Greedy threshold / lowest-quantile pruning with a connectivity guard
#'
#' Candidate edges (weight <= cut-off) are visited in column-major order and
#' removed one at a time only when removal keeps the network weakly connected,
#' exactly matching \code{tna:::prune_default()}.
#' @noRd
.prune_default <- function(W, method, threshold, lowest, labels) {
  pos      <- W > 0
  n_edges  <- sum(pos)
  pos_vals <- W[pos]                                   # diagonal counts equally
  cut_off  <- switch(method,
                     threshold = threshold,
                     lowest    = stats::quantile(pos_vals, probs = lowest))
  # Self-loops (diagonal) are observed data: they inform the cut-off but are
  # never removed, so only off-diagonal edges are removal candidates.
  off_diag <- row(W) != col(W)
  cand <- which(pos & W <= cut_off & off_diag, arr.ind = TRUE)  # column-major

  state <- Reduce(function(s, j) {
    row <- cand[j, 1L]
    col <- cand[j, 2L]
    tmp <- s$W
    tmp[row, col] <- 0
    if (.is_weakly_connected(tmp)) {
      s$W <- tmp
      s$removed[j] <- TRUE
    }
    s
  }, seq_len(nrow(cand)), list(W = W, removed = logical(nrow(cand))))

  removed <- data.frame(
    from   = labels[cand[state$removed, 1L]],
    to     = labels[cand[state$removed, 2L]],
    weight = W[cand[state$removed, , drop = FALSE]],
    stringsAsFactors = FALSE
  )
  n_removed <- sum(state$removed)
  list(weights = state$W, method = method, threshold = threshold,
       lowest = lowest, cut_off = cut_off, removed = removed,
       num_removed = n_removed, num_retained = n_edges - n_removed)
}

#' Serrano disparity filter (backbone significance) -- matches tna
#' @noRd
.disparity_filter <- function(W, level) {
  d          <- ncol(W)
  idx_mat    <- 1 * (W > 0)
  out_edges  <- W / .rowSums(W, m = d, n = d)
  out_degree <- .rowSums(idx_mat, m = d, n = d)
  out_p      <- (1 - out_edges)^(out_degree - 1)
  in_edges   <- t(W) / .colSums(W, m = d, n = d)
  in_degree  <- .colSums(idx_mat, m = d, n = d)
  in_p       <- t((1 - in_edges)^(in_degree - 1))
  p_values   <- pmin(out_p, in_p)
  sig        <- 1 * (p_values < level)
  diag(sig)  <- 0
  sig
}

#' @noRd
.prune_disparity <- function(W, level, labels) {
  n_edges <- sum(W > 0)
  pruned  <- .disparity_filter(W, level) * W
  dimnames(pruned) <- dimnames(W)
  diag(pruned) <- diag(W)                 # self-loops are data: never removed
  rem_idx <- which(pruned == 0 & W != 0, arr.ind = TRUE)
  removed <- data.frame(
    from   = labels[rem_idx[, 1L]],
    to     = labels[rem_idx[, 2L]],
    weight = W[rem_idx],
    stringsAsFactors = FALSE
  )
  list(weights = pruned, method = "disparity", level = level,
       removed = removed, num_removed = nrow(rem_idx),
       num_retained = n_edges - nrow(rem_idx))
}

#' @noRd
.prune_bootstrap <- function(x, boot, labels, dots) {
  if (is.null(boot)) {
    boot <- do.call(bootstrap_network, c(list(x), dots))
  }
  if (!inherits(boot, "net_bootstrap")) {
    stop("'boot' must be a 'net_bootstrap' from bootstrap_network().",
         call. = FALSE)
  }
  W       <- x$weights
  n_edges <- sum(W > 0)
  pruned  <- boot$significant            # significant edges keep weight, rest 0
  dimnames(pruned) <- dimnames(W)
  diag(pruned) <- diag(W)                # self-loops are data: never removed
  rem_idx <- which(W > 0 & pruned == 0, arr.ind = TRUE)
  removed <- data.frame(
    from   = labels[rem_idx[, 1L]],
    to     = labels[rem_idx[, 2L]],
    weight = W[rem_idx],
    stringsAsFactors = FALSE
  )
  list(weights = pruned, method = "bootstrap", removed = removed,
       num_removed = nrow(rem_idx), num_retained = n_edges - nrow(rem_idx))
}

#' Re-sync a netobject's edge table to its (possibly pruned) weights
#' @noRd
.sync_netobject_edges <- function(x) {
  x$edges   <- .extract_edges_from_matrix(x$weights,
                                          directed = isTRUE(x$directed))
  x$n_edges <- nrow(x$edges)
  x
}


# ---- Exported verbs ----

#' Prune a Network's Edges
#'
#' Removes weak or non-significant edges from a network, keeping a record so
#' the operation can be reversed. This is Nestimate's counterpart of
#' \code{tna::prune()}; the \code{net_} prefix avoids a name clash with
#' \code{tna::prune()}.
#'
#' Pruning is non-destructive: the pruned network carries a \code{"pruning"}
#' attribute holding the original weights, the pruned weights, the parameters
#' used, and a tidy table of removed edges. Use \code{\link{net_deprune}} to
#' restore the original weights and \code{\link{net_reprune}} to re-apply the
#' pruning, both without recomputation. \code{\link{net_pruning_details}}
#' reports what was removed.
#'
#' Methods:
#' \describe{
#'   \item{\code{"threshold"}}{Remove edges with weight \eqn{\le}
#'     \code{threshold}.}
#'   \item{\code{"lowest"}}{Remove the lowest \code{lowest} quantile of
#'     non-zero edges.}
#'   \item{\code{"disparity"}}{Serrano disparity-filter backbone at
#'     significance \code{level}.}
#'   \item{\code{"bootstrap"}}{Remove edges deemed non-significant by
#'     \code{\link{bootstrap_network}} (pass a precomputed result via
#'     \code{boot}, or extra bootstrap arguments via \code{...}).}
#' }
#' For \code{"threshold"}, \code{"lowest"}, and \code{"disparity"} an edge is
#' dropped only when its removal leaves the network weakly connected.
#'
#' Diagonal self-loops (self-transitions) are observed data: they are counted
#' equally when computing the cut-off but are never removed by any method.
#' (This is a deliberate divergence from \code{tna::prune()}, which prunes
#' self-loops like any other edge.)
#'
#' @param x A \code{netobject} or \code{netobject_group}.
#' @param method One of \code{"threshold"}, \code{"lowest"},
#'   \code{"disparity"}, \code{"bootstrap"}. Default \code{"threshold"}.
#' @param threshold Numeric cut-off for \code{method = "threshold"}.
#'   Default \code{0.1}.
#' @param lowest Quantile (0-1) for \code{method = "lowest"}. Default
#'   \code{0.05}.
#' @param level Significance level (0-1) for \code{method = "disparity"}.
#'   Default \code{0.5}.
#' @param boot Optional precomputed \code{net_bootstrap} for
#'   \code{method = "bootstrap"}.
#' @param ... Passed to \code{\link{bootstrap_network}} when
#'   \code{method = "bootstrap"}.
#'
#' @return The input network (or group) with pruned \code{$weights} and a
#'   \code{"pruning"} attribute. Class is unchanged.
#'
#' @seealso \code{\link{net_deprune}}, \code{\link{net_reprune}},
#'   \code{\link{net_pruning_details}}
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
#'   V3 = c("C","A","C","B","A"))
#' net    <- build_network(seqs, method = "relative")
#' pruned <- net_prune(net, method = "threshold", threshold = 0.2)
#' net_pruning_details(pruned)
#'
#' @export
net_prune <- function(x, method = "threshold", threshold = 0.1,
                      lowest = 0.05, level = 0.5, boot = NULL, ...) {
  UseMethod("net_prune")
}

#' @rdname net_prune
#' @export
net_prune.netobject <- function(x, method = "threshold", threshold = 0.1,
                                lowest = 0.05, level = 0.5, boot = NULL, ...) {
  method <- match.arg(method,
                      c("threshold", "lowest", "disparity", "bootstrap"))
  if (!is.null(attr(x, "pruning"))) {
    stop("The network has already been pruned. ",
         "Call net_deprune() first to undo it.", call. = FALSE)
  }
  stopifnot(is.numeric(threshold), length(threshold) == 1L)
  if (lowest < 0 || lowest > 1) stop("'lowest' must be in [0, 1].", call. = FALSE)
  if (level  < 0 || level  > 1) stop("'level' must be in [0, 1].",  call. = FALSE)

  labels <- x$nodes$label %||% as.character(seq_len(nrow(x$weights)))
  tmp <- switch(method,
                bootstrap = .prune_bootstrap(x, boot, labels, list(...)),
                disparity = .prune_disparity(x$weights, level, labels),
                .prune_default(x$weights, method, threshold, lowest, labels))

  tmp$original  <- x$weights
  tmp$active    <- TRUE
  x$weights     <- tmp$weights
  attr(x, "pruning") <- tmp
  .sync_netobject_edges(x)
}

#' @rdname net_prune
#' @export
net_prune.netobject_group <- function(x, method = "threshold", threshold = 0.1,
                                      lowest = 0.05, level = 0.5, boot = NULL,
                                      ...) {
  result <- lapply(x, net_prune.netobject, method = method,
                   threshold = threshold, lowest = lowest, level = level,
                   boot = boot, ...)
  names(result) <- names(x)
  class(result) <- class(x)
  result
}

#' @rdname net_prune
#' @export
net_prune.default <- function(x, method = "threshold", threshold = 0.1,
                              lowest = 0.05, level = 0.5, boot = NULL, ...) {
  stop("net_prune() requires a 'netobject' or 'netobject_group'; got '",
       class(x)[1L], "'.", call. = FALSE)
}


#' Undo Network Pruning
#'
#' Restores the original (pre-pruning) weights of a network pruned by
#' \code{\link{net_prune}}, without recomputation. The pruning record is kept,
#' so \code{\link{net_reprune}} can re-apply it.
#'
#' @param x A pruned \code{netobject} or \code{netobject_group}.
#' @param ... Ignored.
#' @return The network (or group) with original weights restored and its
#'   pruning marked inactive.
#' @seealso \code{\link{net_prune}}, \code{\link{net_reprune}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"),
#'                    V3 = c("C","A","B"))
#' net    <- build_network(seqs, method = "relative")
#' pruned <- net_prune(net, threshold = 0.2)
#' net_deprune(pruned)
#' @export
net_deprune <- function(x, ...) UseMethod("net_deprune")

#' @rdname net_deprune
#' @export
net_deprune.netobject <- function(x, ...) {
  tmp <- attr(x, "pruning")
  if (is.null(tmp)) stop("'x' must have been pruned.", call. = FALSE)
  if (!isTRUE(tmp$active)) stop("Pruning is already inactive for 'x'.",
                                call. = FALSE)
  tmp$active <- FALSE
  x$weights  <- tmp$original
  attr(x, "pruning") <- tmp
  .sync_netobject_edges(x)
}

#' @rdname net_deprune
#' @export
net_deprune.netobject_group <- function(x, ...) {
  result <- lapply(x, net_deprune.netobject)
  names(result) <- names(x)
  class(result) <- class(x)
  result
}

#' @rdname net_deprune
#' @export
net_deprune.default <- function(x, ...) {
  stop("net_deprune() requires a pruned 'netobject' or 'netobject_group'.",
       call. = FALSE)
}


#' Re-apply Network Pruning
#'
#' Re-applies a previously computed pruning that was undone by
#' \code{\link{net_deprune}}, without recomputation.
#'
#' @param x A depruned \code{netobject} or \code{netobject_group}.
#' @param ... Ignored.
#' @return The network (or group) with pruned weights re-applied and its
#'   pruning marked active.
#' @seealso \code{\link{net_prune}}, \code{\link{net_deprune}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"),
#'                    V3 = c("C","A","B"))
#' net    <- build_network(seqs, method = "relative")
#' pruned <- net_prune(net, threshold = 0.2)
#' undone <- net_deprune(pruned)
#' net_reprune(undone)
#' @export
net_reprune <- function(x, ...) UseMethod("net_reprune")

#' @rdname net_reprune
#' @export
net_reprune.netobject <- function(x, ...) {
  tmp <- attr(x, "pruning")
  if (is.null(tmp)) stop("'x' must have been pruned.", call. = FALSE)
  if (isTRUE(tmp$active)) stop("Pruning is already active for 'x'.",
                               call. = FALSE)
  tmp$active <- TRUE
  x$weights  <- tmp$weights
  attr(x, "pruning") <- tmp
  .sync_netobject_edges(x)
}

#' @rdname net_reprune
#' @export
net_reprune.netobject_group <- function(x, ...) {
  result <- lapply(x, net_reprune.netobject)
  names(result) <- names(x)
  class(result) <- class(x)
  result
}

#' @rdname net_reprune
#' @export
net_reprune.default <- function(x, ...) {
  stop("net_reprune() requires a depruned 'netobject' or 'netobject_group'.",
       call. = FALSE)
}


#' Report Network Pruning Details
#'
#' Returns the edges removed by \code{\link{net_prune}} as a tidy
#' one-row-per-edge data frame, with the method, cut-off, and retained/removed
#' counts attached as attributes and shown by its print method.
#'
#' @param x A pruned \code{netobject} or \code{netobject_group}.
#' @param ... Ignored.
#' @return For a \code{netobject}: a \code{net_pruning_details} data frame
#'   (columns \code{from}, \code{to}, \code{weight}) of removed edges. For a
#'   \code{netobject_group}: a named list of such data frames.
#' @seealso \code{\link{net_prune}}
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
#'   V3 = c("C","A","C","B","A"))
#' net <- build_network(seqs, method = "relative")
#' net_pruning_details(net_prune(net, threshold = 0.2))
#' @export
net_pruning_details <- function(x, ...) UseMethod("net_pruning_details")

#' @rdname net_pruning_details
#' @export
net_pruning_details.netobject <- function(x, ...) {
  pruning <- attr(x, "pruning")
  if (is.null(pruning)) stop("'x' must have been pruned.", call. = FALSE)
  removed <- pruning$removed
  rownames(removed) <- NULL
  structure(removed,
            class = c("net_pruning_details", "data.frame"),
            method = pruning$method, threshold = pruning$threshold,
            lowest = pruning$lowest, level = pruning$level,
            cut_off = pruning$cut_off, num_removed = pruning$num_removed,
            num_retained = pruning$num_retained)
}

#' @rdname net_pruning_details
#' @export
net_pruning_details.netobject_group <- function(x, ...) {
  result <- lapply(x, net_pruning_details.netobject)
  names(result) <- names(x)
  result
}

#' @rdname net_pruning_details
#' @export
net_pruning_details.default <- function(x, ...) {
  stop("net_pruning_details() requires a pruned 'netobject' or ",
       "'netobject_group'.", call. = FALSE)
}

#' Print method for pruning details
#' @param x A \code{net_pruning_details} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.net_pruning_details <- function(x, ...) {
  method_txt <- switch(attr(x, "method"),
    threshold = paste0("user-specified threshold (", attr(x, "threshold"), ")"),
    lowest    = paste0("lowest ", attr(x, "lowest") * 100,
                       "% of non-zero edges"),
    disparity = paste0("disparity filter (sig. level = ", attr(x, "level"),
                       ")"),
    bootstrap = "bootstrap significance")
  cat("Pruning details\n")
  cat("  Method:  ", method_txt, "\n", sep = "")
  cat("  Removed: ", attr(x, "num_removed"), " edges\n", sep = "")
  cat("  Retained:", attr(x, "num_retained"), "edges\n\n")
  if (nrow(x) == 0L) {
    cat("(no edges removed)\n")
  } else {
    print(as.data.frame(x), row.names = FALSE)
  }
  invisible(x)
}
