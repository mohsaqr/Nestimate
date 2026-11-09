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

#' The removal record shared by every pruning method
#'
#' One entry per edge: a cell of a directed network, an unordered pair of an
#' undirected one (upper triangle, self-loops included). Edges are the
#' non-zero cells, negative ones included, so the retained count is the
#' non-zero edges left after pruning.
#' @noRd
.prune_record <- function(W, pruned, labels, directed) {
  scope <- if (directed) matrix(TRUE, nrow(W), ncol(W)) else
    upper.tri(W, diag = TRUE)
  rem_idx <- which(scope & W != 0 & pruned == 0, arr.ind = TRUE)
  n_edges <- sum(scope & W != 0)
  list(
    removed = data.frame(
      from   = labels[rem_idx[, 1L]],
      to     = labels[rem_idx[, 2L]],
      weight = W[rem_idx],
      stringsAsFactors = FALSE
    ),
    num_removed  = nrow(rem_idx),
    num_retained = n_edges - nrow(rem_idx)
  )
}

#' Greedy threshold / lowest-quantile pruning with a connectivity guard
#'
#' Candidate edges (weight <= cut-off) are visited in column-major order and
#' removed one at a time only when removal keeps the network weakly connected,
#' exactly matching \code{tna:::prune_default()} on a directed network. On an
#' undirected network an edge is a pair: the candidates are the upper
#' triangle, a removal zeroes both cells, and the quantile of `"lowest"`
#' counts each pair once, so the matrix stays symmetric.
#' @noRd
.prune_default <- function(W, method, threshold, lowest, labels,
                           directed = TRUE) {
  scope    <- if (directed) matrix(TRUE, nrow(W), ncol(W)) else
    upper.tri(W, diag = TRUE)
  pos      <- W > 0 & scope
  pos_vals <- W[pos]                                   # diagonal counts equally
  cut_off  <- switch(method,
                     threshold = threshold,
                     lowest    = stats::quantile(pos_vals, probs = lowest))
  # Self-loops (diagonal) are observed data: they inform the cut-off but are
  # never removed, so only off-diagonal edges are removal candidates.
  off_diag <- row(W) != col(W)
  cand <- which(pos & W <= cut_off & off_diag, arr.ind = TRUE)  # column-major

  pruned <- Reduce(function(cur, j) {
    row <- cand[j, 1L]
    col <- cand[j, 2L]
    tmp <- cur
    tmp[row, col] <- 0
    if (!directed) tmp[col, row] <- 0
    if (.is_weakly_connected(tmp)) tmp else cur
  }, seq_len(nrow(cand)), W)

  c(list(weights = pruned, method = method, threshold = threshold,
         lowest = lowest, cut_off = cut_off),
    .prune_record(W, pruned, labels, directed))
}

#' Serrano disparity filter (backbone significance) -- matches tna
#' @noRd
.disparity_filter <- function(W, level) {
  d          <- ncol(W)
  idx_mat    <- 1 * (W > 0)
  # A state with no outgoing (incoming) edges has nothing to normalize;
  # dividing by 1 keeps its absent edges at 0 instead of 0 / 0 = NaN, which
  # until 1.0.5 left NA weights around a source or sink state.
  out_sums   <- .rowSums(W, m = d, n = d)
  in_sums    <- .colSums(W, m = d, n = d)
  out_sums[!(out_sums > 0)] <- 1
  in_sums[!(in_sums > 0)]   <- 1
  out_edges  <- W / out_sums
  out_degree <- .rowSums(idx_mat, m = d, n = d)
  out_p      <- (1 - out_edges)^(out_degree - 1)
  in_edges   <- t(W) / in_sums
  in_degree  <- .colSums(idx_mat, m = d, n = d)
  in_p       <- t((1 - in_edges)^(in_degree - 1))
  p_values   <- pmin(out_p, in_p)
  sig        <- 1 * (p_values < level)
  diag(sig)  <- 0
  sig
}

#' @noRd
.prune_disparity <- function(W, level, labels, directed = TRUE) {
  pruned  <- .disparity_filter(W, level) * W
  dimnames(pruned) <- dimnames(W)
  diag(pruned) <- diag(W)                 # self-loops are data: never removed
  c(list(weights = pruned, method = "disparity", level = level),
    .prune_record(W, pruned, labels, directed))
}

#' @noRd
.prune_bootstrap <- function(x, boot, labels, dots) {
  if (is.null(boot)) {
    boot <- do.call(bootstrap_network, c(list(x), dots))
  }
  if (!inherits(boot, "net_bootstrap")) {
    stop(errorCondition(
      "'boot' must be a 'net_bootstrap' from bootstrap_network().",
      class = "nestimate_bad_input", call = NULL))
  }
  W <- x$weights
  # The significant matrix replaces the weights, so it must come from this
  # network: a bootstrap of another one would paste in its weights.
  boot_w <- boot$original$weights
  if (is.null(boot_w) || !identical(dim(boot_w), dim(W)) ||
      !isTRUE(all.equal(unname(boot_w), unname(W)))) {
    stop(errorCondition(paste0(
      "'boot' was not computed on this network: its weights differ from ",
      "the network's. Pass bootstrap_network() of this network, or leave ",
      "`boot = NULL`."), class = "nestimate_bad_input", call = NULL))
  }
  pruned  <- boot$significant            # significant edges keep weight, rest 0
  dimnames(pruned) <- dimnames(W)
  diag(pruned) <- diag(W)                # self-loops are data: never removed
  c(list(weights = pruned, method = "bootstrap"),
    .prune_record(W, pruned, labels, !isFALSE(x$directed)))
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
#' @description
#' `net_prune()` removes weak or non-significant edges from a network and
#' returns the pruned network with a record of what it removed, so the
#' pruning can be undone and re-applied. It is Nestimate's counterpart of
#' `tna::prune()`; the `net_` prefix avoids a name clash with it. It takes a
#' `netobject` or a `netobject_group`, whose networks are pruned one by one
#' with the same settings, and offers four methods (see *Methods*).
#'
#' Pruning is non-destructive: the network carries a `"pruning"` attribute
#' with the original weights, the pruned weights, the settings used and a
#' table of the removed edges. [net_deprune()] restores the original weights
#' and [net_reprune()] re-applies the pruning, both without recomputing
#' anything; [net_pruning_details()] lists what was removed.
#'
#' `"threshold"` and `"lowest"` remove small weights; they do not show that
#' those edges are absent, and the cut-off is the analyst's choice.
#' `"disparity"` and `"bootstrap"` keep an edge on a statistical criterion.
#' Before reading a pruned network, see which edges went with
#' [net_pruning_details()] and whether the conclusions hold at another
#' cut-off; [bootstrap_network()] gives the uncertainty of each edge.
#'
#' @section Methods:
#' \describe{
#'   \item{`"threshold"`}{Removes positive edges with weight at most
#'     `threshold`.}
#'   \item{`"lowest"`}{Removes positive edges with weight at most the
#'     `lowest` quantile of the positive weights (self-loops included in
#'     the quantile; each pair counted once on an undirected network).}
#'   \item{`"disparity"`}{Keeps the backbone found by the disparity filter
#'     of Serrano, Boguna & Vespignani (2009): an edge stays when, at its
#'     source's outgoing or its target's incoming edges, its share of the
#'     strength is significant at `level` against a uniform split of that
#'     strength; the filter is defined for non-negative weights.}
#'   \item{`"bootstrap"`}{Removes the edges that [bootstrap_network()] finds
#'     non-significant: pass a computed result as `boot`, or bootstrap
#'     arguments through `...`.}
#' }
#' For `"threshold"` and `"lowest"`, the candidate edges are visited one at
#' a time, column by column and down each column, and each is removed
#' only when the network
#' stays weakly connected without it, so a network that is not weakly
#' connected to begin with loses no edge to these two methods.
#' `"disparity"` and `"bootstrap"` apply no connectivity check. Negative
#' edges are never candidates for `"threshold"` and `"lowest"`, whatever
#' their size. On a signed network, `"disparity"` treats negative edges as
#' not significant and removes them at the default `level`; at a `level`
#' above 0.5, a negative edge leaving a state with no positive outgoing edge
#' (or entering one with no positive incoming edge) can stay. `"bootstrap"`
#' removes the non-significant edges of either sign.
#'
#' On an undirected network an edge is a pair of cells: it is removed in
#' both directions or kept in both, so the pruned matrix stays symmetric,
#' and the record lists and counts it once.
#'
#' Self-loops (the diagonal) are observed data: they count towards the
#' `"lowest"` quantile but are never removed by any method. This is a
#' deliberate difference from `tna::prune()`, which prunes self-loops like
#' any other edge.
#'
#' @param x The network: a `netobject` or a `netobject_group`. Another class
#'   raises an error. A network whose pruning is active raises an error of
#'   class `nestimate_bad_input`; after [net_deprune()] it can be pruned
#'   again with other settings, which replaces its pruning record.
#' @param method The pruning method: `"threshold"` (default), `"lowest"`,
#'   `"disparity"` or `"bootstrap"` (see *Methods*).
#' @param threshold For `method = "threshold"`: the largest weight that is
#'   removed. Default `0.1`. Ignored by the other methods, but must still
#'   be a single number.
#' @param lowest For `method = "lowest"`: the quantile, between 0 and 1, of
#'   the positive weights that sets the cut-off. Default `0.05`. Ignored by
#'   the other methods, but checked whatever the method: a value outside
#'   `[0, 1]` raises an error.
#' @param level For `method = "disparity"`: the significance level, between
#'   0 and 1; an edge is kept when its p-value is below it. Default `0.5`.
#'   Ignored by the other methods, but checked whatever the method: a value
#'   outside `[0, 1]` raises an error.
#' @param boot For `method = "bootstrap"`: a `net_bootstrap` from
#'   [bootstrap_network()] of this network. With the default `NULL`, the
#'   bootstrap is run here on `x`. For a `netobject_group`, pass the
#'   `net_bootstrap_group` of [bootstrap_network()] on the group, matched
#'   to the networks by name, or `NULL` to bootstrap each network on its
#'   own. Another object, a bootstrap of another network (its weights
#'   differ), a single `net_bootstrap` for a group, or a group bootstrap
#'   missing one of the networks raises an error of class
#'   `nestimate_bad_input`.
#' @param ... For `method = "bootstrap"` with `boot = NULL`: arguments
#'   passed to [bootstrap_network()], such as `iter`. Ignored otherwise.
#'
#' @return The input, of the same class, with the pruned weights, its edge
#'   table and edge count updated to them, and the `"pruning"` attribute
#'   described above. [net_pruning_details()] reads the record,
#'   [net_deprune()] and [net_reprune()] switch between the original and
#'   the pruned weights. For a `netobject_group`, each network is pruned
#'   and carries its own record.
#'
#' @references
#' Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the
#' multiscale backbone of complex weighted networks. \emph{Proceedings of
#' the National Academy of Sciences}, 106(16), 6483--6488.
#' \doi{10.1073/pnas.0808904106}
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
  # An active pruning must be undone first; a depruned network carries its
  # original weights again and is pruned afresh, replacing the record.
  if (isTRUE(attr(x, "pruning")$active)) {
    stop(errorCondition(paste0(
      "The network has already been pruned. Call net_deprune() first to ",
      "undo it, then net_prune() with the new settings."),
      class = "nestimate_bad_input", call = NULL))
  }
  .assert(is.numeric(threshold), length(threshold) == 1L)
  if (lowest < 0 || lowest > 1) stop(errorCondition(
    "'lowest' must be in [0, 1].",
    class = "nestimate_bad_input", call = NULL))
  if (level  < 0 || level  > 1) stop(errorCondition(
    "'level' must be in [0, 1].",
    class = "nestimate_bad_input", call = NULL))

  labels <- x$nodes$label %||% as.character(seq_len(nrow(x$weights)))
  directed <- !isFALSE(x$directed)
  tmp <- switch(method,
                bootstrap = .prune_bootstrap(x, boot, labels, list(...)),
                disparity = .prune_disparity(x$weights, level, labels,
                                             directed),
                .prune_default(x$weights, method, threshold, lowest, labels,
                               directed))

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
  boots <- .group_boots(x, boot)
  result <- lapply(names(x), function(nm) {
    net_prune.netobject(x[[nm]], method = method, threshold = threshold,
                        lowest = lowest, level = level, boot = boots[[nm]],
                        ...)
  })
  names(result) <- names(x)
  class(result) <- class(x)
  result
}

#' One bootstrap per network of a group, matched by name
#'
#' A single `net_bootstrap` belongs to one network; applying it to every
#' network of a group would replace their weights with its own, so only a
#' `net_bootstrap_group` covering every network is accepted.
#' @noRd
.group_boots <- function(x, boot) {
  if (is.null(boot)) return(stats::setNames(vector("list", length(x)), names(x)))
  if (!inherits(boot, "net_bootstrap_group")) {
    stop(errorCondition(paste0(
      "For a netobject_group, `boot` must be the net_bootstrap_group of ",
      "bootstrap_network(x), one bootstrap per network, or NULL; a single ",
      "net_bootstrap belongs to one network."),
      class = "nestimate_bad_input", call = NULL))
  }
  missing_nm <- setdiff(names(x), names(boot))
  if (length(missing_nm) > 0L) {
    stop(errorCondition(paste0(
      "`boot` has no bootstrap for: ", paste(missing_nm, collapse = ", "),
      "."), class = "nestimate_bad_input", call = NULL))
  }
  lapply(stats::setNames(names(x), names(x)), function(nm) boot[[nm]])
}

#' @rdname net_prune
#' @export
net_prune.default <- function(x, method = "threshold", threshold = 0.1,
                              lowest = 0.05, level = 0.5, boot = NULL, ...) {
  stop(errorCondition(
    paste0(
      "net_prune() requires a 'netobject' or 'netobject_group'; got '",
      class(x)[1L],
      "'."),
    class = "nestimate_bad_input", call = NULL))
}


#' Undo Network Pruning
#'
#' @description
#' `net_deprune()` restores the weights a network had before
#' [net_prune()], without recomputing anything, and marks its pruning
#' inactive. The pruning record stays on the network, so [net_reprune()]
#' can re-apply it and [net_pruning_details()] still lists what the
#' pruning removed. A `netobject_group` is restored network by network.
#'
#' @param x A pruned `netobject` or `netobject_group`. A network that was
#'   never pruned, or whose pruning is already inactive, raises an error,
#'   and so does another class.
#' @param ... Ignored.
#' @return The input, of the same class, with the original weights, its
#'   edge table and edge count updated to them, and its pruning marked
#'   inactive.
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
  if (is.null(tmp)) stop(errorCondition("'x' must have been pruned.",
    class = "nestimate_bad_input", call = NULL))
  if (!isTRUE(tmp$active)) stop(errorCondition(
    "Pruning is already inactive for 'x'.",
    class = "nestimate_bad_input", call = NULL))
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
  stop(errorCondition(
    "net_deprune() requires a pruned 'netobject' or 'netobject_group'.",
    class = "nestimate_bad_input", call = NULL))
}


#' Re-apply Network Pruning
#'
#' @description
#' `net_reprune()` re-applies a pruning that [net_deprune()] undid, without
#' recomputing it: the pruned weights stored by [net_prune()] come back and
#' the pruning is marked active again. A `netobject_group` is handled
#' network by network.
#'
#' @param x A depruned `netobject` or `netobject_group`. A network that was
#'   never pruned, or whose pruning is already active, raises an error, and
#'   so does another class.
#' @param ... Ignored.
#' @return The input, of the same class, with the pruned weights, its edge
#'   table and edge count updated to them, and its pruning marked active.
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
  if (is.null(tmp)) stop(errorCondition("'x' must have been pruned.",
    class = "nestimate_bad_input", call = NULL))
  if (isTRUE(tmp$active)) stop(errorCondition(
    "Pruning is already active for 'x'.",
    class = "nestimate_bad_input", call = NULL))
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
  stop(errorCondition(
    "net_reprune() requires a depruned 'netobject' or 'netobject_group'.",
    class = "nestimate_bad_input", call = NULL))
}


#' Report Network Pruning Details
#'
#' @description
#' `net_pruning_details()` lists the edges that [net_prune()] removed from a
#' network, one row per edge, and prints the method with its setting and
#' the numbers of edges removed and retained above the table. It reads the
#' pruning record, so it works on a network whose pruning was undone by
#' [net_deprune()] as well.
#'
#' @param x A pruned `netobject` or `netobject_group`; a network that was
#'   never pruned, or another class, raises an error. For the `print()`
#'   method: a `net_pruning_details`.
#' @param ... Ignored.
#' @return For a `netobject`: a `net_pruning_details` data frame of the
#'   removed edges, with columns `from`, `to` and `weight` (the weight
#'   before pruning). The method, its setting, the cut-off and the counts
#'   are attributes; `print()` shows the method, its setting and the
#'   counts. The counts are edges: non-zero cells of either sign
#'   (self-loops included), one per pair on an undirected network. The
#'   retained count is the non-zero edges left after pruning. For a
#'   `netobject_group`: a named list of such data frames, one per network.
#'   The `print()` method returns its input invisibly.
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
  if (is.null(pruning)) stop(errorCondition("'x' must have been pruned.",
    class = "nestimate_bad_input", call = NULL))
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
  stop(errorCondition(
    paste0(
      "net_pruning_details() requires a pruned 'netobject' or ",
      "'netobject_group'."),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname net_pruning_details
#' @export
print.net_pruning_details <- function(x, ...) {
  method_txt <- switch(attr(x, "method"),
    threshold = paste0("user-specified threshold (", attr(x, "threshold"), ")"),
    lowest    = paste0("lowest ", attr(x, "lowest") * 100,
                       "% of positive edges"),
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
