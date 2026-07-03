#' Subtract one network from another
#'
#' Returns \code{x - y} as a \code{netdifference} object: the element-wise
#' difference of the two weight matrices. Works on any pair of networks; for an
#' edge-betweenness difference, subtract two \code{\link{net_edge_betweenness}}
#' results. Draw the signed difference network with \code{cograph::splot(d)}
#' or \code{cograph::plot_difference(d)}; cograph handles the colouring and
#' node palette.
#'
#' @param x,y A \code{netobject}, \code{cograph_network}, or numeric square
#'   matrix. Both must share the same nodes in the same order.
#'
#' @return A \code{netdifference} object: a \code{netobject} whose
#'   \code{$weights} and \code{$difference_matrix} are \code{x - y}, carrying
#'   the source matrices \code{$x} and \code{$y}.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B","A"), V2 = c("B","C","B","A","C","B"),
#'   V3 = c("C","A","C","B","A","C"))
#' a <- build_network(seqs, method = "relative")
#' b <- build_network(seqs[1:4, ], method = "relative")
#' subtract_networks(a, b)
#' # edge-betweenness difference:
#' subtract_networks(net_edge_betweenness(a), net_edge_betweenness(b))
#' @export
subtract_networks <- function(x, y) {
  wx <- .weights_of(x)
  wy <- .weights_of(y)
  if (!identical(dim(wx), dim(wy))) {
    stop("x and y must have the same number of nodes.", call. = FALSE)
  }
  lx <- rownames(wx)
  ly <- rownames(wy)
  if (!is.null(lx) && !is.null(ly) && !identical(lx, ly)) {
    stop("x and y must have the same node labels and order.", call. = FALSE)
  }
  directed <- if (inherits(x, "netobject") && !is.null(x$directed)) {
    isTRUE(x$directed)
  } else {
    !isSymmetric(unname(wx))
  }
  out <- .wrap_netobject(wx - wy, method = "difference", directed = directed)
  out$x <- wx
  out$y <- wy
  out$difference_matrix <- out$weights
  # meta$splot contract (cograph renders by this hint, no class branch needed)
  out$meta$splot <- list(renderer = "difference", defaults = list(minimum = 0))
  class(out) <- c("netdifference", class(out))
  out
}

#' Coerce an inferential comparison to a network difference
#'
#' @param x An object with network-difference fields.
#' @param ... Additional arguments passed to methods.
#' @return A \code{netdifference} object suitable for \code{cograph::splot()}.
#' @export
as_netdifference <- function(x, ...) {
  UseMethod("as_netdifference")
}

#' @rdname as_netdifference
#' @param significant_only Logical. For inferential objects, keep only
#'   supported differences in the plotted weight matrix while retaining the
#'   full difference and interval matrices. Default \code{TRUE}.
#' @export
as_netdifference.net_bayes <- function(x, significant_only = TRUE, ...) {
  # A PURE netdifference (no net_bayes / net_permutation classes): the
  # coercion exists so cograph's difference renderer takes over. Keeping the
  # net_permutation class would send splot() to its permutation renderer
  # instead. $weights holds the display matrix (credible-only by default),
  # $difference_matrix the full posterior mean difference.
  display <- if (isTRUE(significant_only)) x$diff_sig else x$diff
  directed <- if (!is.null(x$x$directed)) isTRUE(x$x$directed) else TRUE
  out <- .wrap_netobject(display, method = "difference", directed = directed)
  out$x <- x$prob_x
  out$y <- x$prob_y
  out$difference_matrix <- x$diff
  out$ci_lower <- x$ci_lower
  out$ci_upper <- x$ci_upper
  out$p_values <- x$p_values
  out$p_difference <- x$p_difference
  out$sig <- x$sig
  out$netdifference_type <- "bayes_compare"
  out$significant_only <- isTRUE(significant_only)
  out$meta$splot <- list(renderer = "difference", defaults = list(minimum = 0))
  class(out) <- c("netdifference", class(out))
  out
}

#' @rdname as_netdifference
#' @export
as_netdifference.netdifference <- function(x, ...) {
  x
}

#' @rdname as_netdifference
#' @export
as_netdifference.default <- function(x, ...) {
  stop("Cannot coerce object of class '", class(x)[1L],
       "' to netdifference.", call. = FALSE)
}

#' @rdname subtract_networks
#' @param max_print Integer. Rows to show in \code{print()}. Default \code{12}.
#' @param ... Ignored.
#' @export
print.netdifference <- function(x, max_print = 12L, ...) {
  d <- x$weights
  labs <- rownames(d)
  directed <- isTRUE(x$directed)
  keep <- if (directed) {
    (x$x != 0 | x$y != 0)
  } else {
    (x$x != 0 | x$y != 0) & (row(d) <= col(d))
  }
  idx <- which(keep, arr.ind = TRUE)
  tab <- data.frame(
    from       = labs[idx[, 1L]],
    to         = labs[idx[, 2L]],
    x          = x$x[idx],
    y          = x$y[idx],
    difference = d[idx],
    stringsAsFactors = FALSE
  )
  tab <- tab[order(abs(tab$difference), decreasing = TRUE), , drop = FALSE]
  cat(sprintf("Network difference (x - y): %d nodes, %d differing edges\n",
              nrow(d), sum(tab$difference != 0)))
  cat("Plot: cograph::splot(d) or cograph::plot_difference(d)\n\n")
  print(utils::head(tab, max_print), row.names = FALSE)
  if (nrow(tab) > max_print) {
    cat(sprintf("\n... %d more edges.\n", nrow(tab) - max_print))
  }
  invisible(x)
}
