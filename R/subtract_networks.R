#' Subtract one network from another
#'
#' @description
#' `subtract_networks()` subtracts one network's weights from another's,
#' cell by cell, and returns the difference `x - y` as a `netdifference`
#' object: a network whose edge weights are the signed differences. It takes
#' any pair of networks on the same nodes (`netobject`s, `cograph_network`s
#' or square numeric matrices); for an edge-betweenness difference, subtract
#' two \code{\link{net_edge_betweenness}} results. `cograph::splot(d)` or
#' `cograph::plot_difference(d)` draws the signed difference network, with
#' cograph's colouring and node palette.
#'
#' The difference is taken between two point estimates and carries no
#' uncertainty, so a large difference on a rarely observed edge can be
#' noise. To see which edges differ beyond sampling variation, compare the
#' two networks with [permutation()] or [bayes_compare()];
#' [as_netdifference()] turns a [bayes_compare()] result into the same kind
#' of object, keeping, by default, only the credible differences.
#'
#' @param x,y The two networks: `netobject`s, `cograph_network`s or square
#'   numeric matrices. They must have the same number of nodes, with the
#'   same labels in the same order, or an error is raised; any other kind of
#'   object raises an error too. A matrix without dimnames is labelled `V1`,
#'   `V2`, ..., so it matches only a network with those labels. The
#'   difference is directed when `x` is a directed
#'   `netobject`; when `x` is not a `netobject`, it is directed when `x`'s
#'   weights are not symmetric. For the `print()` method: `x` is the
#'   `netdifference` object.
#'
#' @return A `netdifference` object, which is also a `netobject` (method
#'   `"difference"`). `print()` lists the edges present in either network,
#'   largest absolute difference first, with both weights and their
#'   difference (for an undirected difference, each pair once);
#'   `cograph::splot()` draws it. Beside the usual `netobject` fields, its
#'   `weights` and `difference_matrix` both hold `x - y`, and `x` and `y`
#'   hold the two weight matrices that were subtracted.
#'
#' @examples
#' early <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"))
#' late <- data.frame(
#'   V1 = c("B","A","C","B"), V2 = c("C","B","A","C"),
#'   V3 = c("A","C","B","A"))
#' a <- build_network(early, method = "relative")
#' b <- build_network(late, method = "relative")
#' subtract_networks(a, b)
#' # edge-betweenness difference:
#' subtract_networks(net_edge_betweenness(a), net_edge_betweenness(b))
#' @export
subtract_networks <- function(x, y) {
  wx <- .weights_of(x)
  wy <- .weights_of(y)
  if (!identical(dim(wx), dim(wy))) {
    stop(errorCondition("x and y must have the same number of nodes.",
      class = "nestimate_bad_input", call = NULL))
  }
  lx <- rownames(wx)
  ly <- rownames(wy)
  if (!is.null(lx) && !is.null(ly) && !identical(lx, ly)) {
    stop(errorCondition("x and y must have the same node labels and order.",
      class = "nestimate_bad_input", call = NULL))
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
#' @description
#' `as_netdifference()` turns the result of a network comparison into a
#' plain `netdifference` object, the signed difference network that
#' [subtract_networks()] returns, so that cograph's difference renderer
#' draws it. It reads a [bayes_compare()] result (`net_bayes`), and returns
#' a `netdifference` unchanged. The drawn weights are the posterior mean
#' differences, by default only those of the credible edges, while the full
#' difference, its credible intervals and its probabilities of direction
#' are kept beside them.
#'
#' @param x The comparison: a `net_bayes` from [bayes_compare()], or a
#'   `netdifference`, which is returned as it is. Any other object raises an
#'   error.
#' @param ... Not used.
#' @return A `netdifference` object, which is also a `netobject` (method
#'   `"difference"`); `print()` lists its edges and `cograph::splot()` draws
#'   it. From a `net_bayes`, its `weights` are the posterior mean
#'   differences (all of them, or only the credible ones, by
#'   `significant_only`), and it keeps the full posterior mean difference
#'   (`difference_matrix`), the two posterior mean matrices (`x`, `y`), the
#'   interval bounds (`ci_lower`, `ci_upper`), the Bayesian p-equivalents
#'   (`p_values`), the probabilities of direction (`p_difference`), the
#'   credibility matrix (`sig`), `significant_only`, and
#'   `netdifference_type = "bayes_compare"`.
#' @examples
#' s1 <- data.frame(V1 = c("A", "B", "C"), V2 = c("B", "C", "A"))
#' s2 <- data.frame(V1 = c("A", "C", "B"), V2 = c("C", "B", "A"))
#' b <- bayes_compare(build_network(s1, method = "relative"),
#'                    build_network(s2, method = "relative"),
#'                    draws = 500, seed = 1)
#' as_netdifference(b, significant_only = FALSE)
#' @export
as_netdifference <- function(x, ...) {
  UseMethod("as_netdifference")
}

#' @rdname as_netdifference
#' @param significant_only For a `net_bayes`: whether the drawn weights keep
#'   only the credible differences (those [bayes_compare()] marks in `sig`),
#'   with every other edge set to 0. The full difference and interval
#'   matrices are kept either way. Default `TRUE`.
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
  stop(errorCondition(paste0("Cannot coerce object of class '", class(x)[1L],
       "' to netdifference."),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname subtract_networks
#' @param max_print For the `print()` method: the most edges to list; the
#'   number left out is reported. Default `12`.
#' @param ... For the `print()` method: ignored.
#' @return The `print()` method returns its input invisibly.
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
