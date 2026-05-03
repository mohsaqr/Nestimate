# ==============================================================================
# cluster_choice() -- k + dissimilarity + method sweep for distance clustering.
#
# Parallel to compare_mmm() in R/mmm.R: one call, table back, print/summary/plot
# methods on the table. Where compare_mmm sweeps k for one MMM family,
# cluster_choice sweeps the cartesian product of (k, dissimilarity, method)
# for distance clustering. "all" sentinels expand to .clustering_metrics /
# .clustering_methods so users don't have to type out every option.
#
# Reuses build_clusters() per row -- no new metric formulas. Silhouette is
# already computed inside build_clusters; mean within-cluster distance is the
# same lower.tri block used in summary.net_clustering.
# ==============================================================================

#' Cluster Choice -- sweep k, dissimilarity and method
#'
#' One-call sweep across any combination of k, dissimilarity metric, and
#' clustering algorithm for distance-based sequence clustering. Mirrors
#' \code{\link{compare_mmm}} for model-based clustering: returns a data
#' frame with one row per swept configuration, a \code{best} marker on
#' the silhouette-max row in the print method, and a \code{plot()} that
#' adapts to the swept axes.
#'
#' @param data Sequence data (data frame or matrix) -- forwarded to
#'   \code{\link{build_clusters}}.
#' @param k Integer vector of cluster counts to sweep. Default
#'   \code{2:5}. Each value must be >= 2 and <= n - 1.
#' @param dissimilarity Character vector of dissimilarity metrics. Use
#'   \code{"all"} to expand to every supported metric:
#'   \code{c("hamming", "osa", "lv", "dl", "lcs", "qgram", "cosine",}
#'   \code{"jaccard", "jw")}. Default \code{"hamming"}.
#' @param method Character vector of clustering algorithms. Use
#'   \code{"all"} to expand to every supported method:
#'   \code{c("pam", "ward.D2", "ward.D", "complete", "average", "single",}
#'   \code{"mcquitty", "median", "centroid")}. Default \code{"ward.D2"}.
#' @param ... Other arguments forwarded to
#'   \code{\link{build_clusters}} (\code{weighted}, \code{lambda},
#'   \code{q}, \code{p}, \code{seed}, \code{na_syms}, \code{covariates}).
#'   Note: \code{weighted = TRUE} only works with
#'   \code{dissimilarity = "hamming"} and is rejected up-front when
#'   sweeping mixed dissimilarities.
#' @return A \code{cluster_choice} object (a data.frame subclass) with
#'   one row per (k, dissimilarity, method) combination and columns:
#'   \describe{
#'     \item{k, dissimilarity, method}{The configuration for that row.}
#'     \item{silhouette}{Overall average silhouette width (from
#'       \code{cluster::silhouette}, computed inside
#'       \code{\link{build_clusters}}).}
#'     \item{mean_within_dist}{Size-weighted mean of within-cluster
#'       distances, in the units of the row's dissimilarity.}
#'     \item{min_size, max_size, size_ratio}{Cluster-size balance bounds
#'       and their ratio (\code{max / min}).}
#'   }
#' @seealso \code{\link{build_clusters}}, \code{\link{compare_mmm}} for
#'   the model-based equivalent, \code{\link{cluster_diagnostics}} for
#'   the post-fit diagnostic surface on a single clustering.
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 40, TRUE),
#'                    V2 = sample(c("A","B","C"), 40, TRUE))
#' cluster_choice(seqs, k = 2:4)
#' \donttest{
#' # Sweep dissimilarities at fixed k
#' cluster_choice(seqs, k = 3, dissimilarity = c("hamming", "lcs", "jaccard"))
#'
#' # Full grid of k x dissimilarity
#' cluster_choice(seqs, k = 2:4, dissimilarity = c("hamming", "lcs"))
#'
#' # "all" sentinel
#' cluster_choice(seqs, k = 3, dissimilarity = "all")
#' }
#' @export
cluster_choice <- function(data,
                            k             = 2:5,
                            dissimilarity = "hamming",
                            method        = "ward.D2",
                            ...) {
  # ---- Expand "all" sentinels ----------------------------------------
  dissimilarity <- .expand_all(dissimilarity, .clustering_metrics,
                                "dissimilarity")
  method        <- .expand_all(method,        .clustering_methods,
                                "method")
  k             <- as.integer(k)

  # ---- Pre-validate weighted/dissimilarity combination ---------------
  dots <- list(...)
  if (isTRUE(dots$weighted) &&
      any(dissimilarity != "hamming")) {
    stop("weighted = TRUE requires dissimilarity = \"hamming\". ",
         "Got: ", paste(dissimilarity, collapse = ", "), ".",
         call. = FALSE)
  }

  # ---- Build the cartesian-product grid ------------------------------
  # Order: k inside dissimilarity inside method (so the printed table
  # reads like rows of a method block with dissim sub-blocks).
  grid <- expand.grid(
    k             = k,
    dissimilarity = dissimilarity,
    method        = method,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS  = FALSE
  )

  # ---- Sweep ---------------------------------------------------------
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fit <- do.call(build_clusters, c(
      list(data          = data,
           k             = grid$k[[i]],
           dissimilarity = grid$dissimilarity[[i]],
           method        = grid$method[[i]]),
      dots))
    rows[[i]] <- .cluster_choice_row(fit)
  }
  out <- cbind(grid, do.call(rbind, rows))
  rownames(out) <- NULL
  class(out) <- c("cluster_choice", "data.frame")
  attr(out, "swept") <- list(k = length(k) > 1L,
                              dissimilarity = length(dissimilarity) > 1L,
                              method = length(method) > 1L)
  out
}

# Expand a user-supplied character vector that may contain "all" into the
# canonical full set. Deduplicates while preserving order. Empty / NA
# inputs fall through to build_clusters() which will raise the existing
# input-validation error.
#' @noRd
.expand_all <- function(values, full, name) {
  if (!is.character(values)) {
    values <- as.character(values)
  }
  if (any(values == "all")) {
    values <- unique(c(values[values != "all"], full))
  }
  if (!length(values)) {
    stop("'", name, "' is empty.", call. = FALSE)
  }
  values
}

# Per-fit row: read pre-computed silhouette, compute size-weighted mean
# within-cluster distance the same way summary.net_clustering does
# (R/cluster_data.R:891-905), and read cluster sizes for balance.
#' @noRd
.cluster_choice_row <- function(fit) {
  sizes <- as.integer(fit$sizes %||% tabulate(fit$assignments))
  n     <- sum(sizes)
  k_    <- length(sizes)

  mean_within <- NA_real_
  if (!is.null(fit$distance) && n > 0L) {
    dmat <- as.matrix(fit$distance)
    per_cluster <- vapply(seq_len(k_), function(cl) {
      members <- which(fit$assignments == cl)
      if (length(members) > 1L) {
        sub <- dmat[members, members]
        mean(sub[lower.tri(sub)])
      } else {
        0
      }
    }, numeric(1L))
    mean_within <- stats::weighted.mean(per_cluster, sizes)
  }

  data.frame(
    silhouette       = as.numeric(fit$silhouette),
    mean_within_dist = mean_within,
    min_size         = as.integer(min(sizes)),
    max_size         = as.integer(max(sizes)),
    size_ratio       = if (min(sizes) > 0L) max(sizes) / min(sizes) else
                         NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for cluster_choice
#'
#' @param x A \code{cluster_choice} object.
#' @param digits Integer. Decimal places for floating-point columns.
#'   Default \code{3L}.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.cluster_choice <- function(x, digits = 3L, ...) {
  digits <- as.integer(digits)
  swept  <- attr(x, "swept") %||%
            list(k = TRUE, dissimilarity = TRUE, method = TRUE)
  swept_names <- names(swept)[vapply(swept, isTRUE, logical(1L))]

  if (length(swept_names) == 0L) {
    only_method <- unique(x$method)[1L]
    only_dissim <- unique(x$dissimilarity)[1L]
    cat(sprintf("Cluster Choice [%s / %s]\n\n", only_method, only_dissim))
  } else {
    cat(sprintf("Cluster Choice (sweep: %s)\n\n",
                paste(swept_names, collapse = " x ")))
  }

  # Build a tidy display copy: round numerics, compress min/max into one
  # column, drop constant axis columns, prepend the silhouette-best
  # marker so it sits adjacent to silhouette.
  disp <- data.frame(
    k             = x$k,
    dissimilarity = x$dissimilarity,
    method        = x$method,
    silhouette    = round(x$silhouette,       digits),
    within_dist   = round(x$mean_within_dist, digits),
    sizes         = sprintf("[%d, %d]", x$min_size, x$max_size),
    ratio         = round(x$size_ratio,       digits),
    stringsAsFactors = FALSE
  )

  best_idx <- which.max(x$silhouette)
  best <- rep("", nrow(disp))
  best[best_idx] <- "<-- best"
  disp$best <- best

  # Drop axis columns that don't vary so a single-dissimilarity sweep
  # doesn't carry a redundant column through the printout.
  if (!isTRUE(swept$k))             disp$k             <- NULL
  if (!isTRUE(swept$dissimilarity)) disp$dissimilarity <- NULL
  if (!isTRUE(swept$method))        disp$method        <- NULL

  print.data.frame(disp, row.names = FALSE, right = FALSE)

  if (length(unique(x$dissimilarity)) > 1L) {
    cat("\nNote: within_dist is in the units of each row's ",
        "dissimilarity; cross-row comparisons are only meaningful for ",
        "silhouette.\n", sep = "")
  }

  invisible(x)
}

#' Summary Method for cluster_choice
#'
#' @param object A \code{cluster_choice} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with the swept configurations, all metrics, and
#'   a \code{best} character column flagging the silhouette-max row.
#' @export
summary.cluster_choice <- function(object, ...) {
  best_idx <- which.max(object$silhouette)
  best <- rep("", nrow(object))
  best[best_idx] <- "silhouette"
  out <- as.data.frame(object)
  out$best <- best
  rownames(out) <- NULL
  out
}

#' Plot Method for cluster_choice
#'
#' Adapts to the swept axes: a single-axis sweep produces the same shape
#' as \code{\link{plot.mmm_compare}} (silhouette + within-distance vs k);
#' multi-axis sweeps colour by the secondary axis and facet by any
#' tertiary axis.
#'
#' @param x A \code{cluster_choice} object.
#' @param ... Additional arguments (ignored).
#' @return A \code{ggplot} object, invisibly.
#' @export
plot.cluster_choice <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) { # nocov start
    stop("Package 'ggplot2' required.", call. = FALSE)
  } # nocov end

  swept <- attr(x, "swept") %||%
           list(k = TRUE, dissimilarity = TRUE, method = TRUE)

  # Single-axis sweep: replicate the mmm_compare shape -- silhouette and
  # mean_within_dist as two y-series. Long-format the two metrics so a
  # single ggplot draws both with a colour aesthetic.
  if (sum(unlist(swept)) <= 1L) {
    if (isTRUE(swept$k)) {
      df <- data.frame(
        k         = rep(x$k, 2),
        value     = c(x$silhouette, x$mean_within_dist),
        metric    = rep(c("silhouette", "mean_within_dist"),
                        each = nrow(x)),
        stringsAsFactors = FALSE
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$k, y = .data$value,
                                              color = .data$metric)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_x_continuous(breaks = unique(x$k)) +
        ggplot2::labs(x = "k", y = "value",
                      title = "Cluster Choice", color = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
    } else {
      # dissimilarity-only or method-only sweep -> bar chart of silhouette.
      axis <- if (isTRUE(swept$dissimilarity)) "dissimilarity" else "method"
      df <- data.frame(
        x_axis     = x[[axis]],
        silhouette = x$silhouette,
        stringsAsFactors = FALSE
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(
                                              .data$x_axis,
                                              -.data$silhouette),
                                              y = .data$silhouette)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::labs(x = axis, y = "silhouette",
                      title = "Cluster Choice") +
        ggplot2::theme_minimal()
    }
    print(p)
    return(invisible(p))
  }

  # Multi-axis sweep: silhouette vs k, colour by dissimilarity, facet
  # by method. Falls back gracefully when k isn't swept.
  if (isTRUE(swept$k)) {
    color_axis <- if (isTRUE(swept$dissimilarity)) "dissimilarity" else
                  if (isTRUE(swept$method))        "method"        else NA
    facet_axis <- if (isTRUE(swept$dissimilarity) && isTRUE(swept$method))
                    "method" else NA

    aes_color <- if (!is.na(color_axis))
      ggplot2::aes(x = .data$k, y = .data$silhouette,
                   color = .data[[color_axis]]) else
      ggplot2::aes(x = .data$k, y = .data$silhouette)

    p <- ggplot2::ggplot(x, aes_color) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_x_continuous(breaks = unique(x$k)) +
      ggplot2::labs(x = "k", y = "silhouette",
                    title = "Cluster Choice",
                    color = if (!is.na(color_axis)) color_axis else NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    if (!is.na(facet_axis)) {
      p <- p + ggplot2::facet_wrap(stats::as.formula(
                paste0("~ ", facet_axis)))
    }
  } else {
    # No k axis: treat dissimilarity x method as a heatmap on silhouette.
    p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$dissimilarity,
                                          y = .data$method,
                                          fill = .data$silhouette)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f",
                                                       .data$silhouette))) +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
      ggplot2::labs(title = "Cluster Choice (silhouette)",
                    x = "dissimilarity", y = "method") +
      ggplot2::theme_minimal()
  }

  print(p)
  invisible(p)
}
