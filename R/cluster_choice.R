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

#' Compare Clusterings Across k, Dissimilarity and Method
#'
#' @description
#' `cluster_choice()` clusters the same sequences under every combination of
#' the cluster counts, dissimilarities and algorithms given, each with
#' [build_clusters()], and returns one row per combination with its average
#' silhouette width, mean within-cluster distance and cluster-size balance.
#' It is the distance-based counterpart of [compare_mmm()], which compares
#' mixtures of Markov chains across k. `"all"` expands `dissimilarity` or
#' `method` to every supported option.
#'
#' The silhouette of a sequence (Rousseeuw, 1987) compares its mean distance
#' to its own cluster with that to the nearest other cluster: near 1 it sits
#' well inside its cluster, near 0 between two, and below 0 it is closer to
#' another cluster. Rows that share a dissimilarity compare k and algorithms
#' on one distance scale; across dissimilarities the silhouettes and
#' within-cluster distances rest on different distances and compare less
#' directly. The row with the highest silhouette is marked, but that is not
#' a decision: also weigh the cluster sizes (`min_size`, `size_ratio`), as
#' a high silhouette can come with a cluster of very few sequences; then
#' fit the chosen configuration with [build_clusters()] and read it with its
#' `summary()`, `plot()` and [cluster_diagnostics()].
#'
#' @param data The sequences, passed to [build_clusters()] unchanged: any
#'   input it takes, such as a wide data frame of sequences or a network
#'   built from sequences.
#' @param k The cluster counts to try. Default `2:5`. Each must be a whole
#'   finite number, or an error is raised; [build_clusters()] also requires
#'   each to be at least 2 and at most the number of sequences minus 1.
#' @param dissimilarity The dissimilarity measures to try. Default
#'   `"hamming"`. `"all"` expands to every supported measure: `"hamming"`,
#'   `"osa"`, `"lv"`, `"dl"`, `"lcs"`, `"qgram"`, `"cosine"`, `"jaccard"`
#'   and `"jw"`, and may stand beside named ones (each is fitted once).
#' @param method The clustering algorithms to try. Default `"ward.D2"`
#'   (where [build_clusters()] defaults to `"pam"`). `"all"` expands to every
#'   supported algorithm: `"pam"`, `"ward.D2"`, `"ward.D"`, `"complete"`,
#'   `"average"`, `"single"`, `"mcquitty"`, `"median"` and `"centroid"`.
#' @param ... Further arguments of [build_clusters()], used for every
#'   combination: `weighted`, `lambda`, `q`, `p`, `seed`, `na_syms`,
#'   `covariates` and `estimator`; any other name raises an error.
#'   `weighted = TRUE` works only with `dissimilarity = "hamming"` and is
#'   refused before any fit when the sweep holds another measure. The
#'   `print()`, `summary()` and `plot()` methods take no further arguments
#'   and raise an error on any.
#' @return A `cluster_choice` (a data frame subclass) with one row per
#'   combination, `k` varying fastest, then `dissimilarity`, then `method`.
#'   `print()` shows the table, with only the settings that vary as columns,
#'   and flags the row with the highest silhouette in a `best` column;
#'   `summary()` returns the full table as a plain data frame with a
#'   character `best` column that holds `"silhouette"` on that row;
#'   `plot()` draws the silhouette across the swept settings (see *Plot
#'   types*) and returns the `ggplot` invisibly (with `type = "facet"` and
#'   `combined = FALSE`, a named list of `ggplot`s, returned invisibly and
#'   not drawn). The `print()` method returns its input invisibly. The
#'   columns:
#'   \describe{
#'     \item{k, dissimilarity, method}{The configuration of the row.}
#'     \item{silhouette}{The average silhouette width over all sequences,
#'       computed by [build_clusters()] (from `cluster::pam()` for
#'       `"pam"`, with `cluster::silhouette()` otherwise).}
#'     \item{mean_within_dist}{The mean distance between members of the
#'       same cluster, averaged over clusters weighted by their size (a
#'       one-member cluster counts as 0), in the units of the row's
#'       dissimilarity.}
#'     \item{min_size, max_size, size_ratio}{The sizes of the smallest and
#'       the largest cluster, and their ratio `max_size / min_size`.}
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
#' @references
#' Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the
#' interpretation and validation of cluster analysis. \emph{Journal of
#' Computational and Applied Mathematics}, 20, 53--65.
#' \doi{10.1016/0377-0427(87)90125-7}
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
  if (!is.numeric(k) || !length(k) || any(!is.finite(k)) ||
      any(k != floor(k))) {
    stop(errorCondition(
      "'k' must be a non-empty vector of whole finite numbers.",
      class = "nestimate_bad_input", call = NULL))
  }
  k             <- as.integer(k)

  # ---- Pre-validate weighted/dissimilarity combination ---------------
  dots <- list(...)
  if (isTRUE(dots$weighted) &&
      any(dissimilarity != "hamming")) {
    stop(errorCondition(
      paste0(
        "weighted = TRUE requires dissimilarity = \"hamming\". ",
        "Got: ",
        paste(dissimilarity, collapse = ", "),
        "."),
      class = "nestimate_bad_input", call = NULL))
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
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    fit <- do.call(build_clusters, c(
      list(data          = data,
           k             = grid$k[[i]],
           dissimilarity = grid$dissimilarity[[i]],
           method        = grid$method[[i]]),
      dots))
    .cluster_choice_row(fit)
  })
  out <- cbind(grid, do.call(rbind, rows))
  rownames(out) <- NULL
  class(out) <- c("cluster_choice", "data.frame")
  swept_axes <- character()
  if (length(k)             > 1L) swept_axes <- c(swept_axes, "k")
  if (length(dissimilarity) > 1L) swept_axes <- c(swept_axes, "dissimilarity")
  if (length(method)        > 1L) swept_axes <- c(swept_axes, "method")
  attr(out, "swept") <- swept_axes
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
    stop(errorCondition(
      paste0("'", name, "' is empty."),
      class = "nestimate_bad_input", call = NULL))
  }
  values
}

#' @noRd
.cluster_choice_row <- function(fit) {
  sizes <- as.integer(fit$sizes %||% tabulate(fit$assignments))
  n     <- sum(sizes)
  k_    <- length(sizes)

  mean_within <- NA_real_
  if (!is.null(fit$distance) && n > 0L) {
    per_cluster <- .per_cluster_within_dist(fit$distance, fit$assignments,
                                              k_)
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

#' @rdname cluster_choice
#' @param x For the `print()` and `plot()` methods: a `cluster_choice`.
#' @param digits For the `print()` method: the decimal places of the
#'   silhouette, within-cluster distance and size ratio. Default `3`. Must
#'   be a single non-negative whole number.
#' @param object For the `summary()` method: a `cluster_choice`.
#' @param type For the `plot()` method: the kind of plot, one of `"auto"`
#'   (default), `"lines"`, `"bars"`, `"heatmap"`, `"tradeoff"` and
#'   `"facet"` (see *Plot types*). A type the sweep cannot support raises
#'   an error that names the alternatives, and so does a table with a
#'   single configuration, which has nothing to compare.
#' @param abbrev For the `plot()` method: whether to shorten the
#'   dissimilarity and method names (e.g. `"hamming"` to `"ham"`,
#'   `"ward.D2"` to `"wD2"`) on the axes, point labels and legend. Default
#'   `FALSE`.
#' @param combined For `plot(type = "facet")`: `TRUE` (default) draws one
#'   panel per method in a single plot; `FALSE` returns a named list of
#'   one-method plots instead. Ignored by the other types, but must still
#'   be a single logical value.
#' @export
print.cluster_choice <- function(x, digits = 3L, ...) {
  .cluster_choice_check_unused_dots("print.cluster_choice", ...)
  .cluster_choice_check_digits(digits)
  digits <- as.integer(digits)
  swept  <- attr(x, "swept") %||%
            c("k", "dissimilarity", "method")

  if (length(swept) == 0L) {
    only_method <- unique(x$method)[1L]
    only_dissim <- unique(x$dissimilarity)[1L]
    cat(sprintf("Cluster Choice [%s / %s]\n\n", only_method, only_dissim))
  } else {
    cat(sprintf("Cluster Choice (sweep: %s)\n\n",
                paste(swept, collapse = " x ")))
  }

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

  disp[setdiff(c("k", "dissimilarity", "method"), swept)] <- NULL

  print.data.frame(disp, row.names = FALSE, right = FALSE)

  invisible(x)
}

#' @rdname cluster_choice
#' @export
summary.cluster_choice <- function(object, ...) {
  .cluster_choice_check_unused_dots("summary.cluster_choice", ...)
  best_idx <- which.max(object$silhouette)
  best <- rep("", nrow(object))
  best[best_idx] <- "silhouette"
  out <- as.data.frame(object)
  out$best <- best
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Abbreviation helpers (display-only). The underlying $dissimilarity and
# $method columns always carry the canonical names; abbreviation is a
# rendering concern triggered by `abbrev = TRUE` on plot().
# ---------------------------------------------------------------------------

.dissimilarity_abbrev <- c(
  hamming = "ham", osa = "osa", lv  = "lv",  dl  = "dl",  lcs = "lcs",
  qgram   = "qgr", cosine = "cos", jaccard = "jac", jw  = "jw"
)

.method_abbrev <- c(
  pam = "pam", ward.D2 = "wD2", ward.D = "wD", complete = "cmp",
  average = "avg", single = "sng", mcquitty = "mcq", median = "mdn",
  centroid = "cen"
)

# Map canonical names -> abbreviations; pass through any unknown value.
#' @noRd
.abbrev_dissimilarity <- function(x) {
  out <- unname(.dissimilarity_abbrev[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

#' @noRd
.abbrev_method <- function(x) {
  out <- unname(.method_abbrev[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

#' @rdname cluster_choice
#' @section Plot types:
#' The plots show the values only: they carry no "best" marker and no
#' recommendation. `plot()` needs at least one swept setting.
#' \describe{
#'   \item{`"auto"`}{The default. Picks a type from the swept settings:
#'     `"lines"` when k is swept with at most one other setting, `"bars"`
#'     when only `dissimilarity` or only `method` is swept, `"heatmap"`
#'     when both are swept at a single k, and `"facet"` when all three are
#'     swept.}
#'   \item{`"lines"`}{The silhouette across k, one line per value of the
#'     other swept setting. When only k is swept, `mean_within_dist` is
#'     drawn as a second line on the same axis. With all three swept, the
#'     lines follow `dissimilarity` and the methods share them; use
#'     `"facet"` instead. Needs k swept.}
#'   \item{`"bars"`}{Horizontal bars of the silhouette per value of the
#'     swept setting, or per combination when several are swept, sorted by
#'     silhouette, with the values printed.}
#'   \item{`"heatmap"`}{Tiles of the silhouette over `dissimilarity` and
#'     `method`, with the values printed. Needs both swept, and is meant
#'     for a single k.}
#'   \item{`"tradeoff"`}{The silhouette against `size_ratio`, each point
#'     labelled with its swept settings.}
#'   \item{`"facet"`}{The silhouette across k, one line per dissimilarity
#'     and one panel per method. Needs all three swept.}
#' }
#' @export
plot.cluster_choice <- function(x,
                                 type   = c("auto", "lines", "bars",
                                             "heatmap", "tradeoff",
                                             "facet"),
                                 abbrev = FALSE,
                                 combined = TRUE,
                                 ...) {
  .cluster_choice_check_unused_dots("plot.cluster_choice", ...)
  if (!is.logical(abbrev) || length(abbrev) != 1L || is.na(abbrev)) {
    stop(errorCondition(
      "'abbrev' must be TRUE or FALSE.",
      class = "nestimate_bad_input", call = NULL))
  }
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))
  if (!requireNamespace("ggplot2", quietly = TRUE)) { # nocov start
    stop(errorCondition("Package 'ggplot2' required.",
                        class = "nestimate_missing_package", call = NULL))
  } # nocov end
  type  <- match.arg(type)
  swept <- attr(x, "swept") %||% c("k", "dissimilarity", "method")

  if (!length(swept)) {
    stop(errorCondition(
      paste0("Nothing to plot: this cluster_choice() holds one configuration. ",
             "Sweep at least one of 'k', 'dissimilarity' or 'method' ",
             "(for example k = 2:5)."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (type == "auto") type <- .auto_choice_type(swept)
  .require_type_supported(type, swept)

  df <- as.data.frame(x)
  if (abbrev) {
    df$dissimilarity <- .abbrev_dissimilarity(df$dissimilarity)
    df$method        <- .abbrev_method(df$method)
  }

  if (type == "facet" && !combined) {
    methods <- unique(as.character(df$method))
    plots <- lapply(methods, function(m) {
      sub <- df[df$method == m, , drop = FALSE]
      .plot_choice_facet(sub, swept) +
        ggplot2::labs(title = sprintf("Cluster Choice (method = %s)", m))
    })
    names(plots) <- methods
    return(invisible(plots))
  }

  p <- switch(type,
              lines    = .plot_choice_lines(df, swept),
              bars     = .plot_choice_bars(df, swept),
              heatmap  = .plot_choice_heatmap(df),
              tradeoff = .plot_choice_tradeoff(df, swept),
              facet    = .plot_choice_facet(df, swept))

  print(p)
  invisible(p)
}

#' @noRd
.cluster_choice_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(paste0(
    method, "() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", ")),
    class = "nestimate_bad_input", call = NULL)
  )
}

#' @noRd
.cluster_choice_check_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L ||
      !is.finite(digits) || digits != floor(digits) || digits < 0L) {
    stop(errorCondition(
      "'digits' must be a single non-negative whole finite number.",
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

#' @noRd
.auto_choice_type <- function(swept) {
  n_axes <- length(swept)
  if (n_axes <= 1L) {
    if ("k" %in% swept) "lines" else "bars"
  } else if (all(c("k", "dissimilarity", "method") %in% swept)) {
    "facet"
  } else if ("k" %in% swept) {
    "lines"
  } else {
    "heatmap"
  }
}

#' @noRd
.require_type_supported <- function(type, swept) {
  msg <- function(...) stop(errorCondition(
    paste0(...),
    class = "nestimate_bad_input", call = NULL))
  switch(type,
         lines = if (!("k" %in% swept))
           msg("type = \"lines\" requires k to be swept (length > 1). ",
               "Use type = \"bars\" or type = \"heatmap\" for a fixed-k ",
               "sweep."),
         heatmap = if (!all(c("dissimilarity", "method") %in% swept))
           msg("type = \"heatmap\" requires both dissimilarity and method ",
               "to be swept. Use type = \"bars\" for one categorical axis ",
               "or type = \"lines\" for a k sweep."),
         facet = if (!("k" %in% swept) || length(swept) < 3L)
           msg("type = \"facet\" requires k plus two categorical axes ",
               "(dissimilarity and method) to be swept. Use type = ",
               "\"lines\" for k + one categorical, or type = \"heatmap\" ",
               "for two categoricals at fixed k."),
         bars     = invisible(),
         tradeoff = invisible())
}

# ---------------------------------------------------------------------------
# Per-type plot builders (each returns a ggplot, no editorialising)
# ---------------------------------------------------------------------------

#' @noRd
.plot_choice_lines <- function(df, swept) {
  # k-only sweep: silhouette and within_dist as two y-series.
  if (identical(swept, "k")) {
    long <- data.frame(
      k = rep(df$k, 2L),
      value = c(df$silhouette, df$mean_within_dist),
      metric = rep(c("silhouette", "mean_within_dist"),
                   each = nrow(df)),
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(long, ggplot2::aes(x = .data$k, y = .data$value,
                                           colour = .data$metric,
                                           group  = .data$metric)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2.8) +
        ggplot2::scale_x_continuous(breaks = unique(df$k)) +
        ggplot2::labs(title = "Cluster Choice", x = "k", y = NULL,
                      colour = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "bottom")
    )
  }

  # k + one categorical axis: silhouette vs k, one line per axis level.
  colour_axis <- if ("dissimilarity" %in% swept) "dissimilarity" else
                 if ("method"        %in% swept) "method"        else NA

  aes_obj <- if (!is.na(colour_axis))
    ggplot2::aes(x = .data$k, y = .data$silhouette,
                 colour = .data[[colour_axis]],
                 group  = .data[[colour_axis]])
  else
    ggplot2::aes(x = .data$k, y = .data$silhouette)

  ggplot2::ggplot(df, aes_obj) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_x_continuous(breaks = unique(df$k)) +
    ggplot2::scale_colour_viridis_d(option = "D", end = 0.85) +
    ggplot2::labs(title = "Cluster Choice", x = "k", y = "silhouette",
                  colour = colour_axis) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

#' @noRd
.plot_choice_bars <- function(df, swept) {
  # Pick the categorical axis to bar against. If multiple axes vary,
  # flatten them into one composite "k / dissim / method" label.
  axis <- if (identical(swept, "dissimilarity")) "dissimilarity"
          else if (identical(swept, "method"))   "method"
          else if (identical(swept, "k"))        "k"
          else NA

  if (is.na(axis)) {
    label_parts <- list()
    if ("k"             %in% swept) label_parts$k    <- as.character(df$k)
    if ("dissimilarity" %in% swept) label_parts$diss <- df$dissimilarity
    if ("method"        %in% swept) label_parts$meth <- df$method
    df$.row <- do.call(paste, c(label_parts, sep = " / "))
    axis <- ".row"
  }

  ggplot2::ggplot(df,
                   ggplot2::aes(x = stats::reorder(.data[[axis]],
                                                     .data$silhouette),
                                  y = .data$silhouette)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f",
                                                      .data$silhouette)),
                        hjust = -0.15, size = 3.2, colour = "grey25") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult =
                                                              c(0, 0.18))) +
    ggplot2::labs(title = "Cluster Choice",
                  x = if (axis == ".row") NULL else axis,
                  y = "silhouette") +
    ggplot2::theme_minimal(base_size = 12)
}

#' @noRd
.plot_choice_heatmap <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$dissimilarity,
                                     y = .data$method,
                                     fill = .data$silhouette)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f",
                                                      .data$silhouette)),
                        size = 3.4, colour = "grey15") +
    ggplot2::scale_fill_viridis_c(option = "D", end = 0.95) +
    ggplot2::labs(title = "Cluster Choice", x = "dissimilarity",
                  y = "method", fill = "silhouette") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' @noRd
.plot_choice_tradeoff <- function(df, swept) {
  # Build a label that names whichever axes vary, so the scatter is
  # readable without a legend.
  parts <- list()
  if ("k"             %in% swept) parts$k    <- sprintf("k=%d", df$k)
  if ("dissimilarity" %in% swept) parts$diss <- df$dissimilarity
  if ("method"        %in% swept) parts$meth <- df$method
  df$.label <- do.call(paste, c(parts, sep = "/"))

  colour_axis <- if ("dissimilarity" %in% swept) "dissimilarity" else
                 if ("method"        %in% swept) "method"        else
                 if ("k"             %in% swept) "k"             else NA

  aes_obj <- if (!is.na(colour_axis))
    ggplot2::aes(x = .data$size_ratio, y = .data$silhouette,
                 colour = .data[[colour_axis]])
  else
    ggplot2::aes(x = .data$size_ratio, y = .data$silhouette)

  p <- ggplot2::ggplot(df, aes_obj) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = .data$.label),
                        hjust = -0.15, size = 3, colour = "grey25") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult =
                                                              c(0.05, 0.18))) +
    ggplot2::labs(title = "Cluster Choice",
                  x = "size_ratio", y = "silhouette",
                  colour = if (!is.na(colour_axis)) colour_axis else NULL) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.na(colour_axis) &&
      colour_axis %in% c("dissimilarity", "method")) {
    p <- p + ggplot2::scale_colour_viridis_d(option = "D", end = 0.85)
  }
  p
}

#' @noRd
.plot_choice_facet <- function(df, swept) {
  # k + dissim + method: lines vs k, colour = dissim, facet = method.
  ggplot2::ggplot(df, ggplot2::aes(x = .data$k, y = .data$silhouette,
                                     colour = .data$dissimilarity,
                                     group  = .data$dissimilarity)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = unique(df$k)) +
    ggplot2::scale_colour_viridis_d(option = "D", end = 0.85) +
    ggplot2::facet_wrap(~ method) +
    ggplot2::labs(title = "Cluster Choice", x = "k", y = "silhouette",
                  colour = "dissimilarity") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"))
}
