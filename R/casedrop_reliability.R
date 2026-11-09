# ---- Edge-weight case-dropping stability --------------------------------
# Mirrors centrality_stability() but the stability target is the vector of
# edge weights (flattened adjacency) rather than per-node centralities.
# Answers: "How many cases can be dropped before the full set of edge
# weights loses rank correlation with the original?"

#' Edge-weight Case-dropping Stability
#'
#' @description
#' `casedrop_reliability()` measures how stable a network's edge weights are
#' when cases are left out. It re-estimates the network on random subsets
#' that drop each proportion in `drop_prop` of the cases, compares each
#' re-estimated edge vector with the original, and reports a CS-coefficient
#' for the edge weights: the largest proportion of cases that can be dropped
#' while the correlation between the original and the re-estimated edge
#' weights stays at or above `threshold` in at least `certainty` of the
#' iterations (the case-dropping CS-coefficient of Epskamp, Borsboom & Fried,
#' 2018, applied to the edge weights). Beside the correlation it records the
#' mean, median and largest absolute change of the edge weights. The cases
#' are the sequences of a transition network and the rows of the data
#' otherwise. It takes a single network, or a `netobject_group` or an `mcml`,
#' whose networks are assessed one by one.
#'
#' It complements [centrality_stability()]: that asks whether centrality
#' *rankings* are stable, this asks whether the *edge-weight structure
#' itself* is. For an `mcml` built from an edge list, where each row of the
#' data is one transition, this is case-dropping of **edges**.
#'
#' The CS-coefficient is the largest proportion in `drop_prop` that passes,
#' so it can only be one of the values tried. A smaller proportion need not
#' have passed as well, so look at the correlation at every proportion, not
#' only at the CS-coefficient: [summary()] gives its mean, median and
#' quantiles (`q025`, `q975`) per proportion, and `plot()` draws it. A correlation near 1 can hide a
#' large change in a few edges, which the largest absolute change
#' (`max_abs_dev`) shows. Stability under case-dropping says the estimate
#' does not rest on a few cases; it gives no interval for a single edge,
#' which [bootstrap_network()] does.
#'
#' @param x The network: a `netobject` or `cograph_network` that carries its
#'   data, or a `netobject_group` or `mcml`, whose networks (for an `mcml`,
#'   the macro network and each cluster's) are assessed one by one. A
#'   network without data raises an error, and one that rests on a single
#'   sequence, or one-hot data with no actor or session, has no cases to
#'   drop and raises an error of class `nestimate_single_sequence`. For the
#'   `print()` and `plot()` methods: the result of
#'   `casedrop_reliability()` (for `print()`, also the `summary()` of a
#'   group result).
#' @param iter The number of random subsets drawn for each drop proportion,
#'   at least 2. Default `1000`.
#' @param drop_prop The proportions of cases to drop, each strictly between
#'   0 and 1, in any order. Default `seq(0.1, 0.9, by = 0.1)`.
#'   A proportion that drops no case of a small sample (`n * p < 1`) is
#'   not run (without a warning; its metrics are `NA`) and counts as not
#'   passing. For the `summary()` method of a group result: the one
#'   proportion to report, which must be one of those the result was built
#'   with; a value outside them raises an error listing the available ones.
#'   Default `NULL` there, the middle proportion of the grid (the lower of
#'   the two middle ones for an even number).
#' @param threshold The lowest correlation between the original and the
#'   re-estimated edge weights for an iteration to count as stable, between
#'   0 and 1. Default `0.7`.
#' @param certainty The share of iterations that must reach `threshold` for
#'   a drop proportion to pass, between 0 and 1. Default `0.95`. It is the
#'   share of all iterations: one in which the correlation is undefined (the
#'   edge vector is constant in that subset), or in which the network could
#'   not be re-estimated, counts as not reaching `threshold`.
#' @param method The correlation: `"spearman"` (default) compares ranks, so
#'   it is not affected by the scale of the weights, which often span
#'   several orders of magnitude; `"pearson"` compares the weight
#'   magnitudes; `"kendall"` is Kendall's tau.
#' @param include_diag Whether self-loops are part of the edge vector.
#'   Default `FALSE`.
#' @param seed A number that makes the subsets reproducible. The session's
#'   own random stream is left as it was. For a `netobject_group` or `mcml`, each network
#'   is assessed from the same seed. Default `NULL`.
#' @param target The cells compared. `"existing"` compares only the cells
#'   that are non-zero in the original network; `"all"` compares every cell,
#'   zeros included, so edges that appear once cases are dropped count.
#'   `"auto"` (default) uses `"existing"` for unscaled and unthresholded
#'   transition networks, where an absent transition cannot appear in a
#'   subset, and `"all"` otherwise. `include_diag` applies to either.
#'
#' @return A `net_casedrop_reliability` object. `print()` shows the
#'   CS-coefficient and the mean and standard deviation of each metric at
#'   each drop proportion; [summary()] returns that table in full, one row
#'   per metric and drop proportion; `plot()` draws the metrics across the
#'   drop proportions. A `netobject_group` or `mcml` gives a
#'   `net_casedrop_reliability_group`, a named list of one result per
#'   network: its `print()` gives each network's number of cases, number of
#'   edges and CS-coefficient, its `summary()` one row per network at one
#'   drop proportion, and its `plot()` one curve per network. The fields of
#'   a single result:
#' \describe{
#'   \item{`cs`}{The CS-coefficient: the largest drop proportion at which
#'     the correlation reaches `threshold` in at least `certainty` of the
#'     iterations. `0` when no proportion passes; `NA` when the original
#'     edge weights do not vary (constant up to rounding error, relative to
#'     their size), so the correlation is undefined.}
#'   \item{`summary`}{A data frame with one row per metric and drop
#'     proportion: `metric` (`"mean_abs_dev"`, `"median_abs_dev"`,
#'     `"correlation"`, `"max_abs_dev"`), `drop_prop`, and the `mean`, `sd`,
#'     `median`, `mad`, `q025` and `q975` of the metric over the
#'     iterations.}
#'   \item{`metrics`}{The per-iteration values: a named list of four
#'     matrices, one per metric, with one row per iteration and one column
#'     per drop proportion.}
#'   \item{`correlations`}{The same matrix as the `correlation` entry of
#'     `metrics`.}
#'   \item{`drop_prop`, `threshold`, `certainty`, `iter`, `method`,
#'     `include_diag`}{The settings used.}
#'   \item{`n_cases`}{The number of cases resampled from: sequences for a
#'     transition network, rows of the data otherwise.}
#'   \item{`n_edges`}{The number of cells compared, set by `target` and
#'     `include_diag`.}
#' }
#' The resolved `target` is kept as the result's `target` attribute.
#'
#' When the original edge weights do not vary, a warning is given and the
#' result has `cs = NA`, a `summary` whose statistics are `NA`, and metric
#' matrices of `NA`.
#'
#' @details
#' For each drop proportion \eqn{p} and each iteration, a subset of
#' \eqn{n - \lfloor np \rfloor}{n - floor(np)} of the \eqn{n} cases is
#' drawn **without replacement**, the network is re-estimated with the
#' input's method, scaling and threshold, and the compared cells of the new
#' weight matrix are set against the same cells of the original. For
#' transition networks the per-sequence transition counts are computed once
#' and summed over each subset; other networks are re-estimated from the
#' subset of rows. A re-estimation that fails, or that loses a node (a
#' variable constant in the subset), leaves the iteration without values;
#' such iterations are counted in a warning of class
#' `nestimate_resample_failed` and count as below `threshold`.
#'
#' With `target = "auto"`, unscaled and unthresholded transition networks
#' compare only the existing edges: an absent transition cannot appear in a
#' subset, and the fixed zeros would inflate the correlation. Association
#' networks and scaled or thresholded networks compare all cells, because a
#' zero original weight can become non-zero once cases are dropped; set
#' `target = "existing"` to compare only the original edges there.
#'
#' Case-dropping drops cases at random, so it treats them as independent
#' units. When the cases are nested in actors (several sessions per
#' person), the remaining cases of an actor stand in for the dropped ones,
#' and the CS-coefficient overstates the stability; drop whole actors
#' instead (for example with `tna::estimate_cs(actor = )`).
#'
#' @seealso [centrality_stability()], [bootstrap_network()].
#'
#' @examples
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE),
#'   V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' es  <- casedrop_reliability(net, iter = 50, drop_prop = c(0.1, 0.3, 0.5),
#'                       seed = 1)
#' print(es)
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods} 50(1), 195-212.
#' \doi{10.3758/s13428-017-0862-1}
#'
#' @export
casedrop_reliability <- function(x,
                           iter        = 1000L,
                           drop_prop   = seq(0.1, 0.9, by = 0.1),
                           threshold   = 0.7,
                           certainty   = 0.95,
                           method      = c("spearman", "pearson", "kendall"),
                           include_diag = FALSE,
                           seed        = NULL,
                           target      = c("auto", "existing", "all")) {
  target <- match.arg(target)

  # ---- Dispatch for mcml / group inputs ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "netobject_group")) {
    .stop_single_unit_group(x)
    out <- lapply(x, function(net) {
      casedrop_reliability(net, iter = iter, drop_prop = drop_prop,
                     threshold = threshold, certainty = certainty,
                     method = method, include_diag = include_diag,
                     seed = seed, target = target)
    })
    class(out) <- c("net_casedrop_reliability_group", "list")
    return(out)
  }
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
  }
  if (!inherits(x, "netobject")) {
    stop(errorCondition(
      paste0(
        "'x' must be a netobject from build_network() (or an mcml / ",
        "netobject_group that contains one)."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.null(x$data)) {
    stop(errorCondition(
      "netobject does not contain $data. Rebuild with build_network().",
      class = "nestimate_bad_input", call = NULL))
  }
  .stop_single_unit(x)

  # ---- Input validation ----
  .assert(
    is.numeric(iter), length(iter) == 1L, iter >= 2,
    is.numeric(drop_prop), all(drop_prop > 0), all(drop_prop < 1),
    is.numeric(threshold), length(threshold) == 1L,
    threshold >= 0, threshold <= 1,
    is.numeric(certainty), length(certainty) == 1L,
    certainty >= 0, certainty <= 1,
    is.logical(include_diag), length(include_diag) == 1L
  )
  iter <- as.integer(iter)
  method <- match.arg(method)
  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1L)
    .local_seed(seed)
  }

  # ---- Resolve estimator + reuse bootstrap precompute for transition methods ----
  net_method <- .resampling_method(.resolve_method_alias(x$method), x$params)
  states     <- x$nodes$label
  n_states   <- length(states)
  is_transition <- net_method %in% c("relative", "frequency", "co_occurrence")
  is_relative   <- .counts_are_relative(net_method, x$params)
  scaling <- x$scaling
  thresh  <- x$threshold

  # Original edge vector (flattened)
  idx_mask <- if (include_diag) {
    matrix(TRUE, n_states, n_states)
  } else {
    diag_mask <- matrix(TRUE, n_states, n_states)
    diag(diag_mask) <- FALSE
    diag_mask
  }
  if (target == "auto") {
    target <- if (is_transition && is.null(scaling) && isTRUE(thresh == 0)) {
      "existing"
    } else {
      "all"
    }
  }
  idx_mask <- idx_mask & !is.na(x$weights)
  if (target == "existing") idx_mask <- idx_mask & x$weights != 0
  orig_edges <- as.vector(x$weights[idx_mask])
  if (.is_constant_vector(orig_edges)) {
    warning("Original edge vector has zero variance; CS undefined.",
            call. = FALSE)
    empty <- matrix(NA_real_, iter, length(drop_prop),
                    dimnames = list(NULL, paste0("p", drop_prop)))
    metrics <- list(mean_abs_dev = empty, median_abs_dev = empty,
                    correlation = empty,  max_abs_dev = empty)
    result <- list(
      cs            = NA_real_,
      # Same columns as a defined result, so print/summary/plot still work
      summary       = .casedrop_summary(metrics, drop_prop),
      metrics       = metrics,
      correlations  = empty,
      drop_prop     = drop_prop,
      threshold     = threshold,
      certainty     = certainty,
      iter          = iter,
      method        = method,
      include_diag  = include_diag,
      n_cases       = nrow(x$data),
      n_edges       = length(orig_edges)
    )
    class(result) <- "net_casedrop_reliability"
    attr(result, "target") <- target
    return(result)
  }

  # ---- Build-matrix closure: transition fast path or association slow path ----
  if (is_transition) {
    transition_data <- .resampling_transition_data(x$data, x$metadata, x$params)
    trans_2d <- .precompute_per_sequence(transition_data, net_method, x$params,
                                         states)
    n_cases  <- nrow(trans_2d)

    build_matrix <- function(idx) {
      counts <- colSums(trans_2d[idx, , drop = FALSE])
      mat <- matrix(counts, n_states, n_states, byrow = TRUE)
      if (is_relative) {
        rs <- rowSums(mat)
        nz <- rs > 0
        mat[nz, ] <- mat[nz, ] / rs[nz]
      }
      if (!is.null(scaling)) {
        mat <- .apply_scaling(mat, scaling,
                              include_zeros = .scaling_includes_zeros(net_method))
      }
      if (thresh > 0) mat[abs(mat) < thresh] <- 0
      mat
    }
  } else {
    data      <- x$data
    n_cases   <- nrow(data)
    estimator <- get_estimator(net_method)
    params    <- x$params
    level     <- x$level
    id_col    <- .param_get(x$params, "id") %||%
      .param_get(x$params, "id_col")

    build_matrix <- function(idx) {
      sub <- data[idx, , drop = FALSE]
      if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
        sub <- tryCatch(
          .decompose_multilevel(sub, id_col = id_col, level = level),
          error = function(e) NULL
        )
        if (is.null(sub)) return(NULL)
      }
      est <- tryCatch(
        do.call(estimator$fn, c(list(data = sub), params)),
        error = function(e) NULL
      )
      if (is.null(est)) return(NULL)
      # a state lost to the subset (a column constant in it) leaves no network
      # to compare: an unbuilt draw, which counts as below the threshold
      if (!all(states %in% rownames(est$matrix))) return(NULL)
      mat <- est$matrix[states, states]
      if (!is.null(scaling)) {
        mat <- .apply_scaling(mat, scaling,
                              include_zeros = .scaling_includes_zeros(net_method))
      }
      if (thresh > 0) mat[abs(mat) < thresh] <- 0
      mat
    }
  }

  # ---- Main case-dropping loop ----
  n_prop <- length(drop_prop)
  # Per-iteration, per-drop_prop storage for the four model-level metrics
  metric_names <- c("mean_abs_dev", "median_abs_dev", "correlation",
                    "max_abs_dev")
  metrics <- lapply(metric_names, function(nm) {
    matrix(NA_real_, nrow = iter, ncol = n_prop,
           dimnames = list(NULL, paste0("p", drop_prop)))
  })
  names(metrics) <- metric_names

  case_seq <- seq_len(n_cases)
  n_failed <- 0L
  n_tried <- 0L
  for (p_idx in seq_len(n_prop)) {
    n_drop <- floor(n_cases * drop_prop[p_idx])
    if (n_drop == 0L || n_drop >= n_cases) next
    keep_n <- n_cases - n_drop

    cor_metric <- match.arg(method, c("pearson", "spearman", "kendall"))
    for (it in seq_len(iter)) {
      idx <- sample(case_seq, keep_n, replace = FALSE)
      mat <- build_matrix(idx)
      n_tried <- n_tried + 1L
      if (is.null(mat)) {
        n_failed <- n_failed + 1L
        next
      }
      sub_edges <- as.vector(mat[idx_mask])
      sim <- .network_similarity(orig_edges, sub_edges,
                                 metrics = c("mean_abs_diff",
                                             "median_abs_diff",
                                             cor_metric,
                                             "max_abs_diff"))
      metrics$mean_abs_dev  [it, p_idx] <- sim[["mean_abs_diff"]]
      metrics$median_abs_dev[it, p_idx] <- sim[["median_abs_diff"]]
      metrics$correlation   [it, p_idx] <- sim[[cor_metric]]
      metrics$max_abs_dev   [it, p_idx] <- sim[["max_abs_diff"]]
    }
  }

  .warn_failed_resamples(n_failed, n_tried, "casedrop_reliability()",
                         below_threshold = TRUE)

  # ---- Per-drop-proportion summary across iterations ----
  summary_df <- .casedrop_summary(metrics, drop_prop)

  # ---- CS-coefficient: max drop_prop where corr certainty-fraction >= threshold ----
  cs <- .calculate_cs(metrics$correlation, threshold, certainty, drop_prop)

  structure(
    list(
      cs            = cs,
      summary       = summary_df,
      metrics       = metrics,
      correlations  = metrics$correlation,
      drop_prop     = drop_prop,
      threshold     = threshold,
      certainty     = certainty,
      iter          = iter,
      method        = method,
      include_diag  = include_diag,
      n_cases       = n_cases,
      n_edges       = length(orig_edges)
    ),
    class = "net_casedrop_reliability", target = target
  )
}

# One row per metric and drop proportion. Failed resamples are NA in
# `metrics` and are left out; a proportion without any value (for example
# when the CS-coefficient is undefined) gets NA statistics.
#' @noRd
.casedrop_summary <- function(metrics, drop_prop) {
  rows <- lapply(names(metrics), function(nm) {
    m <- metrics[[nm]]
    means <- colMeans(m, na.rm = TRUE)
    means[is.nan(means)] <- NA_real_
    data.frame(
      metric    = nm,
      drop_prop = drop_prop,
      mean      = means,
      sd        = apply(m, 2L, stats::sd,       na.rm = TRUE),
      median    = apply(m, 2L, stats::median,   na.rm = TRUE),
      mad       = apply(m, 2L, stats::mad,      na.rm = TRUE),
      q025      = apply(m, 2L, stats::quantile,
                        probs = 0.025, na.rm = TRUE),
      q975      = apply(m, 2L, stats::quantile,
                        probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

#' @rdname casedrop_reliability
#' @param digits For the `print()` method of a single result: the decimal
#'   places of the metric table. Default `3`.
#' @param ... For the `print()`, `summary()` and `plot()` methods: ignored.
#' @param object For the `summary()` methods: the result of
#'   `casedrop_reliability()`.
#' @param combined For the `plot()` method of a single result: `TRUE`
#'   (default) draws the four metrics as panels of one plot; `FALSE` returns
#'   a named list of four single-panel plots, one per metric.
#' @param metric For the `plot()` method of a group result: the metric
#'   drawn, `"correlation"` (default), `"mean_abs_dev"`, `"median_abs_dev"`
#'   or `"max_abs_dev"`.
#' @return The `print()` methods return their input invisibly. The
#'   `summary()` method of a single result returns its `summary` field,
#'   described above. The `summary()` method of
#'   a group result returns a `summary.net_casedrop_reliability_group` data
#'   frame with one row per network (named in the row names): `n_edges` and
#'   the four metrics (`cor`, `mean_abs_dev`, `median_abs_dev`,
#'   `max_abs_dev`) at the chosen drop proportion, each written as
#'   `"mean +/- sd"`; its `print()` names that proportion. The `plot()`
#'   methods return a `ggplot`; with `combined = FALSE`, the named list of
#'   plots is returned invisibly. A single result is drawn as one panel per
#'   metric, the mean across iterations as a line with a band of plus and
#'   minus one standard deviation; a group result draws the chosen metric
#'   with one line and band per network. On the correlation, a dashed line
#'   marks `threshold`. A result whose CS-coefficient is undefined is drawn as a
#'   panel that says so, and in a group such networks are left out and
#'   named in the caption.
#' @export
print.net_casedrop_reliability <- function(x, digits = 3, ...) {
  cat(sprintf("Edge-weight Case-dropping Stability\n"))
  cat(sprintf("  Cases (rows of $data) : %d\n", x$n_cases))
  cat(sprintf("  Edges assessed        : %d%s\n", x$n_edges,
              if (!x$include_diag) " (diagonal excluded)" else ""))
  cat(sprintf("  Iterations / prop     : %d\n", x$iter))
  cat(sprintf("  Correlation method    : %s\n", x$method))
  cat(sprintf("  CS-coefficient (r)    : %.2f  (threshold=%.2f, certainty=%.2f)\n",
              x$cs, x$threshold, x$certainty))
  cat("\nModel-level reliability across iterations (mean +/- sd per drop):\n")
  # Pivot summary into a compact matrix for display
  props  <- x$drop_prop
  fmt <- function(v) formatC(v, digits = digits, format = "f", flag = " ")
  show <- function(metric_name, display) {
    rows <- x$summary[x$summary$metric == metric_name, ]
    cat(sprintf("  %-14s ", display))
    cat(paste(sprintf("%s+-%s",
                      fmt(rows$mean), fmt(rows$sd)),
              collapse = "  "))
    cat("\n")
  }
  cat(sprintf("  %-14s %s\n", "drop_prop",
              paste(sprintf("%-11s", sprintf("p=%.1f", props)),
                    collapse = "  ")))
  show("mean_abs_dev",   "mean|diff|")
  show("median_abs_dev", "MAD")
  show("correlation",    "cor")
  show("max_abs_dev",    "max|diff|")
  invisible(x)
}

#' @rdname casedrop_reliability
#' @export
summary.net_casedrop_reliability <- function(object, ...) {
  object$summary
}

#' @rdname casedrop_reliability
#' @export
print.net_casedrop_reliability_group <- function(x, ...) {
  cs        <- vapply(x, function(e) e$cs,       numeric(1))
  n_edges   <- vapply(x, function(e) e$n_edges,  integer(1))
  n_cases   <- vapply(x, function(e) e$n_cases,  integer(1))
  out <- data.frame(
    n_cases = n_cases,
    n_edges = n_edges,
    CS      = round(cs, 2),
    row.names = names(x)
  )
  cat(sprintf("Edge-weight Case-dropping Stability (%d networks, threshold = %.2f)\n",
              length(x), x[[1]]$threshold))
  print(out)
  invisible(x)
}

#' @rdname casedrop_reliability
#' @export
summary.net_casedrop_reliability_group <- function(object,
                                                   drop_prop = NULL, ...) {
  grid <- sort(object[[1L]]$drop_prop)
  if (is.null(drop_prop)) {
    # Default to the object's actual middle grid point (no assumed 0.7).
    drop_prop <- grid[ceiling(length(grid) / 2)]
  }
  .assert(is.numeric(drop_prop), length(drop_prop) == 1L)
  if (!any(abs(grid - drop_prop) < 1e-9)) {
    stop(errorCondition(
      paste0(
        "drop_prop = ",
        drop_prop,
        " is not in the object's grid. ",
        "Available proportions: ",
        paste(grid, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  fmt <- function(m, s) sprintf("%.3f +/- %.3f", m, s)
  one_row <- function(e) {
    s <- e$summary
    r <- s[abs(s$drop_prop - drop_prop) < 1e-9, , drop = FALSE]
    stats::setNames(
      list(fmt(r$mean[r$metric == "correlation"],    r$sd[r$metric == "correlation"]),
           fmt(r$mean[r$metric == "mean_abs_dev"],   r$sd[r$metric == "mean_abs_dev"]),
           fmt(r$mean[r$metric == "median_abs_dev"], r$sd[r$metric == "median_abs_dev"]),
           fmt(r$mean[r$metric == "max_abs_dev"],    r$sd[r$metric == "max_abs_dev"])),
      c("cor", "mean_abs_dev", "median_abs_dev", "max_abs_dev")
    )
  }
  rows <- lapply(object, one_row)
  out  <- do.call(rbind.data.frame, rows)
  out  <- cbind(
    n_edges = vapply(object, function(e) e$n_edges, integer(1)),
    out
  )
  rownames(out) <- names(object)
  structure(out, class = c("summary.net_casedrop_reliability_group", "data.frame"),
            drop_prop = drop_prop)
}

#' @rdname casedrop_reliability
#' @export
print.summary.net_casedrop_reliability_group <- function(x, ...) {
  cat(sprintf("Edge-weight reliability at drop = %.2f (mean +/- sd over iterations)\n",
              attr(x, "drop_prop")))
  print.data.frame(x)
  invisible(x)
}

#' @rdname casedrop_reliability
#' @export
plot.net_casedrop_reliability <- function(x, combined = TRUE, ...) {
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))
  df <- x$summary
  df$metric <- factor(df$metric,
                      levels = c("correlation", "mean_abs_dev",
                                 "median_abs_dev", "max_abs_dev"),
                      labels = c("Correlation",
                                 "Mean |diff|",
                                 "Median |diff| (MAD)",
                                 "Max |diff|"))

  base_p <- function(d, ttl, sub_title, hline_df) {
    p <- ggplot2::ggplot(d,
                          ggplot2::aes(x = .data$drop_prop,
                                       y = .data$mean)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$mean - .data$sd,
                     ymax = .data$mean + .data$sd),
        fill = "#2196F3", alpha = 0.2
      ) +
      ggplot2::geom_line(color = "#0D47A1", linewidth = 0.8) +
      ggplot2::geom_point(color = "#0D47A1", size = 1.5) +
      ggplot2::labs(
        x = "Proportion of cases dropped",
        y = "Mean \u00b1 SD across iterations",
        title = ttl,
        subtitle = sub_title
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(color = "grey40"),
        strip.text    = ggplot2::element_text(face = "bold")
      )
    if (NROW(hline_df) > 0L) {
      p <- p + ggplot2::geom_hline(
        data = hline_df,
        ggplot2::aes(yintercept = .data$yint),
        linetype = "dashed", color = "grey40"
      )
    }
    p
  }

  ttl <- sprintf("Edge-weight reliability (CS = %.2f, %d iter, %d edges)",
                 x$cs, x$iter, x$n_edges)
  sub_title <- sprintf("Case-dropping on %d cases, %s correlation",
                       x$n_cases, x$method)
  if (!any(is.finite(df$mean))) {
    # Nothing to draw: say why instead of an empty or failing panel
    return(.casedrop_undefined_plot(ttl))
  }

  if (!combined) {
    levs <- levels(df$metric)
    plots <- lapply(levs, function(lv) {
      sub <- df[df$metric == lv, , drop = FALSE]
      hline <- if (lv == "Correlation")
        data.frame(yint = x$threshold) else data.frame()
      base_p(sub, sprintf("%s -- %s", ttl, lv), sub_title, hline)
    })
    names(plots) <- levs
    return(invisible(plots))
  }

  hline_df <- data.frame(
    metric = factor("Correlation", levels = levels(df$metric)),
    yint = x$threshold
  )
  base_p(df, ttl, sub_title, hline_df) +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y")
}

#' @rdname casedrop_reliability
#' @export
plot.net_casedrop_reliability_group <- function(x,
                                          metric = c("correlation",
                                                     "mean_abs_dev",
                                                     "median_abs_dev",
                                                     "max_abs_dev"),
                                          ...) {
  metric <- match.arg(metric)

  rows <- lapply(names(x), function(nm) {
    s <- x[[nm]]$summary
    s <- s[s$metric == metric, , drop = FALSE]
    s$network <- nm
    s$cs <- x[[nm]]$cs
    s
  })
  df <- do.call(rbind, rows)
  df$network <- factor(df$network, levels = names(x))
  # Networks whose CS is undefined have no values to draw
  undefined <- names(x)[!vapply(
    names(x),
    function(nm) any(is.finite(df$mean[df$network == nm])),
    logical(1L)
  )]
  df <- df[is.finite(df$mean), , drop = FALSE]
  if (!nrow(df)) {
    return(.casedrop_undefined_plot("Edge-weight case-dropping reliability"))
  }
  undefined_note <- if (length(undefined)) {
    paste0("Not shown (edge weights do not vary, CS undefined): ",
           paste(undefined, collapse = ", "))
  }

  # Use a discrete palette, recycling if there are many networks
  palette <- grDevices::hcl.colors(length(levels(df$network)), "Dynamic")

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$drop_prop, y = .data$mean,
                                    colour = .data$network,
                                    fill   = .data$network)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$mean - .data$sd,
                   ymax = .data$mean + .data$sd),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::labs(
      x = "Proportion of cases dropped",
      y = switch(metric,
                 correlation    = "Edge-vector correlation",
                 mean_abs_dev   = "Mean |diff|",
                 median_abs_dev = "Median |diff| (MAD)",
                 max_abs_dev    = "Max |diff|"),
      colour = "Network", fill = "Network",
      title = "Edge-weight case-dropping reliability",
      subtitle = sprintf("Metric: %s; mean \u00b1 sd across iterations",
                         metric),
      caption = undefined_note
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40"),
      legend.position = "right"
    )

  if (metric == "correlation") {
    threshold <- x[[1]]$threshold
    p <- p + ggplot2::geom_hline(yintercept = threshold,
                                 linetype = "dashed", color = "grey40")
  }

  p
}

# A panel that states why there is nothing to draw
#' @noRd
.casedrop_undefined_plot <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0, size = 4,
      label = paste0("CS undefined: the edge weights of the original ",
                     "network do not vary")
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
