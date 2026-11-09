# ---- Split-Half Reliability ----

#' Split-Half Reliability for Network Estimates
#'
#' @description
#' `network_reliability()` measures how much a network's weights depend on
#' the particular sample: it splits the sequences (for an association
#' network, the rows of its data) at random into two halves, estimates the
#' network on each half with the original network's method and settings,
#' compares the two weight matrices, and repeats this `iter` times. Each
#' comparison gives four numbers: the mean, median and largest absolute
#' difference between the halves' weights and their Pearson correlation,
#' taken over all cells of the weight matrix. It takes
#' one network or several, so the reliability of different methods or
#' groups can be set side by side; a `netobject_group` contributes one
#' network per group.
#'
#' Each half holds only part of the sequences, so the agreement is that of
#' networks estimated from a sample of that size. Differences are on the
#' scale of the weights: counts of a `"frequency"` network and
#' probabilities of a `"relative"` one are not comparable without `scale`.
#' To see how far single edges move, use [bootstrap_network()]; to check
#' centrality rankings, [centrality_stability()].
#'
#' @details
#' The first half holds `round(split * n)` of the `n` sequences (at least
#' one), the second the rest; a `split` that leaves the second half empty
#' (few sequences, `split` far from 0.5) raises an error of class
#' `nestimate_bad_input`. For the transition methods (`"relative"`,
#' `"frequency"`, `"co_occurrence"`, and one-hot `wtna` networks, which are
#' counted the same way), each sequence's transition counts are computed
#' once and each half is summed from them, as in [bootstrap_network()].
#' Other methods are re-estimated on each half of the rows. The network's
#' own scaling and threshold are applied to both halves. A half on which
#' the estimator fails, or that loses a variable (a column constant in
#' that half), gives `NA` for that iteration, counted in a warning of class
#' `nestimate_resample_failed`. The summary leaves out missing values, such
#' as a correlation that is undefined because one half's weights are all
#' equal.
#'
#' @param ... The networks: one or more `netobject`s from
#'   [build_network()]. A `netobject_group` is split into its networks,
#'   named by group; an `mcml` is converted with
#'   \code{\link[=as_tna]{as_tna()}} first, and a bare `cograph_network`
#'   is read as a `netobject`. Each network is named by its argument name
#'   or, failing that, its method; repeated names get a suffix (`"_1"`,
#'   `"_2"`, ...). A network without its sequence data (which
#'   [build_network()] keeps) raises an error, a
#'   network resting on one sequence an error of class
#'   `nestimate_single_sequence`, and networks of different methods with
#'   `scale = "none"` a warning. For the `print()`, `summary()` and
#'   `plot()` methods: ignored.
#' @param iter The number of split-half iterations, at least 2. Default
#'   `1000`.
#' @param split The share of the sequences that goes to the first half, the
#'   rest forming the second, strictly between 0 and 1. Default `0.5`.
#' @param scale The rescaling applied to both halves' weight matrices
#'   before they are compared, so that networks on different scales can be
#'   set side by side: `"none"` (default); `"minmax"`, which maps each
#'   matrix to 0--1; `"standardize"`, which subtracts the mean of its cells
#'   and divides by their SD; or `"proportion"`, which divides by the sum
#'   of its cells. A matrix whose cells are all equal, or sum to 0, is left
#'   as it is.
#' @param seed A seed for the random splits, so the result repeats exactly.
#'   The session's own random stream is left as it was. Default `NULL`.
#'
#' @return A `net_reliability` object. `print()` shows the mean and SD of
#'   each metric over the iterations, per network; [summary()] returns
#'   them as a data frame with one row per network and metric and the
#'   columns `model`, `metric` (`mean_dev`, `median_dev`, `cor`,
#'   `max_dev`), `mean` and `sd`; `plot()` draws their distributions. The
#'   object's fields:
#' \describe{
#'   \item{iterations}{A data frame with one row per iteration and
#'     network, and the columns `model`, `mean_dev`, `median_dev`, `cor`
#'     and `max_dev`.}
#'   \item{summary}{The table `summary()` returns.}
#'   \item{models}{The networks compared, as a named list.}
#'   \item{iter, split, scale}{The settings of the call.}
#' }
#' `print()` returns its input invisibly. `plot()` draws the histograms
#' (see *Methods*) and returns the `ggplot` invisibly; with
#' `combined = FALSE` it returns, without drawing, a list of four
#' `ggplot`s named by the metric labels (`"Max Abs. Diff."`,
#' `"Mean Abs. Diff."`, `"Median Abs. Diff."`, `"Pearson"`).
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' rel <- network_reliability(net, iter = 10)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' rel <- network_reliability(net, iter = 100, seed = 42)
#' print(rel)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_network}}
#'
#' @importFrom stats cor median sd
#' @export
network_reliability <- function(..., iter = 1000L, split = 0.5,
                        scale = "none", seed = NULL) {

  dots <- list(...)

  # ---- Input validation ----
  .assert(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(split), length(split) == 1, split > 0, split < 1,
    is.character(scale), length(scale) == 1
  )
  iter <- as.integer(iter)
  scale <- match.arg(scale, c("none", "minmax", "standardize", "proportion"))

  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
  }

  # ---- Flatten netobject_group elements ----
  objs <- list()
  labels <- character(0)
  for (i in seq_along(dots)) {
    obj <- dots[[i]]
    nm <- names(dots)[i]
    if (inherits(obj, "mcml")) obj <- as_tna(obj)
    if (inherits(obj, "cograph_network")) obj <- .as_netobject(obj)
    if (inherits(obj, "netobject_group")) {
      for (g in names(obj)) {
        objs <- c(objs, list(obj[[g]]))
        labels <- c(labels, g)
      }
    } else if (inherits(obj, "netobject")) {
      label <- if (!is.null(nm) && nzchar(nm)) nm else obj$method
      objs <- c(objs, list(obj))
      labels <- c(labels, label)
    } else {
      stop(errorCondition(
        "All arguments must be netobject or netobject_group objects.",
        class = "nestimate_bad_input", call = NULL))
    }
  }

  if (length(objs) == 0L) {
    stop(errorCondition("At least one netobject is required.",
                        class = "nestimate_bad_input", call = NULL))
  }
  invisible(Map(.stop_single_unit, objs, label = labels))

  # Deduplicate names
  labels <- make.unique(labels, sep = "_")
  model_list <- stats::setNames(objs, labels)

  # ---- Warn if different methods without scaling ----
  methods <- vapply(model_list, function(m) m$method, character(1))
  if (length(unique(methods)) > 1L && scale == "none") {
    warning(
      "Models use different methods (",
      paste(unique(methods), collapse = ", "),
      "). Consider setting scale = 'minmax' for comparable results.",
      call. = FALSE
    )
  }

  # ---- Run split-half per model ----
  all_iters <- lapply(names(model_list), function(model_name) {
    net <- model_list[[model_name]]

    if (is.null(net$data)) {
      stop(errorCondition(paste0("netobject '", model_name,
           "' does not contain $data. Rebuild with build_network()."),
           class = "nestimate_bad_input", call = NULL))
    }

    method <- .resampling_method(.resolve_method_alias(net$method), net$params)
    is_transition <- method %in% c("relative", "frequency", "co_occurrence")

    if (is_transition) {
      .reliability_transition(net, method, iter, split, scale)
    } else {
      .reliability_association(net, method, iter, split, scale)
    }
  })
  names(all_iters) <- names(model_list)

  # ---- Assemble iterations data frame ----
  iterations <- do.call(rbind, lapply(names(all_iters), function(nm) {
    df <- all_iters[[nm]]
    df$model <- nm
    df[, c("model", "mean_dev", "median_dev", "cor", "max_dev")]
  }))
  rownames(iterations) <- NULL

  # ---- Compute summary ----
  metric_names <- c("mean_dev", "median_dev", "cor", "max_dev")
  summary_rows <- do.call(rbind, lapply(names(model_list), function(nm) {
    sub <- iterations[iterations$model == nm, , drop = FALSE]
    do.call(rbind, lapply(metric_names, function(met) {
      vals <- sub[[met]]
      data.frame(
        model = nm,
        metric = met,
        mean = mean(vals, na.rm = TRUE),
        sd = sd(vals, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary_rows) <- NULL

  result <- list(
    iterations = iterations,
    summary = summary_rows,
    models = model_list,
    iter = iter,
    split = split,
    scale = scale
  )
  class(result) <- "net_reliability"
  result
}


# ---- Transition fast path ----

#' Split-half reliability for transition methods
#' @noRd
.reliability_transition <- function(net, method, iter, split, scale) {
  states <- net$nodes$label
  n_states <- length(states)
  nbins <- n_states * n_states
  is_relative <- .counts_are_relative(method, net$params)
  scaling <- net$scaling
  threshold <- net$threshold

  # Pre-compute per-sequence count matrix (reuse bootstrap infrastructure)
  transition_data <- .resampling_transition_data(net$data, net$metadata,
                                                 net$params)
  trans_2d <- .precompute_per_sequence(transition_data, method, net$params,
                                       states)
  n_seq <- nrow(trans_2d)
  n_half <- max(1L, round(n_seq * split))
  .reliability_check_halves(n_half, n_seq, split)

  # Post-process raw counts into the appropriate matrix form
  postprocess <- function(counts) {
    mat <- matrix(counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) { # nocov start
      mat <- .apply_scaling(mat, scaling,
                            include_zeros = .scaling_includes_zeros(method))
    } # nocov end
    if (threshold > 0) mat[abs(mat) < threshold] <- 0 # nocov
    mat
  }

  # Run iterations
  results <- vapply(seq_len(iter), function(i) {
    idx_a <- sample.int(n_seq, n_half, replace = FALSE)
    idx_b <- setdiff(seq_len(n_seq), idx_a)

    counts_a <- colSums(trans_2d[idx_a, , drop = FALSE])
    counts_b <- colSums(trans_2d[idx_b, , drop = FALSE])

    mat_a <- postprocess(counts_a)
    mat_b <- postprocess(counts_b)

    # Apply reliability scaling
    mat_a <- .scale_matrix(mat_a, scale)
    mat_b <- .scale_matrix(mat_b, scale)

    .split_half_metrics(mat_a, mat_b)
  }, numeric(4))

  data.frame(
    mean_dev = results[1, ],
    median_dev = results[2, ],
    cor = results[3, ],
    max_dev = results[4, ],
    stringsAsFactors = FALSE
  )
}


# ---- Association path ----

#' Split-half reliability for association methods
#' @noRd
.reliability_association <- function(net, method, iter, split, scale) {
  data <- net$data
  states <- net$nodes$label
  n_states <- length(states)
  nbins <- n_states * n_states
  scaling <- net$scaling
  threshold <- net$threshold
  params <- net$params
  level <- net$level
  id_col <- .param_get(params, "id") %||% .param_get(params, "id_col")

  estimator <- get_estimator(method)
  n <- nrow(data)
  n_half <- max(1L, round(n * split))
  .reliability_check_halves(n_half, n, split)

  postprocess <- function(mat) {
    if (!is.null(scaling)) { # nocov start
      mat <- .apply_scaling(mat, scaling,
                            include_zeros = .scaling_includes_zeros(method))
    }
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    mat # nocov end
  }

  results <- vapply(seq_len(iter), function(i) {
    idx_a <- sample.int(n, n_half, replace = FALSE)
    idx_b <- setdiff(seq_len(n), idx_a)

    build_half <- function(sub_data) {
      if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
        sub_data <- tryCatch( # nocov start
          .decompose_multilevel(sub_data, id_col = id_col, level = level),
          error = function(e) NULL
        )
        if (is.null(sub_data)) return(NULL) # nocov end
      }
      tryCatch(
        do.call(estimator$fn, c(list(data = sub_data), params)),
        error = function(e) NULL
      )
    }

    est_a <- build_half(data[idx_a, , drop = FALSE])
    est_b <- build_half(data[idx_b, , drop = FALSE])

    if (is.null(est_a) || is.null(est_b)) return(rep(NA_real_, 4))
    # a state lost to one half (a column constant in it) leaves no pair to compare
    if (!all(states %in% rownames(est_a$matrix)) ||
        !all(states %in% rownames(est_b$matrix))) return(rep(NA_real_, 4))

    mat_a <- est_a$matrix[states, states] # nocov start
    mat_b <- est_b$matrix[states, states]

    mat_a <- postprocess(mat_a)
    mat_b <- postprocess(mat_b)

    mat_a <- .scale_matrix(mat_a, scale)
    mat_b <- .scale_matrix(mat_b, scale)

    .split_half_metrics(mat_a, mat_b) # nocov end
  }, numeric(4))
  .warn_failed_resamples(sum(is.na(results[1, ])), iter,
                         "network_reliability()")

  data.frame(
    mean_dev = results[1, ],
    median_dev = results[2, ],
    cor = results[3, ],
    max_dev = results[4, ],
    stringsAsFactors = FALSE
  )
}


# ---- Helpers ----

# Both halves need at least one unit: with few units a split far from 0.5
# leaves the second half empty, and comparing it with zeros means nothing.
.reliability_check_halves <- function(n_half, n, split) {
  if (n_half >= n) {
    stop(errorCondition(sprintf(paste0(
      "`split` = %g puts all %d units in the first half, so the second half ",
      "is empty. Use a split nearer 0.5 or more units."), split, n),
      class = "nestimate_bad_input", call = NULL))
  }
}

#' Compute split-half metrics between two matrices
#'
#' Thin wrapper around `.network_similarity()` (R/compare_network.R) that
#' returns the four metrics the reliability summary table reports, in the
#' fixed order `c(mean_dev, median_dev, cor, max_dev)`.
#' @noRd
.split_half_metrics <- function(mat_a, mat_b) {
  unname(.network_similarity(mat_a, mat_b,
                             metrics = c("mean_abs_diff", "median_abs_diff",
                                         "pearson", "max_abs_diff")))
}


#' Scale a matrix for cross-method comparison
#' @noRd
.scale_matrix <- function(mat, method) {
  switch(method,
    none = mat,
    minmax = {
      rng <- range(mat)
      if (rng[1] == rng[2]) mat
      else (mat - rng[1]) / (rng[2] - rng[1])
    },
    standardize = {
      s <- sd(as.vector(mat))
      if (s == 0) mat
      else (mat - mean(mat)) / s
    },
    proportion = {
      total <- sum(mat)
      if (total == 0) mat
      else mat / total
    }
  )
}


# ---- S3 Methods ----

#' @rdname network_reliability
#' @param x For the `print()` and `plot()` methods: a `net_reliability`
#'   object.
#' @param object For the `summary()` method: a `net_reliability` object.
#' @param bins For the `plot()` method: the number of histogram bins per
#'   panel. Default `60`.
#' @param combined For the `plot()` method: `TRUE` (default) draws the four
#'   metrics as panels of one plot; `FALSE` returns a named list of four
#'   single-panel plots, one per metric, without drawing them.
#' @export
print.net_reliability <- function(x, ...) {
  cat(sprintf("Split-Half Reliability (%d iterations, split = %.0f%%",
              x$iter, x$split * 100))
  if (x$scale != "none") {
    cat(sprintf(", scale = %s", x$scale))
  }
  cat(")\n")

  models <- unique(x$summary$model)
  metric_labels <- c(
    mean_dev = "Mean Abs. Diff.",
    median_dev = "Median Abs. Diff.",
    cor = "Pearson",
    max_dev = "Max Abs. Diff."
  )

  for (m in models) {
    if (length(models) > 1L) cat(sprintf("\n  %s:\n", m))
    sub <- x$summary[x$summary$model == m, , drop = FALSE]
    for (i in seq_len(nrow(sub))) {
      label <- metric_labels[sub$metric[i]]
      prefix <- if (length(models) > 1L) "    " else "  "
      cat(sprintf("%s%-18s  mean = %.4f  sd = %.4f\n",
                  prefix, label, sub$mean[i], sub$sd[i]))
    }
  }

  invisible(x)
}


#' @rdname network_reliability
#' @export
summary.net_reliability <- function(object, ...) {
  object$summary
}


#' @rdname network_reliability
#' @section Methods:
#' * `plot()`: histograms of the split-half metrics over the iterations, one
#'   panel per metric, the mean marked by a dashed line and a label. With
#'   several networks, their histograms are overlaid, one colour per
#'   network, with a legend below.
#' @export
plot.net_reliability <- function(x, bins = 60L, combined = TRUE, ...) {
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))
  iters <- x$iterations
  models <- unique(iters$model)
  multi <- length(models) > 1L

  metric_labels <- c(
    max_dev = "Max Abs. Diff.",
    mean_dev = "Mean Abs. Diff.",
    median_dev = "Median Abs. Diff.",
    cor = "Pearson"
  )

  metric_cols <- names(metric_labels)
  long <- do.call(rbind, lapply(metric_cols, function(met) {
    data.frame(
      model = iters$model,
      metric = metric_labels[met],
      value = iters[[met]],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  long$metric <- factor(long$metric, levels = metric_labels)

  means <- aggregate(value ~ model + metric, data = long, FUN = mean)
  means$label <- sprintf("%.2f", means$value)

  palette <- c("#4DA167", "#E89AB4", "#2B6CB0", "#DD8452",
               "#8172B2", "#937860", "#DA8BC3", "#8C8C8C")
  fill_vals <- palette[seq_along(models)]
  names(fill_vals) <- models

  base_p <- function(d, m_subset, ttl) {
    ggplot2::ggplot(d, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = .data$model),
        bins = bins, alpha = 0.65, position = "identity",
        color = NA
      ) +
      ggplot2::geom_vline(
        data = m_subset,
        ggplot2::aes(xintercept = .data$value, color = .data$model),
        linetype = "dashed", linewidth = 0.6,
        show.legend = FALSE
      ) +
      ggplot2::geom_label(
        data = m_subset,
        ggplot2::aes(x = .data$value, y = Inf, label = .data$label,
                     color = .data$model),
        vjust = 1.2, size = 3.3, label.size = 0.4,
        fill = "white", show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = fill_vals, name = "Model") +
      ggplot2::scale_color_manual(values = fill_vals, guide = "none") +
      ggplot2::labs(x = "Metric Value", y = "Frequency", title = ttl) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = if (multi) "bottom" else "none",
        strip.text = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  if (!combined) {
    levs <- levels(long$metric)
    plots <- lapply(levs, function(lv) {
      sub_long  <- long[long$metric == lv, , drop = FALSE]
      sub_means <- means[means$metric == lv, , drop = FALSE]
      base_p(sub_long, sub_means, lv)
    })
    names(plots) <- levs
    return(invisible(plots))
  }

  p <- base_p(long, means,
              if (multi) "Split-Half Reliability" else NULL) +
    ggplot2::facet_wrap(~ metric, scales = "free", ncol = 2L)

  print(p)
  invisible(p)
}
