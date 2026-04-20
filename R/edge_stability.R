# ---- Edge-weight case-dropping stability --------------------------------
# Mirrors centrality_stability() but the stability target is the vector of
# edge weights (flattened adjacency) rather than per-node centralities.
# Answers: "How many cases can be dropped before the full set of edge
# weights loses rank correlation with the original?"

#' Edge-weight Case-dropping Stability
#'
#' Computes a **CS-coefficient for the edge-weight vector** of a network:
#' the maximum proportion of cases (rows of `x$data`) that can be dropped
#' while the flattened edge-weight vector of the re-estimated network
#' still correlates with the original above `threshold` in at least
#' `certainty` of iterations.
#'
#' Complements [centrality_stability()]: that function asks whether
#' centrality *rankings* are stable; this one asks whether the *edge-weight
#' structure itself* is stable. For MCML-derived networks where each row
#' of `$data` is one transition, this is case-dropping of **edges**.
#'
#' @param x A `netobject`, `cograph_network`, `netobject_group`, or `mcml`.
#'   For the group types this function iterates over each element and
#'   returns a named list.
#' @param iter Integer. Iterations per drop proportion. Default `1000`.
#' @param drop_prop Numeric vector of proportions to evaluate. Each entry
#'   must lie strictly between 0 and 1. Default `seq(0.1, 0.9, by = 0.1)`.
#' @param threshold Numeric in `[0, 1]`. Minimum edge-vector correlation
#'   for an iteration to count as stable. Default `0.7`.
#' @param certainty Numeric in `[0, 1]`. Required fraction of iterations
#'   whose correlation must exceed `threshold` for a drop proportion to
#'   qualify. Default `0.95`.
#' @param method Correlation method: `"pearson"` (weight magnitudes),
#'   `"spearman"` (ranks, robust to scale), or `"kendall"`. Default
#'   `"spearman"` because edge weights often span several orders of
#'   magnitude and rank stability is the typical target.
#' @param include_diag Logical. Include diagonal (self-loop) edges in the
#'   edge vector. Default `FALSE`.
#' @param seed Optional integer for reproducibility.
#'
#' @return An object of class `net_edge_stability` with:
#' \describe{
#'   \item{`cs`}{Scalar CS-coefficient — the maximum drop proportion for
#'     which the edge-vector correlation remains >= `threshold` in at
#'     least `certainty` of iterations. Zero if no proportion qualifies.}
#'   \item{`correlations`}{`iter` x `length(drop_prop)` matrix of per-
#'     iteration correlations.}
#'   \item{`drop_prop`, `threshold`, `certainty`, `iter`, `method`}{Inputs.}
#' }
#'
#' @details
#' For each `drop_prop` p and each iteration, a size `n_cases * (1 - p)`
#' subset of `$data` rows is selected **without replacement**, the network
#' is re-estimated using the same method/scaling/threshold as the input,
#' and the upper/lower-triangle (directed: all off-diagonal entries) of
#' the new weight matrix is flattened and correlated with the
#' corresponding vector of the original matrix. The correlation method
#' defaults to Spearman for robustness to the wide dynamic range of
#' transition probabilities.
#'
#' Unlike bootstrap CIs, case-dropping does not estimate sampling variance
#' and so does not rely on the i.i.d. assumption. This makes it the
#' appropriate robustness check for **edgelist-derived** networks (where
#' rows of `$data` lack actor grouping), since dropping rows at random is
#' a well-posed operation regardless of within-actor correlation.
#'
#' @seealso [centrality_stability()], [bootstrap_network()].
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE),
#'   V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' es  <- edge_stability(net, iter = 50, drop_prop = c(0.1, 0.3, 0.5),
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
edge_stability <- function(x,
                           iter        = 1000L,
                           drop_prop   = seq(0.1, 0.9, by = 0.1),
                           threshold   = 0.7,
                           certainty   = 0.95,
                           method      = c("spearman", "pearson", "kendall"),
                           include_diag = FALSE,
                           seed        = NULL) {

  # ---- Dispatch for mcml / group inputs ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      edge_stability(net, iter = iter, drop_prop = drop_prop,
                     threshold = threshold, certainty = certainty,
                     method = method, include_diag = include_diag,
                     seed = seed)
    })
    class(out) <- c("net_edge_stability_group", "list")
    return(out)
  }
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
  }
  if (!inherits(x, "netobject")) {
    stop("'x' must be a netobject from build_network() (or an mcml / ",
         "netobject_group that contains one).", call. = FALSE)
  }
  if (is.null(x$data)) {
    stop("netobject does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  # ---- Input validation ----
  stopifnot(
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
    stopifnot(is.numeric(seed), length(seed) == 1L)
    set.seed(seed)
  }

  # ---- Resolve estimator + reuse bootstrap precompute for transition methods ----
  net_method <- .resolve_method_alias(x$method)
  states     <- x$nodes$label
  n_states   <- length(states)
  is_transition <- net_method %in% c("relative", "frequency", "co_occurrence")
  is_relative   <- net_method == "relative"
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
  orig_edges <- as.vector(x$weights[idx_mask])
  if (stats::sd(orig_edges, na.rm = TRUE) == 0) {
    warning("Original edge vector has zero variance; CS undefined.",
            call. = FALSE)
    result <- list(
      cs = 0,
      correlations = matrix(NA_real_, iter, length(drop_prop)),
      drop_prop = drop_prop, threshold = threshold,
      certainty = certainty, iter = iter, method = method
    )
    class(result) <- "net_edge_stability"
    return(result)
  }

  # ---- Build-matrix closure: transition fast path or association slow path ----
  if (is_transition) {
    trans_2d <- .precompute_per_sequence(x$data, net_method, x$params, states)
    n_cases  <- nrow(trans_2d)

    build_matrix <- function(idx) {
      counts <- colSums(trans_2d[idx, , drop = FALSE])
      mat <- matrix(counts, n_states, n_states, byrow = TRUE)
      if (is_relative) {
        rs <- rowSums(mat)
        nz <- rs > 0
        mat[nz, ] <- mat[nz, ] / rs[nz]
      }
      if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
      if (thresh > 0) mat[abs(mat) < thresh] <- 0
      mat
    }
  } else {
    data      <- x$data
    n_cases   <- nrow(data)
    estimator <- get_estimator(net_method)
    params    <- x$params
    level     <- x$level
    id_col    <- x$params$id %||% x$params$id_col

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
      mat <- est$matrix[states, states]
      if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
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
  for (p_idx in seq_len(n_prop)) {
    n_drop <- floor(n_cases * drop_prop[p_idx])
    if (n_drop == 0L || n_drop >= n_cases) next
    keep_n <- n_cases - n_drop

    for (it in seq_len(iter)) {
      idx <- sample(case_seq, keep_n, replace = FALSE)
      mat <- build_matrix(idx)
      if (is.null(mat)) next
      sub_edges <- as.vector(mat[idx_mask])
      diffs <- abs(orig_edges - sub_edges)
      r <- if (stats::sd(sub_edges, na.rm = TRUE) > 0) {
        stats::cor(orig_edges, sub_edges, method = method,
                   use = "complete.obs")
      } else NA_real_
      metrics$mean_abs_dev  [it, p_idx] <- mean(diffs,   na.rm = TRUE)
      metrics$median_abs_dev[it, p_idx] <- stats::median(diffs, na.rm = TRUE)
      metrics$correlation   [it, p_idx] <- r
      metrics$max_abs_dev   [it, p_idx] <- max(diffs,    na.rm = TRUE)
    }
  }

  # ---- Per-drop-proportion summary across iterations ----
  summary_rows <- lapply(metric_names, function(nm) {
    m <- metrics[[nm]]
    data.frame(
      metric    = nm,
      drop_prop = drop_prop,
      mean      = colMeans(m, na.rm = TRUE),
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
  summary_df <- do.call(rbind, summary_rows)

  # ---- CS-coefficient: max drop_prop where corr certainty-fraction >= threshold ----
  prop_above <- apply(metrics$correlation, 2L, function(c) {
    mean(c >= threshold, na.rm = TRUE)
  })
  qualifying <- which(prop_above >= certainty)
  cs <- if (length(qualifying) > 0L) max(drop_prop[qualifying]) else 0

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
    class = "net_edge_stability"
  )
}

#' @param x An edge-stability object.
#' @param digits Digits to display. Default `3`.
#' @param ... Additional arguments (ignored).
#' @return The input `x` invisibly.
#' @rdname edge_stability
#' @export
print.net_edge_stability <- function(x, digits = 3, ...) {
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

#' @rdname edge_stability
#' @export
print.net_edge_stability_group <- function(x, ...) {
  cat("Edge-weight Case-dropping Stability (grouped)\n")
  cat(sprintf("  %-12s  %s\n", "Network", "CS-coefficient"))
  for (nm in names(x)) {
    cat(sprintf("  %-12s  %.2f\n", nm, x[[nm]]$cs))
  }
  invisible(x)
}
