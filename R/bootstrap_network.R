# ---- Bootstrap Network Estimation ----

#' Bootstrap a Network Estimate
#'
#' @description
#' Non-parametric bootstrap for any network estimated by
#' \code{\link{build_network}}. Works with all built-in methods
#' (transition and association) as well as custom registered estimators.
#'
#' For transition methods (\code{"relative"}, \code{"frequency"},
#' \code{"co_occurrence"}), uses a fast pre-computation strategy:
#' per-sequence count matrices are computed once, and each bootstrap
#' iteration only resamples sequences via \code{colSums} (C-level)
#' plus lightweight post-processing. Data must be in wide format for
#' transition bootstrap; use \code{\link{convert_sequence_format}} to
#' convert long-format data first.
#'
#' For association methods (\code{"cor"}, \code{"pcor"}, \code{"glasso"},
#' and custom estimators), the full estimator is called on resampled rows
#' each iteration.
#'
#' @param x A \code{netobject} from \code{\link{build_network}}.
#'   The data, method, params, scaling, threshold, and level are all
#'   extracted from this object.
#' @param iter Integer. Number of bootstrap iterations (default: 1000).
#' @param ci_level Numeric. Significance level for CIs and p-values
#'   (default: 0.05).
#' @param inference Character. \code{"stability"} (default) tests whether
#'   bootstrap replicates fall within a multiplicative consistency range
#'   around the original weight. \code{"threshold"} tests whether
#'   replicates exceed a fixed edge threshold.
#' @param consistency_range Numeric vector of length 2. Multiplicative
#'   bounds for stability inference (default: \code{c(0.75, 1.25)}).
#' @param edge_threshold Numeric or NULL. Fixed threshold for
#'   \code{inference = "threshold"}. If NULL, defaults to the 10th
#'   percentile of absolute original edge weights.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"net_bootstrap"} containing:
#' \describe{
#'   \item{original}{The original \code{netobject}.}
#'   \item{mean}{Bootstrap mean weight matrix.}
#'   \item{sd}{Bootstrap SD matrix.}
#'   \item{p_values}{P-value matrix.}
#'   \item{significant}{Original weights where p < ci_level, else 0.}
#'   \item{ci_lower}{Lower CI bound matrix.}
#'   \item{ci_upper}{Upper CI bound matrix.}
#'   \item{cr_lower}{Consistency range lower bound (stability only).}
#'   \item{cr_upper}{Consistency range upper bound (stability only).}
#'   \item{summary}{Long-format data frame of edge-level statistics.}
#'   \item{model}{Pruned \code{netobject} (non-significant edges zeroed).}
#'   \item{method, params, iter, ci_level, inference}{Bootstrap config.}
#'   \item{consistency_range, edge_threshold}{Inference parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Build network, then bootstrap it
#' net <- build_network(tna::group_regulation, method = "relative")
#' boot <- bootstrap_network(net, iter = 500)
#' print(boot)
#' summary(boot)
#'
#' # One-hot data
#' net_wtna <- build_network(onehot_df, method = "wtna", actor = "student")
#' boot_wtna <- bootstrap_network(net_wtna, iter = 200)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{print.net_bootstrap}},
#'   \code{\link{summary.net_bootstrap}}, \code{\link{plot.net_bootstrap}}
#'
#' @importFrom stats quantile sd
#' @export
bootstrap_network <- function(x,
                              iter = 1000L,
                              ci_level = 0.05,
                              inference = "stability",
                              consistency_range = c(0.75, 1.25),
                              edge_threshold = NULL,
                              seed = NULL) {

  # ---- netobject_group dispatch: bootstrap each element ----
  if (inherits(x, "netobject_group")) {
    results <- lapply(x, function(net) {
      bootstrap_network(net, iter = iter, ci_level = ci_level,
                        inference = inference,
                        consistency_range = consistency_range,
                        edge_threshold = edge_threshold, seed = seed)
    })
    return(results)
  }

  # Accept netobject or cograph_network
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (!inherits(x, "netobject")) {
    stop("'x' must be a netobject from build_network().", call. = FALSE)
  }
  if (is.null(x$data)) {
    stop("netobject does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  data <- x$data
  method <- x$method
  params <- x$params
  scaling <- x$scaling
  threshold <- x$threshold
  level <- x$level

  # ---- Input validation ----
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
    is.character(inference), length(inference) == 1,
    is.numeric(consistency_range), length(consistency_range) == 2
  )
  iter <- as.integer(iter)
  inference <- match.arg(inference, c("stability", "threshold"))

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Resolve method and get estimator ----
  method <- .resolve_method_alias(method)
  estimator <- get_estimator(method)
  directed <- estimator$directed

  # ---- Original network is x itself ----
  original <- x
  states <- original$nodes
  n_states <- length(states)

  # ---- Dispatch bootstrap ----
  # Fast precomputed path only works when data has usable columns
  has_data <- is.data.frame(data) && ncol(data) > 0L && nrow(data) > 0L
  if (method %in% c("relative", "frequency", "co_occurrence")) {
    if (!has_data) {
      stop("Bootstrap requires the original data stored in the netobject. ",
           "For wtna/cna networks, use wtna() directly instead of ",
           "build_network(method='cna').", call. = FALSE)
    }
    boot_matrices <- .bootstrap_transition(
      data = data, method = method, params = params, states = states,
      scaling = scaling, threshold = threshold, iter = iter
    )
  } else {
    boot_matrices <- .bootstrap_association(
      data = data, estimator = estimator, params = params, states = states,
      scaling = scaling, threshold = threshold, iter = iter,
      level = level, id_col = params$id
    )
  }

  # ---- Auto edge_threshold for threshold inference ----
  if (inference == "threshold" && is.null(edge_threshold)) {
    abs_weights <- abs(as.vector(original$matrix))
    nz_weights <- abs_weights[abs_weights > 0]
    edge_threshold <- if (length(nz_weights) > 0) {
      quantile(nz_weights, probs = 0.10)
    } else {
      0
    }
  }

  # ---- Compute statistics ----
  stats <- .compute_bootstrap_stats(
    boot_matrices = boot_matrices,
    original_matrix = original$matrix,
    states = states,
    directed = directed,
    iter = iter,
    ci_level = ci_level,
    inference = inference,
    consistency_range = consistency_range,
    edge_threshold = edge_threshold
  )

  # ---- Build summary data frame ----
  summary_df <- .build_bootstrap_summary(
    stats = stats,
    original_matrix = original$matrix,
    states = states,
    directed = directed,
    ci_level = ci_level,
    inference = inference
  )

  # ---- Build pruned model ----
  pruned_matrix <- stats$significant
  pruned_edges <- .extract_edges_from_matrix(pruned_matrix, directed = directed)

  model <- list(
    matrix = pruned_matrix,
    nodes = states,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = n_states,
    n_edges = nrow(pruned_edges),
    edges = pruned_edges,
    level = level
  )
  class(model) <- "netobject"

  # ---- Assemble result ----
  result <- list(
    original          = original,
    mean              = stats$mean,
    sd                = stats$sd,
    p_values          = stats$p_values,
    significant       = stats$significant,
    ci_lower          = stats$ci_lower,
    ci_upper          = stats$ci_upper,
    cr_lower          = stats$cr_lower,
    cr_upper          = stats$cr_upper,
    summary           = summary_df,
    model             = model,
    method            = method,
    params            = params,
    iter              = iter,
    ci_level          = ci_level,
    inference         = inference,
    consistency_range = consistency_range,
    edge_threshold    = edge_threshold
  )
  class(result) <- "net_bootstrap"
  result
}


# ---- Transition fast path ----

#' Bootstrap transition networks via pre-computed per-sequence counts
#' @noRd
.bootstrap_transition <- function(data, method, params, states,
                                  scaling, threshold, iter) {
  n_states <- length(states)
  nbins <- n_states * n_states
  is_relative <- method == "relative"

  # Pre-compute per-sequence count matrix ONCE
  trans_2d <- .precompute_per_sequence(data, method, params, states)
  n_seq <- nrow(trans_2d)

  # vapply: each iteration resamples sequences + sums + post-processes
  boot_flat <- vapply(seq_len(iter), function(i) {
    idx <- sample.int(n_seq, n_seq, replace = TRUE)
    boot_counts <- colSums(trans_2d[idx, , drop = FALSE])
    # Inline post-processing (no dimnames — we flatten immediately)
    mat <- matrix(boot_counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    as.vector(mat)
  }, numeric(nbins))

  t(boot_flat)
}


#' Pre-compute per-sequence count matrices
#'
#' Returns an n_seq x (n_states^2) matrix where each row holds the
#' flattened transition (or co-occurrence) counts for one sequence.
#'
#' Only wide-format data is supported for the bootstrap fast path.
#' For long-format data, convert first using
#' \code{\link{convert_sequence_format}}.
#' @noRd
.precompute_per_sequence <- function(data, method, params, states) {
  # Determine format (mirrors estimator dispatch)
  format <- params$format %||% "auto"
  action <- params$action %||% "Action"

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  if (format == "long") {
    stop(
      "Bootstrap fast path requires wide-format data. ",
      "Convert with convert_sequence_format() first.",
      call. = FALSE
    )
  }

  id_col <- params$id %||% params$id_col
  cols <- params$cols
  .precompute_per_sequence_wide(data, method, cols, id_col, states)
}


#' Pre-compute per-sequence counts from wide format (transition + co-occurrence)
#'
#' Single function handling both consecutive-pair transitions
#' (relative/frequency) and all-column-pair co-occurrences.
#' @noRd
.precompute_per_sequence_wide <- function(data, method, cols, id_col, states) {
  state_cols <- .select_state_cols(data, id_col, cols)
  mat <- as.matrix(data[, state_cols, drop = FALSE])
  nc <- ncol(mat)
  nr <- nrow(mat)
  n_states <- length(states)
  nbins <- n_states * n_states

  # Integer-encode the entire matrix once
  int_mat <- matrix(match(mat, states), nrow = nr, ncol = nc)

  if (method %in% c("relative", "frequency")) {
    # Consecutive column pairs: (col[i], col[i+1])
    from_mat <- int_mat[, -nc, drop = FALSE]
    to_mat <- int_mat[, -1L, drop = FALSE]

    row_ids <- rep(seq_len(nr), times = nc - 1L)
    from_vec <- as.vector(from_mat)
    to_vec <- as.vector(to_mat)

    valid <- !is.na(from_vec) & !is.na(to_vec)
    row_ids <- row_ids[valid]
    from_vec <- from_vec[valid]
    to_vec <- to_vec[valid]

    pair_idx <- (from_vec - 1L) * n_states + to_vec
    combined_idx <- (row_ids - 1L) * nbins + pair_idx

    counts <- tabulate(combined_idx, nbins = nr * nbins)
    matrix(as.numeric(counts), nrow = nr, ncol = nbins, byrow = TRUE)

  } else {
    # Co-occurrence: all column pairs (i, j) where i < j
    result <- matrix(0, nrow = nr, ncol = nbins)

    for (i in seq_len(nc - 1L)) {
      col_i <- int_mat[, i]
      for (j in seq(i + 1L, nc)) {
        col_j <- int_mat[, j]
        valid <- !is.na(col_i) & !is.na(col_j)
        fi <- col_i[valid]
        tj <- col_j[valid]
        rows_valid <- which(valid)

        # Bidirectional: forward + reverse
        idx_fwd <- (fi - 1L) * n_states + tj
        idx_rev <- (tj - 1L) * n_states + fi
        combined <- c(
          (rows_valid - 1L) * nbins + idx_fwd,
          (rows_valid - 1L) * nbins + idx_rev
        )

        tab <- tabulate(combined, nbins = nr * nbins)
        result <- result + matrix(tab, nrow = nr, ncol = nbins, byrow = TRUE)
      }
    }

    # Self-pairs are double-counted by bidirectional approach
    diag_indices <- (seq_len(n_states) - 1L) * n_states + seq_len(n_states)
    result[, diag_indices] <- result[, diag_indices] / 2

    result
  }
}


# ---- Association path ----

#' Bootstrap association networks via full estimator calls
#' @noRd
.bootstrap_association <- function(data, estimator, params, states,
                                   scaling, threshold, iter,
                                   level, id_col) {
  n <- nrow(data)
  n_states <- length(states)
  nbins <- n_states * n_states

  boot_flat <- vapply(seq_len(iter), function(i) {
    boot_data <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]

    # For multilevel: apply decomposition to bootstrap sample
    if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
      boot_data <- tryCatch(
        .decompose_multilevel(boot_data, id_col = id_col, level = level),
        error = function(e) NULL
      )
      if (is.null(boot_data)) return(rep(NA_real_, nbins))
    }

    est <- tryCatch(
      do.call(estimator$fn, c(list(data = boot_data), params)),
      error = function(e) NULL
    )
    if (is.null(est)) return(rep(NA_real_, nbins))

    mat <- est$matrix
    # Align to expected states order
    if (!identical(rownames(mat), states)) {
      common <- intersect(states, rownames(mat))
      if (length(common) < n_states) return(rep(NA_real_, nbins))
      mat <- mat[states, states]
    }

    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    as.vector(mat)
  }, numeric(nbins))

  t(boot_flat)
}


# ---- Vectorized statistics ----

#' Compute bootstrap statistics from the iter x nbins matrix
#' @noRd
.compute_bootstrap_stats <- function(boot_matrices, original_matrix, states,
                                     directed, iter, ci_level, inference,
                                     consistency_range, edge_threshold) {
  n_states <- length(states)
  orig_flat <- as.vector(original_matrix)

  # Mean and SD (vectorized)
  weights_mean <- colMeans(boot_matrices, na.rm = TRUE)
  weights_sd <- apply(boot_matrices, 2, sd, na.rm = TRUE)

  # Percentile CIs (vectorized via apply)
  ci_lower <- apply(boot_matrices, 2, quantile,
                    probs = ci_level / 2, na.rm = TRUE)
  ci_upper <- apply(boot_matrices, 2, quantile,
                    probs = 1 - ci_level / 2, na.rm = TRUE)

  # p-values via vectorized sweep
  valid_rows <- rowSums(is.na(boot_matrices)) == 0
  bm <- boot_matrices[valid_rows, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 1) {
    p_values <- rep(1, length(orig_flat))
  } else if (inference == "stability") {
    cr_low <- pmin(orig_flat * consistency_range[1],
                   orig_flat * consistency_range[2])
    cr_high <- pmax(orig_flat * consistency_range[1],
                    orig_flat * consistency_range[2])
    below <- sweep(bm, 2, cr_low, "<")
    above <- sweep(bm, 2, cr_high, ">")
    p_counts <- colSums(below | above)
    p_values <- (p_counts + 1) / (n_valid + 1)
  } else {
    below_thresh <- sweep(abs(bm), 2, edge_threshold, "<")
    p_counts <- colSums(below_thresh)
    p_values <- (p_counts + 1) / (n_valid + 1)
  }

  # Reshape all to n_states x n_states matrices
  to_mat <- function(x) {
    matrix(x, n_states, n_states, dimnames = list(states, states))
  }

  p_mat <- to_mat(p_values)
  sig_mask <- (p_mat < ci_level) * 1
  sig_mat <- original_matrix * sig_mask

  list(
    mean        = to_mat(weights_mean),
    sd          = to_mat(weights_sd),
    ci_lower    = to_mat(ci_lower),
    ci_upper    = to_mat(ci_upper),
    p_values    = p_mat,
    significant = sig_mat,
    cr_lower    = to_mat(pmin(orig_flat * consistency_range[1],
                              orig_flat * consistency_range[2])),
    cr_upper    = to_mat(pmax(orig_flat * consistency_range[1],
                              orig_flat * consistency_range[2]))
  )
}


#' Build long-format summary data frame from bootstrap stats
#' @noRd
.build_bootstrap_summary <- function(stats, original_matrix, states, directed,
                                     ci_level, inference) {
  n <- length(states)
  dt <- data.table::data.table(
    from     = rep(states, each = n),
    to       = rep(states, times = n),
    weight   = as.vector(t(original_matrix)),
    mean     = as.vector(t(stats$mean)),
    sd       = as.vector(t(stats$sd)),
    p_value  = as.vector(t(stats$p_values)),
    sig      = as.vector(t(stats$p_values)) < ci_level,
    ci_lower = as.vector(t(stats$ci_lower)),
    ci_upper = as.vector(t(stats$ci_upper))
  )
  if (inference == "stability") {
    dt[, cr_lower := as.vector(t(stats$cr_lower))]
    dt[, cr_upper := as.vector(t(stats$cr_upper))]
  }

  # Filter: keep only non-zero original edges (including self-loops)
  if (directed) {
    dt <- dt[weight != 0]
  } else {
    dt <- dt[weight != 0 & from <= to]
  }

  as.data.frame(dt)
}


# ---- S3 Methods ----

#' Print Method for net_bootstrap
#'
#' @param x A \code{net_bootstrap} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.net_bootstrap <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network",
    attention     = "Attention Network (decay-weighted transitions)",
    wtna          = "Window TNA (transitions)"
  )
  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label <- if (x$original$directed) " [directed]" else " [undirected]"

  cat("Bootstrap:", label, dir_label, "\n", sep = "")
  cat(sprintf("  Iterations: %d  |  CI level: %.2f  |  Inference: %s\n",
              x$iter, x$ci_level, x$inference))
  cat(sprintf("  Nodes: %d  |  Original edges: %d  |  Significant edges: %d\n",
              x$original$n_nodes, x$original$n_edges, x$model$n_edges))

  if (x$inference == "stability") {
    cat(sprintf("  Consistency range: [%.2f, %.2f]\n",
                x$consistency_range[1], x$consistency_range[2]))
  } else {
    cat(sprintf("  Edge threshold: %g\n", x$edge_threshold))
  }

  invisible(x)
}


#' Summary Method for net_bootstrap
#'
#' @param object A \code{net_bootstrap} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with edge-level bootstrap statistics.
#'
#' @export
summary.net_bootstrap <- function(object, ...) {
  object$summary
}


#' Plot Method for net_bootstrap
#'
#' @description
#' Plots the original and bootstrapped (significant-only) networks
#' side by side using \code{cograph::splot()}.
#'
#' @param x A \code{net_bootstrap} object.
#' @param ... Additional arguments passed to \code{cograph::splot()}.
#'
#' @importFrom graphics par
#' @export
plot.net_bootstrap <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }

  old_par <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old_par))

  node_cols <- .node_colors(x$original$n_nodes)

  dots <- list(
    directed = x$original$directed,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    ...
  )

  dots_orig <- c(list(x = x$original$matrix, title = "Original"), dots)
  dots_sig <- c(list(x = x$significant, title = "Significant"), dots)

  do.call(cograph::splot, dots_orig)
  do.call(cograph::splot, dots_sig)
}
