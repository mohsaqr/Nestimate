# ---- Vertex Bootstrap (Snijders & Borgatti 1999) ----

#' Vertex Bootstrap for Network-Level Statistics
#'
#' @description
#' Non-parametric vertex bootstrap of a single observed network
#' (Snijders & Borgatti 1999). Vertices are resampled with replacement and
#' the weight matrix is rebuilt from the original entries of the resampled
#' vertex pairs; network-level statistics computed on each replicate give
#' bootstrap distributions, standard errors, and confidence intervals.
#'
#' Unlike \code{\link{bootstrap_network}}, which resamples the underlying
#' cases (sequences or rows) and therefore requires the raw data stored in
#' the netobject, the vertex bootstrap needs \strong{only the weight
#' matrix}. It works on any \code{netobject} — including data-less ones
#' such as \code{\link{build_mlvar}} constituents or \code{as_tna(mcml)}
#' elements — and on plain weight matrices. The two procedures answer
#' different questions: the case bootstrap quantifies sampling-of-subjects
#' uncertainty in the edge weights; the vertex bootstrap quantifies
#' structural uncertainty of whole-network descriptives given the one
#' network you observed.
#'
#' @details
#' Each replicate draws \code{n} vertex indices with replacement and sets
#' \code{W_b[i, j] = W[idx_i, idx_j]}. When the same original vertex is
#' drawn for two different positions, the off-diagonal cell would be a
#' structural self-pair; following Snijders & Borgatti, such cells are
#' filled with the weight of a randomly chosen pair of distinct original
#' vertices. Diagonal entries carry the original self-weight of the
#' resampled vertex (\code{W[idx_i, idx_i]}) — self-loops are meaningful
#' in transition networks and are never altered. For undirected networks
#' the substitution is applied symmetrically so replicates stay symmetric.
#'
#' Built-in statistics (all computed on the off-diagonal part of the
#' weight matrix):
#' \describe{
#'   \item{\code{density}}{Proportion of non-zero off-diagonal cells.}
#'   \item{\code{mean_weight}}{Mean of the non-zero off-diagonal weights.}
#'   \item{\code{centralization}}{Freeman-type strength centralization:
#'     \code{sum(max(s) - s) / ((n - 1) * max(s))} where \code{s} is total
#'     node strength on absolute weights. 0 when all nodes have equal
#'     strength, approaching 1 for a star.}
#'   \item{\code{reciprocity}}{Directed networks only. Weighted
#'     reciprocity \code{sum(pmin(|W|, |t(W)|)) / sum(|W|)} over
#'     off-diagonal cells: the proportion of total weight that is
#'     reciprocated.}
#' }
#'
#' @param x A \code{netobject} (from \code{\link{build_network}} or any
#'   builder), a \code{cograph_network}, or a square numeric weight matrix.
#' @param iter Integer. Number of bootstrap replicates (default 1000).
#' @param ci_level Numeric. Significance level for the confidence
#'   intervals (default 0.05 for 95% CIs).
#' @param ci_method Character. \code{"percentile"} (default) for empirical
#'   quantile intervals, or \code{"basic"} for intervals reflected around
#'   the observed value (Davison & Hinkley 1997, eq. 5.6).
#' @param statistics Character vector selecting built-in statistics
#'   (see Details). Default: all applicable to the network's directedness.
#' @param statistic_fn Optional named list of functions, each taking the
#'   weight matrix and returning a single numeric value. Computed alongside
#'   the built-ins.
#' @param directed Logical or NULL. Directedness of the network. NULL
#'   (default) reads \code{x$directed} when available, otherwise falls
#'   back to matrix symmetry.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"net_vertex_bootstrap"} containing:
#' \describe{
#'   \item{summary}{Tidy data frame, one row per statistic: \code{statistic},
#'     \code{observed}, \code{boot_mean}, \code{boot_sd}, \code{bias},
#'     \code{ci_lower}, \code{ci_upper}.}
#'   \item{boot_stats}{\code{iter} x n_statistics matrix of replicate values.}
#'   \item{observed}{Named vector of observed statistics.}
#'   \item{iter, ci_level, ci_method, directed, n_nodes}{Configuration.}
#' }
#'
#' @references
#' Snijders, T. A. B., & Borgatti, S. P. (1999). Non-parametric standard
#' errors and tests for network statistics. \emph{Connections}, 22(2),
#' 161-170.
#'
#' Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and
#' their Application}. Cambridge University Press.
#'
#' @seealso \code{\link{bootstrap_network}} for case-resampling edge-weight
#'   inference, \code{\link{centrality_stability}} for case-dropping
#'   centrality stability.
#'
#' @examples
#' seqs <- data.frame(
#'   T1 = c("plan", "code", "debug", "plan", "test", "code"),
#'   T2 = c("code", "debug", "code", "plan", "code", "test"),
#'   T3 = c("debug", "code", "plan", "code", "debug", "plan"),
#'   T4 = c("test", "plan", "test", "debug", "plan", "code")
#' )
#' net <- build_network(seqs, method = "relative")
#' vb <- vertex_bootstrap(net, iter = 100, seed = 1)
#' vb$summary
#' \donttest{
#' plot(vb)
#' }
#'
#' @export
vertex_bootstrap <- function(x,
                             iter = 1000L,
                             ci_level = 0.05,
                             ci_method = c("percentile", "basic"),
                             statistics = NULL,
                             statistic_fn = NULL,
                             directed = NULL,
                             seed = NULL) {
  ci_method <- match.arg(ci_method)
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
    is.null(directed) ||
      (is.logical(directed) && length(directed) == 1L && !is.na(directed))
  )
  iter <- as.integer(iter)

  # ---- Extract weight matrix and directedness ----
  if (is.matrix(x)) {
    W <- x
    if (is.null(directed)) directed <- !isSymmetric(unname(W))
  } else if (inherits(x, c("netobject", "cograph_network"))) {
    W <- x$weights
    if (is.null(directed)) {
      directed <- if (is.logical(x$directed) && length(x$directed) == 1L &&
                      !is.na(x$directed)) {
        x$directed
      } else {
        !isSymmetric(unname(W))
      }
    }
  } else {
    stop("'x' must be a netobject, cograph_network, or square numeric ",
         "weight matrix.", call. = FALSE)
  }
  if (!is.matrix(W) || !is.numeric(W) || nrow(W) != ncol(W)) {
    stop("Weight matrix must be square and numeric.", call. = FALSE)
  }
  n <- nrow(W)
  if (n < 3) {
    stop("Vertex bootstrap requires at least 3 nodes.", call. = FALSE)
  }

  # ---- Resolve statistics ----
  builtin <- .vb_builtin_statistics(directed)
  if (is.null(statistics)) {
    statistics <- names(builtin)
  } else {
    stopifnot(is.character(statistics), length(statistics) >= 1)
    bad <- setdiff(statistics, names(builtin))
    if (length(bad) > 0) {
      stop("Unknown statistics: ", paste(bad, collapse = ", "),
           ". Available: ", paste(names(builtin), collapse = ", "),
           call. = FALSE)
    }
  }
  stat_fns <- builtin[statistics]
  if (!is.null(statistic_fn)) {
    stopifnot(
      is.list(statistic_fn), length(statistic_fn) >= 1,
      !is.null(names(statistic_fn)), all(nzchar(names(statistic_fn))),
      all(vapply(statistic_fn, is.function, logical(1)))
    )
    stat_fns <- c(stat_fns, statistic_fn)
  }
  stat_names <- names(stat_fns)
  n_stats <- length(stat_fns)

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Observed values ----
  observed <- vapply(stat_fns, function(f) as.numeric(f(W)), numeric(1))
  names(observed) <- stat_names

  # ---- Bootstrap replicates ----
  boot_raw <- vapply(seq_len(iter), function(b) {
    Wb <- .vertex_resample(W, directed)
    vapply(stat_fns, function(f) as.numeric(f(Wb)), numeric(1))
  }, numeric(n_stats))
  # vapply returns a vector for a single statistic, a matrix otherwise
  boot_stats <- if (n_stats == 1L) matrix(boot_raw, ncol = 1L) else t(boot_raw)
  colnames(boot_stats) <- stat_names

  # ---- CIs ----
  probs <- c(ci_level / 2, 1 - ci_level / 2)
  ci <- apply(boot_stats, 2, quantile, probs = probs, na.rm = TRUE)
  ci_lower <- ci[1, ]
  ci_upper <- ci[2, ]
  if (identical(ci_method, "basic")) {
    basic_lower <- 2 * observed - ci_upper
    ci_upper <- 2 * observed - ci_lower
    ci_lower <- basic_lower
  }

  boot_mean <- colMeans(boot_stats, na.rm = TRUE)
  boot_sd <- apply(boot_stats, 2, sd, na.rm = TRUE)

  summary_df <- data.frame(
    statistic = stat_names,
    observed = unname(observed),
    boot_mean = unname(boot_mean),
    boot_sd = unname(boot_sd),
    bias = unname(boot_mean - observed),
    ci_lower = unname(ci_lower),
    ci_upper = unname(ci_upper),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  result <- list(
    summary = summary_df,
    boot_stats = boot_stats,
    observed = observed,
    iter = iter,
    ci_level = ci_level,
    ci_method = ci_method,
    directed = directed,
    n_nodes = n
  )
  class(result) <- "net_vertex_bootstrap"
  result
}


#' Built-in network-level statistics for vertex_bootstrap
#' @noRd
.vb_builtin_statistics <- function(directed) {
  off <- function(W) W[row(W) != col(W)]
  fns <- list(
    density = function(W) {
      mean(off(W) != 0)
    },
    mean_weight = function(W) {
      w <- off(W)
      w <- w[w != 0]
      if (length(w) == 0) return(NA_real_)
      mean(w)
    },
    centralization = function(W) {
      A <- abs(W)
      s <- rowSums(A) + colSums(A) - 2 * diag(A)
      s_max <- max(s)
      if (s_max == 0) return(NA_real_)
      sum(s_max - s) / ((length(s) - 1) * s_max)
    }
  )
  if (directed) {
    fns$reciprocity <- function(W) {
      A <- abs(W)
      keep <- row(A) != col(A)
      total <- sum(A[keep])
      if (total == 0) return(NA_real_)
      sum(pmin(A, t(A))[keep]) / total
    }
  }
  fns
}


#' One vertex-bootstrap replicate of a weight matrix
#'
#' Resamples vertex indices with replacement and rebuilds the matrix from
#' original entries. Off-diagonal cells whose two positions drew the same
#' original vertex are filled with the weight of a random distinct pair
#' (Snijders & Borgatti 1999). The diagonal carries the original
#' self-weights of the resampled vertices.
#' @noRd
.vertex_resample <- function(W, directed) {
  n <- nrow(W)
  idx <- sample.int(n, n, replace = TRUE)
  Wb <- W[idx, idx, drop = FALSE]

  off_diag <- row(Wb) != col(Wb)
  collision <- outer(idx, idx, "==") & off_diag

  if (directed) {
    n_coll <- sum(collision)
    if (n_coll > 0) {
      Wb[collision] <- W[.vb_random_distinct_pairs(n, n_coll)]
    }
  } else {
    upper_coll <- collision & upper.tri(Wb)
    n_coll <- sum(upper_coll)
    if (n_coll > 0) {
      Wb[upper_coll] <- W[.vb_random_distinct_pairs(n, n_coll)]
      # Mirror so the replicate stays symmetric
      Wb[lower.tri(Wb)] <- t(Wb)[lower.tri(Wb)]
    }
  }

  diag(Wb) <- diag(W)[idx]
  Wb
}


#' k random ordered pairs of distinct vertex indices
#' @noRd
.vb_random_distinct_pairs <- function(n, k) {
  a <- sample.int(n, k, replace = TRUE)
  b <- sample.int(n, k, replace = TRUE)
  eq <- a == b
  while (any(eq)) {
    b[eq] <- sample.int(n, sum(eq), replace = TRUE)
    eq <- a == b
  }
  cbind(a, b)
}


#' Print a Vertex Bootstrap Result
#'
#' @param x A \code{net_vertex_bootstrap} object.
#' @param digits Number of digits to display (default 3).
#' @param ... Additional arguments (ignored).
#' @return \code{x}, invisibly.
#' @export
print.net_vertex_bootstrap <- function(x, digits = 3, ...) {
  cat("Vertex Bootstrap (Snijders & Borgatti)\n")
  cat(sprintf("  Nodes: %d | Directed: %s | Replicates: %d\n",
              x$n_nodes, x$directed, x$iter))
  cat(sprintf("  %.0f%% CIs (%s method)\n\n",
              100 * (1 - x$ci_level), x$ci_method))
  df <- x$summary
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  print(df, row.names = FALSE)
  invisible(x)
}


#' Summarize a Vertex Bootstrap Result
#'
#' @param object A \code{net_vertex_bootstrap} object.
#' @param ... Additional arguments (ignored).
#' @return The tidy summary data frame (one row per statistic).
#' @export
summary.net_vertex_bootstrap <- function(object, ...) {
  object$summary
}


#' Plot Vertex Bootstrap Distributions
#'
#' Histogram of the bootstrap distribution per statistic, with the observed
#' value (solid line) and confidence bounds (dashed lines).
#'
#' @param x A \code{net_vertex_bootstrap} object.
#' @param bins Number of histogram bins (default 30).
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.net_vertex_bootstrap <- function(x, bins = 30, ...) {
  long <- data.frame(
    statistic = rep(colnames(x$boot_stats), each = nrow(x$boot_stats)),
    value = as.vector(x$boot_stats),
    stringsAsFactors = FALSE
  )
  refs <- x$summary

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#4A6FE3",
                            color = "white", linewidth = 0.2) +
    ggplot2::geom_vline(data = refs,
                        ggplot2::aes(xintercept = .data$observed),
                        color = "#D33F6A", linewidth = 0.7) +
    ggplot2::geom_vline(data = refs,
                        ggplot2::aes(xintercept = .data$ci_lower),
                        linetype = "dashed", color = "gray30",
                        linewidth = 0.4) +
    ggplot2::geom_vline(data = refs,
                        ggplot2::aes(xintercept = .data$ci_upper),
                        linetype = "dashed", color = "gray30",
                        linewidth = 0.4) +
    ggplot2::facet_wrap(~statistic, scales = "free") +
    ggplot2::labs(
      x = "Bootstrap value", y = "Count",
      title = "Vertex bootstrap distributions",
      subtitle = sprintf(
        "%d replicates | observed (solid), %.0f%% CI (dashed, %s)",
        x$iter, 100 * (1 - x$ci_level), x$ci_method)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


# ---- Two-Network Vertex-Bootstrap Comparison ----

#' Compare Network-Level Statistics of Two Networks
#'
#' @description
#' Snijders & Borgatti (1999) two-network test: each network's statistics
#' get vertex-bootstrap standard errors, and each difference is tested with
#' \deqn{z = (\hat{\theta}_x - \hat{\theta}_y) /
#'   \sqrt{SE_x^2 + SE_y^2}}
#' against a standard normal reference. This is the comparison the vertex
#' bootstrap was originally proposed for: deciding whether two observed
#' networks differ in density, centralization, reciprocity, or any other
#' whole-network descriptive.
#'
#' @param x,y The two networks: \code{netobject}s, \code{cograph_network}s,
#'   square weight matrices, or precomputed \code{net_vertex_bootstrap}
#'   objects (then \code{iter}, \code{statistics}, \code{statistic_fn},
#'   \code{directed}, and \code{seed} are ignored for that argument).
#'   The two sides must cover exactly the same statistics; a mismatch
#'   (e.g., a directed network's \code{reciprocity} against an undirected
#'   one, or precomputed objects built with different \code{statistics}
#'   selections) is an error, never a silent subset.
#' @param labels Character vector of length 2 naming the networks in the
#'   output (default \code{c("x", "y")}).
#' @inheritParams vertex_bootstrap
#'
#' @return An object of class \code{"net_vertex_comparison"} containing:
#' \describe{
#'   \item{summary}{Tidy data frame, one row per statistic:
#'     \code{statistic}, the two observed values, \code{diff},
#'     \code{se_diff}, \code{z}, \code{p_value}, and a normal-approximation
#'     confidence interval for the difference.}
#'   \item{x, y}{The two \code{net_vertex_bootstrap} results.}
#'   \item{labels, ci_level}{Configuration.}
#' }
#' When both bootstrap SEs are zero (a statistic with no resampling
#' variation in either network) \code{z} and \code{p_value} are \code{NA}.
#'
#' @references
#' Snijders, T. A. B., & Borgatti, S. P. (1999). Non-parametric standard
#' errors and tests for network statistics. \emph{Connections}, 22(2),
#' 161-170.
#'
#' @seealso \code{\link{vertex_bootstrap}}, \code{\link{nct}} for the
#'   permutation-based comparison of edge-level structure when raw data
#'   are available, \code{\link{permutation}}.
#'
#' @examples
#' states <- c("plan", "code", "debug", "test")
#' s1 <- data.frame(
#'   T1 = rep(states, 5), T2 = rep(rev(states), 5),
#'   T3 = rep(states[c(2, 3, 4, 1)], 5)
#' )
#' s2 <- data.frame(
#'   T1 = rep(states[c(3, 1, 4, 2)], 5), T2 = rep(states, 5),
#'   T3 = rep(states[c(4, 3, 1, 2)], 5)
#' )
#' net1 <- build_network(s1, method = "relative")
#' net2 <- build_network(s2, method = "relative")
#' cmp <- vertex_compare(net1, net2, iter = 100, seed = 1)
#' cmp$summary
#'
#' @export
vertex_compare <- function(x, y,
                           iter = 1000L,
                           ci_level = 0.05,
                           statistics = NULL,
                           statistic_fn = NULL,
                           directed = NULL,
                           seed = NULL,
                           labels = c("x", "y")) {
  stopifnot(is.character(labels), length(labels) == 2, all(nzchar(labels)))

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  as_vb <- function(obj) {
    if (inherits(obj, "net_vertex_bootstrap")) return(obj)
    vertex_bootstrap(obj, iter = iter, ci_level = ci_level,
                     statistics = statistics, statistic_fn = statistic_fn,
                     directed = directed)
  }
  vb_x <- as_vb(x)
  vb_y <- as_vb(y)

  stats_x <- vb_x$summary$statistic
  stats_y <- vb_y$summary$statistic
  only_one_side <- c(setdiff(stats_x, stats_y), setdiff(stats_y, stats_x))
  if (length(only_one_side) > 0) {
    shared <- intersect(stats_x, stats_y)
    hint <- if (length(shared) > 0) {
      paste0(" Compute both with the same `statistics =` selection, e.g. ",
             "statistics = c(\"", paste(shared, collapse = "\", \""), "\").")
    } else {
      " The two cover no common statistics at all."
    }
    stop("The two networks cover different statistics (",
         paste(only_one_side, collapse = ", "),
         " present on only one side).", hint, call. = FALSE)
  }
  common <- stats_x
  sx <- vb_x$summary
  sy <- vb_y$summary[match(common, stats_y), ]

  diff <- sx$observed - sy$observed
  se_diff <- sqrt(sx$boot_sd^2 + sy$boot_sd^2)
  z <- ifelse(se_diff > 0, diff / se_diff, NA_real_)
  p_value <- 2 * pnorm(-abs(z))
  z_crit <- qnorm(1 - ci_level / 2)

  summary_df <- data.frame(
    statistic = common,
    observed_x = sx$observed,
    observed_y = sy$observed,
    diff = diff,
    se_diff = se_diff,
    z = z,
    p_value = p_value,
    ci_lower = diff - z_crit * se_diff,
    ci_upper = diff + z_crit * se_diff,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  names(summary_df)[2:3] <- paste0("observed_", make.names(labels))

  result <- list(
    summary = summary_df,
    x = vb_x,
    y = vb_y,
    labels = labels,
    ci_level = ci_level
  )
  class(result) <- "net_vertex_comparison"
  result
}


#' Print a Two-Network Vertex-Bootstrap Comparison
#'
#' @param x A \code{net_vertex_comparison} object.
#' @param digits Number of digits to display (default 3).
#' @param ... Additional arguments (ignored).
#' @return \code{x}, invisibly.
#' @export
print.net_vertex_comparison <- function(x, digits = 3, ...) {
  cat("Two-Network Vertex Bootstrap Comparison (Snijders & Borgatti)\n")
  cat(sprintf("  %s (%d nodes) vs %s (%d nodes) | %d replicates each\n\n",
              x$labels[1], x$x$n_nodes, x$labels[2], x$y$n_nodes, x$x$iter))
  df <- x$summary
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  stars <- ifelse(is.na(x$summary$p_value), "",
           ifelse(x$summary$p_value < 0.001, "***",
           ifelse(x$summary$p_value < 0.01, "**",
           ifelse(x$summary$p_value < 0.05, "*", ""))))
  df$sig <- stars
  print(df, row.names = FALSE)
  cat("---\nSignif. codes: *** p<0.001, ** p<0.01, * p<0.05\n")
  invisible(x)
}


#' Summarize a Two-Network Vertex-Bootstrap Comparison
#'
#' @param object A \code{net_vertex_comparison} object.
#' @param ... Additional arguments (ignored).
#' @return The tidy summary data frame (one row per statistic).
#' @export
summary.net_vertex_comparison <- function(object, ...) {
  object$summary
}


#' Plot a Two-Network Vertex-Bootstrap Comparison
#'
#' Forest plot of the statistic differences with normal-approximation
#' confidence intervals; differences whose interval excludes zero are the
#' statistically distinguishable ones.
#'
#' @param x A \code{net_vertex_comparison} object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.net_vertex_comparison <- function(x, ...) {
  df <- x$summary
  df$excludes_zero <- df$ci_lower > 0 | df$ci_upper < 0
  ggplot2::ggplot(df, ggplot2::aes(x = .data$diff, y = .data$statistic)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray40") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper,
                   color = .data$excludes_zero),
      width = 0.25, linewidth = 0.6
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$excludes_zero),
                        size = 2.2) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#4A6FE3", `TRUE` = "#D33F6A"),
      guide = "none"
    ) +
    ggplot2::labs(
      x = sprintf("Difference (%s - %s)", x$labels[1], x$labels[2]),
      y = NULL,
      title = "Network-level differences with vertex-bootstrap CIs",
      subtitle = sprintf(
        "%.0f%% normal-approximation CIs; red = interval excludes 0",
        100 * (1 - x$ci_level))
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
