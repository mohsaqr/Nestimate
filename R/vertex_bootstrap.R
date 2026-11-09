# ---- Vertex Bootstrap (Snijders & Borgatti 1999) ----

#' Vertex Bootstrap for Network-Level Statistics
#'
#' @description
#' `vertex_bootstrap()` resamples the nodes of one observed network and
#' returns bootstrap standard errors, bias and confidence intervals for
#' network-level statistics: density, mean weight, strength centralisation
#' and, for a directed network, reciprocity, plus any statistic of your own
#' (Snijders & Borgatti 1999). It needs only the weight matrix, so it works
#' on any `netobject` or `cograph_network`, including those that carry no
#' data, such as the constituents of [build_mlvar()] or the networks
#' [as_tna()] returns for an `mcml` object built from a weight matrix, and
#' on a plain square weight matrix.
#'
#' The intervals describe how much a whole-network statistic depends on
#' which nodes happen to make up the network, with the edge weights taken
#' as given. They say nothing about how the weights themselves vary from
#' sample to sample; when the sequences or rows behind a network are
#' available, [bootstrap_network()] resamples those instead. Two networks
#' are compared on these statistics with [vertex_compare()].
#'
#' @details
#' Each replicate draws `n` node indices with replacement and sets
#' `W_b[i, j] = W[idx_i, idx_j]`. When the same original node is drawn for
#' two different positions, the off-diagonal cell between them would be a
#' self-pair; following Snijders & Borgatti, it is filled with the weight of
#' a randomly chosen pair of distinct original nodes. Diagonal cells keep the
#' original self-weight of the resampled node (`W[idx_i, idx_i]`), since
#' self-loops are meaningful in transition networks. For an undirected
#' network the substitution is made in the upper triangle and mirrored, so
#' every replicate stays symmetric.
#'
#' The built-in statistics are computed on the off-diagonal cells:
#' \describe{
#'   \item{\code{density}}{The share of non-zero off-diagonal cells.}
#'   \item{\code{mean_weight}}{The mean of the non-zero off-diagonal
#'     weights.}
#'   \item{\code{centralization}}{Freeman-type strength centralisation,
#'     \code{sum(max(s) - s) / ((n - 1) * max(s))}, where \code{s} is each
#'     node's total (in plus out) absolute strength without its self-loop:
#'     0 when every node has the same strength, approaching 1 for a star.}
#'   \item{\code{reciprocity}}{Directed networks only. Weighted reciprocity,
#'     \code{sum(pmin(|W|, |t(W)|)) / sum(|W|)} over the off-diagonal cells:
#'     the share of the total weight that is reciprocated.}
#' }
#' A statistic is `NA` when it is undefined: `mean_weight` with no non-zero
#' cell, `centralization` when every strength is 0, `reciprocity` when the
#' total weight is 0. `NA` replicates are left out of the bootstrap mean,
#' standard deviation and quantiles.
#'
#' @param x The network: a `netobject` (from [build_network()] or any
#'   builder), a `cograph_network`, or a square numeric weight matrix, with
#'   at least 3 nodes. Other input raises an error. For the `print()` and
#'   `plot()` methods: a `net_vertex_bootstrap` result.
#' @param iter The number of bootstrap replicates, at least 2. Default
#'   `1000`.
#' @param ci_level The significance level of the confidence intervals, not
#'   their coverage: `0.05` (default) gives 95% intervals. Must lie strictly
#'   between 0 and 1.
#' @param ci_method How the intervals are formed. `"percentile"` (default)
#'   takes the `ci_level / 2` and `1 - ci_level / 2` quantiles of the
#'   replicates; `"basic"` reflects them around the observed value,
#'   `2 * observed - upper` and `2 * observed - lower` (Davison & Hinkley
#'   1997, eq. 5.6).
#' @param statistics The built-in statistics to compute (see *Details*).
#'   Default `NULL`: all of them, with `reciprocity` only for a directed
#'   network. A name that is not available for the network, such as
#'   `reciprocity` on an undirected one, raises an error.
#' @param statistic_fn Statistics of your own, as a named list of functions,
#'   each taking the weight matrix and returning a single number. They are
#'   computed beside the built-ins on the whole matrix, diagonal included.
#'   Default `NULL`.
#' @param directed Whether the network is directed, which decides how
#'   collisions are filled and whether `reciprocity` is available. Default
#'   `NULL` reads the network's own directedness, and for a matrix (or a
#'   network that records none) treats an asymmetric matrix as directed.
#' @param seed A seed passed to [set.seed()] before resampling, so the run
#'   can be repeated; the session's random-number stream then continues from
#'   it. Default `NULL` uses the current stream.
#'
#' @return A `net_vertex_bootstrap` object. Printing it shows the number of
#'   nodes, the directedness, the number of replicates and the interval
#'   level and method, then the summary table; `summary()` returns that
#'   table as a data frame with one row per statistic and the columns
#'   `statistic`, `observed`, `boot_mean`, `boot_sd`, `bias`
#'   (`boot_mean - observed`), `ci_lower` and `ci_upper`; `plot()` draws the
#'   bootstrap distributions. The object holds:
#' \describe{
#'   \item{summary}{The summary table described above.}
#'   \item{boot_stats}{A matrix of the replicate values, `iter` rows by one
#'     column per statistic.}
#'   \item{observed}{The observed statistics, as a named numeric vector.}
#'   \item{iter, ci_level, ci_method, directed, n_nodes}{The settings used
#'     and the number of nodes.}
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
#' summary(vb)
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
  .assert(
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
    stop(errorCondition(
      paste0(
        "'x' must be a netobject, cograph_network, or square numeric ",
        "weight matrix."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.matrix(W) || !is.numeric(W) || nrow(W) != ncol(W)) {
    stop(errorCondition("Weight matrix must be square and numeric.",
      class = "nestimate_bad_input", call = NULL))
  }
  n <- nrow(W)
  if (n < 3) {
    stop(errorCondition("Vertex bootstrap requires at least 3 nodes.",
      class = "nestimate_bad_input", call = NULL))
  }

  # ---- Resolve statistics ----
  builtin <- .vb_builtin_statistics(directed)
  if (is.null(statistics)) {
    statistics <- names(builtin)
  } else {
    .assert(is.character(statistics), length(statistics) >= 1)
    bad <- setdiff(statistics, names(builtin))
    if (length(bad) > 0) {
      stop(errorCondition(
        paste0(
          "Unknown statistics: ",
          paste(bad, collapse = ", "),
          ". Available: ",
          paste(names(builtin), collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
  }
  stat_fns <- builtin[statistics]
  if (!is.null(statistic_fn)) {
    .assert(
      is.list(statistic_fn), length(statistic_fn) >= 1,
      !is.null(names(statistic_fn)), all(nzchar(names(statistic_fn))),
      all(vapply(statistic_fn, is.function, logical(1)))
    )
    stat_fns <- c(stat_fns, statistic_fn)
  }
  stat_names <- names(stat_fns)
  n_stats <- length(stat_fns)

  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
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


#' @rdname vertex_bootstrap
#' @param digits For the `print()` method: the number of decimals shown.
#'   Default `3`.
#' @param ... Not used.
#' @param object For the `summary()` method: a `net_vertex_bootstrap`
#'   result.
#' @param bins For the `plot()` method: the number of histogram bins.
#'   Default `30`.
#' @return The `print()` method returns its input invisibly, `summary()` the
#'   summary data frame and `plot()` a `ggplot` object.
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


#' @rdname vertex_bootstrap
#' @export
summary.net_vertex_bootstrap <- function(object, ...) {
  object$summary
}


#' @rdname vertex_bootstrap
#' @section Methods:
#' `plot()` draws a histogram of the replicate values for each statistic,
#' one panel per statistic with its own scales, a solid line at the observed
#' value and dashed lines at the interval bounds.
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
#' `vertex_compare()` tests whether two observed networks differ in their
#' network-level statistics, the two-network test Snijders & Borgatti
#' (1999) proposed the vertex bootstrap for. Each network's statistics get
#' vertex-bootstrap standard errors from [vertex_bootstrap()], and each
#' difference is tested with
#' \deqn{z = (\hat{\theta}_x - \hat{\theta}_y) /
#'   \sqrt{SE_x^2 + SE_y^2}}
#' against a standard normal distribution, with a normal-approximation
#' confidence interval for the difference. It covers density, mean weight,
#' centralisation, reciprocity and any statistic of your own.
#'
#' The test is about whole-network descriptives, and its standard errors
#' reflect the resampling of nodes, not of the cases behind the weights.
#' The p-values are not adjusted for the number of statistics tested; pass
#' them to [stats::p.adjust()] when several are read together. When the raw
#' data of both networks are available, [permutation()] and [nct()] compare
#' them edge by edge.
#'
#' @param x,y The two networks: `netobject`s, `cograph_network`s, square
#'   weight matrices, or `net_vertex_bootstrap` results computed beforehand.
#'   A precomputed result is used as it is, so `iter`, `statistics`,
#'   `statistic_fn` and `directed` do not apply to it. The two must cover
#'   the same statistics: a mismatch, such as a directed network's
#'   `reciprocity` against an undirected network, or two precomputed results
#'   with different `statistics`, raises an error rather than comparing a
#'   subset. For the `print()` and `plot()` methods, `x` is a
#'   `net_vertex_comparison` result.
#' @param ci_level The significance level of the normal-approximation
#'   intervals for the differences: `0.05` (default) gives 95% intervals. It
#'   is also passed to the bootstraps run here; with a precomputed result,
#'   the differences still use this level.
#' @param statistics The built-in statistics to compare, as named for
#'   [vertex_bootstrap()]: `density`, `mean_weight`, `centralization` and,
#'   for directed networks, `reciprocity`. Default `NULL`: all that apply.
#'   A name not available for a network raises an error.
#' @param seed A seed passed to [set.seed()] once, before both bootstraps
#'   run, so the comparison can be repeated; the session's random-number
#'   stream then continues from it. Default `NULL` uses the current stream.
#' @param labels Two names for the networks, used in the print header, the
#'   plot axis and the column names `observed_<label>` (made syntactic with
#'   [make.names()]). Default `c("x", "y")`.
#' @param directed Whether the networks are directed, which decides how
#'   collisions are filled and whether `reciprocity` is available; one value
#'   applies to both networks, not to a precomputed result. Default `NULL`
#'   reads each network's own directedness, and for a matrix (or a network
#'   that records none) treats an asymmetric matrix as directed.
#' @inheritParams vertex_bootstrap
#'
#' @return A `net_vertex_comparison` object. Printing it shows the two
#'   networks' sizes and the first network's number of replicates, then the
#'   summary table
#'   with significance stars; `summary()` returns that table as a data frame
#'   with one row per statistic and the columns `statistic`, the two
#'   observed values (`observed_<label>`), `diff` (first minus second),
#'   `se_diff`, `z`, `p_value` (two-sided), `ci_lower` and `ci_upper`;
#'   `plot()` draws the differences with their intervals. When a statistic
#'   has no bootstrap variation in either network (both standard errors 0),
#'   its `z` and `p_value` are `NA`. The object holds:
#' \describe{
#'   \item{summary}{The summary table described above.}
#'   \item{x, y}{The two `net_vertex_bootstrap` results.}
#'   \item{labels, ci_level}{The settings used.}
#' }
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
#' summary(cmp)
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
  .assert(is.character(labels), length(labels) == 2, all(nzchar(labels)))

  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
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
    stop(errorCondition(paste0("The two networks cover different statistics (",
         paste(only_one_side, collapse = ", "),
         " present on only one side).", hint),
      class = "nestimate_bad_input", call = NULL))
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


#' @rdname vertex_compare
#' @param digits For the `print()` method: the number of decimals shown.
#'   Default `3`.
#' @param ... Not used.
#' @param object For the `summary()` method: a `net_vertex_comparison`
#'   result.
#' @return The `print()` method returns its input invisibly, `summary()` the
#'   summary data frame and `plot()` a `ggplot` object.
#' @export
print.net_vertex_comparison <- function(x, digits = 3, ...) {
  cat("Two-Network Vertex Bootstrap Comparison (Snijders & Borgatti)\n")
  # Two precomputed bootstraps may have used different numbers of replicates.
  replicates <- if (isTRUE(x$x$iter == x$y$iter)) {
    sprintf("%d replicates each", x$x$iter)
  } else {
    sprintf("%d and %d replicates", x$x$iter, x$y$iter)
  }
  cat(sprintf("  %s (%d nodes) vs %s (%d nodes) | %s\n\n",
              x$labels[1], x$x$n_nodes, x$labels[2], x$y$n_nodes, replicates))
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


#' @rdname vertex_compare
#' @export
summary.net_vertex_comparison <- function(object, ...) {
  object$summary
}


#' @rdname vertex_compare
#' @section Methods:
#' `plot()` draws each statistic's difference as a point with its
#' normal-approximation interval, beside a dashed line at zero; intervals
#' that exclude zero are drawn in red, the others in blue.
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
