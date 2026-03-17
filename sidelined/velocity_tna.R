#' Velocity TNA: Derivative Analysis for Transition Networks
#'
#' @description
#' Computes the rate of change (velocity) of transition probabilities over time.
#' For each edge P(i -> j) in a series of transition matrices, the function
#' estimates how quickly that transition probability is changing. The result is a
#' velocity matrix (same n_states x n_states shape) showing which transitions
#' are strengthening or weakening over time.
#'
#' Three estimation methods are available:
#' \describe{
#'   \item{\code{"regression"}}{Fits OLS or beta regression per edge across
#'     all time points. Most robust to noise. Returns slopes, standard errors,
#'     p-values, and R-squared per edge. Default method.}
#'   \item{\code{"glla"}}{Generalized Local Linear Approximation (dynEGA-style).
#'     Local polynomial derivative estimation with implicit smoothing.}
#'   \item{\code{"finite_difference"}}{Direct numerical differentiation.
#'     Simplest but most sensitive to noise.}
#' }
#'
#' @param data Input data in one of three formats:
#'   \enumerate{
#'     \item Output from \code{cograph::tna_windows()} (detected by
#'       \code{$windows} field where each element has \code{$weights}).
#'     \item A list of square numeric matrices (transition matrices over time).
#'     \item A data frame of long-format sequence data (requires
#'       \code{time_col}).
#'   }
#' @param method Character. Estimation method: \code{"regression"} (default),
#'   \code{"glla"}, or \code{"finite_difference"}.
#' @param regression_type Character. For \code{method = "regression"}: type of
#'   regression. \code{"ols"} (default), \code{"beta"} (requires the
#'   \pkg{betareg} package), or \code{"logistic"} (per-edge logistic
#'   regression on individual transitions --- requires wide-format data,
#'   bypasses matrix aggregation, uses all raw transitions for maximum
#'   statistical power).
#' @param weights Numeric vector or NULL. Regression weights per time point
#'   (e.g., number of transitions per matrix). If \code{NULL} (default),
#'   weights are auto-computed from the matrix element sums, giving more
#'   influence to matrices based on larger samples. Set to \code{rep(1, T)}
#'   for unweighted regression. Only used when \code{method = "regression"}.
#' @param n_embed Integer. GLLA embedding dimension (must be odd, >= 3).
#'   Default: \code{5L}.
#' @param delta Numeric. Time interval between consecutive observations.
#'   Default: \code{1}.
#' @param order Integer. Polynomial/derivative order: 1 = linear velocity,
#'   2 = quadratic (also returns acceleration). Default: \code{1L}.
#' @param diff_type Character. For \code{method = "finite_difference"}: type of
#'   differencing. One of \code{"central"}, \code{"forward"}, \code{"backward"}.
#'   Default: \code{"central"}.
#' @param time_col Character or NULL. Column name for the time window
#'   (required for data frame input). Default: \code{NULL}.
#' @param id_col Character or NULL. Column name for the sequence/individual ID
#'   (data frame input). Default: \code{NULL}.
#' @param action_col Character or NULL. Column name for the action/state
#'   (data frame input). Default: \code{NULL}.
#' @param scaling Character vector or NULL. Optional post-estimation scaling
#'   to apply when building matrices from raw data. Default: \code{NULL}.
#'
#' @return An object of class \code{"tna_velocity"} containing:
#' \describe{
#'   \item{velocity_matrix}{Velocity matrix (n_states x n_states). For
#'     regression: slope coefficients. For GLLA/FD: mean derivative.}
#'   \item{velocity_series}{List of velocity matrices per time point. For
#'     regression: predicted velocity at each time point (constant for
#'     order=1, varying for order=2).}
#'   \item{acceleration_matrix}{Acceleration matrix (NULL if order < 2).}
#'   \item{acceleration_series}{List of acceleration matrices (NULL if
#'     order < 2).}
#'   \item{smoothed_matrices}{Fitted/smoothed matrices at each time point.}
#'   \item{original_matrices}{Input transition matrices.}
#'   \item{edge_stats}{For regression: data frame with per-edge statistics
#'     (from, to, slope, se, t_value, p_value, r_squared). NULL for other
#'     methods.}
#'   \item{nodes}{Character vector of state names.}
#'   \item{n_timepoints}{Number of original time points.}
#'   \item{n_embedded}{Number of output time points.}
#'   \item{method}{Method used.}
#'   \item{n_embed}{Embedding dimension (GLLA only).}
#'   \item{delta}{Time interval used.}
#'   \item{order}{Derivative/polynomial order.}
#'   \item{directed}{\code{TRUE} (velocity networks are always directed).}
#'   \item{n_nodes}{Number of nodes.}
#'   \item{n_edges}{Number of non-zero velocity edges.}
#'   \item{edges}{Data frame with \code{from}, \code{to}, \code{weight}.}
#' }
#'
#' @examples
#' \dontrun{
#' # From a list of transition matrices
#' mats <- lapply(1:10, function(t) {
#'   m <- matrix(runif(9), 3, 3)
#'   m / rowSums(m)
#' })
#' vel <- velocity_tna(mats)
#' print(vel)
#' summary(vel)
#' plot(vel)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_network}}
#'
#' @importFrom stats embed lm glm binomial coef predict poly
#' @importFrom graphics matplot abline legend image axis
#' @importFrom grDevices colorRampPalette
#' @export
velocity_tna <- function(data,
                         method = "regression",
                         regression_type = "ols",
                         weights = NULL,
                         n_embed = 5L,
                         delta = 1,
                         order = 1L,
                         diff_type = "central",
                         time_col = NULL,
                         id_col = NULL,
                         action_col = NULL,
                         scaling = NULL) {
  # ---- Input validation ----
  method <- match.arg(method, c("regression", "glla", "finite_difference"))
  regression_type <- match.arg(regression_type, c("ols", "beta", "logistic"))
  diff_type <- match.arg(diff_type, c("central", "forward", "backward"))
  stopifnot(
    is.numeric(delta), length(delta) == 1, delta > 0,
    is.numeric(order), length(order) == 1, order >= 1L
  )
  order <- as.integer(order)
  n_embed <- as.integer(n_embed)

  # ---- Logistic regression: works on raw wide data directly ----
  if (method == "regression" && regression_type == "logistic") {
    if (!is.data.frame(data) &&
        !(is.list(data) && !is.null(data$windows))) {
      # Accept wide data frames or try to detect columns
      stop("regression_type = 'logistic' requires a wide-format data frame ",
           "(rows = sequences, columns = time points).", call. = FALSE)
    }
    result_logistic <- .velocity_logistic(data, delta, order)
    return(result_logistic)
  }

  # ---- Extract matrices from input ----
  mat_list <- .detect_and_extract_matrices(
    data,
    time_col = time_col,
    id_col = id_col,
    action_col = action_col,
    scaling = scaling
  )

  # ---- Validate and harmonize ----
  mat_list <- .validate_matrix_list(mat_list)
  n_states <- nrow(mat_list[[1]])
  nodes <- rownames(mat_list[[1]])
  T_len <- length(mat_list)

  if (T_len < 2L) {
    stop("At least 2 time points are required.", call. = FALSE)
  }

  # ---- Weights for regression ----
  if (method == "regression") {
    if (is.null(weights)) {
      weights <- rep(1, T_len)
    } else {
      stopifnot(is.numeric(weights), length(weights) == T_len)
      weights <- weights / mean(weights)  # normalize to mean = 1
    }
  }

  # ---- Method dispatch ----
  if (method == "regression") {
    if (regression_type == "beta" &&
        !requireNamespace("betareg", quietly = TRUE)) {
      stop("Package 'betareg' is required for beta regression. ",
           "Install it with: install.packages('betareg')",
           call. = FALSE)
    }
    result <- .velocity_regression(mat_list, delta, order, regression_type,
                                   weights)

  } else if (method == "glla") {
    if (n_embed < 3L) {
      stop("n_embed must be >= 3 for GLLA.", call. = FALSE)
    }
    if (n_embed %% 2L == 0L) {
      n_embed <- n_embed + 1L
      warning("n_embed must be odd; increased to ", n_embed, ".",
              call. = FALSE)
    }
    if (n_embed > T_len) {
      old_embed <- n_embed
      n_embed <- if (T_len %% 2L == 1L) T_len else T_len - 1L
      if (n_embed < 3L) {
        stop("Too few time points (", T_len,
             ") for GLLA with minimum n_embed=3.", call. = FALSE)
      }
      warning("n_embed (", old_embed, ") > number of time points (", T_len,
              "); reduced to ", n_embed, ".", call. = FALSE)
    }
    if (order >= n_embed) {
      stop("order (", order, ") must be < n_embed (", n_embed, ").",
           call. = FALSE)
    }
    result <- .velocity_glla(mat_list, n_embed, delta, order)

  } else {
    # finite_difference
    if (diff_type == "central" && T_len < 3L) {
      warning("Central differences require >= 3 time points; ",
              "falling back to forward differences.", call. = FALSE)
      diff_type <- "forward"
    }
    result <- .velocity_finite_diff(mat_list, delta, diff_type)
  }

  # ---- Build output ----
  vel_mat <- result$velocity_mean
  vel_series <- result$velocity_series

  acc_mat <- result$acceleration_mean
  acc_series <- result$acceleration_series
  smoothed <- result$smoothed_matrices

  # Extract edges from velocity matrix
  edges <- .extract_edges_from_matrix(vel_mat, directed = TRUE)

  structure(
    list(
      velocity_matrix     = vel_mat,
      velocity_series     = vel_series,
      acceleration_matrix = acc_mat,
      acceleration_series = acc_series,
      smoothed_matrices   = smoothed,
      original_matrices   = mat_list,
      edge_stats          = result$edge_stats,
      nodes               = nodes,
      n_timepoints        = T_len,
      n_embedded          = length(vel_series),
      method              = method,
      regression_type     = if (method == "regression") regression_type
                            else NULL,
      n_embed             = if (method == "glla") n_embed else NULL,
      delta               = delta,
      order               = order,
      directed            = TRUE,
      n_nodes             = n_states,
      n_edges             = nrow(edges),
      edges               = edges
    ),
    class = "tna_velocity"
  )
}


# ---- Regression method ----

#' Regression-based velocity estimation
#'
#' For each edge P(i,j), fits a polynomial regression across all time points.
#' OLS or beta regression. Returns slope (velocity), standard errors, p-values,
#' R-squared, and fitted values.
#'
#' @param mat_list List of square numeric matrices.
#' @param delta Numeric. Time interval.
#' @param order Integer. Polynomial order (1 = linear, 2 = quadratic).
#' @param regression_type Character. "ols" or "beta".
#'
#' @return List with velocity_mean, velocity_series, acceleration_mean,
#'   acceleration_series, smoothed_matrices, edge_stats.
#' @noRd
.velocity_regression <- function(mat_list, delta, order, regression_type,
                                weights = NULL) {
  n_states <- nrow(mat_list[[1]])
  T_len <- length(mat_list)
  n_edges_total <- n_states * n_states
  nodes <- rownames(mat_list[[1]])

  # Time vector

  t_vec <- seq_len(T_len) * delta

  # Stack all matrices: T x n_edges
  all_edges <- do.call(rbind, lapply(mat_list, as.vector))

  # Pre-allocate output
  vel_mat <- matrix(0, n_states, n_states, dimnames = list(nodes, nodes))
  acc_mat <- matrix(0, n_states, n_states, dimnames = list(nodes, nodes))
  fitted_arr <- matrix(0, T_len, n_edges_total)
  vel_arr <- matrix(0, T_len, n_edges_total)
  acc_arr <- matrix(0, T_len, n_edges_total)

  # Edge stats storage
  edge_from <- character(0)
  edge_to <- character(0)
  edge_slope <- numeric(0)
  edge_se <- numeric(0)
  edge_t <- numeric(0)
  edge_p <- numeric(0)
  edge_r2 <- numeric(0)
  edge_baseline <- numeric(0)
  edge_sd_y <- numeric(0)

  sd_t <- stats::sd(t_vec)
  total_span <- max(t_vec) - min(t_vec)

  for (j in seq_len(n_edges_total)) {
    y <- all_edges[, j]
    # as.vector(matrix) is column-major: col varies slowest
    row_idx <- ((j - 1L) %% n_states) + 1L
    col_idx <- ((j - 1L) %/% n_states) + 1L

    if (regression_type == "beta") {
      fit_result <- .fit_beta_edge(y, t_vec, order, weights)
    } else {
      fit_result <- .fit_ols_edge(y, t_vec, order, weights)
    }

    vel_mat[row_idx, col_idx] <- fit_result$slope
    fitted_arr[, j] <- fit_result$fitted
    vel_arr[, j] <- fit_result$velocity_at_t
    acc_arr[, j] <- fit_result$acceleration_at_t

    if (order >= 2L) {
      acc_mat[row_idx, col_idx] <- fit_result$acceleration
    }

    edge_from <- c(edge_from, nodes[row_idx])
    edge_to <- c(edge_to, nodes[col_idx])
    edge_slope <- c(edge_slope, fit_result$slope)
    edge_se <- c(edge_se, fit_result$se)
    edge_t <- c(edge_t, fit_result$t_value)
    edge_p <- c(edge_p, fit_result$p_value)
    edge_r2 <- c(edge_r2, fit_result$r_squared)
    edge_baseline <- c(edge_baseline, mean(y))
    edge_sd_y <- c(edge_sd_y, stats::sd(y))
  }

  # Reshape to lists of matrices
  .arr_to_mat_list <- function(arr) {
    lapply(seq_len(nrow(arr)), function(i) {
      m <- matrix(arr[i, ], n_states, n_states, dimnames = list(nodes, nodes))
      m
    })
  }

  smoothed <- .arr_to_mat_list(fitted_arr)
  vel_series <- .arr_to_mat_list(vel_arr)
  acc_series <- if (order >= 2L) .arr_to_mat_list(acc_arr) else NULL

  # Derived effect-size metrics
  # Standardized coefficient: beta* = slope * sd(t) / sd(y)
  standardized <- ifelse(edge_sd_y > .Machine$double.eps,
                         edge_slope * sd_t / edge_sd_y, 0)
  # Percent change per step: slope / baseline * 100
  pct_change <- ifelse(abs(edge_baseline) > .Machine$double.eps,
                       edge_slope / edge_baseline * 100, 0)
  # Total predicted change over full observation window
  total_change <- edge_slope * total_span

  # Build edge_stats data frame
  edge_stats <- data.frame(
    from = edge_from,
    to = edge_to,
    slope = edge_slope,
    se = edge_se,
    t_value = edge_t,
    p_value = edge_p,
    r_squared = edge_r2,
    standardized = standardized,
    pct_change = pct_change,
    total_change = total_change,
    stringsAsFactors = FALSE
  )
  # Sort by absolute slope
  edge_stats <- edge_stats[order(-abs(edge_stats$slope)), ]
  rownames(edge_stats) <- NULL

  list(
    velocity_mean       = vel_mat,
    velocity_series     = vel_series,
    acceleration_mean   = if (order >= 2L) acc_mat else NULL,
    acceleration_series = acc_series,
    smoothed_matrices   = smoothed,
    edge_stats          = edge_stats
  )
}


#' Logistic regression on individual transitions (no matrix aggregation)
#'
#' For each edge (i -> j), fits: glm(is_j ~ time_step, family = binomial,
#' subset = from == i) across all individual transitions. Uses all raw data
#' for maximum statistical power.
#'
#' @param data Wide-format data frame (rows = sequences, cols = time points).
#' @param delta Numeric. Time interval.
#' @param order Integer. Polynomial order.
#'
#' @return A tna_velocity object.
#' @noRd
.velocity_logistic <- function(data, delta, order) {
  # Detect state columns (T1, T2, ... or V1, V2, ...)
  if (is.data.frame(data)) {
    state_cols <- grep("^[TV][0-9]+$", names(data), value = TRUE)
    if (length(state_cols) < 2L) {
      state_cols <- names(data)[vapply(data, function(x) {
        is.character(x) || is.factor(x)
      }, logical(1))]
    }
    if (length(state_cols) < 2L) {
      stop("Cannot detect time-point columns. Ensure columns are named ",
           "T1, T2, ... or V1, V2, ...", call. = FALSE)
    }
    wide <- data[, state_cols, drop = FALSE]
  } else {
    stop("regression_type = 'logistic' requires a wide-format data frame.",
         call. = FALSE)
  }

  n_cols <- ncol(wide)
  n_steps <- n_cols - 1

  # Build flat transition table: from, to, time_step
  from_vec <- character(0)
  to_vec <- character(0)
  time_vec <- integer(0)

  for (s in seq_len(n_steps)) {
    fr <- wide[[s]]
    tr <- wide[[s + 1L]]
    valid <- !is.na(fr) & !is.na(tr)
    from_vec <- c(from_vec, as.character(fr[valid]))
    to_vec <- c(to_vec, as.character(tr[valid]))
    time_vec <- c(time_vec, rep(s, sum(valid)))
  }

  all_states <- sort(unique(c(from_vec, to_vec)))
  n_states <- length(all_states)
  t_scaled <- time_vec * delta

  cat(sprintf("  Logistic regression on %d individual transitions (%d states)\n",
              length(from_vec), n_states))

  # Per-edge logistic regression
  vel_mat <- matrix(0, n_states, n_states,
                    dimnames = list(all_states, all_states))
  acc_mat <- matrix(0, n_states, n_states,
                    dimnames = list(all_states, all_states))

  edge_from <- character(0)
  edge_to <- character(0)
  edge_slope <- numeric(0)
  edge_se <- numeric(0)
  edge_z <- numeric(0)
  edge_p <- numeric(0)
  edge_n <- integer(0)
  edge_or <- numeric(0)

  for (i_state in all_states) {
    mask_from <- from_vec == i_state
    n_from <- sum(mask_from)
    if (n_from < 5L) {
      # Too few transitions from this state
      for (j_state in all_states) {
        edge_from <- c(edge_from, i_state)
        edge_to <- c(edge_to, j_state)
        edge_slope <- c(edge_slope, 0)
        edge_se <- c(edge_se, 0)
        edge_z <- c(edge_z, 0)
        edge_p <- c(edge_p, 1)
        edge_n <- c(edge_n, 0L)
        edge_or <- c(edge_or, 1)
      }
      next
    }

    t_sub <- t_scaled[mask_from]
    to_sub <- to_vec[mask_from]

    for (j_state in all_states) {
      y_binary <- as.integer(to_sub == j_state)

      # Skip if no variation in outcome
      if (all(y_binary == 0) || all(y_binary == 1)) {
        edge_from <- c(edge_from, i_state)
        edge_to <- c(edge_to, j_state)
        edge_slope <- c(edge_slope, 0)
        edge_se <- c(edge_se, 0)
        edge_z <- c(edge_z, 0)
        edge_p <- c(edge_p, 1)
        edge_n <- c(edge_n, n_from)
        edge_or <- c(edge_or, 1)
        next
      }

      fit <- tryCatch({
        if (order == 1L) {
          stats::glm(y_binary ~ t_sub, family = stats::binomial())
        } else {
          stats::glm(y_binary ~ stats::poly(t_sub, degree = order,
                                             raw = TRUE),
                     family = stats::binomial())
        }
      }, warning = function(w) {
        suppressWarnings(
          if (order == 1L) {
            stats::glm(y_binary ~ t_sub, family = stats::binomial())
          } else {
            stats::glm(y_binary ~ stats::poly(t_sub, degree = order,
                                               raw = TRUE),
                       family = stats::binomial())
          }
        )
      }, error = function(e) NULL)

      if (is.null(fit)) {
        edge_from <- c(edge_from, i_state)
        edge_to <- c(edge_to, j_state)
        edge_slope <- c(edge_slope, 0)
        edge_se <- c(edge_se, 0)
        edge_z <- c(edge_z, 0)
        edge_p <- c(edge_p, 1)
        edge_n <- c(edge_n, n_from)
        edge_or <- c(edge_or, 1)
        next
      }

      sm <- suppressWarnings(summary(fit))
      coefs <- sm$coefficients

      b1 <- coefs[2, 1]      # log-odds slope
      se1 <- coefs[2, 2]
      z1 <- coefs[2, 3]
      p1 <- coefs[2, 4]

      # Average marginal effect: mean(p_hat * (1 - p_hat)) * b1
      p_hat <- stats::predict(fit, type = "response")
      ame <- mean(p_hat * (1 - p_hat)) * b1

      vel_mat[i_state, j_state] <- ame
      if (order >= 2L && nrow(coefs) >= 3) {
        b2 <- coefs[3, 1]
        acc_mat[i_state, j_state] <- mean(p_hat * (1 - p_hat)) * 2 * b2
      }

      edge_from <- c(edge_from, i_state)
      edge_to <- c(edge_to, j_state)
      edge_slope <- c(edge_slope, ame)
      edge_se <- c(edge_se, se1)
      edge_z <- c(edge_z, z1)
      edge_p <- c(edge_p, p1)
      edge_n <- c(edge_n, n_from)
      edge_or <- c(edge_or, exp(b1))
    }
  }

  # Build edge_stats
  edge_stats <- data.frame(
    from = edge_from,
    to = edge_to,
    slope = edge_slope,
    se = edge_se,
    z_value = edge_z,
    p_value = edge_p,
    odds_ratio = edge_or,
    n_transitions = edge_n,
    stringsAsFactors = FALSE
  )

  # Add standardized / pct_change / total_change
  # For logistic, baseline = mean predicted probability
  edge_baseline <- vapply(seq_len(nrow(edge_stats)), function(k) {
    i_st <- edge_stats$from[k]
    j_st <- edge_stats$to[k]
    mask <- from_vec == i_st
    if (sum(mask) == 0) return(0)
    mean(to_vec[mask] == j_st)
  }, numeric(1))

  total_span <- n_steps * delta
  edge_stats$standardized <- ifelse(
    edge_baseline > 0 & edge_baseline < 1,
    edge_stats$slope / sqrt(edge_baseline * (1 - edge_baseline)), 0)
  edge_stats$pct_change <- ifelse(
    abs(edge_baseline) > .Machine$double.eps,
    edge_stats$slope / edge_baseline * 100, 0)
  edge_stats$total_change <- edge_stats$slope * total_span
  edge_stats$r_squared <- NA_real_  # not meaningful for logistic

  edge_stats <- edge_stats[order(-abs(edge_stats$slope)), ]
  rownames(edge_stats) <- NULL

  # Build velocity series (predicted probabilities at each time step)
  vel_series <- lapply(seq_len(n_steps), function(s) {
    vel_mat  # constant for linear model
  })

  # Build smoothed matrices (predicted P at each time step from logistic)
  smoothed <- lapply(seq_len(n_steps), function(s) {
    # For simplicity, use the per-step observed frequencies
    mask <- time_vec == s
    if (sum(mask) == 0) return(matrix(0, n_states, n_states,
                                       dimnames = list(all_states, all_states)))
    m <- matrix(0, n_states, n_states,
                dimnames = list(all_states, all_states))
    fr <- from_vec[mask]
    tr <- to_vec[mask]
    for (st in all_states) {
      st_mask <- fr == st
      if (sum(st_mask) > 0) {
        tab <- table(factor(tr[st_mask], levels = all_states))
        m[st, ] <- tab / sum(tab)
      }
    }
    m
  })

  edges <- .extract_edges_from_matrix(vel_mat, directed = TRUE)

  structure(
    list(
      velocity_matrix     = vel_mat,
      velocity_series     = vel_series,
      acceleration_matrix = if (order >= 2L) acc_mat else NULL,
      acceleration_series = NULL,
      smoothed_matrices   = smoothed,
      original_matrices   = smoothed,
      edge_stats          = edge_stats,
      nodes               = all_states,
      n_timepoints        = n_steps,
      n_embedded          = n_steps,
      method              = "regression",
      regression_type     = "logistic",
      n_embed             = NULL,
      delta               = delta,
      order               = order,
      directed            = TRUE,
      n_nodes             = n_states,
      n_edges             = nrow(edges),
      edges               = edges
    ),
    class = "tna_velocity"
  )
}


#' Fit OLS regression for a single edge
#' @noRd
.fit_ols_edge <- function(y, t_vec, order, weights = NULL) {
  # Handle constant edges
  if (stats::var(y) < .Machine$double.eps) {
    n <- length(y)
    return(list(
      slope = 0, se = 0, t_value = 0, p_value = 1, r_squared = 0,
      acceleration = 0,
      fitted = rep(mean(y), n),
      velocity_at_t = rep(0, n),
      acceleration_at_t = rep(0, n)
    ))
  }

  if (order == 1L) {
    fit <- stats::lm(y ~ t_vec, weights = weights)
    slope <- stats::coef(fit)[2]
    sm <- suppressWarnings(summary(fit))
    se <- sm$coefficients[2, 2]
    t_value <- sm$coefficients[2, 3]
    p_value <- sm$coefficients[2, 4]
    r_squared <- sm$r.squared
    fitted_vals <- stats::predict(fit)
    vel_at_t <- rep(slope, length(t_vec))
    acc_at_t <- rep(0, length(t_vec))
    acceleration <- 0
  } else {
    fit <- stats::lm(y ~ stats::poly(t_vec, degree = order, raw = TRUE),
                     weights = weights)
    sm <- suppressWarnings(summary(fit))
    coeffs <- stats::coef(fit)
    # slope = b1 (linear coefficient)
    slope <- coeffs[2]
    se <- sm$coefficients[2, 2]
    t_value <- sm$coefficients[2, 3]
    p_value <- sm$coefficients[2, 4]
    r_squared <- sm$r.squared
    fitted_vals <- stats::predict(fit)

    # velocity at each t: dy/dt = b1 + 2*b2*t (+ 3*b3*t^2 ...)
    vel_at_t <- rep(coeffs[2], length(t_vec))
    if (order >= 2L && length(coeffs) >= 3) {
      vel_at_t <- vel_at_t + 2 * coeffs[3] * t_vec
    }
    # acceleration = 2*b2 (constant for quadratic)
    acceleration <- if (length(coeffs) >= 3) 2 * coeffs[3] else 0
    acc_at_t <- rep(acceleration, length(t_vec))
  }

  list(
    slope = unname(slope),
    se = unname(se),
    t_value = unname(t_value),
    p_value = unname(p_value),
    r_squared = unname(r_squared),
    acceleration = unname(acceleration),
    fitted = unname(fitted_vals),
    velocity_at_t = unname(vel_at_t),
    acceleration_at_t = unname(acc_at_t)
  )
}


#' Fit beta regression for a single edge
#'
#' Beta regression is appropriate for bounded (0 to 1) data. Values at exactly
#' 0 or 1 are squeezed using the Smithson-Verkuilen transformation.
#'
#' @noRd
.fit_beta_edge <- function(y, t_vec, order, weights = NULL) {
  n <- length(y)

  # Handle constant or all-zero edges
  if (stats::var(y) < .Machine$double.eps || all(y == 0)) {
    return(list(
      slope = 0, se = 0, t_value = 0, p_value = 1, r_squared = 0,
      acceleration = 0,
      fitted = rep(mean(y), n),
      velocity_at_t = rep(0, n),
      acceleration_at_t = rep(0, n)
    ))
  }

  # Smithson-Verkuilen squeeze: (y*(n-1) + 0.5) / n
  y_sq <- (y * (n - 1) + 0.5) / n

  if (order == 1L) {
    formula <- y_sq ~ t_vec
  } else {
    formula <- y_sq ~ stats::poly(t_vec, degree = order, raw = TRUE)
  }

  fit <- tryCatch(
    betareg::betareg(formula, weights = weights),
    error = function(e) NULL
  )

  # Fall back to OLS if beta regression fails
  if (is.null(fit)) {
    return(.fit_ols_edge(y, t_vec, order, weights))
  }

  sm <- summary(fit)
  coeffs <- stats::coef(fit)

  # Mean model coefficients (before precision model)
  mean_coeffs <- coeffs[seq_len(order + 1L)]

  slope <- mean_coeffs[2]
  se <- sm$coefficients$mean[2, 2]
  t_value <- sm$coefficients$mean[2, 3]
  p_value <- sm$coefficients$mean[2, 4]

  # Pseudo R-squared: squared correlation between fitted and observed
  fitted_link <- stats::predict(fit, type = "response")
  r_squared <- stats::cor(y, fitted_link)^2

  # Velocity on response scale (marginal effect via chain rule)
  # For logit link: dy/dt = mu*(1-mu) * d(eta)/dt
  fitted_mu <- stats::predict(fit, type = "response")
  deta_dt <- rep(mean_coeffs[2], n)
  if (order >= 2L && length(mean_coeffs) >= 3) {
    deta_dt <- deta_dt + 2 * mean_coeffs[3] * t_vec
  }
  vel_at_t <- fitted_mu * (1 - fitted_mu) * deta_dt

  # Mean velocity = mean of marginal effects
  mean_vel <- mean(vel_at_t)

  # Acceleration
  if (order >= 2L && length(mean_coeffs) >= 3) {
    # d2y/dt2 on response scale (approximate)
    d2eta_dt2 <- 2 * mean_coeffs[3]
    acc_at_t <- fitted_mu * (1 - fitted_mu) * d2eta_dt2 +
      (1 - 2 * fitted_mu) * fitted_mu * (1 - fitted_mu) * deta_dt^2
    acceleration <- mean(acc_at_t)
  } else {
    acc_at_t <- rep(0, n)
    acceleration <- 0
  }

  list(
    slope = unname(mean_vel),
    se = unname(se),
    t_value = unname(t_value),
    p_value = unname(p_value),
    r_squared = unname(r_squared),
    acceleration = unname(acceleration),
    fitted = unname(fitted_mu),
    velocity_at_t = unname(vel_at_t),
    acceleration_at_t = unname(acc_at_t)
  )
}


# ---- GLLA method ----

#' Compute GLLA weight matrix
#'
#' @param n_embed Integer. Embedding dimension.
#' @param delta Numeric. Time interval.
#' @param order Integer. Maximum derivative order.
#'
#' @return Matrix of size (order+1) x n_embed.
#' @noRd
.glla_weight_matrix <- function(n_embed, delta, order) {
  t_j <- (seq_len(n_embed) - (n_embed + 1) / 2) * delta
  L <- outer(t_j, seq(0, order), "^")
  W <- solve(crossprod(L), t(L))
  # Scale by k! to get actual derivatives
  for (k in seq(0, order)) {
    W[k + 1L, ] <- W[k + 1L, ] * factorial(k)
  }
  W
}


#' GLLA velocity computation
#' @noRd
.velocity_glla <- function(mat_list, n_embed, delta, order) {
  n_states <- nrow(mat_list[[1]])
  T_len <- length(mat_list)
  n_edges <- n_states * n_states

  all_edges <- do.call(rbind, lapply(mat_list, as.vector))
  W <- .glla_weight_matrix(n_embed, delta, order)
  n_emb <- T_len - n_embed + 1L

  deriv_arrays <- lapply(seq(0, order), function(o) {
    matrix(0, nrow = n_emb, ncol = n_edges)
  })

  for (j in seq_len(n_edges)) {
    emb <- stats::embed(all_edges[, j], n_embed)
    emb <- emb[, rev(seq_len(n_embed)), drop = FALSE]
    derivs <- emb %*% t(W)
    for (o in seq(0, order)) {
      deriv_arrays[[o + 1L]][, j] <- derivs[, o + 1L]
    }
  }

  .reshape_to_mat_list <- function(arr, n_states) {
    lapply(seq_len(nrow(arr)), function(i) {
      m <- matrix(arr[i, ], nrow = n_states, ncol = n_states)
      dimnames(m) <- dimnames(mat_list[[1]])
      m
    })
  }

  smoothed <- .reshape_to_mat_list(deriv_arrays[[1]], n_states)
  vel_series <- .reshape_to_mat_list(deriv_arrays[[2]], n_states)
  vel_mean <- Reduce("+", vel_series) / length(vel_series)

  acc_mean <- NULL
  acc_series <- NULL
  if (order >= 2L) {
    acc_series <- .reshape_to_mat_list(deriv_arrays[[3]], n_states)
    acc_mean <- Reduce("+", acc_series) / length(acc_series)
  }

  list(
    velocity_mean       = vel_mean,
    velocity_series     = vel_series,
    acceleration_mean   = acc_mean,
    acceleration_series = acc_series,
    smoothed_matrices   = smoothed,
    edge_stats          = NULL
  )
}


# ---- Finite difference method ----

#' Finite difference velocity computation
#' @noRd
.velocity_finite_diff <- function(mat_list, delta, diff_type) {
  n_states <- nrow(mat_list[[1]])
  T_len <- length(mat_list)

  all_edges <- do.call(rbind, lapply(mat_list, as.vector))

  if (diff_type == "central") {
    vel_arr <- (all_edges[3:T_len, , drop = FALSE] -
                  all_edges[1:(T_len - 2L), , drop = FALSE]) / (2 * delta)
  } else if (diff_type == "forward") {
    vel_arr <- diff(all_edges) / delta
  } else {
    vel_arr <- diff(all_edges) / delta
  }

  vel_series <- lapply(seq_len(nrow(vel_arr)), function(i) {
    m <- matrix(vel_arr[i, ], nrow = n_states, ncol = n_states)
    dimnames(m) <- dimnames(mat_list[[1]])
    m
  })

  vel_mean <- Reduce("+", vel_series) / length(vel_series)

  list(
    velocity_mean       = vel_mean,
    velocity_series     = vel_series,
    acceleration_mean   = NULL,
    acceleration_series = NULL,
    smoothed_matrices   = NULL,
    edge_stats          = NULL
  )
}


# ---- Input detection & validation helpers ----

#' Detect input format and extract transition matrices
#' @noRd
.detect_and_extract_matrices <- function(data,
                                         time_col = NULL,
                                         id_col = NULL,
                                         action_col = NULL,
                                         scaling = NULL) {
  if (is.list(data) && !is.data.frame(data) && !is.null(data$windows)) {
    if (length(data$windows) > 0 && !is.null(data$windows[[1]]$weights)) {
      return(lapply(data$windows, function(w) w$weights))
    }
  }

  if (is.list(data) && !is.data.frame(data) && length(data) > 0 &&
      is.matrix(data[[1]])) {
    return(data)
  }

  if (is.data.frame(data)) {
    if (is.null(time_col)) {
      stop("'time_col' is required when 'data' is a data frame.", call. = FALSE)
    }
    return(.matrices_from_raw_data(
      data,
      time_col = time_col,
      id_col = id_col,
      action_col = action_col,
      scaling = scaling
    ))
  }

  stop("Unrecognised input format. Provide a list of matrices, ",
       "cograph::tna_windows() output, or a data frame with 'time_col'.",
       call. = FALSE)
}


#' Validate and harmonize a list of matrices
#' @noRd
.validate_matrix_list <- function(mat_list) {
  if (!is.list(mat_list) || length(mat_list) < 1L) {
    stop("mat_list must be a non-empty list.", call. = FALSE)
  }

  is_mat <- vapply(mat_list, is.matrix, logical(1))
  if (!all(is_mat)) {
    stop("All elements must be matrices.", call. = FALSE)
  }

  dims <- vapply(mat_list, function(m) dim(m), integer(2))
  if (!all(dims[1, ] == dims[2, ])) {
    stop("All matrices must be square.", call. = FALSE)
  }
  if (length(unique(dims[1, ])) > 1L) {
    stop("All matrices must have the same dimensions.", call. = FALSE)
  }

  is_num <- vapply(mat_list, is.numeric, logical(1))
  if (!all(is_num)) {
    stop("All matrices must be numeric.", call. = FALSE)
  }

  n <- nrow(mat_list[[1]])
  default_names <- paste0("S", seq_len(n))

  mat_list <- lapply(mat_list, function(m) {
    if (is.null(rownames(m))) rownames(m) <- default_names
    if (is.null(colnames(m))) colnames(m) <- default_names
    m
  })

  mat_list
}


#' Build transition matrices from raw long-format data
#' @noRd
.matrices_from_raw_data <- function(data,
                                    time_col,
                                    id_col = NULL,
                                    action_col = NULL,
                                    scaling = NULL) {
  stopifnot(is.data.frame(data))

  if (!time_col %in% names(data)) {
    stop("time_col '", time_col, "' not found in data.", call. = FALSE)
  }

  if (is.null(action_col)) action_col <- "Action"
  if (!action_col %in% names(data)) {
    stop("action_col '", action_col, "' not found in data.", call. = FALSE)
  }

  time_vals <- data[[time_col]]
  time_levels <- sort(unique(time_vals))

  if (length(time_levels) < 2L) {
    stop("At least 2 unique time windows are required.", call. = FALSE)
  }

  mat_list <- lapply(time_levels, function(tw) {
    sub <- data[time_vals == tw, , drop = FALSE]
    .count_transitions_long(sub, action = action_col, id = id_col,
                            time = time_col)
  })

  all_states <- sort(unique(unlist(lapply(mat_list, rownames))))
  n_states <- length(all_states)

  mat_list <- lapply(mat_list, function(m) {
    full <- matrix(0, nrow = n_states, ncol = n_states,
                   dimnames = list(all_states, all_states))
    shared_rows <- intersect(rownames(m), all_states)
    shared_cols <- intersect(colnames(m), all_states)
    full[shared_rows, shared_cols] <- m[shared_rows, shared_cols]
    full
  })

  mat_list <- lapply(mat_list, function(m) {
    rs <- rowSums(m)
    nonzero <- rs > 0
    m[nonzero, ] <- m[nonzero, ] / rs[nonzero]
    storage.mode(m) <- "double"
    m
  })

  if (!is.null(scaling)) {
    mat_list <- lapply(mat_list, function(m) .apply_scaling(m, scaling))
  }

  mat_list
}


# ---- S3 methods ----

#' Print Method for Velocity TNA Object
#'
#' @param x A \code{tna_velocity} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.tna_velocity <- function(x, ...) {
  cat("Transition Network Velocity [directed]\n")

  if (x$method == "regression") {
    reg_label <- if (!is.null(x$regression_type) &&
                     x$regression_type == "beta") "Beta"
                 else if (!is.null(x$regression_type) &&
                     x$regression_type == "logistic") "Logistic"
                 else "OLS"
    cat(sprintf("  Method: %s Regression (order=%d)\n", reg_label, x$order))
  } else if (x$method == "glla") {
    cat(sprintf("  Method: GLLA (n_embed=%d, delta=%g, order=%d)\n",
                x$n_embed, x$delta, x$order))
  } else {
    cat(sprintf("  Method: Finite Difference (delta=%g)\n", x$delta))
  }

  cat(sprintf("  Time points: %d\n", x$n_timepoints))
  cat(sprintf("  Nodes: %d  |  Velocity edges: %d\n",
              x$n_nodes, x$n_edges))

  # Top edges with effect sizes for regression
  if (x$method == "regression" && !is.null(x$edge_stats) &&
      "standardized" %in% names(x$edge_stats)) {
    es <- x$edge_stats[x$edge_stats$from != x$edge_stats$to, ]
    if (nrow(es) > 0) {
      top_pos <- es[which.max(es$slope), ]
      top_neg <- es[which.min(es$slope), ]
      if (top_pos$slope > 0) {
        cat(sprintf(
          "\n  Top strengthening:  %s -> %s: %+.4f/t (beta=%.3f, %+.1f%%/t, total=%+.3f)\n",
          top_pos$from, top_pos$to, top_pos$slope,
          top_pos$standardized, top_pos$pct_change, top_pos$total_change))
      }
      if (top_neg$slope < 0) {
        cat(sprintf(
          "  Top weakening:      %s -> %s: %+.4f/t (beta=%.3f, %+.1f%%/t, total=%+.3f)\n",
          top_neg$from, top_neg$to, top_neg$slope,
          top_neg$standardized, top_neg$pct_change, top_neg$total_change))
      }
    }
    sig <- es[es$p_value < 0.05, ]
    cat(sprintf("\n  Significant edges (p < .05): %d\n", nrow(sig)))
  } else {
    # Non-regression methods
    vel <- x$velocity_matrix
    diag(vel) <- 0
    if (any(vel != 0)) {
      edges_df <- .extract_edges_from_matrix(vel, directed = TRUE)
      if (nrow(edges_df) > 0) {
        top_pos <- edges_df[which.max(edges_df$weight), ]
        top_neg <- edges_df[which.min(edges_df$weight), ]
        if (top_pos$weight > 0) {
          cat(sprintf("\n  Top strengthening:  %s -> %s: %+.4f/t\n",
                      top_pos$from, top_pos$to, top_pos$weight))
        }
        if (top_neg$weight < 0) {
          cat(sprintf("  Top weakening:      %s -> %s: %+.4f/t\n",
                      top_neg$from, top_neg$to, top_neg$weight))
        }
      }
    }
  }

  invisible(x)
}


#' Summary Method for Velocity TNA Object
#'
#' @param object A \code{tna_velocity} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns: \code{from}, \code{to}, \code{velocity},
#'   \code{abs_velocity}, \code{direction}, \code{consistency}. For regression,
#'   also includes \code{se}, \code{p_value}, \code{r_squared}.
#'
#' @export
summary.tna_velocity <- function(object, ...) {
  vel <- object$velocity_matrix
  n <- nrow(vel)
  nodes <- rownames(vel)

  # All directed edges (excluding self-loops)
  idx <- which(row(vel) != col(vel), arr.ind = TRUE)

  from <- nodes[idx[, 1]]
  to <- nodes[idx[, 2]]
  velocity <- vel[idx]
  abs_velocity <- abs(velocity)

  direction <- ifelse(
    abs_velocity < 1e-10, "stable",
    ifelse(velocity > 0, "strengthening", "weakening")
  )

  # Consistency
  consistency <- vapply(seq_len(nrow(idx)), function(k) {
    i <- idx[k, 1]
    j <- idx[k, 2]
    mean_vel <- vel[i, j]
    if (abs(mean_vel) < 1e-10) return(1)
    signs <- vapply(object$velocity_series, function(m) {
      sign(m[i, j])
    }, numeric(1))
    mean(signs == sign(mean_vel))
  }, numeric(1))

  result <- data.frame(
    from = from,
    to = to,
    velocity = velocity,
    abs_velocity = abs_velocity,
    direction = direction,
    consistency = consistency,
    stringsAsFactors = FALSE
  )

  # Merge regression stats if available
  if (object$method == "regression" && !is.null(object$edge_stats)) {
    es <- object$edge_stats
    # Exclude self-loops from edge_stats for merge
    es <- es[es$from != es$to, ]
    merge_cols <- c("from", "to", "se", "p_value", "r_squared",
                    "standardized", "pct_change", "total_change")
    merge_cols <- intersect(merge_cols, names(es))
    result <- merge(result, es[, c("from", "to", merge_cols), drop = FALSE],
                    by = c("from", "to"), all.x = TRUE)
  }

  result <- result[order(-result$abs_velocity), ]
  rownames(result) <- NULL
  result
}


#' Plot Method for Velocity TNA Object
#'
#' @param x A \code{tna_velocity} object.
#' @param type Character. Plot type: \code{"network"} (default),
#'   \code{"series"}, or \code{"heatmap"}.
#' @param time Integer or NULL. For \code{type = "network"}: which time point
#'   to plot (NULL = mean velocity).
#' @param n_top Integer. For \code{type = "series"}: number of top edges to
#'   show. Default: \code{10}.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @importFrom graphics matplot abline legend image axis mtext
#' @importFrom grDevices colorRampPalette
#' @export
plot.tna_velocity <- function(x, type = "network", time = NULL,
                              n_top = 10, ...) {
  type <- match.arg(type, c("network", "series", "heatmap"))

  if (type == "network") {
    .plot_velocity_network(x, time = time, ...)
  } else if (type == "series") {
    .plot_velocity_series(x, n_top = n_top, ...)
  } else {
    .plot_velocity_heatmap(x, time = time, ...)
  }
}


#' @noRd
.plot_velocity_network <- function(x, time = NULL, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Package 'cograph' is required for network plots. ",
         "Install it with: install.packages('cograph')", call. = FALSE)
  }

  if (is.null(time)) {
    mat <- x$velocity_matrix
    title_str <- "Mean Velocity Network"
  } else {
    stopifnot(is.numeric(time), length(time) == 1,
              time >= 1, time <= x$n_embedded)
    mat <- x$velocity_series[[time]]
    title_str <- sprintf("Velocity Network (t=%d)", time)
  }

  node_cols <- .node_colors(x$n_nodes)

  dots <- list(
    x = mat,
    directed = TRUE,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    edge_positive_color = "#2E7D32",
    edge_negative_color = "#C62828",
    title = title_str,
    ...
  )

  do.call(cograph::splot, dots)
}


#' @noRd
.plot_velocity_series <- function(x, n_top = 10, ...) {
  n_emb <- x$n_embedded
  nodes <- x$nodes

  idx <- which(row(x$velocity_matrix) != col(x$velocity_matrix),
               arr.ind = TRUE)
  edge_labels <- paste0(nodes[idx[, 1]], "->", nodes[idx[, 2]])

  # For regression, plot observed edge values + fitted line
  if (x$method == "regression" && !is.null(x$original_matrices)) {
    # Plot observed edge time series with regression fit
    obs_mat <- do.call(rbind, lapply(x$original_matrices, function(m) m[idx]))
    fit_mat <- do.call(rbind, lapply(x$smoothed_matrices, function(m) m[idx]))
    colnames(obs_mat) <- edge_labels
    colnames(fit_mat) <- edge_labels

    # Select top N by absolute slope
    mean_abs_vel <- abs(x$velocity_matrix[idx])
    top_idx <- utils::head(order(-mean_abs_vel), min(n_top, length(mean_abs_vel)))

    obs_top <- obs_mat[, top_idx, drop = FALSE]
    fit_top <- fit_mat[, top_idx, drop = FALSE]

    mean_vel <- x$velocity_matrix[idx[top_idx, , drop = FALSE]]
    colors <- ifelse(mean_vel > 0, "#2E7D32", "#C62828")

    ylim <- range(c(obs_top, fit_top), na.rm = TRUE)

    graphics::matplot(
      seq_len(n_emb), obs_top,
      type = "p", pch = 16, col = adjustcolor(colors, 0.4), cex = 0.8,
      xlab = "Time point", ylab = "Transition probability",
      main = "Edge Trajectories with Regression Fit",
      ylim = ylim,
      ...
    )
    graphics::matlines(seq_len(n_emb), fit_top,
                       lty = 1, col = colors, lwd = 2)
    graphics::legend(
      "topright",
      legend = colnames(obs_top),
      col = colors, lty = 1, lwd = 2,
      cex = 0.65, bg = "white"
    )
  } else {
    series_mat <- do.call(rbind, lapply(x$velocity_series, function(m) m[idx]))
    colnames(series_mat) <- edge_labels

    mean_abs <- colMeans(abs(series_mat))
    top_idx <- utils::head(order(-mean_abs), min(n_top, ncol(series_mat)))
    series_top <- series_mat[, top_idx, drop = FALSE]

    mean_vel <- colMeans(series_top)
    colors <- ifelse(mean_vel > 0, "#2E7D32", "#C62828")

    graphics::matplot(
      seq_len(n_emb), series_top,
      type = "l", lty = 1, col = colors, lwd = 1.5,
      xlab = "Embedded time point", ylab = "Velocity",
      main = "Velocity Time Series (Top Edges)",
      ...
    )
    graphics::abline(h = 0, lty = 2, col = "gray50")
    graphics::legend(
      "topright",
      legend = colnames(series_top),
      col = colors, lty = 1, lwd = 1.5,
      cex = 0.7, bg = "white"
    )
  }
}


#' @noRd
.plot_velocity_heatmap <- function(x, time = NULL, ...) {
  if (is.null(time)) {
    mat <- x$velocity_matrix
    title_str <- "Mean Velocity Heatmap"
  } else {
    stopifnot(is.numeric(time), length(time) == 1,
              time >= 1, time <= x$n_embedded)
    mat <- x$velocity_series[[time]]
    title_str <- sprintf("Velocity Heatmap (t=%d)", time)
  }

  n <- nrow(mat)
  nodes <- rownames(mat)

  pal <- grDevices::colorRampPalette(
    c("#2166AC", "#F7F7F7", "#B2182B")
  )(101)

  max_abs <- max(abs(mat))
  if (max_abs == 0) max_abs <- 1

  graphics::image(
    seq_len(n), seq_len(n),
    t(mat[rev(seq_len(n)), ]),
    col = pal,
    zlim = c(-max_abs, max_abs),
    axes = FALSE,
    xlab = "To", ylab = "From",
    main = title_str,
    ...
  )
  graphics::axis(1, at = seq_len(n), labels = nodes, las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(n), labels = rev(nodes), las = 2,
                 cex.axis = 0.8)
}
