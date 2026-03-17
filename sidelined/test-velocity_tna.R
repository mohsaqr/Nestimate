# ---- Synthetic data generators ----

.make_linear_matrices <- function(n_states = 3, T_len = 10, slope = 0.01,
                                  seed = 42) {
  set.seed(seed)
  base <- matrix(1 / n_states, n_states, n_states)
  states <- paste0("S", seq_len(n_states))
  dimnames(base) <- list(states, states)

  lapply(seq_len(T_len), function(t) {
    m <- base
    m[1, 2] <- base[1, 2] + slope * (t - 1)
    m[1, 1] <- base[1, 1] - slope * (t - 1)
    m
  })
}

.make_constant_matrices <- function(n_states = 3, T_len = 10, seed = 42) {
  set.seed(seed)
  base <- matrix(runif(n_states * n_states), n_states, n_states)
  rs <- rowSums(base)
  base <- base / rs
  states <- paste0("S", seq_len(n_states))
  dimnames(base) <- list(states, states)
  rep(list(base), T_len)
}

.make_quadratic_matrices <- function(n_states = 3, T_len = 15, accel = 0.001,
                                     seed = 42) {
  set.seed(seed)
  base <- matrix(1 / n_states, n_states, n_states)
  states <- paste0("S", seq_len(n_states))
  dimnames(base) <- list(states, states)

  lapply(seq_len(T_len), function(t) {
    m <- base
    m[1, 2] <- base[1, 2] + 0.5 * accel * (t - 1)^2
    m[1, 1] <- base[1, 1] - 0.5 * accel * (t - 1)^2
    m
  })
}


# ---- Input validation tests ----

test_that("velocity_tna errors on fewer than 2 time points", {
  mats <- list(matrix(c(0.5, 0.5, 0.3, 0.7), 2, 2))
  expect_error(velocity_tna(mats), "At least 2 time points")
})

test_that("velocity_tna errors on non-matrix list elements", {
  bad <- list(data.frame(a = 1), data.frame(b = 2))
  expect_error(velocity_tna(bad), "matrices")
})

test_that("velocity_tna errors on non-square matrices", {
  mats <- list(matrix(1:6, 2, 3), matrix(1:6, 2, 3))
  expect_error(velocity_tna(mats), "square")
})

test_that("velocity_tna errors on mismatched dimensions", {
  m1 <- matrix(1, 2, 2)
  m2 <- matrix(1, 3, 3)
  expect_error(velocity_tna(list(m1, m2)), "same dimensions")
})

test_that("velocity_tna errors on n_embed < 3 for GLLA", {
  mats <- .make_constant_matrices(T_len = 5)
  expect_error(velocity_tna(mats, method = "glla", n_embed = 2L),
               "n_embed must be >= 3")
})

test_that("velocity_tna errors on order >= n_embed", {
  mats <- .make_constant_matrices(T_len = 10)
  expect_error(velocity_tna(mats, method = "glla", n_embed = 3L, order = 3L),
               "order.*must be < n_embed")
})

test_that("velocity_tna warns and adjusts even n_embed", {
  mats <- .make_constant_matrices(T_len = 10)
  expect_warning(
    vel <- velocity_tna(mats, method = "glla", n_embed = 4L),
    "n_embed must be odd"
  )
  expect_equal(vel$n_embed, 5L)
})

test_that("velocity_tna warns and reduces n_embed > T", {
  mats <- .make_constant_matrices(T_len = 5)
  expect_warning(
    vel <- velocity_tna(mats, method = "glla", n_embed = 11L),
    "reduced to"
  )
  expect_equal(vel$n_embed, 5L)
})

test_that("velocity_tna errors on data.frame without time_col", {
  df <- data.frame(Action = letters[1:10], Time = 1:10)
  expect_error(velocity_tna(df), "time_col.*required")
})

test_that("central diff falls back to forward with T=2", {
  m1 <- matrix(c(0.5, 0.5, 0.3, 0.7), 2, 2,
               dimnames = list(c("A", "B"), c("A", "B")))
  m2 <- matrix(c(0.4, 0.6, 0.2, 0.8), 2, 2,
               dimnames = list(c("A", "B"), c("A", "B")))
  expect_warning(
    vel <- velocity_tna(list(m1, m2), method = "finite_difference",
                        diff_type = "central"),
    "falling back to forward"
  )
  expect_equal(vel$n_embedded, 1L)
})


# ---- Regression tests ----

test_that("OLS regression recovers linear slope exactly", {
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "regression", order = 1L)

  expect_s3_class(vel, "tna_velocity")
  expect_equal(vel$method, "regression")
  expect_equal(vel$velocity_matrix[1, 2], slope, tolerance = 1e-6)
  expect_equal(vel$velocity_matrix[1, 1], -slope, tolerance = 1e-6)
  # Other edges should be 0
  expect_true(all(abs(vel$velocity_matrix[2:3, ]) < 1e-10))
})

test_that("OLS regression returns edge_stats with correct columns", {
  mats <- .make_linear_matrices(T_len = 10, slope = 0.01)
  vel <- velocity_tna(mats, method = "regression")

  expect_false(is.null(vel$edge_stats))
  expect_true(is.data.frame(vel$edge_stats))
  expect_true(all(c("from", "to", "slope", "se", "t_value",
                     "p_value", "r_squared", "standardized",
                     "pct_change", "total_change") %in% names(vel$edge_stats)))
  # 3x3 = 9 edges (including diagonal)
  expect_equal(nrow(vel$edge_stats), 9)
})

test_that("OLS regression p-values: trending edge significant, constant not", {
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  # S1->S2 should be highly significant
  p_12 <- es[es$from == "S1" & es$to == "S2", "p_value"]
  expect_true(p_12 < 0.001)

  # S2->S3 (constant) should NOT be significant
  p_23 <- es[es$from == "S2" & es$to == "S3", "p_value"]
  expect_true(p_23 > 0.05 || abs(es[es$from == "S2" & es$to == "S3", "slope"]) < 1e-10)
})

test_that("OLS regression R-squared: perfect linear trend = 1.0", {
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  r2_12 <- es[es$from == "S1" & es$to == "S2", "r_squared"]
  expect_equal(r2_12, 1.0, tolerance = 1e-8)
})

test_that("OLS regression on constant series gives velocity = 0", {
  mats <- .make_constant_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "regression")
  expect_true(all(abs(vel$velocity_matrix) < 1e-10))
})

test_that("OLS regression returns smoothed_matrices (fitted values)", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "regression")

  expect_false(is.null(vel$smoothed_matrices))
  expect_equal(length(vel$smoothed_matrices), 10)
  expect_equal(dim(vel$smoothed_matrices[[1]]), c(3, 3))
})

test_that("OLS regression order=2 recovers quadratic acceleration", {
  accel <- 0.002
  mats <- .make_quadratic_matrices(n_states = 3, T_len = 15, accel = accel)
  vel <- velocity_tna(mats, method = "regression", order = 2L)

  expect_false(is.null(vel$acceleration_matrix))
  expect_equal(vel$acceleration_matrix[1, 2], accel, tolerance = 1e-4)
  expect_true(all(abs(vel$acceleration_matrix[2:3, ]) < 1e-10))
})

test_that("Regression summary includes all effect-size columns", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "regression")
  s <- summary(vel)

  expect_true("p_value" %in% names(s))
  expect_true("r_squared" %in% names(s))
  expect_true("se" %in% names(s))
  expect_true("standardized" %in% names(s))
  expect_true("pct_change" %in% names(s))
  expect_true("total_change" %in% names(s))
})

test_that("Regression print shows significant edge count", {
  mats <- .make_linear_matrices(T_len = 10, slope = 0.02)
  vel <- velocity_tna(mats, method = "regression")
  expect_output(print(vel), "Significant edges")
  expect_output(print(vel), "Regression")
})

test_that("Regression series plot shows observed + fitted", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "regression")
  expect_no_error(plot(vel, type = "series"))
})


# ---- Effect-size metrics tests ----

test_that("Standardized coefficient = 1.0 for perfect linear trend", {
  # For a perfect linear trend, cor(y, t) = +/-1, so |beta*| = 1
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  beta_12 <- es[es$from == "S1" & es$to == "S2", "standardized"]
  expect_equal(abs(beta_12), 1.0, tolerance = 1e-8)

  # Constant edges should have standardized = 0
  beta_23 <- es[es$from == "S2" & es$to == "S3", "standardized"]
  expect_equal(beta_23, 0)
})

test_that("Standardized coefficient sign matches slope sign", {
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  # S1->S2 has positive slope -> positive beta*
  beta_12 <- es[es$from == "S1" & es$to == "S2", "standardized"]
  expect_true(beta_12 > 0)
  # S1->S1 has negative slope -> negative beta*
  beta_11 <- es[es$from == "S1" & es$to == "S1", "standardized"]
  expect_true(beta_11 < 0)
})

test_that("Standardized coefficient formula: slope * sd(t) / sd(y)", {
  slope <- 0.015
  T_len <- 12
  mats <- .make_linear_matrices(n_states = 3, T_len = T_len, slope = slope)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  beta_12 <- es[es$from == "S1" & es$to == "S2", "standardized"]

  # Manual calculation
  t_vec <- seq_len(T_len) * 1  # delta = 1
  y <- vapply(mats, function(m) m[1, 2], numeric(1))
  expected_beta <- slope * sd(t_vec) / sd(y)

  expect_equal(beta_12, expected_beta, tolerance = 1e-8)
})

test_that("Percent change per step: slope / baseline * 100", {
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  pct_12 <- es[es$from == "S1" & es$to == "S2", "pct_change"]

  # Baseline for S1->S2: mean across time of the edge values
  y <- vapply(mats, function(m) m[1, 2], numeric(1))
  expected_pct <- slope / mean(y) * 100

  expect_equal(pct_12, expected_pct, tolerance = 1e-6)
})

test_that("Percent change is 0 for constant edges", {
  mats <- .make_linear_matrices(n_states = 3, T_len = 10, slope = 0.01)
  vel <- velocity_tna(mats, method = "regression")

  es <- vel$edge_stats
  # S2->S3 is constant (slope=0), so pct_change = 0
  pct_23 <- es[es$from == "S2" & es$to == "S3", "pct_change"]
  expect_equal(pct_23, 0)
})

test_that("Total change = slope * time span", {
  slope <- 0.02
  T_len <- 20
  delta <- 1
  mats <- .make_linear_matrices(n_states = 3, T_len = T_len, slope = slope)
  vel <- velocity_tna(mats, method = "regression", delta = delta)

  es <- vel$edge_stats
  total_12 <- es[es$from == "S1" & es$to == "S2", "total_change"]

  # Total span = max(t_vec) - min(t_vec) = T_len*delta - 1*delta = (T_len-1)*delta
  expected_total <- slope * (T_len - 1) * delta
  expect_equal(total_12, expected_total, tolerance = 1e-8)
})

test_that("Total change with delta != 1", {
  # With delta=0.5, t_vec = 0.5, 1.0, ..., 5.0
  # Regression slope = slope/delta (steeper), span = (T-1)*delta (shorter)
  # Total = (slope/delta) * (T-1)*delta = slope * (T-1), same as delta=1
  slope <- 0.02
  T_len <- 10
  delta <- 0.5
  mats <- .make_linear_matrices(n_states = 3, T_len = T_len, slope = slope)
  vel <- velocity_tna(mats, method = "regression", delta = delta)

  es <- vel$edge_stats
  total_12 <- es[es$from == "S1" & es$to == "S2", "total_change"]
  # Total change in y is independent of time scaling
  expected_total <- slope * (T_len - 1)
  expect_equal(total_12, expected_total, tolerance = 1e-8)

  # But reported slope should be slope/delta
  slope_12 <- es[es$from == "S1" & es$to == "S2", "slope"]
  expect_equal(slope_12, slope / delta, tolerance = 1e-8)
})

test_that("Summary includes effect-size columns for regression", {
  mats <- .make_linear_matrices(T_len = 10, slope = 0.01)
  vel <- velocity_tna(mats, method = "regression")
  s <- summary(vel)

  expect_true(all(c("standardized", "pct_change", "total_change") %in% names(s)))
  # S1->S2 should have |standardized| = 1
  s12 <- s[s$from == "S1" & s$to == "S2", ]
  expect_equal(abs(s12$standardized), 1.0, tolerance = 1e-8)
})

test_that("Print shows beta and pct for regression", {
  mats <- .make_linear_matrices(T_len = 10, slope = 0.02)
  vel <- velocity_tna(mats, method = "regression")
  expect_output(print(vel), "beta=")
  expect_output(print(vel), "%/t")
  expect_output(print(vel), "total=")
})

test_that("Effect-size metrics are absent for non-regression methods", {
  mats <- .make_linear_matrices(T_len = 10)
  vel_glla <- velocity_tna(mats, method = "glla", n_embed = 5L)
  expect_null(vel_glla$edge_stats)

  s <- summary(vel_glla)
  expect_false("standardized" %in% names(s))
  expect_false("pct_change" %in% names(s))
  expect_false("total_change" %in% names(s))
})

test_that("All three methods agree on linear data velocity", {
  slope <- 0.01
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)

  vel_reg <- velocity_tna(mats, method = "regression")
  vel_glla <- velocity_tna(mats, method = "glla", n_embed = 5L)
  vel_fd <- velocity_tna(mats, method = "finite_difference",
                         diff_type = "central")

  expect_equal(vel_reg$velocity_matrix[1, 2], slope, tolerance = 1e-6)
  expect_equal(vel_glla$velocity_matrix[1, 2], slope, tolerance = 1e-6)
  expect_equal(vel_fd$velocity_matrix[1, 2], slope, tolerance = 1e-6)
})


# ---- GLLA weight matrix tests ----

test_that(".glla_weight_matrix has correct dimensions", {
  W <- .glla_weight_matrix(5L, 1, 1L)
  expect_equal(dim(W), c(2, 5))

  W2 <- .glla_weight_matrix(7L, 1, 2L)
  expect_equal(dim(W2), c(3, 7))
})

test_that("GLLA on constant series gives velocity ~ 0", {
  mats <- .make_constant_matrices(n_states = 3, T_len = 10)
  vel <- velocity_tna(mats, method = "glla", n_embed = 5L)
  expect_true(all(abs(vel$velocity_matrix) < 1e-10))
})

test_that("GLLA on linear series recovers correct slope", {
  slope <- 0.02
  mats <- .make_linear_matrices(n_states = 3, T_len = 20, slope = slope)
  vel <- velocity_tna(mats, method = "glla", n_embed = 5L)

  expect_equal(vel$velocity_matrix[1, 2], slope, tolerance = 1e-6)
  expect_equal(vel$velocity_matrix[1, 1], -slope, tolerance = 1e-6)
  expect_true(all(abs(vel$velocity_matrix[2:3, ]) < 1e-10))
})

test_that("velocity_tna GLLA returns correct structure", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "glla", n_embed = 5L)

  expect_s3_class(vel, "tna_velocity")
  expect_equal(vel$method, "glla")
  expect_equal(vel$n_embed, 5L)
  expect_equal(vel$n_timepoints, 10)
  expect_equal(vel$n_embedded, 6)
  expect_equal(vel$n_nodes, 3)
  expect_true(vel$directed)
  expect_equal(dim(vel$velocity_matrix), c(3, 3))
  expect_equal(length(vel$velocity_series), 6)
  expect_null(vel$acceleration_matrix)
  expect_equal(length(vel$smoothed_matrices), 6)
  expect_equal(length(vel$original_matrices), 10)
})

test_that("GLLA weight matrix rows sum correctly for constant input", {
  W <- .glla_weight_matrix(5L, 1, 1L)
  const <- rep(1, 5)
  result <- W %*% const
  expect_equal(result[1, 1], 1, tolerance = 1e-10)
  expect_equal(result[2, 1], 0, tolerance = 1e-10)
})

test_that("GLLA weight matrix recovers linear function derivative", {
  W <- .glla_weight_matrix(5L, 1, 1L)
  t_j <- c(-2, -1, 0, 1, 2)
  y <- 2 * t_j
  result <- W %*% y
  expect_equal(result[1, 1], 0, tolerance = 1e-10)
  expect_equal(result[2, 1], 2, tolerance = 1e-10)
})


# ---- Finite difference tests ----

test_that("Forward finite difference on linear data matches slope", {
  slope <- 0.015
  mats <- .make_linear_matrices(n_states = 3, T_len = 10, slope = slope)
  vel <- velocity_tna(mats, method = "finite_difference",
                      diff_type = "forward")

  expect_equal(vel$velocity_matrix[1, 2], slope, tolerance = 1e-10)
  expect_equal(vel$velocity_matrix[1, 1], -slope, tolerance = 1e-10)
  expect_equal(vel$n_embedded, 9)
  expect_null(vel$smoothed_matrices)
})

test_that("Central finite difference on linear data matches slope", {
  slope <- 0.015
  mats <- .make_linear_matrices(n_states = 3, T_len = 10, slope = slope)
  vel <- velocity_tna(mats, method = "finite_difference",
                      diff_type = "central")

  expect_equal(vel$velocity_matrix[1, 2], slope, tolerance = 1e-10)
  expect_equal(vel$n_embedded, 8)
})

test_that("Finite difference constant series gives velocity ~ 0", {
  mats <- .make_constant_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "finite_difference",
                      diff_type = "forward")
  expect_true(all(abs(vel$velocity_matrix) < 1e-10))
})


# ---- Acceleration tests ----

test_that("GLLA on quadratic data recovers constant acceleration", {
  accel <- 0.002
  mats <- .make_quadratic_matrices(n_states = 3, T_len = 15, accel = accel)
  vel <- velocity_tna(mats, method = "glla", n_embed = 5L, order = 2L)

  expect_equal(vel$acceleration_matrix[1, 2], accel, tolerance = 1e-4)
  expect_true(all(abs(vel$acceleration_matrix[2:3, ]) < 1e-10))
})

test_that("Acceleration series all agree for constant acceleration", {
  accel <- 0.002
  mats <- .make_quadratic_matrices(n_states = 3, T_len = 15, accel = accel)
  vel <- velocity_tna(mats, method = "glla", n_embed = 5L, order = 2L)

  acc_vals <- vapply(vel$acceleration_series, function(m) m[1, 2], numeric(1))
  expect_true(all(abs(acc_vals - accel) < 1e-4))
})


# ---- Input format tests ----

test_that("velocity_tna works with list of matrices (default regression)", {
  mats <- .make_linear_matrices(T_len = 8)
  vel <- velocity_tna(mats)
  expect_s3_class(vel, "tna_velocity")
  expect_equal(vel$method, "regression")
  expect_equal(vel$n_timepoints, 8)
})

test_that("velocity_tna works with tna_windows mock", {
  mats <- .make_linear_matrices(T_len = 6)
  mock_windows <- list(
    windows = lapply(mats, function(m) list(weights = m))
  )
  vel <- velocity_tna(mock_windows)
  expect_s3_class(vel, "tna_velocity")
  expect_equal(vel$n_timepoints, 6)
})

test_that("velocity_tna works with raw data.frame", {
  set.seed(123)
  n_per_window <- 50
  actions <- c("A", "B", "C")
  df_list <- lapply(1:6, function(t) {
    data.frame(
      Time = t,
      ID = rep(1:10, each = n_per_window / 10),
      Action = sample(actions, n_per_window, replace = TRUE),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df_list)

  vel <- velocity_tna(df, time_col = "Time", id_col = "ID",
                      action_col = "Action")
  expect_s3_class(vel, "tna_velocity")
  expect_equal(vel$n_timepoints, 6)
  expect_equal(vel$n_nodes, 3)
})

test_that("velocity_tna assigns default dimnames when missing", {
  mats <- lapply(1:5, function(t) {
    m <- matrix(runif(4), 2, 2)
    m / rowSums(m)
  })
  vel <- velocity_tna(mats)
  expect_equal(vel$nodes, c("S1", "S2"))
})

test_that("Missing states in raw data windows are expanded", {
  df <- data.frame(
    Time = c(rep(1, 20), rep(2, 30)),
    ID = c(rep(1:4, each = 5), rep(1:6, each = 5)),
    Action = c(
      sample(c("A", "B"), 20, replace = TRUE),
      sample(c("A", "B", "C"), 30, replace = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  set.seed(42)
  vel <- velocity_tna(df, method = "finite_difference",
                      diff_type = "forward",
                      time_col = "Time", id_col = "ID",
                      action_col = "Action")
  expect_equal(vel$n_nodes, 3)
  expect_true(all(c("A", "B", "C") %in% vel$nodes))
})


# ---- Logistic regression tests ----

test_that("Logistic regression works on wide-format data", {
  set.seed(42)
  n_seq <- 200
  n_cols <- 10
  states <- c("A", "B", "C")
  wide <- as.data.frame(matrix(
    sample(states, n_seq * n_cols, replace = TRUE),
    nrow = n_seq, ncol = n_cols,
    dimnames = list(NULL, paste0("T", seq_len(n_cols)))
  ), stringsAsFactors = FALSE)

  vel <- velocity_tna(wide, method = "regression", regression_type = "logistic")
  expect_s3_class(vel, "tna_velocity")
  expect_equal(vel$regression_type, "logistic")
  expect_equal(vel$n_nodes, 3)
  expect_equal(sort(vel$nodes), c("A", "B", "C"))
})

test_that("Logistic regression returns edge_stats with correct columns", {
  set.seed(42)
  n_seq <- 200
  states <- c("A", "B")
  wide <- as.data.frame(matrix(
    sample(states, n_seq * 8, replace = TRUE),
    nrow = n_seq, ncol = 8,
    dimnames = list(NULL, paste0("T", 1:8))
  ), stringsAsFactors = FALSE)

  vel <- velocity_tna(wide, method = "regression", regression_type = "logistic")
  es <- vel$edge_stats
  expect_true(all(c("from", "to", "slope", "se", "z_value",
                     "p_value", "odds_ratio", "n_transitions",
                     "standardized", "pct_change", "total_change") %in% names(es)))
})

test_that("Logistic regression detects known trend in synthetic data", {
  # Create data where A->B becomes more likely over time
  set.seed(123)
  n_seq <- 500
  n_cols <- 15
  wide <- matrix(NA_character_, nrow = n_seq, ncol = n_cols)
  for (i in seq_len(n_seq)) {
    for (t in seq_len(n_cols)) {
      if (t == 1) {
        wide[i, t] <- sample(c("A", "B"), 1)
      } else {
        prev <- wide[i, t - 1]
        if (prev == "A") {
          p_b <- 0.3 + 0.02 * t  # increasing P(A->B)
          wide[i, t] <- sample(c("A", "B"), 1, prob = c(1 - p_b, p_b))
        } else {
          wide[i, t] <- sample(c("A", "B"), 1, prob = c(0.4, 0.6))
        }
      }
    }
  }
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(n_cols))

  vel <- velocity_tna(wide, method = "regression", regression_type = "logistic")
  es <- vel$edge_stats
  # A->B should have positive slope and be significant
  ab <- es[es$from == "A" & es$to == "B", ]
  expect_true(ab$slope > 0)
  expect_true(ab$p_value < 0.05)
})

test_that("Logistic regression errors on non-dataframe input", {
  mats <- .make_linear_matrices(T_len = 5)
  expect_error(
    velocity_tna(mats, method = "regression", regression_type = "logistic"),
    "wide-format data frame"
  )
})

test_that("Logistic print shows Logistic label", {
  set.seed(42)
  wide <- as.data.frame(matrix(
    sample(c("A", "B"), 100 * 6, replace = TRUE),
    nrow = 100, ncol = 6,
    dimnames = list(NULL, paste0("T", 1:6))
  ), stringsAsFactors = FALSE)
  vel <- velocity_tna(wide, method = "regression", regression_type = "logistic")
  expect_output(print(vel), "Logistic")
})


# ---- S3 method tests ----

test_that("print.tna_velocity works for regression", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "regression")
  expect_output(print(vel), "Transition Network Velocity")
  expect_output(print(vel), "Regression")
})

test_that("print.tna_velocity works for GLLA", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "glla", n_embed = 5L)
  expect_output(print(vel), "GLLA")
})

test_that("print.tna_velocity works for finite_difference", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "finite_difference")
  expect_output(print(vel), "Finite Difference")
})

test_that("summary.tna_velocity returns correct data.frame", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats, method = "regression")
  s <- summary(vel)

  expect_true(is.data.frame(s))
  expect_true(all(c("from", "to", "velocity", "abs_velocity",
                     "direction", "consistency") %in% names(s)))
  expect_equal(nrow(s), 6)
  expect_true(all(diff(s$abs_velocity) <= 0))
})

test_that("summary direction labels correct for regression", {
  mats <- .make_linear_matrices(T_len = 10, slope = 0.01)
  vel <- velocity_tna(mats, method = "regression")
  s <- summary(vel)

  s12 <- s[s$from == "S1" & s$to == "S2", ]
  expect_equal(s12$direction, "strengthening")

  s_other <- s[!(s$from == "S1" & s$to == "S2"), ]
  expect_true(all(s_other$direction == "stable"))
})

test_that("plot.tna_velocity network type works", {
  skip_if_not_installed("cograph")
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats)
  expect_no_error(plot(vel, type = "network"))
})

test_that("plot.tna_velocity series type works", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats)
  expect_no_error(plot(vel, type = "series"))
})

test_that("plot.tna_velocity heatmap type works", {
  mats <- .make_linear_matrices(T_len = 10)
  vel <- velocity_tna(mats)
  expect_no_error(plot(vel, type = "heatmap"))
})


# ---- Edge cases ----

test_that("velocity_tna works with 2x2 matrices", {
  mats <- lapply(1:6, function(t) {
    p <- 0.3 + 0.02 * t
    matrix(c(1 - p, p, 0.4, 0.6), 2, 2, byrow = TRUE,
           dimnames = list(c("A", "B"), c("A", "B")))
  })
  vel <- velocity_tna(mats)
  expect_equal(vel$n_nodes, 2)
  expect_equal(vel$nodes, c("A", "B"))
})

test_that("velocity_tna handles matrices with zero rows", {
  mats <- lapply(1:6, function(t) {
    m <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"),
                                          c("A", "B", "C")))
    m[1, 2] <- 0.5 + 0.01 * t
    m[1, 1] <- 0.5 - 0.01 * t
    m[2, 3] <- 1
    m
  })
  vel <- velocity_tna(mats)
  expect_s3_class(vel, "tna_velocity")
  expect_true(all(vel$velocity_matrix[3, ] == 0))
})
