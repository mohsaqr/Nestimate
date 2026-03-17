# ---- estimate_network() Tests ----
# estimate_network() is deprecated; these tests verify backward compat.

# Helper: generate reproducible frequency-like data
.make_assoc_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}

# Helper: generate wide sequence data
.make_wide_seq <- function(n = 50, t = 10, states = c("A", "B", "C"),
                           seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}


# ---- Deprecation ----

test_that("estimate_network emits deprecation warning", {
  wide <- .make_wide_seq()
  expect_warning(
    estimate_network(wide, method = "relative"),
    "deprecated"
  )
})


# ---- Transition methods ----

test_that("estimate_network works with method='relative'", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "relative"))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "relative")
  expect_true(net$directed)
  expect_true(is.matrix(net$weights))
  expect_equal(nrow(net$weights), 3)
  expect_equal(ncol(net$weights), 3)

  # Rows should sum to approximately 1
  row_sums <- rowSums(net$weights)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  expect_equal(net$n_nodes, 3)
  expect_true(net$n_edges > 0)
  expect_true(is.data.frame(net$edges))
})

test_that("estimate_network works with method='frequency'", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "frequency"))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "frequency")
  expect_true(net$directed)
  expect_true(is.integer(net$weights))
  expect_true(all(net$weights >= 0))

  # frequency_matrix should be available
  expect_true(!is.null(net$frequency_matrix))
})

test_that("estimate_network works with method='co_occurrence'", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "co_occurrence"))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "co_occurrence")
  expect_false(net$directed)
  # Symmetric
  expect_equal(net$weights, t(net$weights))
  # Diagonal contains self-co-occurrence counts (non-negative)
  expect_true(all(diag(net$weights) >= 0))
})


# ---- Association methods ----

test_that("estimate_network works with method='glasso'", {
  df <- .make_assoc_data(n = 80, p = 6)
  net <- suppressWarnings(
    estimate_network(df, method = "glasso", params = list(nlambda = 20L))
  )

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "glasso")
  expect_false(net$directed)
  expect_true(is.matrix(net$weights))
  expect_equal(nrow(net$weights), 6)
  expect_true(all(diag(net$weights) == 0))
  # Symmetric
  expect_equal(net$weights, t(net$weights))
  # Method-specific extras
  expect_true(!is.null(net$precision_matrix))
  expect_true(!is.null(net$gamma))
  expect_true(!is.null(net$lambda_selected))
  expect_true(!is.null(net$cor_matrix))
})

test_that("estimate_network works with method='pcor'", {
  df <- .make_assoc_data(n = 80, p = 5)
  net <- suppressWarnings(estimate_network(df, method = "pcor"))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_false(net$directed)
  expect_true(all(diag(net$weights) == 0))
  expect_equal(net$weights, t(net$weights))
  expect_true(!is.null(net$precision_matrix))
})

test_that("estimate_network works with method='cor'", {
  df <- .make_assoc_data(n = 80, p = 5)
  net <- suppressWarnings(estimate_network(df, method = "cor"))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "cor")
  expect_false(net$directed)
  expect_true(all(diag(net$weights) == 0))
  expect_equal(net$weights, t(net$weights))
})


# ---- Method aliases ----

test_that("method aliases resolve correctly", {
  df <- .make_assoc_data(n = 80, p = 4)
  wide <- .make_wide_seq()

  net1 <- suppressWarnings(
    estimate_network(df, method = "ebicglasso", params = list(nlambda = 20L))
  )
  expect_equal(net1$method, "glasso")

  net2 <- suppressWarnings(
    estimate_network(df, method = "regularized", params = list(nlambda = 20L))
  )
  expect_equal(net2$method, "glasso")

  net3 <- suppressWarnings(estimate_network(df, method = "partial"))
  expect_equal(net3$method, "pcor")

  net4 <- suppressWarnings(estimate_network(df, method = "correlation"))
  expect_equal(net4$method, "cor")

  net5 <- suppressWarnings(estimate_network(wide, method = "transition"))
  expect_equal(net5$method, "relative")

  net6 <- suppressWarnings(estimate_network(wide, method = "counts"))
  expect_equal(net6$method, "frequency")
})


# ---- Scaling ----

test_that("scaling='minmax' works", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", scaling = "minmax")
  )

  nz <- net$weights[net$weights != 0]
  if (length(nz) > 1) {
    expect_true(min(nz) >= 0)
    expect_true(max(nz) <= 1 + 1e-10)
  }
})

test_that("scaling='max' normalizes by max absolute value", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", scaling = "max")
  )
  expect_true(max(abs(net$weights)) <= 1 + 1e-10)
})

test_that("scaling='rank' replaces values with ranks", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", scaling = "rank")
  )
  nz <- net$weights[net$weights != 0]
  if (length(nz) > 0) {
    # Ranks should be positive integers or half-integers (from ties)
    expect_true(all(nz > 0))
    expect_true(all(nz <= length(nz)))
  }
})

test_that("scaling='normalize' makes rows sum to ~1 (absolute values)", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", scaling = "normalize")
  )
  rs <- rowSums(abs(net$weights))
  nonzero <- rs > 0
  expect_true(all(abs(rs[nonzero] - 1) < 1e-10))
})

test_that("combined scaling works", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative",
                     scaling = c("rank", "minmax"))
  )
  expect_s3_class(net, "netobject")
  nz <- net$weights[net$weights != 0]
  if (length(nz) > 1) {
    expect_true(min(nz) >= 0)
    expect_true(max(nz) <= 1 + 1e-10)
  }
})

test_that("invalid scaling errors", {
  wide <- .make_wide_seq()
  expect_error(
    suppressWarnings(
      estimate_network(wide, method = "relative", scaling = "invalid")
    ),
    "Unknown scaling"
  )
})


# ---- Threshold ----

test_that("threshold filters weak edges", {
  df <- .make_assoc_data(n = 100, p = 5)
  net_low <- suppressWarnings(
    estimate_network(df, method = "cor", threshold = 0.01)
  )
  net_high <- suppressWarnings(
    estimate_network(df, method = "cor", threshold = 0.3)
  )

  expect_true(net_high$n_edges <= net_low$n_edges)
  # All non-zero values should be >= threshold
  nz <- abs(net_high$weights[net_high$weights != 0])
  if (length(nz) > 0) {
    expect_true(all(nz >= 0.3))
  }
})


# ---- Directed vs Undirected edge extraction ----

test_that("directed network has directional edges", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "relative"))

  # For a 3-state directed network, can have up to 3*2=6 directed edges
  # (excluding self-loops on diagonal)
  expect_true(net$directed)
  # Edge from A->B and B->A should be separate
  if (nrow(net$edges) > 0) {
    expect_true(all(c("from", "to", "weight") %in% names(net$edges)))
  }
})

test_that("undirected network uses upper triangle only", {
  df <- .make_assoc_data(n = 80, p = 5)
  net <- suppressWarnings(estimate_network(df, method = "cor"))

  expect_false(net$directed)
  # n_edges should match upper triangle non-zeros
  mat <- net$weights
  upper_nz <- sum(upper.tri(mat) & mat != 0)
  expect_equal(net$n_edges, upper_nz)
})


# ---- Multilevel ----

test_that("multilevel level='between' works with association methods", {
  set.seed(42)
  n_persons <- 30
  obs_per <- 5
  df <- data.frame(
    person = rep(seq_len(n_persons), each = obs_per),
    state_1 = rpois(n_persons * obs_per, 10),
    state_2 = rpois(n_persons * obs_per, 10),
    state_3 = rpois(n_persons * obs_per, 10)
  )

  net <- suppressWarnings(
    estimate_network(df, method = "cor", level = "between",
                     id_col = "person")
  )

  expect_s3_class(net, "netobject")
  expect_equal(net$level, "between")
  expect_equal(net$n, 30)
})

test_that("multilevel level='within' works", {
  set.seed(42)
  n_persons <- 30
  obs_per <- 5
  df <- data.frame(
    person = rep(seq_len(n_persons), each = obs_per),
    state_1 = rpois(n_persons * obs_per, 10),
    state_2 = rpois(n_persons * obs_per, 10),
    state_3 = rpois(n_persons * obs_per, 10)
  )

  net <- suppressWarnings(
    estimate_network(df, method = "cor", level = "within",
                     id_col = "person")
  )

  expect_s3_class(net, "netobject")
  expect_equal(net$level, "within")
  expect_equal(net$n, 150)
})

test_that("multilevel level='both' returns netobject_ml", {
  set.seed(42)
  n_persons <- 30
  obs_per <- 5
  df <- data.frame(
    person = rep(seq_len(n_persons), each = obs_per),
    state_1 = rpois(n_persons * obs_per, 10),
    state_2 = rpois(n_persons * obs_per, 10),
    state_3 = rpois(n_persons * obs_per, 10)
  )

  net <- suppressWarnings(
    estimate_network(df, method = "cor", level = "both",
                     id_col = "person")
  )

  expect_s3_class(net, "netobject_ml")
  expect_s3_class(net$between, "netobject")
  expect_s3_class(net$within, "netobject")
  expect_equal(net$method, "cor")
})

test_that("level requires id_col", {
  df <- .make_assoc_data()
  expect_error(
    suppressWarnings(
      estimate_network(df, method = "cor", level = "between")
    ),
    "id.*required"
  )
})

test_that("level requires data frame", {
  m <- diag(5)
  expect_error(
    suppressWarnings(
      estimate_network(m, method = "cor", level = "between",
                       id_col = "person", params = list(n = 50))
    ),
    "data frame"
  )
})


# ---- Params composability ----

test_that("params list is stored and reusable", {
  wide <- .make_wide_seq()
  params <- list(format = "wide")
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", params = params)
  )

  expect_identical(net$params, params)

  # Can be replayed
  net2 <- build_network(wide, method = net$method, params = net$params)
  expect_equal(net$weights, net2$weights)
})


# ---- Print methods ----

test_that("print.netobject works for transition", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "relative"))
  out <- capture.output(print(net))
  expect_true(any(grepl("Transition Network", out)))
  expect_true(any(grepl("directed", out)))
  expect_true(any(grepl("Nodes", out)))
  expect_true(any(grepl("Edges", out)))
})

test_that("print.netobject works for association", {
  df <- .make_assoc_data(n = 80, p = 5)
  net <- suppressWarnings(
    estimate_network(df, method = "glasso", params = list(nlambda = 20L))
  )
  out <- capture.output(print(net))
  expect_true(any(grepl("EBICglasso", out)))
  expect_true(any(grepl("undirected", out)))
  expect_true(any(grepl("Gamma:", out)))
})

test_that("print.netobject shows scaling", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", scaling = "minmax")
  )
  out <- capture.output(print(net))
  expect_true(any(grepl("Scaling:", out)))
})

test_that("print.netobject shows threshold", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(
    estimate_network(wide, method = "relative", threshold = 0.1)
  )
  out <- capture.output(print(net))
  expect_true(any(grepl("Threshold:", out)))
})

test_that("print.netobject returns invisible(x)", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "relative"))
  ret <- capture.output(result <- print(net))
  expect_identical(result, net)
})

test_that("print.netobject_ml works", {
  set.seed(42)
  df <- data.frame(
    person = rep(1:20, each = 5),
    state_1 = rpois(100, 10),
    state_2 = rpois(100, 10),
    state_3 = rpois(100, 10)
  )
  net <- suppressWarnings(
    estimate_network(df, method = "cor", level = "both",
                     id_col = "person")
  )
  out <- capture.output(print(net))
  expect_true(any(grepl("Multilevel", out)))
  expect_true(any(grepl("Between-person", out)))
  expect_true(any(grepl("Within-person", out)))
})


# ---- Edge data frame correctness ----

test_that("edges match non-zero entries in matrix (directed)", {
  wide <- .make_wide_seq()
  net <- suppressWarnings(estimate_network(wide, method = "relative"))

  mat <- net$weights
  # All non-diagonal non-zero entries
  expected_n <- sum(mat != 0 & row(mat) != col(mat))
  expect_equal(nrow(net$edges), expected_n)
})

test_that("edges match upper triangle (undirected)", {
  df <- .make_assoc_data(n = 80, p = 5)
  net <- suppressWarnings(estimate_network(df, method = "cor"))

  mat <- net$weights
  expected_n <- sum(upper.tri(mat) & mat != 0)
  expect_equal(nrow(net$edges), expected_n)
})


# ---- tna package integration ----

test_that("estimate_network works with tna::group_regulation (wide)", {
  skip_if_not_installed("tna")

  net <- suppressWarnings(
    estimate_network(tna::group_regulation, method = "relative")
  )

  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  expect_equal(net$n_nodes, 9)
  # Rows sum to ~1
  row_sums <- rowSums(net$weights)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("estimate_network frequency matches frequencies() output", {
  skip_if_not_installed("tna")

  net <- suppressWarnings(
    estimate_network(tna::group_regulation, method = "frequency")
  )
  freq_mat <- frequencies(tna::group_regulation, format = "wide")

  expect_equal(net$weights, freq_mat)
})
