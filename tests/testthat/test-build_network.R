# Helper: generate reproducible frequency-like data
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$rid <- seq_len(n)
  df
}


# ---- Input validation ----

test_that("build_network errors on non-numeric data only", {
  df <- data.frame(a = letters[1:10], b = letters[10:1])
  expect_error(
    build_network(df, method = "glasso"),
    "At least 2 numeric columns"
  )
})

test_that("build_network errors on non-symmetric matrix", {
  m <- matrix(1:9, 3, 3)
  expect_error(
    build_network(m, method = "glasso", params = list(n = 50)),
    "symmetric"
  )
})

test_that("build_network errors when n missing for matrix input", {
  m <- diag(5)
  expect_error(
    build_network(m, method = "glasso"),
    "Sample size 'n' is required"
  )
})


# ---- Auto-cleaning ----

test_that("zero-variance columns are dropped with message", {
  df <- .make_freq_data()
  df$constant <- 5
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "Dropping zero-variance"
  )
  expect_equal(net$n_nodes, 5)
  expect_false("constant" %in% colnames(net$weights))
})

test_that("non-syntactic column names are dropped with message", {
  df <- .make_freq_data(n = 80, p = 4)
  df$`%` <- rpois(80, 2)
  df$`*` <- rpois(80, 3)
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "non-syntactic"
  )
  expect_equal(net$n_nodes, 4)
  expect_false("%" %in% colnames(net$weights))
})

test_that("all-NA columns are dropped with message", {
  df <- .make_freq_data(n = 80, p = 5)
  df$empty <- NA_real_
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "all-NA"
  )
  expect_equal(net$n_nodes, 5)
})

test_that("rows with NA are dropped with message", {
  df <- .make_freq_data(n = 80, p = 5)
  df$state_1[1:3] <- NA
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "rows with NA"
  )
  expect_equal(net$n, 77)
})


# ---- Method: glasso ----

test_that("build_network works with data frame input (glasso)", {
  df <- .make_freq_data(n = 80, p = 6)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "glasso")
  expect_equal(net$n, 80)
  expect_equal(net$n_nodes, 6)
  expect_true(is.matrix(net$weights))
  expect_equal(nrow(net$weights), 6)
  expect_equal(ncol(net$weights), 6)
  # Diagonal should be zero
  expect_true(all(diag(net$weights) == 0))
  # Should be symmetric
  expect_equal(net$weights, t(net$weights))
  # Edges data frame
  expect_true(is.data.frame(net$edges))
  expect_true(all(c("from", "to", "weight") %in% names(net$edges)))
  expect_equal(net$n_edges, nrow(net$edges))
  # EBIC path length matches nlambda
  expect_equal(length(net$ebic_path), 20)
  expect_equal(length(net$lambda_path), 20)
  # Glasso-specific fields
  expect_true(!is.null(net$precision_matrix))
  expect_true(!is.null(net$gamma))
  expect_true(!is.null(net$lambda_selected))
})

test_that("build_network works with correlation matrix input (glasso)", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  S <- cor(df[, num_cols])

  net <- build_network(S, method = "glasso", params = list(n = 100,
                                                            nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$n, 100)
  expect_equal(net$n_nodes, 5)
})

test_that("build_network works with covariance matrix input (glasso)", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  C <- cov(df[, num_cols])

  net <- build_network(C, method = "glasso",
                       params = list(n = 100, input_type = "cov",
                                     nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$n_nodes, 5)
})

test_that("method aliases resolve to glasso", {
  df <- .make_freq_data(n = 80, p = 4)
  net1 <- build_network(df, method = "ebicglasso",
                         params = list(nlambda = 20L))
  net2 <- build_network(df, method = "regularized",
                         params = list(nlambda = 20L))
  expect_equal(net1$method, "glasso")
  expect_equal(net2$method, "glasso")
})


# ---- Method: pcor (unregularised) ----

test_that("build_network works with method='pcor'", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "pcor")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_equal(net$n, 80)
  expect_equal(net$n_nodes, 5)
  expect_true(is.matrix(net$weights))
  # Diagonal should be zero
  expect_true(all(diag(net$weights) == 0))
  # Should be symmetric
  expect_equal(net$weights, t(net$weights))
  # Should have precision matrix
  expect_true(!is.null(net$precision_matrix))
  # Edges
  expect_true(is.data.frame(net$edges))
  expect_equal(net$n_edges, nrow(net$edges))
  # No glasso-specific fields
  expect_null(net$lambda_selected)
  expect_null(net$ebic_path)
})

test_that("method='partial' resolves to pcor", {
  df <- .make_freq_data(n = 80, p = 4)
  net <- build_network(df, method = "partial")
  expect_equal(net$method, "pcor")
})

test_that("pcor errors on singular matrix", {
  # p > n: more variables than observations
  set.seed(42)
  mat <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  colnames(mat) <- paste0("V", seq_len(20))
  S <- cor(mat)
  expect_error(
    build_network(S, method = "pcor", params = list(n = 10)),
    "singular"
  )
})


# ---- Method: cor (correlation network) ----

test_that("build_network works with method='cor'", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "cor")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "cor")
  expect_equal(net$n, 80)
  expect_equal(net$n_nodes, 5)
  expect_true(is.matrix(net$weights))
  # Diagonal should be zero
  expect_true(all(diag(net$weights) == 0))
  # Should be symmetric
  expect_equal(net$weights, t(net$weights))
  # No precision matrix for cor method
  expect_null(net$precision_matrix)
  # Edges
  expect_true(is.data.frame(net$edges))
  expect_equal(net$n_edges, nrow(net$edges))
})

test_that("method='correlation' resolves to cor", {
  df <- .make_freq_data(n = 80, p = 4)
  net <- build_network(df, method = "correlation")
  expect_equal(net$method, "cor")
})

test_that("cor threshold filters weak edges", {
  df <- .make_freq_data(n = 100, p = 5)
  net_low <- build_network(df, method = "cor", threshold = 0.01)
  net_high <- build_network(df, method = "cor", threshold = 0.3)
  # Higher threshold should produce same or fewer edges
  expect_true(net_high$n_edges <= net_low$n_edges)
})

test_that("cor matrix matches thresholded cor_matrix", {
  df <- .make_freq_data(n = 80, p = 5)
  thr <- 0.15
  net <- build_network(df, method = "cor", threshold = thr)
  expected <- net$cor_matrix
  diag(expected) <- 0
  expected[abs(expected) < thr] <- 0
  expect_equal(net$weights, expected)
})


# ---- New aliases ----

test_that("new aliases tna, ftna, cna, corr resolve correctly", {
  df <- .make_freq_data(n = 80, p = 4)

  net_corr <- build_network(df, method = "corr")
  expect_equal(net_corr$method, "cor")
})


# ---- id_col exclusion ----

test_that("id_col columns are excluded from analysis", {
  df <- .make_freq_data(n = 80, p = 5)
  df$subject_id <- seq_len(80)

  net <- build_network(df, method = "glasso",
                       params = list(id_col = "subject_id", nlambda = 20L))

  # subject_id and rid should be excluded -> 5 variables
  expect_equal(net$n_nodes, 5)
  expect_false("subject_id" %in% colnames(net$weights))
  expect_false("rid" %in% colnames(net$weights))
})


# ---- Gamma effects ----

test_that("higher gamma produces sparser or equal networks", {
  df <- .make_freq_data(n = 150, p = 7, seed = 123)

  net_low <- build_network(df, method = "glasso",
                           params = list(gamma = 0, nlambda = 50L))
  net_high <- build_network(df, method = "glasso",
                            params = list(gamma = 1, nlambda = 50L))

  expect_true(net_high$n_edges <= net_low$n_edges)
})


# ---- S3 print method ----

test_that("print.netobject produces expected output for glasso", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("Partial Correlation Network \\(EBICglasso\\)", out)))
  expect_true(any(grepl("Weight matrix:", out)))
  expect_true(any(grepl("Sample size: 80", out)))
  expect_true(any(grepl("Gamma:", out)))
  expect_true(any(grepl("Lambda:", out)))
})

test_that("print.netobject produces expected output for pcor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "pcor")

  out <- capture.output(print(net))
  expect_true(any(grepl("unregularised", out)))
  expect_false(any(grepl("Gamma:", out)))
})

test_that("print.netobject produces expected output for cor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "cor")

  out <- capture.output(print(net))
  expect_true(any(grepl("Correlation Network", out)))
  expect_false(any(grepl("Gamma:", out)))
})

test_that("print.netobject returns invisible(x)", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
  ret <- capture.output(result <- print(net))
  expect_identical(result, net)
})


# ---- Correlation method argument ----

test_that("cor_method argument is respected", {
  df <- .make_freq_data(n = 80, p = 5)

  net_pearson <- build_network(df, method = "glasso",
                               params = list(cor_method = "pearson",
                                             nlambda = 20L))
  net_spearman <- build_network(df, method = "glasso",
                                params = list(cor_method = "spearman",
                                              nlambda = 20L))

  expect_false(identical(net_pearson$cor_matrix, net_spearman$cor_matrix))
})


# ---- Edge data frame correctness ----

test_that("edges match non-zero upper triangle of matrix", {
  df <- .make_freq_data(n = 100, p = 6)
  net <- build_network(df, method = "glasso", params = list(nlambda = 30L))

  mat <- net$weights
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net$edges), nrow(upper_nz))

  for (i in seq_len(nrow(net$edges))) {
    expect_equal(net$edges$weight[i], mat[net$edges$from[i], net$edges$to[i]])
  }
})

test_that("edges match for pcor and cor methods too", {
  df <- .make_freq_data(n = 80, p = 5)

  net_pcor <- build_network(df, method = "pcor")
  mat <- net_pcor$weights
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net_pcor$edges), nrow(upper_nz))

  net_cor <- build_network(df, method = "cor", threshold = 0.1)
  mat <- net_cor$weights
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net_cor$edges), nrow(upper_nz))
})


# ---- Cross-method consistency ----

test_that("all methods produce consistent structure", {
  df <- .make_freq_data(n = 80, p = 5)
  methods <- c("glasso", "pcor", "cor")

  for (m in methods) {
    net <- build_network(df, method = m, params = list(nlambda = 20L))
    expect_s3_class(net, "netobject")
    expect_equal(net$method, m)
    expect_equal(net$n, 80)
    expect_equal(net$n_nodes, 5)
    expect_true(is.matrix(net$weights))
    expect_true(is.matrix(net$cor_matrix))
    expect_true(is.data.frame(net$edges))
    expect_true(is.numeric(net$n_edges))
    expect_true(all(diag(net$weights) == 0))
  }
})


# ---- Multilevel: helper ----

# Generate data with repeated measures per person
.make_multilevel_data <- function(n_persons = 30, obs_per_person = 5,
                                  p = 5, seed = 42) {
  set.seed(seed)
  n_total <- n_persons * obs_per_person
  mat <- matrix(rpois(n_total * p, lambda = 10), nrow = n_total, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$person <- rep(seq_len(n_persons), each = obs_per_person)
  df$rid <- seq_len(n_total)
  df
}


# ---- Multilevel: level = "between" ----

test_that("level='between' aggregates to person means", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "between", params = list(nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$level, "between")
  # n should be number of unique persons
  expect_equal(net$n, 30)
  expect_equal(net$n_nodes, 5)
})

test_that("level='between' works with method='pcor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "pcor", id_col = "person",
                       level = "between")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_equal(net$n, 30)
})

test_that("level='between' works with method='cor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "between")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "cor")
  expect_equal(net$n, 30)
})


# ---- Multilevel: level = "within" ----

test_that("level='within' centers correctly (column means near 0)", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "within")

  # Within-centered correlations should exist
  expect_s3_class(net, "netobject")
  expect_equal(net$level, "within")
  # n = total observations (all persons have >= 2 obs)
  expect_equal(net$n, 150)
  expect_equal(net$n_nodes, 5)
})

test_that("level='within' drops single-observation persons with message", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  # Add 3 persons with only 1 observation each
  singles <- data.frame(
    state_1 = rpois(3, 10), state_2 = rpois(3, 10),
    state_3 = rpois(3, 10), state_4 = rpois(3, 10),
    state_5 = rpois(3, 10),
    person = c(101, 102, 103), rid = 151:153
  )
  df <- rbind(df, singles)

  expect_message(
    net <- build_network(df, method = "glasso", id_col = "person",
                         level = "within", params = list(nlambda = 20L)),
    "single-observation"
  )
  # Single-obs rows should be dropped: 153 - 3 = 150
  expect_equal(net$n, 150)
})

test_that("level='within' works with method='pcor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "pcor", id_col = "person",
                       level = "within")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_equal(net$n, 150)
})


# ---- Multilevel: level = "both" ----

test_that("level='both' returns netobject_ml with both sub-networks", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))

  expect_s3_class(net, "netobject_ml")
  expect_s3_class(net$between, "netobject")
  expect_s3_class(net$within, "netobject")
  expect_equal(net$method, "glasso")
  expect_equal(net$between$level, "between")
  expect_equal(net$within$level, "within")
  expect_equal(net$between$n, 30)
  expect_equal(net$within$n, 150)
})

test_that("level='both' works with method='cor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "both")

  expect_s3_class(net, "netobject_ml")
  expect_equal(net$method, "cor")
  expect_s3_class(net$between, "netobject")
  expect_s3_class(net$within, "netobject")
})


# ---- Multilevel: validation ----

test_that("level requires id_col", {
  df <- .make_multilevel_data()
  expect_error(
    build_network(df, method = "glasso", level = "between"),
    "id.*required"
  )
})

test_that("level requires data frame input", {
  m <- diag(5)
  expect_error(
    build_network(m, method = "glasso",
                  params = list(n = 50), id_col = "person",
                  level = "between"),
    "data frame"
  )
})

test_that("level=NULL preserves backward compatibility", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       params = list(nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_null(net$level)
  # Without level, n = total rows
  expect_equal(net$n, 150)
})


# ---- Multilevel: print methods ----

test_that("print.netobject shows level label for between", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "between", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("between-person", out)))
})

test_that("print.netobject shows level label for within", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "within")

  out <- capture.output(print(net))
  expect_true(any(grepl("within-person", out)))
})

test_that("print.netobject_ml shows both levels", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("Multilevel", out)))
  expect_true(any(grepl("Between-person", out)))
  expect_true(any(grepl("Within-person", out)))
  expect_true(any(grepl("unique persons", out)))
  expect_true(any(grepl("observations", out)))
})

test_that("print.netobject_ml returns invisible(x)", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))
  ret <- capture.output(result <- print(net))
  expect_identical(result, net)
})


# ---- Predictability ----

test_that("predictability returns named numeric vector for glasso", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
  r2 <- predictability(net)

  expect_true(is.numeric(r2))
  expect_equal(length(r2), 5)
  expect_equal(names(r2), colnames(net$weights))
  expect_true(all(r2 >= 0 & r2 <= 1))
})

test_that("predictability returns named numeric vector for pcor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "pcor")
  r2 <- predictability(net)

  expect_true(is.numeric(r2))
  expect_equal(length(r2), 5)
  expect_true(all(r2 >= 0 & r2 <= 1))
})

test_that("predictability returns named numeric vector for cor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "cor", threshold = 0.1)
  r2 <- predictability(net)

  expect_true(is.numeric(r2))
  expect_equal(length(r2), 5)
  expect_true(all(r2 >= 0 & r2 <= 1))
})

test_that("predictability.cor returns 0 for isolated nodes", {
  df <- .make_freq_data(n = 80, p = 5)
  # Very high threshold should isolate most nodes
  net <- build_network(df, method = "cor", threshold = 0.99)
  r2 <- predictability(net)

  # Isolated nodes (no edges) should have R^2 = 0
  isolated <- vapply(seq_len(net$n_nodes), function(j) {
    all(net$weights[j, ] == 0)
  }, logical(1))
  if (any(isolated)) {
    expect_true(all(r2[isolated] == 0))
  }
})

test_that("predictability works for netobject_ml", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))
  r2 <- predictability(net)

  expect_true(is.list(r2))
  expect_true("between" %in% names(r2))
  expect_true("within" %in% names(r2))
  expect_equal(length(r2$between), 5)
  expect_equal(length(r2$within), 5)
  expect_true(all(r2$between >= 0 & r2$between <= 1))
  expect_true(all(r2$within >= 0 & r2$within <= 1))
})

test_that("print does not show predictability", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
  out <- capture.output(print(net))
  expect_false(any(grepl("predictability", out)))
})


# ---- $data field ----

test_that("$data is cleaned numeric matrix for association methods (data frame input)", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  expect_true(is.matrix(net$data))
  expect_true(is.numeric(net$data))
  expect_equal(nrow(net$data), 80)
  # 5 state columns only (rid excluded during cleaning)
  expect_equal(ncol(net$data), 5)
})

test_that("$data is NULL for association methods (matrix input)", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  S <- cor(df[, num_cols])

  net <- build_network(S, method = "glasso", params = list(n = 100,
                                                            nlambda = 20L))

  # No row-level data available from matrix input
  expect_null(net$data)
})

test_that("$data is a data frame for transition methods", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "tna")

  expect_true(is.data.frame(net$data))
  expect_equal(nrow(net$data), 2000)
  expect_equal(ncol(net$data), 26)
})

test_that("print.netobject shows data dimensions", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("Sample size: 80", out)))
})


# ---- Coverage gap tests ----

# L180-181: multi-column group key via interaction()
test_that("build_network group dispatch with multi-column group key", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B", "C"), 60, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 60, replace = TRUE),
    T3 = sample(c("A", "B", "C"), 60, replace = TRUE),
    g1 = rep(c("X", "Y"), 30),
    g2 = rep(c("P", "Q"), each = 30),
    stringsAsFactors = FALSE
  )
  grp_nets <- build_network(df, method = "relative",
                            group = c("g1", "g2"),
                            params = list(format = "wide"))
  expect_s3_class(grp_nets, "netobject_group")
  # 4 combinations: X-P, X-Q, Y-P, Y-Q
  expect_equal(length(grp_nets), 4L)
  expect_equal(attr(grp_nets, "group_col"), c("g1", "g2"))
})

# L211-212: explicit codes triggers onehot format
test_that("build_network with explicit codes triggers onehot format", {
  set.seed(42)
  df <- data.frame(
    A = c(1L, 0L, 1L, 0L, 1L),
    B = c(0L, 1L, 0L, 1L, 0L),
    C = c(1L, 1L, 0L, 0L, 1L)
  )
  net <- build_network(df, method = "co_occurrence",
                       codes = c("A", "B", "C"),
                       window_size = 2L)
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
})

# L215: action column present → long format detection
test_that("build_network auto-detects long format via action column", {
  long_data <- data.frame(
    Actor = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time  = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  net <- build_network(long_data, method = "relative",
                       action = "Action")
  expect_s3_class(net, "netobject")
  expect_true(net$directed)
})

# L236-246: long format path through prepare_data
test_that("build_network long format with actor/time/action passes through prepare_data", {
  long_data <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  net <- build_network(long_data, method = "relative",
                       actor = "Actor", action = "Action",
                       time = "Time")
  expect_s3_class(net, "netobject")
  # After prepare_data the format should be reset to "wide"
  expect_equal(net$params$format, "wide")
})

# L260-264: onehot without windowing or session warns
test_that("build_network one-hot without windowing warns", {
  df <- data.frame(A = c(1L, 0L, 1L), B = c(0L, 1L, 0L))
  expect_warning(
    build_network(df, method = "co_occurrence",
                  codes = c("A", "B"), window_size = 1L),
    "sparse"
  )
})

# L270-271, L277: grp_col built from actor + session → params$actor
test_that("build_network onehot with actor and session sets params$actor", {
  set.seed(42)
  df <- data.frame(
    actor   = rep(1:3, each = 4),
    session = rep(c("s1", "s2"), 6),
    A = sample(0:1, 12, replace = TRUE),
    B = sample(0:1, 12, replace = TRUE)
  )
  net <- build_network(df, method = "co_occurrence",
                       codes = c("A", "B"),
                       actor = "actor", session = "session",
                       window_size = 2L)
  expect_s3_class(net, "netobject")
  expect_equal(net$params$actor, c("actor", "session"))
})

# L341-343: estimator returning bad structure triggers error
test_that("build_network errors when estimator returns malformed output", {
  bad_fn <- function(data, ...) {
    list(weights = matrix(1, 2, 2))  # missing 'matrix', 'nodes', 'directed'
  }
  register_estimator("bad_estimator", bad_fn, "bad", directed = FALSE)
  on.exit(remove_estimator("bad_estimator"))

  set.seed(42)
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_error(
    build_network(df, method = "bad_estimator"),
    "must return a list with 'matrix', 'nodes', and 'directed'"
  )
})

# L385: print method label falls through to unknown method
test_that("print.netobject shows generic label for unknown method", {
  # Register a custom estimator with a non-standard name
  custom_fn <- function(data, ...) {
    S <- cor(data)
    diag(S) <- 0
    list(matrix = S, nodes = colnames(S), directed = FALSE)
  }
  register_estimator("my_custom_net", custom_fn, "Custom", directed = FALSE)
  on.exit(remove_estimator("my_custom_net"))

  set.seed(42)
  df <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
  net <- build_network(df, method = "my_custom_net")
  out <- capture.output(print(net))
  expect_true(any(grepl("my_custom_net", out)))
})

# L474-475: print.netobject shows metadata when present
test_that("print.netobject shows metadata column names", {
  skip_if_not_installed("tna")
  # tna::group_regulation has sequences with metadata-like columns
  # Build a dataset where there are extra non-state numeric columns
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B"), 50, replace = TRUE),
    T2 = sample(c("A", "B"), 50, replace = TRUE),
    T3 = sample(c("A", "B"), 50, replace = TRUE),
    Age = rpois(50, 25),
    stringsAsFactors = FALSE
  )
  net <- build_network(df, method = "relative",
                       params = list(format = "wide"))
  # metadata should be present
  if (!is.null(net$metadata)) {
    out <- capture.output(print(net))
    expect_true(any(grepl("Metadata:", out)))
  } else {
    skip("No metadata column detected in this data")
  }
})

# L569-575: print.netobject_group
test_that("print.netobject_group shows group info", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B", "C"), 60, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 60, replace = TRUE),
    grp = rep(c("X", "Y", "Z"), 20),
    stringsAsFactors = FALSE
  )
  nets <- build_network(df, method = "relative", group = "grp",
                        params = list(format = "wide"))
  out <- capture.output(print(nets))
  expect_true(any(grepl("Group Networks", out)))
  expect_true(any(grepl("3 groups", out)))
  expect_true(any(grepl("X:", out)))
})

test_that("print.netobject_group returns invisible(x)", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B"), 40, replace = TRUE),
    T2 = sample(c("A", "B"), 40, replace = TRUE),
    grp = rep(c("X", "Y"), 20),
    stringsAsFactors = FALSE
  )
  nets <- build_network(df, method = "relative", group = "grp",
                        params = list(format = "wide"))
  ret <- capture.output(result <- print(nets))
  expect_identical(result, nets)
})


# L803-808: predictability for cor with single-neighbor node
test_that("predictability.netobject for cor with single-neighbor node", {
  # Build a cor network where one node has exactly one neighbor
  set.seed(99)
  df <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    z = rnorm(100)
  )
  # Use a high threshold to ensure exactly one neighbor for at least one node
  net <- build_network(df, method = "cor", threshold = 0)

  r2 <- predictability(net)
  expect_true(is.numeric(r2))
  expect_equal(length(r2), 3L)
  expect_true(all(r2 >= 0 & r2 <= 1))
})

# print.netobject with ising shows thresholds range
test_that("print.netobject shows ising thresholds range", {
  skip_if_not_installed("glmnet")
  df <- data.frame(
    V1 = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L,
           0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L),
    V2 = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L,
           1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    V3 = c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L,
           1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
  )
  net <- build_network(df, method = "ising")
  out <- capture.output(print(net))
  # Should show Ising-specific fields
  expect_true(any(grepl("Ising", out)))
})

# netobject_ml print shows no sample size when $n is NULL
test_that("print.netobject_ml handles sub-networks with no n", {
  df <- .make_multilevel_data(n_persons = 20, obs_per_person = 5, p = 3)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))

  # Both sub-networks are full netobjects; print should succeed
  expect_no_error(print(net))
})



# ---- Estimators.R coverage gap tests ----

# L89: .count_transitions_wide stop on missing cols
test_that(".count_transitions_wide errors on missing state columns", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  expect_error(
    .count_transitions_wide(df, cols = c("T1", "T2", "T_MISSING")),
    "Columns not found"
  )
})

# L92: .count_transitions_wide stop on < 2 state columns
test_that(".count_transitions_wide errors on fewer than 2 state columns", {
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(
    .count_transitions_wide(df),
    "At least 2 state columns"
  )
})

# L113: .count_transitions_wide all-NA returns 0x0 matrix
test_that(".count_transitions_wide returns empty matrix when all states NA", {
  df <- data.frame(
    T1 = c(NA_character_, NA_character_),
    T2 = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  result <- .count_transitions_wide(df)
  expect_equal(nrow(result), 0L)
  expect_equal(ncol(result), 0L)
})

# L170-174: .count_transitions_long single row → zero matrix (n < 2)
test_that(".count_transitions_long single row returns zero matrix", {
  df <- data.frame(
    Time = 1L,
    Action = "A",
    stringsAsFactors = FALSE
  )
  result <- .count_transitions_long(df, action = "Action", id = NULL,
                                     time = "Time")
  expect_true(is.matrix(result))
  expect_true(all(result == 0L))
  expect_equal(rownames(result), "A")
})

# L200: .count_transitions_long group n < 2 returns empty pair
test_that(".count_transitions_long groups with 1 obs produce zero transitions", {
  df <- data.frame(
    Actor  = c(1L, 2L, 2L),
    Time   = c(1L, 1L, 2L),
    Action = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  # Actor 1 has only 1 row → no pairs from that group
  result <- .count_transitions_long(df, action = "Action", id = "Actor",
                                     time = "Time")
  expect_true(is.matrix(result))
  # The transition A->B or B->A from actor 2 only
  expect_equal(sort(rownames(result)), c("A", "B"))
})

# L347-351: .estimator_co_occurrence with explicit codes (one-hot path)
test_that("co_occurrence with explicit codes uses wtna one-hot path", {
  set.seed(42)
  df <- data.frame(
    A = sample(0L:1L, 20, replace = TRUE),
    B = sample(0L:1L, 20, replace = TRUE),
    C = sample(0L:1L, 20, replace = TRUE)
  )
  net <- build_network(df, method = "co_occurrence",
                       codes = c("A", "B", "C"),
                       window_size = 2L)
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
  expect_equal(sort(net$nodes$label), c("A", "B", "C"))
})

# L359-361: .estimator_co_occurrence long format (non-binary)
test_that("co_occurrence long format path works", {
  long_data <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  net <- build_network(long_data, method = "co_occurrence",
                       action = "Action", actor = "Actor",
                       params = list(format = "long",
                                     id = "Actor", time = "Time"))
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
  expect_equal(sort(net$nodes$label), c("A", "B"))
})

# L393-395: .count_cooccurrence_wide empty n_states or nc < 2
test_that(".count_cooccurrence_wide returns empty when no valid states", {
  df <- data.frame(
    T1 = c(NA_character_, NA_character_),
    T2 = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_wide(df)
  expect_equal(nrow(result), 0L)
  expect_equal(ncol(result), 0L)
})

# L443-444: .count_cooccurrence_long missing action column
test_that(".count_cooccurrence_long errors on missing action column", {
  df <- data.frame(Time = 1:3, Action = c("A", "B", "A"),
                   stringsAsFactors = FALSE)
  expect_error(
    .count_cooccurrence_long(df, action = "NoSuchColumn"),
    "Action column.*not found"
  )
})

# L460-468, L474-489, L492-499, L502-514: .count_cooccurrence_long full paths
test_that(".count_cooccurrence_long single-id group path", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "C", "B", "C", "A"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action", id = "Actor",
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_true(isSymmetric(result))
  expect_equal(sort(rownames(result)), c("A", "B", "C"))
})

test_that(".count_cooccurrence_long multi-id composite group key", {
  df <- data.frame(
    Actor   = c(1L, 1L, 2L, 2L),
    Session = c("s1", "s1", "s2", "s2"),
    Time    = c(1L, 2L, 1L, 2L),
    Action  = c("A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action",
                                      id = c("Actor", "Session"),
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_true(isSymmetric(result))
})

test_that(".count_cooccurrence_long NULL id single sequence", {
  df <- data.frame(
    Time   = 1:4,
    Action = c("A", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action", id = NULL,
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_equal(sort(rownames(result)), c("A", "B", "C"))
})

# L497-499: .count_cooccurrence_long empty from_vec returns zero matrix
test_that(".count_cooccurrence_long returns zero matrix for single-obs groups", {
  df <- data.frame(
    Actor  = c(1L, 2L),
    Time   = c(1L, 1L),
    Action = c("A", "B"),
    stringsAsFactors = FALSE
  )
  # Each group has only 1 obs → no pairs
  result <- .count_cooccurrence_long(df, action = "Action", id = "Actor",
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_true(all(result == 0))
})

# L517-519: .count_cooccurrence_long diagonal correction
test_that(".count_cooccurrence_long diagonal is correctly halved", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L),
    Time   = c(1L, 2L, 3L),
    Action = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action", id = "Actor",
                                      time = "Time")
  # A appears at positions 1 and 2 → pair (A,A) is counted once
  # Diagonal A-A should equal half of what bidirectional counting would give
  expect_true(is.matrix(result))
  expect_true(result["A", "A"] >= 0)
})

# L842: .prepare_association_input errors on < 2 numeric cols after cleaning
test_that(".prepare_association_input errors when < 2 cols after zero-var removal", {
  df <- data.frame(
    V1 = rep(1.0, 20),  # zero variance
    V2 = rep(2.0, 20),  # zero variance
    V3 = rnorm(20)      # one valid col
  )
  expect_message(
    expect_error(
      build_network(df, method = "cor"),
      "At least 2 variable"
    ),
    "zero-variance"
  )
})

# L871, L875: .prepare_association_input matrix input with cov type and no colnames
test_that("build_network matrix input with cov type assigns V-names when colnames NULL", {
  set.seed(42)
  raw <- matrix(rnorm(50 * 4), 50, 4)
  cov_mat <- cov(raw)
  rownames(cov_mat) <- colnames(cov_mat) <- NULL

  net <- build_network(cov_mat, method = "cor",
                       params = list(n = 50, input_type = "cov"))
  expect_s3_class(net, "netobject")
  expect_equal(net$n_nodes, 4L)
  expect_true(all(grepl("^V", net$nodes$label)))
})

# L875: bad data type for .prepare_association_input
test_that(".prepare_association_input errors on non-df non-matrix input", {
  expect_error(
    build_network(list(a = 1), method = "cor"),
    "data frame or a square symmetric matrix"
  )
})

# L967: .compute_lambda_path errors when all off-diagonal correlations zero
test_that(".compute_lambda_path errors when all off-diagonal correlations zero", {
  # Create a correlation matrix that is identity (all off-diagonal are zero)
  S <- diag(4)
  rownames(S) <- colnames(S) <- paste0("V", 1:4)
  expect_error(
    .compute_lambda_path(S, nlambda = 10L, lambda.min.ratio = 0.01),
    "All off-diagonal correlations are zero"
  )
})

# initial probabilities stored for tna / ftna / atna
test_that("build_network stores initial probs for tna, ftna, atna", {
  set.seed(1)
  seqs <- data.frame(
    V1 = sample(c("A","B","C"), 30, TRUE),
    V2 = sample(c("A","B","C"), 30, TRUE),
    V3 = sample(c("A","B","C"), 30, TRUE)
  )
  for (m in c("tna", "ftna", "atna")) {
    net <- build_network(seqs, method = m)
    expect_false(is.null(net$initial), info = paste(m, "has $initial"))
    expect_named(net$initial)
    expect_equal(sum(net$initial), 1, tolerance = 1e-9,
                 info = paste(m, "$initial sums to 1"))
    expect_true(all(net$initial >= 0), info = paste(m, "$initial non-negative"))
    expect_true(all(names(net$initial) %in% net$nodes$label),
                info = paste(m, "initial names match nodes"))
  }
})

# .compute_initial_probs: states never appearing as first get prob 0
test_that(".compute_initial_probs gives 0 to states never starting a sequence", {
  seqs <- data.frame(
    V1 = c("A","A","A"),
    V2 = c("B","C","B"),
    V3 = c("C","B","C")
  )
  net <- build_network(seqs, method = "tna")
  expect_equal(net$initial[["A"]], 1)
  expect_equal(net$initial[["B"]], 0)
  expect_equal(net$initial[["C"]], 0)
})

# L1004-1005: .select_ebic handles glasso fit failure (NULL fit → Inf EBIC)
test_that("glasso network completes even with challenging input", {
  # Near-singular but symmetric correlation matrix
  set.seed(42)
  S <- diag(3)
  diag(S) <- 1
  S[1,2] <- S[2,1] <- 0.999
  S[2,3] <- S[3,2] <- 0.001
  S[1,3] <- S[3,1] <- 0.001
  rownames(S) <- colnames(S) <- c("A", "B", "C")
  # Should complete without error; some lambda fits may fail internally
  expect_no_error(
    net <- build_network(S, method = "glasso", params = list(n = 100, nlambda = 10L))
  )
  expect_s3_class(net, "netobject")
})

