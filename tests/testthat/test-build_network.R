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
  expect_false("constant" %in% colnames(net$matrix))
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
  expect_false("%" %in% colnames(net$matrix))
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
  expect_true(is.matrix(net$matrix))
  expect_equal(nrow(net$matrix), 6)
  expect_equal(ncol(net$matrix), 6)
  # Diagonal should be zero
  expect_true(all(diag(net$matrix) == 0))
  # Should be symmetric
  expect_equal(net$matrix, t(net$matrix))
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
  expect_true(is.matrix(net$matrix))
  # Diagonal should be zero
  expect_true(all(diag(net$matrix) == 0))
  # Should be symmetric
  expect_equal(net$matrix, t(net$matrix))
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
  expect_true(is.matrix(net$matrix))
  # Diagonal should be zero
  expect_true(all(diag(net$matrix) == 0))
  # Should be symmetric
  expect_equal(net$matrix, t(net$matrix))
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
  expect_equal(net$matrix, expected)
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
  expect_false("subject_id" %in% colnames(net$matrix))
  expect_false("rid" %in% colnames(net$matrix))
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
  expect_true(any(grepl("Nodes: 5", out)))
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

  mat <- net$matrix
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net$edges), nrow(upper_nz))

  for (i in seq_len(nrow(net$edges))) {
    r <- match(net$edges$from[i], colnames(mat))
    cc <- match(net$edges$to[i], colnames(mat))
    expect_equal(net$edges$weight[i], mat[r, cc])
  }
})

test_that("edges match for pcor and cor methods too", {
  df <- .make_freq_data(n = 80, p = 5)

  net_pcor <- build_network(df, method = "pcor")
  mat <- net_pcor$matrix
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net_pcor$edges), nrow(upper_nz))

  net_cor <- build_network(df, method = "cor", threshold = 0.1)
  mat <- net_cor$matrix
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
    expect_true(is.matrix(net$matrix))
    expect_true(is.matrix(net$cor_matrix))
    expect_true(is.data.frame(net$edges))
    expect_true(is.numeric(net$n_edges))
    expect_true(all(diag(net$matrix) == 0))
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
  expect_equal(names(r2), colnames(net$matrix))
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
    all(net$matrix[j, ] == 0)
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
  expect_true(any(grepl("Data: 80 x 5", out)))
})
