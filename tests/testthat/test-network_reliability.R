testthat::skip_on_cran()

# ---- Tests for network_reliability() ----

test_that(".scale_matrix handles edge cases", {
  mat <- matrix(5, 3, 3)  # constant matrix

  # minmax of constant → unchanged
  scaled <- Nestimate:::.scale_matrix(mat, "minmax")
  expect_equal(scaled, mat)

  # standardize of constant → unchanged
  scaled <- Nestimate:::.scale_matrix(mat, "standardize")
  expect_equal(scaled, mat)

  # proportion of zero matrix → unchanged
  zero_mat <- matrix(0, 3, 3)
  scaled <- Nestimate:::.scale_matrix(zero_mat, "proportion")
  expect_equal(scaled, zero_mat)
})


# ---- netobject_group flattening (L206-207) ----

test_that("network_reliability works for cor (association) method", {
  set.seed(42)
  df <- as.data.frame(matrix(rpois(100 * 5, 10), nrow = 100))
  colnames(df) <- paste0("V", 1:5)
  net <- build_network(df, method = "cor")
  rel <- network_reliability(net, iter = 30L, seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_equal(nrow(rel$iterations), 30L)
  expect_true(all(c("mean_dev", "median_dev", "cor", "max_dev") %in%
                    names(rel$iterations)))
})


# ---- association path: pcor method (L244-305) ----

test_that("network_reliability works for pcor (association) method", {
  set.seed(7)
  df <- as.data.frame(matrix(rnorm(80 * 4), nrow = 80))
  colnames(df) <- paste0("V", 1:4)
  net <- build_network(df, method = "pcor")
  rel <- network_reliability(net, iter = 20L, seed = 7)

  expect_s3_class(rel, "net_reliability")
  expect_equal(nrow(rel$iterations), 20L)
  cors <- rel$iterations$cor[!is.na(rel$iterations$cor)]
  expect_true(all(cors >= -1 & cors <= 1))
})


# ---- association path with scale (L293-294) ----

test_that("network_reliability applies scale in association path", {
  set.seed(3)
  df <- as.data.frame(matrix(rpois(80 * 5, 10), nrow = 80))
  colnames(df) <- paste0("V", 1:5)
  net <- build_network(df, method = "cor")
  rel <- network_reliability(net, iter = 20L, scale = "minmax", seed = 3)

  expect_s3_class(rel, "net_reliability")
  expect_equal(rel$scale, "minmax")
})


# ---- association path: summary values (L244-305) ----

test_that("network_reliability summary mean_dev is non-negative for cor method", {
  set.seed(10)
  df <- as.data.frame(matrix(rpois(100 * 4, 10), nrow = 100))
  colnames(df) <- paste0("V", 1:4)
  net <- build_network(df, method = "cor")
  rel <- network_reliability(net, iter = 20L, seed = 10)

  expect_true(all(rel$iterations$mean_dev >= 0, na.rm = TRUE))
  expect_true(all(rel$iterations$max_dev >= 0, na.rm = TRUE))
})
