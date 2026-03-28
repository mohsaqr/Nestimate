# ---- Tests for internal utility functions ----

check_val_in_range <- Nestimate:::check_val_in_range
safe_median <- Nestimate:::safe_median
safe_mean <- Nestimate:::safe_mean
safe_sd <- Nestimate:::safe_sd

test_that("check_val_in_range returns TRUE for values inside range", {
  expect_true(check_val_in_range(5, c(1, 10)))
  expect_true(check_val_in_range(0.5, c(0, 1)))
  expect_true(check_val_in_range(-3, c(-5, 0)))
})

test_that("check_val_in_range returns TRUE at boundaries", {
  expect_true(check_val_in_range(1, c(1, 10)))
  expect_true(check_val_in_range(10, c(1, 10)))
  expect_true(check_val_in_range(0, c(0, 0)))
})

test_that("check_val_in_range returns FALSE for values outside range", {
  expect_false(check_val_in_range(0, c(1, 10)))
  expect_false(check_val_in_range(11, c(1, 10)))
  expect_false(check_val_in_range(-1, c(0, 1)))
})

test_that("check_val_in_range returns TRUE when range_val is NULL", {
  expect_true(check_val_in_range(42, NULL))
  expect_true(check_val_in_range(-999, NULL))
})

test_that("check_val_in_range returns FALSE for NA or non-numeric", {
  expect_false(check_val_in_range(NA, c(1, 10)))
  expect_false(check_val_in_range("a", c(1, 10)))
  expect_false(check_val_in_range(NA_real_, c(0, 1)))
})

test_that("safe_median returns correct median", {
  expect_equal(safe_median(c(1, 2, 3)), 2)
  expect_equal(safe_median(c(1, 2, 3, 4)), 2.5)
  expect_equal(safe_median(c(1, NA, 3)), 2)
})

test_that("safe_median returns NA_real_ for empty vector", {
  expect_identical(safe_median(numeric(0)), NA_real_)
})

test_that("safe_mean returns correct mean", {
  expect_equal(safe_mean(c(2, 4, 6)), 4)
  expect_equal(safe_mean(c(1, NA, 3)), 2)
})

test_that("safe_mean returns NA_real_ for empty vector", {
  expect_identical(safe_mean(numeric(0)), NA_real_)
})

test_that("safe_sd returns correct sd", {
  expect_equal(safe_sd(c(1, 2, 3)), sd(c(1, 2, 3)))
  expect_equal(safe_sd(c(10, NA, 30)), sd(c(10, 30), na.rm = TRUE))
})

test_that("safe_sd returns NA_real_ for single or empty vector", {
  expect_identical(safe_sd(numeric(0)), NA_real_)
  expect_identical(safe_sd(5), NA_real_)
  expect_identical(safe_sd(c(NA, NA)), NA_real_)
})


# ---- .as_netobject() coverage ----

test_that(".as_netobject returns netobject unchanged (L120 short-circuit)", {
  net <- build_network(tna::group_regulation, method = "relative")
  result <- Nestimate:::.as_netobject(net)
  expect_identical(result, net)
})

test_that(".as_netobject errors on non-cograph_network input (L170)", {
  expect_error(
    Nestimate:::.as_netobject(list(x = 1)),
    "Expected a netobject or cograph_network"
  )
})

test_that(".as_netobject infers co_occurrence method from symmetric matrix (L119-120)", {
  # Build a symmetric cograph_network manually (no tna metadata)
  m <- matrix(c(0, 0.5, 0.5, 0), nrow = 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  nodes_df <- data.frame(id = 1:2, label = c("A", "B"), name = c("A", "B"),
                         x = c(0, 1), y = c(0, 0), stringsAsFactors = FALSE)
  cg <- structure(list(
    weights = m, nodes = nodes_df, directed = FALSE,
    meta = list(source = "test", layout = NULL, tna = list(method = NULL)),
    node_groups = NULL, data = NULL
  ), class = c("cograph_network", "list"))

  result <- Nestimate:::.as_netobject(cg)
  expect_s3_class(result, "netobject")
  expect_equal(result$method, "co_occurrence")
})

test_that(".as_netobject infers relative method from asymmetric matrix (L121-122)", {
  m <- matrix(c(0, 0.3, 0.7, 0.4, 0, 0.6, 0.5, 0.5, 0), nrow = 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  nodes_df <- data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"), x = 1:3, y = 1:3,
                         stringsAsFactors = FALSE)
  cg <- structure(list(
    weights = m, nodes = nodes_df, directed = TRUE,
    meta = list(source = "test", layout = NULL, tna = list(method = NULL)),
    node_groups = NULL, data = NULL
  ), class = c("cograph_network", "list"))

  result <- Nestimate:::.as_netobject(cg)
  expect_s3_class(result, "netobject")
  expect_equal(result$method, "relative")
})

test_that(".as_netobject decodes integer-encoded data for sequence methods (L133-140)", {
  # Simulate cograph_network with integer-encoded raw data (sequence method)
  m <- matrix(c(0, 1, 1, 0), nrow = 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  nodes_df <- data.frame(id = 1:2, label = c("A", "B"),
                         name = c("A", "B"), x = 1:2, y = 1:2,
                         stringsAsFactors = FALSE)
  int_data <- data.frame(V1 = c(1L, 2L), V2 = c(2L, 1L))
  cg <- structure(list(
    weights = m, nodes = nodes_df, directed = TRUE,
    meta = list(source = "test", layout = NULL,
                tna = list(method = "relative")),
    node_groups = NULL, data = int_data
  ), class = c("cograph_network", "list"))

  result <- Nestimate:::.as_netobject(cg)
  expect_s3_class(result, "netobject")
  expect_true(is.character(result$data[[1]]))
  expect_true(all(result$data[[1]] %in% c("A", "B")))
})


test_that(".as_netobject converts data.frame weights to matrix (L160)", {
  nodes_df <- data.frame(id = 1:2, label = c("A","B"), name = c("A","B"),
                         x = 1:2, y = 1:2, stringsAsFactors = FALSE)
  w_df <- data.frame(A = c(0, 0.3), B = c(0.7, 0))
  rownames(w_df) <- c("A", "B")
  cg <- structure(list(
    weights = w_df, nodes = nodes_df, directed = TRUE,
    meta = list(source = "test", layout = NULL, tna = list(method = "relative")),
    node_groups = NULL, data = NULL
  ), class = c("cograph_network", "list"))
  result <- Nestimate:::.as_netobject(cg)
  expect_true(is.matrix(result$weights))
  expect_true(is.numeric(result$weights))
})

test_that(".as_netobject converts non-numeric matrix to double (L161)", {
  nodes_df <- data.frame(id = 1:2, label = c("A","B"), name = c("A","B"),
                         x = 1:2, y = 1:2, stringsAsFactors = FALSE)
  w_mat <- matrix(c("0", "0.3", "0.7", "0"), nrow = 2,
                  dimnames = list(c("A","B"), c("A","B")))
  cg <- structure(list(
    weights = w_mat, nodes = nodes_df, directed = TRUE,
    meta = list(source = "test", layout = NULL, tna = list(method = "relative")),
    node_groups = NULL, data = NULL
  ), class = c("cograph_network", "list"))
  result <- Nestimate:::.as_netobject(cg)
  expect_true(is.matrix(result$weights))
  expect_true(is.numeric(result$weights))
})

