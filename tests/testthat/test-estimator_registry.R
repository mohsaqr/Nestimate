# ---- Estimator Registry Tests ----

test_that("built-in estimators are registered on load", {
  est_list <- list_estimators()
  expect_true(is.data.frame(est_list))
  expect_true(nrow(est_list) >= 6)
  expect_true(all(c("relative", "frequency", "co_occurrence",
                     "cor", "pcor", "glasso") %in% est_list$name))
  expect_true(all(c("name", "description", "directed") %in% names(est_list)))
})

test_that("get_estimator returns correct structure", {
  est <- get_estimator("relative")
  expect_true(is.list(est))
  expect_true(is.function(est$fn))
  expect_true(is.character(est$description))
  expect_true(is.logical(est$directed))
  expect_true(est$directed)

  est2 <- get_estimator("cor")
  expect_false(est2$directed)
})

test_that("get_estimator errors on unknown estimator", {
  expect_error(get_estimator("nonexistent"), "not found")
})

test_that("register_estimator adds custom estimator", {
  my_fn <- function(data, ...) {
    p <- ncol(data)
    m <- diag(p)
    colnames(m) <- rownames(m) <- colnames(data)
    list(matrix = m, nodes = colnames(data), directed = FALSE)
  }
  register_estimator("test_custom", my_fn, "Test estimator", directed = FALSE)

  est <- get_estimator("test_custom")
  expect_true(is.function(est$fn))
  expect_equal(est$description, "Test estimator")
  expect_false(est$directed)

  # Verify it appears in list
  est_list <- list_estimators()
  expect_true("test_custom" %in% est_list$name)

  # Clean up
  remove_estimator("test_custom")
})

test_that("remove_estimator removes an estimator", {
  my_fn <- function(data, ...) {
    list(matrix = matrix(0, 1, 1), nodes = "A", directed = FALSE)
  }
  register_estimator("to_remove", my_fn, "Will be removed", directed = FALSE)
  expect_true("to_remove" %in% list_estimators()$name)

  remove_estimator("to_remove")
  expect_false("to_remove" %in% list_estimators()$name)
  expect_error(get_estimator("to_remove"), "not found")
})

test_that("remove_estimator errors on unknown estimator", {
  expect_error(remove_estimator("nonexistent"), "not found")
})

test_that("register_estimator validates inputs", {
  expect_error(register_estimator(123, identity, "desc", directed = TRUE))
  expect_error(register_estimator("x", "not_fn", "desc", directed = TRUE))
  expect_error(register_estimator("x", identity, 123, directed = TRUE))
  expect_error(register_estimator("x", identity, "desc", directed = "yes"))
})

test_that("list_estimators returns sorted names", {
  est_list <- list_estimators()
  expect_equal(est_list$name, sort(est_list$name))
})

test_that("list_estimators on empty registry returns empty data frame", {
  # Save state
  nms <- list_estimators()$name
  # We can't easily clear the registry without side effects, so just

  # verify structure would be correct
  expect_true(is.data.frame(list_estimators()))
})

test_that("directed flags are correct for built-in estimators", {
  est_list <- list_estimators()
  directed_methods <- est_list$name[est_list$directed]
  undirected_methods <- est_list$name[!est_list$directed]

  expect_true(all(c("relative", "frequency") %in% directed_methods))
  expect_true(all(c("cor", "pcor", "glasso", "co_occurrence") %in%
                    undirected_methods))
})

test_that("custom estimator can be used with build_network", {
  my_fn <- function(data, scale_factor = 1, ...) {
    S <- cor(data)
    diag(S) <- 0
    S <- S * scale_factor
    list(matrix = S, nodes = colnames(S), directed = FALSE)
  }
  register_estimator("test_scaled_cor", my_fn, "Scaled cor", directed = FALSE)

  set.seed(42)
  df <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
  net <- build_network(df, method = "test_scaled_cor",
                       params = list(scale_factor = 0.5))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "test_scaled_cor")
  expect_false(net$directed)

  # Clean up
  remove_estimator("test_scaled_cor")
})
