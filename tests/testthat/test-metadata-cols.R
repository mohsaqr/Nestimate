# ---- metadata_cols / state_cols decide what the estimator reads ----

test_that("a declared metadata column in wide data is not a state", {
  d <- data.frame(V1 = rep("A", 8), V2 = rep("B", 8),
                  b = rep(c("s1", "s2"), each = 4))
  net <- build_network(d, method = "relative", metadata_cols = "b")
  expect_identical(net$nodes$label, c("A", "B"))
  expect_identical(names(net$metadata), "b")
  by_state <- build_network(d, method = "relative", state_cols = c("V1", "V2"))
  expect_identical(by_state$nodes$label, c("A", "B"))
  expect_equal(by_state$weights, net$weights)
})

test_that("metadata columns stay out of grouped wide networks", {
  d <- data.frame(V1 = rep("A", 8), V2 = rep("B", 8),
                  b = rep(c("s1", "s2"), each = 4),
                  Opinion = rep(c("High", "Low"), 4))
  nets <- build_network(d, group = "Opinion", method = "relative",
                        metadata_cols = "b")
  expect_identical(nets[[1]]$nodes$label, c("A", "B"))
  expect_identical(names(nets[[1]]$metadata), "b")
  # blocking on the metadata column now tests the real A -> B network
  res <- suppressWarnings(permutation(nets, actor = "b", iter = 19L,
                                      seed = 1))[[1]]
  expect_identical(nrow(res$summary), 1L)
})

test_that("metadata columns are not variables of association networks", {
  set.seed(1)
  a <- data.frame(id = rep(1:10, each = 6), x = rnorm(60), y = rnorm(60),
                  z = rnorm(60))
  expect_identical(build_network(a, method = "glasso",
                                 metadata_cols = "id")$nodes$label,
                   c("x", "y", "z"))
  expect_identical(build_network(a, method = "cor",
                                 metadata_cols = "id")$nodes$label,
                   c("x", "y", "z"))
})
