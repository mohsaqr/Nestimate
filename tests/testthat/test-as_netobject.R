# Boundary tests: psychnet (cograph_network) <-> Nestimate netobject.

# A faithful stand-in for a psychnet EBICglasso result, matching the live
# format: class c("psychnet", "cograph_network"), labelled square $weights,
# id/label/name $nodes, CHARACTER-labelled $edges, plus method-specific extras
# including the $kkt optimality certificate. Built by hand so the test needs no
# psychnet install.
make_psychnet_stub <- function() {
  labs <- c("V1", "V2", "V3")
  w <- matrix(c(0, 0.4, 0,
                0.4, 0, 0.2,
                0, 0.2, 0), 3, 3, dimnames = list(labs, labs))
  structure(
    list(
      weights = w,
      nodes = data.frame(id = 1:3, label = labs, name = labs,
                         stringsAsFactors = FALSE),
      edges = data.frame(from = c("V1", "V2"), to = c("V2", "V3"),
                         weight = c(0.4, 0.2), stringsAsFactors = FALSE),
      directed = FALSE,
      method = "EBICglasso",
      n = 300L,
      precision = solve(diag(3) - w / 2),
      lambda = 0.07,
      kkt = 0
    ),
    class = c("psychnet", "cograph_network")
  )
}

test_that("as_netobject() promotes a psychnet object to a dual-class netobject", {
  net <- as_netobject(make_psychnet_stub())
  expect_true(inherits(net, "netobject"))
  expect_true(inherits(net, "cograph_network"))
})

test_that("as_netobject.psychnet preserves the estimator name (not symmetry-inferred)", {
  net <- as_netobject(make_psychnet_stub())
  expect_identical(net$method, "EBICglasso")
  expect_identical(net$meta$tna$method, "EBICglasso")
  expect_identical(net$meta$source, "psychnet")
})

test_that("as_netobject.psychnet integer-indexes the edges", {
  net <- as_netobject(make_psychnet_stub())
  expect_true(is.integer(net$edges$from) || is.numeric(net$edges$from))
  expect_false(is.character(net$edges$from))
  # endpoints index into $nodes
  expect_true(all(net$edges$from >= 1L & net$edges$from <= nrow(net$nodes)))
})

test_that("as_netobject.psychnet keeps the KKT certificate and extras under $meta", {
  net <- as_netobject(make_psychnet_stub())
  expect_identical(net$meta$psychnet$kkt, 0)
  expect_identical(net$meta$psychnet$lambda, 0.07)
  expect_true(is.matrix(net$meta$psychnet$precision))
})

test_that("as_netobject() is idempotent on a real netobject", {
  net <- build_cor(data.frame(a = rnorm(60), b = rnorm(60), c = rnorm(60)))
  expect_identical(as_netobject(net), net)
})

test_that("as_netobject() errors on an unsupported object", {
  expect_error(as_netobject(list(1, 2, 3)), "psychnet, cograph_network, or netobject")
})

test_that("validate_netobject() passes a real Nestimate netobject", {
  net <- build_cor(data.frame(a = rnorm(60), b = rnorm(60), c = rnorm(60)))
  expect_invisible(validate_netobject(net))
  expect_true(validate_netobject(net))
})

test_that("validate_netobject() passes a psychnet object (character edges allowed)", {
  expect_true(validate_netobject(make_psychnet_stub()))
})

test_that("validate_netobject() and the adapter compose: psychnet -> netobject is valid", {
  expect_true(validate_netobject(as_netobject(make_psychnet_stub())))
})

test_that("validate_netobject() reports concrete violations", {
  bad <- make_psychnet_stub()
  bad$edges <- data.frame(src = "V1", dst = "V2", w = 0.4) # wrong column names
  expect_error(validate_netobject(bad), "\\$edges must be a data.frame")

  bad2 <- make_psychnet_stub()
  bad2$weights <- bad2$weights[, 1:2] # non-square
  expect_error(validate_netobject(bad2), "square")

  bad3 <- make_psychnet_stub()
  class(bad3) <- "list"
  expect_error(validate_netobject(bad3), "cograph_network")
})
