testthat::skip_on_cran()

# Tests for cograph_network → netobject coercion across Nestimate functions

# ---- Helper: build a cograph_network from synthetic data ----
make_cograph_net <- function() {
  skip_if_not_installed("cograph")
  skip_if_not_installed("tna")

  model <- tna::tna(tna::group_regulation)
  cg <- cograph::as_cograph(model)
  skip_if(is.null(cg$weights), "cograph_network has no $weights")
  cg
}

# ---- .as_netobject converter ----

test_that(".as_netobject converts cograph_network to netobject", {
  cg <- make_cograph_net()

  result <- Nestimate:::.as_netobject(cg)

  expect_s3_class(result, "netobject")
  expect_true(is.matrix(result$weights))
  expect_true(is.data.frame(result$nodes))
  expect_equal(nrow(result$weights), nrow(result$nodes))
  expect_equal(ncol(result$weights), nrow(result$nodes))
  expect_true(result$directed)
  expect_equal(result$method, "relative")  # TNA is directed → "relative"
})

test_that(".as_netobject decodes integer-encoded data", {
  cg <- make_cograph_net()

  result <- Nestimate:::.as_netobject(cg)

  if (!is.null(result$data)) {
    # All values should be character state labels, not integers
    vals <- unlist(result$data)
    non_na <- vals[!is.na(vals)]
    expect_true(all(non_na %in% result$nodes$label))
  }
})

test_that(".as_netobject errors on unsupported input", {
  expect_error(
    Nestimate:::.as_netobject("not_a_network"),
    "Expected a netobject or cograph_network"
  )
})

test_that(".as_netobject preserves matrix values from cograph_network", {
  cg <- make_cograph_net()
  result <- Nestimate:::.as_netobject(cg)

  # Matrix should be identical to the weights in the cograph_network
  expect_equal(unname(result$weights), unname(cg$weights))
})

# ---- bootstrap_network ----

test_that("bootstrap_network works with cograph_network", {
  cg <- make_cograph_net()

  boot <- bootstrap_network(cg, iter = 20, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$iter, 20L)
  expect_true(is.matrix(boot$mean))
  expect_true(is.matrix(boot$p_values))
})

# ---- permutation ----

test_that("network_reliability works with cograph_network", {
  cg <- make_cograph_net()

  rel <- network_reliability(cg, iter = 20, seed = 1)

  expect_s3_class(rel, "net_reliability")
  expect_equal(rel$iter, 20L)
  expect_true(nrow(rel$iterations) > 0)
})

# ---- centrality_stability ----

test_that("centrality_stability works with cograph_network", {
  cg <- make_cograph_net()

  cs <- centrality_stability(
    cg, iter = 10, drop_prop = c(0.3, 0.5), seed = 1,
    measures = c("InStrength", "OutStrength")
  )

  expect_s3_class(cs, "net_stability")
  expect_true(length(cs$cs) > 0)
})

# ---- build_mmm ----

test_that("build_mmm works with cograph_network", {
  cg <- make_cograph_net()

  mmm <- build_mmm(cg, k = 2, n_starts = 5, seed = 1)

  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 2L)
  expect_equal(length(mmm$models), 2)
})

# ---- cluster_summary (mcml.R) ----

test_that("net_centrality() works on cograph_network", {
  cg <- make_cograph_net()
  cent <- net_centrality(cg)
  expect_true(is.data.frame(cent))
  expect_true(nrow(cent) > 0)
})


# ---- cluster_network() with cograph_network ----
