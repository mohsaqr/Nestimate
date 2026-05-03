# Equivalence test_that() blocks extracted from
# tests/testthat/test-cograph_network.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

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

test_that(".as_netobject passes through netobject unchanged", {
  skip_if_pkg_broken("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  result <- Nestimate:::.as_netobject(net)
  expect_identical(result, net)
})

test_that("permutation works with cograph_network inputs", {
  skip_if_not_installed("cograph")
  skip_if_not_installed("tna")

  d1 <- tna::group_regulation[1:1000, ]
  d2 <- tna::group_regulation[1001:2000, ]
  net1 <- build_network(d1, method = "relative")
  net2 <- build_network(d2, method = "relative")

  # Convert one to cograph_network
  cg1 <- cograph::as_cograph(tna::tna(d1))
  skip_if(is.null(cg1$weights), "cograph_network has no $weights")

  # Both netobject (baseline)
  perm_base <- permutation(net1, net2, iter = 20, seed = 1)
  expect_s3_class(perm_base, "net_permutation")

  # cograph_network as x
  perm_cg <- permutation(cg1, net2, iter = 20, seed = 1)
  expect_s3_class(perm_cg, "net_permutation")
})

# ---- network_reliability ----

test_that("cluster_summary works with cograph_network", {
  skip_if_not_installed("cograph")
  skip_if_not_installed("tna")

  model <- tna::tna(tna::group_regulation)
  cg <- cograph::as_cograph(model)
  skip_if(is.null(cg$weights), "cograph_network has no $weights")

  states <- cg$nodes$label
  clusters <- list(
    A = states[1:3],
    B = states[4:length(states)]
  )

  cs <- cluster_summary(cg, clusters = clusters)
  expect_s3_class(cs, "mcml")
})

# ---- boot_glasso (requires glasso-compatible cograph_network) ----

test_that("boot_glasso works with cograph_network wrapping glasso netobject", {
  skip_if_not_installed("cograph")
  skip_if_pkg_broken("tna")

  # Build a glasso network, convert to cograph, then feed to boot_glasso
  freq <- convert_sequence_format(tna::group_regulation, format = "frequency")
  net <- build_network(freq, method = "glasso")

  # Create a fake cograph_network that wraps the glasso data
  cg <- structure(list(
    nodes = data.frame(
      id = seq_along(net$nodes$label),
      label = net$nodes$label,
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = integer(0), to = integer(0)),
    directed = FALSE,
    weights = net$weights,
    data = net$data,
    meta = list(source = "test", tna = list(method = "glasso"))
  ), class = c("cograph_network", "list"))

  boot <- boot_glasso(cg, iter = 20, cs_iter = 10, seed = 1,
                      centrality = c("strength", "expected_influence"))
  expect_s3_class(boot, "boot_glasso")
  expect_equal(boot$iter, 20L)
})


# ---- net_centrality() with cograph_network ----

test_that("cluster_network() works on cograph_network", {
  skip_if_not_installed("cograph")
  seqs <- data.frame(
    V1 = sample(LETTERS[1:4], 30, TRUE),
    V2 = sample(LETTERS[1:4], 30, TRUE),
    V3 = sample(LETTERS[1:4], 30, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  cg <- structure(list(
    weights = net$weights,
    nodes = net$nodes,
    edges = net$edges,
    directed = net$directed,
    data = net$data,
    meta = list(source = "test", tna = list(method = "relative"))
  ), class = c("cograph_network", "list"))

  grp <- cluster_network(cg, k = 2, method = "relative")
  expect_s3_class(grp, "netobject_group")
})


# ---- .coerce_sequence_input() with cograph_network ----

test_that(".coerce_sequence_input handles cograph_network", {
  skip_if_not_installed("cograph")
  seqs <- data.frame(
    V1 = sample(LETTERS[1:4], 20, TRUE),
    V2 = sample(LETTERS[1:4], 20, TRUE),
    V3 = sample(LETTERS[1:4], 20, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  cg <- structure(list(
    weights = net$weights,
    nodes = net$nodes,
    edges = net$edges,
    directed = net$directed,
    data = net$data,
    meta = list(source = "test", tna = list(method = "relative"))
  ), class = c("cograph_network", "list"))

  df <- Nestimate:::.coerce_sequence_input(cg)
  expect_true(is.data.frame(df))
  expect_equal(ncol(df), 3)
})