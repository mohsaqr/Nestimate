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

# ---- permutation_test ----

test_that("permutation_test works with cograph_network inputs", {
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
  perm_base <- permutation_test(net1, net2, iter = 20, seed = 1)
  expect_s3_class(perm_base, "net_permutation")

  # cograph_network as x
  perm_cg <- permutation_test(cg1, net2, iter = 20, seed = 1)
  expect_s3_class(perm_cg, "net_permutation")
})

# ---- reliability ----

test_that("reliability works with cograph_network", {
  cg <- make_cograph_net()

  rel <- reliability(cg, iter = 20, seed = 1)

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

test_that("net_centrality() works on cograph_network", {
  cg <- make_cograph_net()
  cent <- net_centrality(cg)
  expect_true(is.data.frame(cent))
  expect_true(nrow(cent) > 0)
})


# ---- cluster_network() with cograph_network ----

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
