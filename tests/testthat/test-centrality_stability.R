# ---- Tests for centrality_stability() ----

test_that("basic structure is correct", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 50,
                             measures = c("InStrength", "OutStrength"),
                             seed = 42)

  expect_s3_class(cs, "net_stability")
  expect_named(cs, c("cs", "correlations", "measures", "drop_prop",
                      "threshold", "certainty", "iter", "method"))
  expect_equal(cs$iter, 50L)
  expect_equal(cs$threshold, 0.7)
  expect_equal(cs$certainty, 0.95)
  expect_equal(cs$method, "pearson")
  expect_equal(cs$measures, c("InStrength", "OutStrength"))
})


test_that("CS coefficients are valid values from drop_prop", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 50,
                             measures = c("InStrength"),
                             seed = 42)

  expect_true(cs$cs["InStrength"] %in% c(0, cs$drop_prop))
})


test_that("correlation matrices have correct dimensions", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  dp <- seq(0.1, 0.5, by = 0.1)
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength", "Betweenness"),
                             drop_prop = dp, seed = 42)

  for (m in cs$measures) {
    expect_equal(nrow(cs$correlations[[m]]), 30)
    expect_equal(ncol(cs$correlations[[m]]), length(dp))
  }
})


test_that("correlation values are between -1 and 1", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 50,
                             measures = c("InStrength"),
                             seed = 42)

  vals <- cs$correlations[["InStrength"]]
  vals <- vals[!is.na(vals)]
  expect_true(all(vals >= -1 & vals <= 1))
})


test_that("seed produces reproducible results", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")

  cs1 <- centrality_stability(net, iter = 50,
                              measures = c("InStrength"), seed = 99)
  cs2 <- centrality_stability(net, iter = 50,
                              measures = c("InStrength"), seed = 99)

  expect_equal(cs1$correlations, cs2$correlations)
  expect_equal(cs1$cs, cs2$cs)
})


test_that("betweenness works with igraph", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("Betweenness"),
                             seed = 42)

  expect_s3_class(cs, "net_stability")
  expect_true("Betweenness" %in% names(cs$cs))
})


test_that("closeness measures work", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InCloseness", "OutCloseness"),
                             seed = 42)

  expect_s3_class(cs, "net_stability")
  expect_length(cs$cs, 2)
})


test_that("frequency method works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "frequency")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength", "OutStrength"),
                             seed = 42)

  expect_s3_class(cs, "net_stability")
})


test_that("loops = TRUE includes diagonal", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")

  cs_no <- centrality_stability(net, iter = 50,
                                measures = c("OutStrength"),
                                loops = FALSE, seed = 42)
  cs_yes <- centrality_stability(net, iter = 50,
                                 measures = c("OutStrength"),
                                 loops = TRUE, seed = 42)

  # With loops = TRUE on relative, OutStrength = 1 always → zero variance
  # → CS should be 0
  expect_equal(cs_yes$cs["OutStrength"], c(OutStrength = 0))
  # Without loops, OutStrength has variance → CS > 0
  expect_true(cs_no$cs["OutStrength"] > 0)
})


test_that("custom drop_prop works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  dp <- c(0.25, 0.5, 0.75)
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength"),
                             drop_prop = dp, seed = 42)

  expect_equal(cs$drop_prop, dp)
  expect_equal(ncol(cs$correlations[["InStrength"]]), 3)
})


test_that("custom threshold and certainty work", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 50,
                             measures = c("InStrength"),
                             threshold = 0.5, certainty = 0.90,
                             seed = 42)

  expect_equal(cs$threshold, 0.5)
  expect_equal(cs$certainty, 0.90)
})


test_that("spearman method works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength"),
                             method = "spearman", seed = 42)

  expect_equal(cs$method, "spearman")
})


test_that("input validation catches bad arguments", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")

  expect_error(centrality_stability("not_a_netobject"))
  expect_error(centrality_stability(net, iter = 0))
  expect_error(centrality_stability(net, threshold = 2))
  expect_error(centrality_stability(net, certainty = -0.1))
  expect_error(centrality_stability(net, measures = c("FakeMeasure")))
  expect_error(centrality_stability(net, method = "invalid"))
})


test_that("engagement dataset matches tna", {
  skip_if_not_installed("tna")
  net_e <- build_network(tna::engagement, method = "relative")
  m_e <- tna::tna(tna::engagement)

  cs <- centrality_stability(net_e, iter = 200,
                             measures = c("InStrength", "OutStrength"),
                             seed = 42)
  tcs <- tna:::estimate_cs(m_e, iter = 200,
                           measures = c("InStrength", "OutStrength"))

  # CS coefficients should match
  expect_equal(cs$cs["InStrength"], tcs$InStrength$cs_coefficient,
               ignore_attr = TRUE)
  expect_equal(cs$cs["OutStrength"], tcs$OutStrength$cs_coefficient,
               ignore_attr = TRUE)
})


test_that("print method works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength"), seed = 42)

  out <- capture.output(print(cs))
  expect_true(any(grepl("Centrality Stability", out)))
  expect_true(any(grepl("InStrength", out)))
})


test_that("summary method returns correct data frame", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength", "OutStrength"),
                             seed = 42)

  s <- summary(cs)
  expect_s3_class(s, "data.frame")
  expect_named(s, c("measure", "drop_prop", "mean_cor", "sd_cor", "prop_above"))
  expect_equal(nrow(s), 2 * length(cs$drop_prop))
})


test_that("plot method returns ggplot", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InStrength", "OutStrength"),
                             seed = 42)

  p <- plot(cs)
  expect_s3_class(p, "ggplot")
})
