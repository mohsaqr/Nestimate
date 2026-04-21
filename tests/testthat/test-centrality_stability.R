testthat::skip_on_cran()

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
                             measures = c("InStrength", "OutStrength"),
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


test_that("betweenness works with centrality_fn", {
  skip_if_not_installed("tna")
  skip_if_not_installed("igraph")
  net <- build_network(tna::group_regulation, method = "relative")
  my_fn <- function(mat) {
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "directed",
                                              weighted = TRUE)
    w_inv <- 1 / igraph::E(g)$weight
    list(Betweenness = igraph::betweenness(g, weights = w_inv))
  }
  cs <- centrality_stability(net, iter = 30,
                             measures = c("Betweenness"),
                             centrality_fn = my_fn, seed = 42)

  expect_s3_class(cs, "net_stability")
  expect_true("Betweenness" %in% names(cs$cs))
})

test_that("betweenness works without centrality_fn (built-in)", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30, measures = c("Betweenness"),
                              seed = 42)
  expect_s3_class(cs, "net_stability")
  expect_true("Betweenness" %in% names(cs$cs))
})

test_that("closeness measures work with centrality_fn", {
  skip_if_not_installed("tna")
  skip_if_not_installed("igraph")
  net <- build_network(tna::group_regulation, method = "relative")
  my_fn <- function(mat) {
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "directed",
                                              weighted = TRUE)
    w_inv <- 1 / igraph::E(g)$weight
    list(
      InCloseness = igraph::closeness(g, mode = "in", weights = w_inv),
      OutCloseness = igraph::closeness(g, mode = "out", weights = w_inv)
    )
  }
  cs <- centrality_stability(net, iter = 30,
                             measures = c("InCloseness", "OutCloseness"),
                             centrality_fn = my_fn, seed = 42)

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

  expect_error(centrality_stability("not_a_netobject"), "must be a netobject")
  expect_error(centrality_stability(net, iter = 0), "iter >= 2")
  expect_error(centrality_stability(net, threshold = 2), "threshold <= 1")
  expect_error(centrality_stability(net, certainty = -0.1), "certainty >= 0")
  expect_error(centrality_stability(net, measures = c("FakeMeasure")), "Unknown measures")
  expect_error(centrality_stability(net, method = "invalid"), "should be one of")
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


# ---- missing $data error (L81-82) ----

test_that("centrality_stability errors when $data is NULL", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  net$data <- NULL
  expect_error(centrality_stability(net, iter = 10L),
               "does not contain \\$data")
})


# ---- zero-variance warning and early return (L124-140) ----

test_that("centrality_stability warns and returns zeros when all measures have zero variance", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  # Force weights to a constant value (all rows identical → all centralities identical)
  n <- nrow(net$weights)
  net$weights[] <- 1 / n

  # OutStrength for relative = 1 always (zero variance); InStrength constant too
  expect_warning(
    cs <- centrality_stability(net, iter = 20L,
                               measures = c("InStrength", "OutStrength"),
                               seed = 1),
    "zero variance"
  )
  expect_s3_class(cs, "net_stability")
  expect_true(all(cs$cs == 0))
  # Correlation matrices should be all NA
  expect_true(all(is.na(cs$correlations[["InStrength"]])))
})


# ---- association path setup (L156-161) ----

test_that("centrality_stability works for cor (association) method", {
  set.seed(5)
  df <- as.data.frame(matrix(rpois(100 * 5, 10), nrow = 100))
  colnames(df) <- paste0("V", 1:5)
  net <- build_network(df, method = "cor")
  cs <- centrality_stability(net, iter = 20L,
                             measures = c("InStrength", "OutStrength"),
                             drop_prop = c(0.2, 0.4),
                             seed = 5)

  expect_s3_class(cs, "net_stability")
  expect_true(all(cs$cs %in% c(0, cs$drop_prop)))
})


# ---- association build_matrix function (L181-197) ----

test_that("centrality_stability association path tolerates estimator errors gracefully", {
  set.seed(9)
  df <- as.data.frame(matrix(rpois(60 * 4, 10), nrow = 60))
  colnames(df) <- paste0("V", 1:4)
  net <- build_network(df, method = "pcor")
  cs <- centrality_stability(net, iter = 20L,
                             measures = c("InStrength"),
                             drop_prop = c(0.3, 0.6),
                             seed = 9)

  expect_s3_class(cs, "net_stability")
  # CS should be 0 or a valid drop_prop value
  expect_true(cs$cs["InStrength"] %in% c(0, cs$drop_prop))
})


# ---- single-measure storage path (L216 and L222) ----

test_that("centrality_stability handles single measure correctly", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cs <- centrality_stability(net, iter = 30L,
                             measures = "InStrength",
                             drop_prop = c(0.1, 0.3, 0.5),
                             seed = 77)

  expect_s3_class(cs, "net_stability")
  expect_equal(cs$measures, "InStrength")
  expect_equal(ncol(cs$correlations[["InStrength"]]), 3L)
  expect_equal(nrow(cs$correlations[["InStrength"]]), 30L)
  # CS value must be 0 or in drop_prop
  expect_true(cs$cs["InStrength"] %in% c(0, cs$drop_prop))
})


# ---- CS-coefficient computation (L296-297) ----

test_that(".calculate_cs returns 0 when no prop_above meets certainty", {
  # Build a correlation matrix where certainty is never met
  iter <- 10L
  n_prop <- 3L
  corr_mat <- matrix(0, nrow = iter, ncol = n_prop)  # all zeros < threshold
  result <- Nestimate:::.calculate_cs(corr_mat, threshold = 0.7, certainty = 0.95,
                                      drop_prop = c(0.1, 0.3, 0.5))
  expect_equal(result, 0)
})

test_that(".calculate_cs returns max valid drop_prop when certainty is met", {
  iter <- 20L
  n_prop <- 3L
  corr_mat <- matrix(1, nrow = iter, ncol = n_prop)  # all ones >= threshold
  result <- Nestimate:::.calculate_cs(corr_mat, threshold = 0.7, certainty = 0.95,
                                      drop_prop = c(0.1, 0.3, 0.5))
  expect_equal(result, 0.5)
})


# ---- cograph_network input (L76) ----

test_that("centrality_stability accepts cograph_network input", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  cograph_net <- net
  class(cograph_net) <- "cograph_network"
  cs <- centrality_stability(cograph_net, iter = 20L,
                             measures = "InStrength",
                             seed = 1)
  expect_s3_class(cs, "net_stability")
})


# ---- net_centrality() (centrality_measures.R) ----

test_that("centrality.netobject returns correct directed defaults (L163-177)", {
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C")
  )
  net <- build_network(seqs, method = "relative")
  c1 <- net_centrality(net)
  expect_true(is.data.frame(c1))
  expect_equal(nrow(c1), 3)
  expect_true(all(c("InStrength", "OutStrength", "Betweenness") %in% names(c1)))
})

test_that("centrality.netobject returns correct undirected defaults (L163-177)", {
  set.seed(42)
  panel <- data.frame(V1 = rnorm(50), V2 = rnorm(50), V3 = rnorm(50))
  net_ud <- build_network(panel, method = "cor")
  c2 <- net_centrality(net_ud)
  expect_true(is.data.frame(c2))
  expect_true(all(c("Closeness", "Betweenness") %in% names(c2)))
})

test_that("centrality.netobject_group returns list of data frames (L185-188)", {
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C"),
    grp = c("X","X","X","Y","Y","Y")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  c3 <- net_centrality(nets)
  expect_true(is.list(c3))
  expect_equal(length(c3), 2)
  expect_true(all(vapply(c3, is.data.frame, logical(1))))
})

test_that(".betweenness returns zeros for n < 3 (L53)", {
  W <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("A","B"), c("A","B")))
  btw <- Nestimate:::.betweenness(W, directed = TRUE)
  expect_equal(unname(btw), c(0, 0))
  expect_equal(names(btw), c("A", "B"))
})

test_that(".compute_centralities handles external centrality_fn (L325-340)", {
  seqs <- data.frame(
    V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
    V3 = c("C","A","C","B")
  )
  net <- build_network(seqs, method = "relative")
  custom_fn <- function(mat) {
    list(MyMeasure = setNames(rowSums(abs(mat)), rownames(mat)))
  }
  c4 <- net_centrality(net, measures = c("InStrength", "MyMeasure"),
                    centrality_fn = custom_fn)
  expect_true("MyMeasure" %in% names(c4))
  expect_true(is.data.frame(c4))
})

test_that(".compute_centralities errors when external measure lacks centrality_fn (L325-329)", {
  seqs <- data.frame(
    V1 = c("A","B","A"), V2 = c("B","C","B"), V3 = c("C","A","C")
  )
  net <- build_network(seqs, method = "relative")
  expect_error(
    net_centrality(net, measures = c("InStrength", "BadMeasure")),
    "centrality_fn is required"
  )
})

# ---- Regression: non-square matrix $data (pre-2026-04-21 bug) ----

test_that("centrality_stability works when $data is a raw numeric matrix", {
  # For association methods (glasso/pcor/cor), build_network() stores $data
  # as a numeric matrix (not data.frame). When centrality_stability resamples
  # rows and re-invokes the estimator, .prepare_association_input was
  # failing the nrow==ncol check and silently producing NULL, which showed
  # up as all-NaN correlations or a zero-variance warning. Fixed by making
  # the matrix branch coerce non-square inputs to data.frame.
  skip_if_not_installed("glasso")
  set.seed(1)
  n <- 150; p <- 6
  Sigma <- diag(p)
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    Sigma[i, j] <- Sigma[j, i] <- 0.4 ^ (j - i)
  }
  L   <- chol(Sigma)
  df  <- as.data.frame(matrix(rnorm(n * p), n, p) %*% L)
  colnames(df) <- paste0("V", seq_len(p))
  net <- build_network(df, method = "glasso", params = list(nlambda = 20))
  # Critical invariant: $data IS a matrix for association methods
  expect_true(is.matrix(net$data))
  # Fix should allow the re-estimation loop to run without all-NaN fallout
  cs <- suppressWarnings(suppressMessages(
    centrality_stability(net, iter = 5L)
  ))
  expect_true(inherits(cs, "net_stability"))
  # At least one correlation must be finite — if the bug regressed, every
  # re-estimation returns NULL and all correlations are NaN.
  corrs <- cs$correlations
  if (is.data.frame(corrs)) corrs <- corrs$correlation
  expect_true(any(is.finite(unlist(corrs))))
})