# Equivalence test_that() blocks extracted from
# tests/testthat/test-boot_glasso.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- Tests for boot_glasso() ----
skip_on_cran()

# ---- Helper: generate structured data with known edges ----
.make_test_data <- function(n = 200, p = 5, seed = 42) {
  set.seed(seed)
  # Create correlated data via Cholesky decomposition
  Sigma <- diag(p)
  # Add some correlations between adjacent variables
  for (i in seq_len(p - 1)) {
    Sigma[i, i + 1] <- 0.4
    Sigma[i + 1, i] <- 0.4
  }
  L <- chol(Sigma)
  mat <- matrix(rnorm(n * p), n, p) %*% L
  df <- as.data.frame(mat)
  names(df) <- paste0("V", seq_len(p))
  df
}

# Small run settings for speed
SMALL_ITER <- 20L
SMALL_CS_ITER <- 10L
SMALL_CS_DROP <- c(0.25, 0.5, 0.75)
FAST_CENT <- c("strength", "expected_influence")


# ========================================================================
# Input validation
# ========================================================================

test_that("boot_glasso original network matches bootnet EBICglasso", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  # Run across 5 seeds
  for (s in c(100, 200, 300, 400, 500)) {
    set.seed(s)
    n <- 200
    x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
    x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
    x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
    df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

    bn_net <- suppressWarnings(suppressMessages(
      bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                                corMethod = "cor", missing = "listwise",
                                verbose = FALSE)
    ))

    bg <- boot_glasso(df, iter = 10, cs_iter = 10, cs_drop = c(0.5),
                       centrality = "strength", gamma = 0.5, seed = 1)

    # Original networks should be very close (same EBIC/glasso)
    max_diff <- max(abs(bn_net$graph - bg$original_pcor))
    expect_true(max_diff < 0.01,
                info = sprintf("seed=%d, max_diff=%.6f", s, max_diff))
  }
})

test_that("boot_glasso edge CIs correlate highly with bootnet", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  bn_net <- suppressWarnings(suppressMessages(
    bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                              corMethod = "cor", missing = "listwise",
                              verbose = FALSE)
  ))

  set.seed(1)
  bn_boot <- suppressWarnings(suppressMessages(
    bootnet::bootnet(bn_net, nBoots = 200, type = "nonparametric",
                      nCores = 1, verbose = FALSE,
                      statistics = c("edge", "strength"))
  ))

  bg <- boot_glasso(df, iter = 200, cs_iter = 10, cs_drop = c(0.5),
                     centrality = "strength", gamma = 0.5, seed = 1)

  # Extract bootnet edge CIs
  bn_edges <- bn_boot$bootTable[bn_boot$bootTable$type == "edge", ]
  edge_ids <- unique(bn_edges$id)
  bn_ci <- do.call(rbind, lapply(edge_ids, function(eid) {
    vals <- bn_edges$value[bn_edges$id == eid]
    data.frame(edge = eid, ci_lower = quantile(vals, 0.025),
               ci_upper = quantile(vals, 0.975), stringsAsFactors = FALSE)
  }))

  bg_ci <- bg$edge_ci
  bg_ci$bn_id <- gsub(" -- ", "--", bg_ci$edge)
  merged <- merge(bn_ci, bg_ci, by.x = "edge", by.y = "bn_id")

  # Edge CIs should correlate highly (r > 0.99)
  edge_ci_r <- cor(merged$ci_lower.x, merged$ci_lower.y)
  expect_true(edge_ci_r > 0.99,
              info = sprintf("Edge CI r = %.4f", edge_ci_r))

  # Inclusion probabilities should correlate highly
  bn_incl <- vapply(edge_ids, function(eid) {
    mean(bn_edges$value[bn_edges$id == eid] != 0)
  }, numeric(1))
  bg_incl <- bg$edge_inclusion[gsub("--", " -- ", names(bn_incl))]
  incl_r <- cor(bn_incl, bg_incl)
  expect_true(incl_r > 0.99,
              info = sprintf("Inclusion r = %.4f", incl_r))
})

test_that("boot_glasso strength CIs correlate highly with bootnet", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  bn_net <- suppressWarnings(suppressMessages(
    bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                              corMethod = "cor", missing = "listwise",
                              verbose = FALSE)
  ))

  set.seed(1)
  bn_boot <- suppressWarnings(suppressMessages(
    bootnet::bootnet(bn_net, nBoots = 200, type = "nonparametric",
                      nCores = 1, verbose = FALSE,
                      statistics = c("edge", "strength"))
  ))

  bg <- boot_glasso(df, iter = 200, cs_iter = 10, cs_drop = c(0.5),
                     centrality = "strength", gamma = 0.5, seed = 1)

  bn_str <- bn_boot$bootTable[bn_boot$bootTable$type == "strength", ]
  bn_str_ci <- do.call(rbind, lapply(unique(bn_str$node1), function(nd) {
    vals <- bn_str$value[bn_str$node1 == nd]
    data.frame(node = nd, ci_lower = quantile(vals, 0.025),
               ci_upper = quantile(vals, 0.975), stringsAsFactors = FALSE)
  }))

  bg_str_ci <- bg$centrality_ci$strength
  str_merged <- merge(bn_str_ci, bg_str_ci, by = "node")
  str_ci_r <- cor(str_merged$ci_lower.x, str_merged$ci_lower.y)
  expect_true(str_ci_r > 0.99,
              info = sprintf("Strength CI r = %.4f", str_ci_r))
})

test_that("boot_glasso CS-coefficient matches bootnet corStability", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  # Run 5 seeds, require 80% agreement within 0.15
  cs_diffs <- numeric(5)

  for (i in seq_len(5)) {
    s <- i * 100
    set.seed(s)
    n <- 200
    x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
    x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
    x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
    df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

    bn_net <- suppressWarnings(suppressMessages(
      bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                                corMethod = "cor", missing = "listwise",
                                verbose = FALSE)
    ))

    set.seed(s + 1000)
    bn_case <- suppressWarnings(suppressMessages(
      bootnet::bootnet(bn_net, nBoots = 500, type = "case", nCores = 1,
                        verbose = FALSE, statistics = c("strength"),
                        caseMin = 0.05, caseMax = 0.75, caseN = 10)
    ))
    bn_cs <- suppressWarnings(suppressMessages(
      bootnet::corStability(bn_case, verbose = FALSE)
    ))

    bg <- boot_glasso(df, iter = 50, cs_iter = 500,
                       cs_drop = seq(0.05, 0.75, length.out = 10),
                       centrality = "strength", gamma = 0.5, seed = s)

    cs_diffs[i] <- abs(as.numeric(bn_cs["strength"]) -
                        as.numeric(bg$cs_coefficient["strength"]))
  }

  # At least 80% should match within 0.15
  pct_match <- mean(cs_diffs <= 0.15)
  expect_true(pct_match >= 0.8,
              info = sprintf("CS match: %.0f%%, diffs: %s",
                             pct_match * 100,
                             paste(round(cs_diffs, 2), collapse = ", ")))

  # Mean diff should be small
  expect_true(mean(cs_diffs) < 0.1,
              info = sprintf("Mean CS diff = %.3f", mean(cs_diffs)))
})

test_that("boot_glasso strength matches bootnet strength definition", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.7 * x1 + rnorm(n, sd = 0.5)
  x3 <- 0.5 * x2 + rnorm(n, sd = 0.5); x4 <- rnorm(n)
  x5 <- 0.6 * x4 + rnorm(n, sd = 0.5)
  df <- data.frame(V1 = x1, V2 = x2, V3 = x3, V4 = x4, V5 = x5)

  bn_net <- suppressWarnings(suppressMessages(
    bootnet::estimateNetwork(df, default = "EBICglasso", tuning = 0.5,
                              corMethod = "cor", missing = "listwise",
                              verbose = FALSE)
  ))

  bg <- boot_glasso(df, iter = 10, cs_iter = 10, cs_drop = c(0.5),
                     centrality = "strength", gamma = 0.5, seed = 1)

  # bootnet uses qgraph::centrality()$InDegree for strength
  bn_str <- qgraph::centrality(bn_net$graph)$InDegree
  bg_str <- bg$original_centrality$strength

  # Should match within tolerance (tiny diff from pcor rounding)
  expect_equal(unname(bg_str), unname(bn_str), tolerance = 0.01)
})


# ========================================================================
# Coverage gap: netobject without $data (L141-142)
# ========================================================================

