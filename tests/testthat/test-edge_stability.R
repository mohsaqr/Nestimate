# ---- edge_stability() tests ----------------------------------------------

testthat::skip_on_cran()

# Helper -------------------------------------------------------------------

.es_seq <- function() {
  set.seed(7)
  data.frame(
    V1 = sample(LETTERS[1:5], 80, TRUE),
    V2 = sample(LETTERS[1:5], 80, TRUE),
    V3 = sample(LETTERS[1:5], 80, TRUE),
    stringsAsFactors = FALSE
  )
}

# Structure ----------------------------------------------------------------

test_that("returns a net_edge_stability with required fields", {
  net <- build_network(.es_seq(), method = "relative")
  es <- edge_stability(net, iter = 30, drop_prop = c(0.1, 0.3, 0.5), seed = 1)
  expect_s3_class(es, "net_edge_stability")
  expect_true(all(c("cs", "summary", "metrics", "correlations",
                    "drop_prop", "threshold", "certainty", "iter",
                    "method", "include_diag", "n_cases", "n_edges")
                  %in% names(es)))
  expect_equal(dim(es$correlations), c(30L, 3L))
  expect_equal(es$n_cases, nrow(.es_seq()))
  expect_true(es$cs >= 0 && es$cs < 1)
})

test_that("metrics contains four per-iteration matrices of correct shape", {
  net <- build_network(.es_seq(), method = "relative")
  es <- edge_stability(net, iter = 25, drop_prop = c(0.2, 0.5), seed = 1)
  expect_named(es$metrics,
               c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev"))
  for (nm in names(es$metrics)) {
    expect_equal(dim(es$metrics[[nm]]), c(25L, 2L))
  }
  # correlations field is an alias for metrics$correlation
  expect_identical(es$correlations, es$metrics$correlation)
})

test_that("summary table has 4 metrics x n_prop rows with expected columns", {
  net <- build_network(.es_seq(), method = "relative")
  es <- edge_stability(net, iter = 20, drop_prop = c(0.2, 0.5, 0.7), seed = 1)
  expect_s3_class(es$summary, "data.frame")
  expect_equal(nrow(es$summary), 4L * 3L)  # 4 metrics x 3 drop_prop
  expect_true(all(c("metric", "drop_prop", "mean", "sd", "median",
                    "mad", "q025", "q975") %in% names(es$summary)))
})

test_that("mean_abs_dev and max_abs_dev are non-negative; correlation in [-1, 1]", {
  net <- build_network(.es_seq(), method = "relative")
  es <- edge_stability(net, iter = 30, drop_prop = 0.3, seed = 1)
  expect_true(all(es$summary$mean[es$summary$metric == "mean_abs_dev"] >= 0))
  expect_true(all(es$summary$mean[es$summary$metric == "max_abs_dev"] >= 0))
  cor_means <- es$summary$mean[es$summary$metric == "correlation"]
  expect_true(all(cor_means >= -1 & cor_means <= 1))
})

# Diagonal handling --------------------------------------------------------

test_that("include_diag controls whether self-loops enter the edge vector", {
  net <- build_network(.es_seq(), method = "relative")
  es_no  <- edge_stability(net, iter = 20, drop_prop = 0.2,
                           include_diag = FALSE, seed = 1)
  es_yes <- edge_stability(net, iter = 20, drop_prop = 0.2,
                           include_diag = TRUE, seed = 1)
  n <- length(net$nodes$label)
  expect_equal(es_no$n_edges, n * (n - 1L))
  expect_equal(es_yes$n_edges, n * n)
})

# Reproducibility ----------------------------------------------------------

test_that("same seed yields identical correlation matrix", {
  net <- build_network(.es_seq(), method = "relative")
  es1 <- edge_stability(net, iter = 15, drop_prop = 0.3, seed = 42)
  es2 <- edge_stability(net, iter = 15, drop_prop = 0.3, seed = 42)
  expect_equal(es1$correlations, es2$correlations)
  expect_equal(es1$cs, es2$cs)
})

# Large drop lowers correlation -------------------------------------------

test_that("mean correlation decreases as drop proportion increases", {
  net <- build_network(.es_seq(), method = "relative")
  es <- edge_stability(net, iter = 50,
                       drop_prop = c(0.1, 0.5, 0.9), seed = 3)
  mean_corr <- colMeans(es$correlations, na.rm = TRUE)
  # Monotone non-increasing (allow small MC wiggle)
  expect_true(mean_corr[1] >= mean_corr[2] - 0.05)
  expect_true(mean_corr[2] >= mean_corr[3] - 0.05)
})

# Group dispatch via mcml --------------------------------------------------

test_that("mcml input dispatches: one CS per sub-network", {
  edges <- data.frame(
    from   = c("a1", "a1", "a2", "b1", "b2", "b1", "a1", "b2", "a2", "b1"),
    to     = c("a2", "b1", "a1", "b2", "b1", "a2", "b2", "a1", "a1", "b1"),
    weight = 1,
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  mc <- build_mcml(edges, clusters)
  es <- edge_stability(mc, iter = 30, drop_prop = c(0.2, 0.5), seed = 2)
  expect_s3_class(es, "net_edge_stability_group")
  expect_true(all(c("macro", "A", "B") %in% names(es)))
  for (nm in names(es)) {
    expect_s3_class(es[[nm]], "net_edge_stability")
  }
})

# Input validation --------------------------------------------------------

test_that("rejects non-netobject input", {
  expect_error(edge_stability(matrix(0, 3, 3)), "netobject")
})

test_that("rejects netobject without \\$data", {
  net <- build_network(.es_seq(), method = "relative")
  net$data <- NULL
  expect_error(edge_stability(net), "does not contain \\$data")
})

test_that("drop_prop must be in (0, 1)", {
  net <- build_network(.es_seq(), method = "relative")
  expect_error(edge_stability(net, drop_prop = 0))
  expect_error(edge_stability(net, drop_prop = 1))
  expect_error(edge_stability(net, drop_prop = 1.2))
})

# Print method -------------------------------------------------------------

test_that("print methods run without error", {
  net <- build_network(.es_seq(), method = "relative")
  es <- edge_stability(net, iter = 10, drop_prop = 0.3, seed = 1)
  expect_invisible(print(es))

  edges <- data.frame(
    from = c("a1", "a2", "b1", "b2"),
    to   = c("a2", "a1", "b2", "b1"),
    stringsAsFactors = FALSE
  )
  mc <- build_mcml(edges, list(A = c("a1","a2"), B = c("b1","b2")))
  esg <- edge_stability(mc, iter = 5, drop_prop = 0.3, seed = 1)
  expect_invisible(print(esg))
})
