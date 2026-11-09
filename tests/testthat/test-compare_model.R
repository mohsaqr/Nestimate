# compare_model() descriptive metrics on hand-computed vectors.

test_that("distance correlation is the square root of the ratio, 1 for identical vectors, NA for a constant", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 1, 4, 3, 6)
  # dCov^2 / sqrt(dVar_x^2 dVar_y^2) = 0.7561906 for these vectors; the
  # distance correlation is its square root (Szekely, Rizzo & Bakirov 2007)
  expect_equal(.distance_correlation(x, y), sqrt(0.7561906), tolerance = 1e-6)
  expect_equal(.distance_correlation(x, x), 1)
  expect_equal(.distance_correlation(x, 2 * x + 1), 1)
  expect_true(is.na(.distance_correlation(x, rep(1, 5))))
  expect_true(.distance_correlation(x, y) >= 0 && .distance_correlation(x, y) <= 1)
})

test_that("Canberra counts an edge present in only one network as 1 and skips edges absent from both", {
  x <- c(1, 0, 2, 0)
  y <- c(1, 3, 0, 0)
  # shared edge 0/2 = 0; one-sided edges 3/3 and 2/2 = 1 each; 0/0 skipped
  expect_equal(unname(.network_similarity(x, y, "canberra")), 2)
  expect_equal(unname(.network_similarity(x, x, "canberra")), 0)
  expect_equal(unname(.network_similarity(c(0, 0), c(0, 0), "canberra")), 0)
})

test_that("compare_model() carries the corrected Canberra and distance correlation", {
  x <- matrix(c(0, 1, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0), 4,
              byrow = TRUE, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  y <- x
  y["A", "C"] <- 3   # present only in y
  y["B", "A"] <- 0   # present only in x
  cm <- compare_model(x, y)
  sm <- cm$summary_metrics
  expect_equal(sm$value[sm$metric == "Canberra"], 2)
  dcor <- sm$value[sm$category == "Correlations" & grepl("Distance", sm$metric)]
  expect_length(dcor, 1L)
  expect_equal(dcor, .distance_correlation(as.vector(x), as.vector(y)))
})

test_that("compare_model metrics take direction from the network, not symmetry", {
  edge_count <- function(res) {
    metrics <- res$network_metrics
    metrics$x[metrics$metric == "Edge Count"]
  }
  density <- function(res) {
    metrics <- res$network_metrics
    metrics$x[metrics$metric == "Network Density"]
  }
  # Four arcs (A<->B, A<->C) make a symmetric matrix on a directed network.
  directed <- build_network(data.frame(T1 = c("A", "B", "A", "C"),
                                       T2 = c("B", "A", "C", "A")),
                            method = "frequency")
  res <- compare_model(directed, directed)
  expect_identical(edge_count(res), 4)
  expect_equal(density(res), 4 / 6)

  # The same pattern as a plain matrix has no flag, so symmetry stands in.
  symmetric <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3,
                      dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  res <- compare_model(symmetric, symmetric)
  expect_identical(edge_count(res), 2)
  expect_equal(density(res), 2 / 3)
  arcs <- symmetric
  arcs[2L, 1L] <- 0
  res <- compare_model(arcs, arcs)
  expect_identical(edge_count(res), 3)
  expect_equal(density(res), 3 / 6)

  # Each side uses its own flag.
  undirected <- directed
  undirected$directed <- FALSE
  res <- compare_model(directed, undirected)
  metrics <- res$network_metrics
  expect_identical(metrics$y[metrics$metric == "Edge Count"], 2)
  expect_identical(edge_count(res), 4)
})

# ---- Review fixes 2026-11-08 ----

test_that("scaling = 'max' divides by the largest absolute weight, keeping signs", {
  # All-negative x: max |w| = 4, so x / 4 keeps every sign (w / max(w)
  # divided by -1 and flipped them).
  x <- matrix(c(-1, -2, -3, -4), 2, 2)
  y <- matrix(c(1, 0.5, -2, 0), 2, 2)
  r <- compare_model(x, y, scaling = "max", network = FALSE)
  expect_equal(unname(r$matrices$x), matrix(c(-0.25, -0.5, -0.75, -1), 2, 2))
  # Mixed signs: y / max |y| = y / 2
  expect_equal(unname(r$matrices$y), matrix(c(0.5, 0.25, -1, 0), 2, 2))
})

test_that("scaling = 'log' refuses weights at or below zero with a classed error", {
  dn <- list(c("A", "B"), c("A", "B"))
  x <- matrix(c(0, 1, 2, 3), 2, 2, dimnames = dn)
  y <- matrix(c(1, 1, 2, 3), 2, 2, dimnames = dn)
  expect_error(compare_model(x, y, scaling = "log", network = FALSE),
               class = "nestimate_bad_input")
  ok <- compare_model(y, y, scaling = "log", network = FALSE)
  expect_equal(unname(ok$matrices$x), log(unname(y)))
})

test_that("heatmap subtitle matches its fill: blue where x > y", {
  skip_if_not_installed("ggplot2")
  dn <- list(c("A", "B"), c("A", "B"))
  x <- matrix(c(0, 1, 0, 0), 2, 2, dimnames = dn)
  y <- matrix(0, 2, 2, dimnames = dn)
  p <- plot(compare_model(x, y, network = FALSE), type = "heatmap")
  fills <- ggplot2::ggplot_build(p)$data[[1]]$fill
  expect_true("#4A6FE3" %in% fills)
  expect_match(p$labels$subtitle, "^Blue: x > y")
})
