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
