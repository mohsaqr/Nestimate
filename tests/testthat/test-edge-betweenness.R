test_that("net_edge_betweenness returns a betweenness netobject", {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B"), V2 = c("B", "C", "B", "A", "C"),
    V3 = c("C", "A", "C", "B", "A"), stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  eb  <- net_edge_betweenness(net)

  expect_true(inherits(eb, "netobject"))
  expect_true(inherits(eb, "cograph_network"))
  expect_identical(eb$method, "edge_betweenness")
  expect_identical(dim(eb$weights), dim(net$weights))
  expect_identical(dimnames(eb$weights), dimnames(net$weights))
  # betweenness is non-negative and only on existing edges
  expect_true(all(eb$weights >= 0))
  expect_true(all(eb$weights[net$weights == 0] == 0))
})

test_that("net_edge_betweenness flips probability weights by default", {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C"), V2 = c("B", "C", "B", "A"),
    V3 = c("C", "A", "C", "B"), stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  # Default must invert (distance = 1/prob), matching tna::betweenness_network.
  expect_equal(net_edge_betweenness(net)$weights,
               net_edge_betweenness(net, invert = TRUE)$weights)

  # On a network with skewed edge strengths, the flip changes which route is
  # the geodesic, so invert = TRUE and invert = FALSE diverge.
  skew <- matrix(c(0, 0.9, 0.1, 0,
                   0, 0,   0.5, 0.5,
                   0, 0.5, 0,   0.5,
                   0.5, 0, 0.5, 0), 4, 4, byrow = TRUE,
                 dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  sn <- build_network(skew, method = "relative")
  expect_false(isTRUE(all.equal(
    net_edge_betweenness(sn, invert = TRUE)$weights,
    net_edge_betweenness(sn, invert = FALSE)$weights
  )))
})

test_that("undirected networks give a symmetric betweenness matrix", {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B"), V2 = c("B", "C", "B", "A", "C"),
    V3 = c("C", "A", "C", "B", "A"), stringsAsFactors = FALSE
  )
  und <- build_network(seqs, method = "co_occurrence")
  eb  <- net_edge_betweenness(und)
  expect_true(isSymmetric(unname(eb$weights)))
  expect_false(eb$directed)
})

test_that("net_edge_betweenness maps over a netobject_group", {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B", "C"), V2 = c("B", "C", "B", "A", "C", "A"),
    V3 = c("C", "A", "C", "B", "A", "B"), stringsAsFactors = FALSE
  )
  n1  <- build_network(seqs[1:3, ], method = "relative")
  n2  <- build_network(seqs[4:6, ], method = "relative")
  grp <- structure(list(`Cluster 1` = n1, `Cluster 2` = n2),
                   class = "netobject_group")
  eb  <- net_edge_betweenness(grp)

  expect_true(inherits(eb, "netobject_group"))
  expect_identical(names(eb), c("Cluster 1", "Cluster 2"))
  expect_true(all(vapply(eb, inherits, logical(1), "netobject")))
})

test_that("net_edge_betweenness rejects unsupported input", {
  expect_error(net_edge_betweenness(1:5), "netobject")
  expect_error(net_edge_betweenness(matrix(0, 2, 2)), "netobject")
})

test_that("net_edge_betweenness carries net_edge_betweenness class without breaking netobject dispatch", {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B"), V2 = c("B", "C", "B", "A", "C"),
    V3 = c("C", "A", "C", "B", "A"), stringsAsFactors = FALSE
  )
  eb <- net_edge_betweenness(build_network(seqs, method = "relative"))
  expect_true(inherits(eb, "net_edge_betweenness"))
  expect_identical(class(eb)[1L], "net_edge_betweenness")
  # netobject accessors still resolve via inheritance
  expect_s3_class(extract_edges(eb), "data.frame")
  expect_silent(invisible(capture.output(print(eb))))
})

test_that("plot.net_edge_betweenness returns a ranked ggplot", {
  skip_if_not_installed("ggplot2")
  seqs <- as.data.frame(matrix(c(
    "A", "B", "C", "D", "E",
    "A", "B", "C", "D", "E",
    "B", "C", "D", "E", "A"), ncol = 5, byrow = TRUE),
    stringsAsFactors = FALSE)
  eb <- net_edge_betweenness(build_network(seqs, method = "relative"))
  p  <- plot(eb)
  expect_s3_class(p, "ggplot")
  # forest and delta styles also return a ggplot
  expect_s3_class(plot(eb, style = "forest"), "ggplot")
  pd <- plot(eb, style = "delta")
  expect_s3_class(pd, "ggplot")
  # edge deviations from the mean sum to ~0
  expect_true(abs(sum(pd$data$delta)) < 1e-8)
  # top_n trims the edge set
  p2 <- plot(eb, top_n = 3)
  expect_s3_class(p2, "ggplot")
  expect_equal(nrow(p2$data), 3L)
  # bars are sorted descending by betweenness
  expect_true(all(diff(p$data$weight) <= 0))
})
