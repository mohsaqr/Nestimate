# subtract_networks() / netdifference class

test_that("subtract_networks returns x - y with source matrices", {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B", "A"), V2 = c("B", "C", "B", "A", "C", "B"),
    V3 = c("C", "A", "C", "B", "A", "C"), stringsAsFactors = FALSE
  )
  a <- build_network(seqs, method = "relative")
  b <- build_network(seqs[1:4, ], method = "relative")
  d <- subtract_networks(a, b)

  expect_s3_class(d, "netdifference")
  expect_true(inherits(d, "netobject"))
  expect_equal(d$weights, a$weights - b$weights)
  expect_equal(d$difference_matrix, d$weights)
  expect_equal(d$x, a$weights)
  expect_equal(d$y, b$weights)
  expect_true(d$directed)
  expect_identical(d$method, "difference")
})

test_that("subtract_networks accepts plain matrices and validates inputs", {
  lab <- c("A", "B")
  m1 <- matrix(c(0, 2, 1, 0), 2, 2, dimnames = list(lab, lab))
  m2 <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(lab, lab))
  d <- subtract_networks(m1, m2)
  expect_equal(unname(d$weights), unname(m1 - m2))

  m3 <- matrix(0, 3, 3)
  expect_error(subtract_networks(m1, m3), "same number of nodes")
  m4 <- m2
  dimnames(m4) <- list(rev(lab), rev(lab))
  expect_error(subtract_networks(m1, m4), "same node labels")
})

test_that("print.netdifference is tidy and respects max_print", {
  lab <- c("A", "B", "C")
  m1 <- matrix(c(0, 3, 1, 2, 0, 1, 0, 4, 0), 3, 3, byrow = TRUE,
               dimnames = list(lab, lab))
  m2 <- matrix(c(0, 1, 1, 3, 0, 2, 0, 1, 0), 3, 3, byrow = TRUE,
               dimnames = list(lab, lab))
  d <- subtract_networks(m1, m2)
  out <- capture.output(print(d, max_print = 2L))
  expect_true(any(grepl("Network difference", out)))
  expect_true(any(grepl("more edges", out)))
})

test_that("as_netdifference dispatch: identity, net_bayes coercion, default error", {
  seqs <- data.frame(V1 = c("A", "B", "C", "A"), V2 = c("B", "C", "A", "B"),
                     stringsAsFactors = FALSE)
  a <- build_network(seqs, method = "relative")
  d <- subtract_networks(a, a)
  expect_identical(as_netdifference(d), d)
  expect_error(as_netdifference(42), "Cannot coerce")

  ebd <- subtract_networks(net_edge_betweenness(a), net_edge_betweenness(a))
  expect_s3_class(ebd, "netdifference")
  expect_true(all(ebd$weights == 0))
})
