test_that("certainty returns a net_bootstrap-format object", {
  set.seed(1)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:4], 40, TRUE), V2 = sample(LETTERS[1:4], 40, TRUE),
    V3 = sample(LETTERS[1:4], 40, TRUE), V4 = sample(LETTERS[1:4], 40, TRUE)
  )
  net <- build_network(seqs, method = "relative")
  cert <- certainty(net)

  expect_true(inherits(cert, "net_certainty"))
  expect_true(inherits(cert, "net_bootstrap"))   # drop-in
  # exact net_bootstrap slot set
  boot_slots <- c("original", "mean", "sd", "p_values", "significant",
                  "ci_lower", "ci_upper", "cr_lower", "cr_upper", "summary",
                  "model", "method", "params", "iter", "ci_level", "inference",
                  "consistency_range", "edge_threshold", "ci_method")
  expect_true(all(boot_slots %in% names(cert)))
  # summary has the same columns as a bootstrap stability summary
  expect_true(all(c("from", "to", "weight", "mean", "sd", "p_value", "sig",
                    "ci_lower", "ci_upper", "cr_lower", "cr_upper")
                  %in% names(summary(cert))))
})

test_that("inherited net_bootstrap summary method works (true drop-in)", {
  seqs <- data.frame(V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
                     V3 = c("C","A","C","B","A"))
  net <- build_network(seqs, method = "relative")
  cert <- certainty(net)
  # summary() dispatches to summary.net_bootstrap and returns the data frame
  expect_identical(summary(cert), cert$summary)
  expect_s3_class(cert$model, "netobject")
})

test_that("posterior mean and ci match the closed-form Beta", {
  seqs <- data.frame(V1 = c("A","A","B","B","C","C"),
                     V2 = c("B","C","A","C","A","B"))
  net <- build_network(seqs, method = "relative")
  cert <- certainty(net, prior = 0.5, ci_level = 0.05)

  nodes <- net$nodes$label
  C <- net$frequency_matrix[nodes, nodes]
  a <- C + 0.5
  rowA <- rowSums(a)
  expected_mean <- a / rowA
  expect_equal(cert$mean, expected_mean, tolerance = 1e-12)

  b <- matrix(rowA, length(nodes), length(nodes)) - a
  expect_equal(cert$ci_lower, matrix(qbeta(0.025, a, b), length(nodes),
               dimnames = dimnames(C)), tolerance = 1e-12)
})

test_that("stability p-value is posterior mass outside the consistency band", {
  seqs <- data.frame(V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
                     V3 = c("C","A","C","B","A"))
  net <- build_network(seqs, method = "relative")
  cert <- certainty(net, consistency_range = c(0.75, 1.25))
  nodes <- net$nodes$label
  C <- net$frequency_matrix[nodes, nodes]; W <- net$weights[nodes, nodes]
  a <- C + 0.5; b <- matrix(rowSums(a), length(nodes), length(nodes)) - a
  cr_lo <- pmin(W * 0.75, W * 1.25); cr_hi <- pmax(W * 0.75, W * 1.25)
  expected_p <- pbeta(cr_lo, a, b) + (1 - pbeta(cr_hi, a, b))
  expect_equal(cert$p_values, matrix(expected_p, length(nodes),
               dimnames = dimnames(C)), tolerance = 1e-12)
})

test_that("certainty rejects non-relative networks", {
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(120), ncol = 4))
  net <- build_network(d, method = "cor")
  expect_error(certainty(net), "transition-probability network")
})

test_that("certainty rejects scaled networks (posterior is on probability scale)", {
  set.seed(1)
  seqs <- data.frame(V1 = sample(LETTERS[1:4], 40, TRUE),
                     V2 = sample(LETTERS[1:4], 40, TRUE),
                     V3 = sample(LETTERS[1:4], 40, TRUE))
  net <- build_network(seqs, method = "relative", scaling = "max")
  expect_error(certainty(net), "unscaled transition network")
  # unscaled is accepted
  expect_s3_class(certainty(build_network(seqs, method = "relative")),
                  "net_certainty")
})

test_that("certainty dispatches over netobject_group", {
  s <- data.frame(V1 = c("A","B","A","C","B","A"),
                  V2 = c("B","C","B","A","C","B"),
                  V3 = c("C","A","C","B","A","C"),
                  grp = c("X","X","X","Y","Y","Y"))
  nets <- build_network(s, method = "relative", group = "grp")
  cg <- certainty(nets)
  expect_true(inherits(cg, "net_certainty_group"))
  expect_true(inherits(cg, "net_bootstrap_group"))
  expect_length(cg, 2L)
  expect_true(inherits(cg[[1]], "net_certainty"))
})
