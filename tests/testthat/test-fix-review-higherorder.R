testthat::skip_on_cran()

test_that("large uniform hyperedges retain normalized tensor eigenvectors", {
  lapply(c(150L, 400L), function(n) {
    hg <- bipartite_groups(data.frame(player = seq_len(n), session = "G"),
                           "player", "session")
    expect_no_warning(result <- hypergraph_centrality(hg, type = c("Z", "H")))
    lapply(result, function(v) {
      expect_equal(unname(v), rep(1 / sqrt(n), n), tolerance = 1e-8)
      expect_equal(sum(v^2), 1, tolerance = 1e-8)
    })
  })
})

test_that("large asymmetric uniform hyperedges solve the analytic sunflower", {
  # Three equal-size petals share only their hub. Z's hub/petal ratio is
  # sqrt(3); H's ratio is 3^(1/k), from their respective eigen-equations.
  k <- 400L
  edges <- lapply(seq_len(3), function(j) {
    c("hub", paste0("p", j, "_", seq_len(k - 1L)))
  })
  hg <- bipartite_groups(data.frame(player = unlist(edges),
                                    session = rep(seq_len(3), lengths(edges))),
                         "player", "session")
  expect_no_warning(result <- hypergraph_centrality(hg, type = c("Z", "H")))
  expect_equal(unname(result$Z["hub"] / result$Z["p1_1"]), sqrt(3), tolerance = 1e-6)
  expect_equal(unname(result$H["hub"] / result$H["p1_1"]), 3^(1/k), tolerance = 1e-6)
  expect_equal(sum(result$Z^2), 1, tolerance = 1e-8)
  expect_equal(sum(result$H^2), 1, tolerance = 1e-8)
})

test_that("MOGen summaries preserve finite path DOFs beyond integer range", {
  set.seed(4)
  seqs <- replicate(100, sample(LETTERS[1:10], 14, TRUE), simplify = FALSE)
  fit <- build_mogen(seqs, max_order = 10)
  expect_gt(max(fit$layer_dof), .Machine$integer.max)
  expect_no_warning(invisible(capture.output(tab <- summary(fit))))
  expect_equal(tab$layer_dof, unname(fit$layer_dof))
  expect_equal(tab$cum_dof, unname(fit$dof))
})

test_that("case-dropping assesses newly appearing association edges by default", {
  set.seed(19)
  d <- as.data.frame(matrix(rnorm(80 * 5), 80, 5))
  d$V2 <- d$V1 + rnorm(80, sd = .7)
  d$V3 <- .4 * d$V1 + rnorm(80)
  net <- build_network(d, method = "cor", threshold = .2)
  set.seed(2)
  keep <- sample(seq_len(nrow(d)), 40)
  sub <- build_network(d[keep, ], method = "cor", threshold = .2)
  expect_true(any(net$weights == 0 & sub$weights != 0))
  all_cells <- casedrop_reliability(net, iter = 2, drop_prop = .5, seed = 2)
  existing <- casedrop_reliability(net, iter = 2, drop_prop = .5, seed = 2,
                                   target = "existing")
  expect_identical(attr(all_cells, "target"), "all")
  expect_identical(attr(existing, "target"), "existing")
  expect_equal(all_cells$n_edges, 20L)
  expect_equal(unname(all_cells$metrics$max_abs_dev[1, 1]),
               max(abs(net$weights - sub$weights)))
  expect_lt(existing$metrics$max_abs_dev[1, 1],
            all_cells$metrics$max_abs_dev[1, 1])
})

test_that("case-dropping target auto preserves only guaranteed structural zeros", {
  seqs <- data.frame(V1 = c("A", "A", "B", "B", "C", "C"),
                     V2 = c("A", "B", "A", "C", "B", "C"),
                     V3 = c("B", "C", "C", "A", "A", "B"))
  net <- build_network(seqs, method = "relative")
  out <- casedrop_reliability(net, iter = 2, drop_prop = .5, seed = 4)
  expect_identical(attr(out, "target"), "existing")
  expect_named(out,c("cs","summary","metrics","correlations","drop_prop",
    "threshold","certainty","iter","method","include_diag","n_cases","n_edges"))
  thresholded <- build_network(seqs, method = "relative", threshold = .1)
  full <- casedrop_reliability(thresholded, iter = 2, drop_prop = .5, seed = 4)
  expect_identical(attr(full, "target"), "all")
  grouped <- build_network(transform(seqs, group = rep(c("x", "y"), each = 3)),
                           group = "group", method = "relative")
  result <- casedrop_reliability(grouped, iter = 2, drop_prop = .5, seed = 4,
                                 target = "all")
  expect_true(all(vapply(result, function(x) attr(x, "target") == "all", logical(1))))
})
