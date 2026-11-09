make_net <- function() {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B", "D"), V2 = c("B", "C", "B", "A", "C", "A"),
    V3 = c("C", "A", "C", "B", "A", "B"), V4 = c("D", "B", "A", "C", "D", "C"),
    stringsAsFactors = FALSE
  )
  build_network(seqs, method = "relative")
}

test_that("net_prune threshold removes weak edges and records them", {
  net    <- make_net()
  pruned <- net_prune(net, method = "threshold", threshold = 0.2)

  expect_true(inherits(pruned, "netobject"))
  pr <- attr(pruned, "pruning")
  expect_false(is.null(pr))
  expect_true(pr$active)
  # every removed edge had weight <= cut-off in the original
  expect_true(all(pruned$weights[net$weights > 0 & net$weights <= 0] == 0))
  # pruned weights are a subset: nothing gained, some lost
  expect_true(all((pruned$weights > 0) <= (net$weights > 0)))
  # edge table re-synced
  expect_identical(pruned$n_edges, sum(pruned$weights != 0))
})

test_that("diagonal self-loops are data and are never removed", {
  # a self-loop weaker than the threshold must still be retained
  W <- matrix(c(0.05, 0.95, 0,
                0.40, 0.00, 0.60,
                0.50, 0.50, 0.00), 3, 3, byrow = TRUE,
              dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  net    <- build_network(W, method = "relative")
  before <- net$weights["A", "A"]
  pruned <- net_prune(net, method = "threshold", threshold = 0.1)
  expect_gt(before, 0)                         # the self-loop exists...
  expect_equal(pruned$weights["A", "A"], before)  # ...and survives pruning
  # no removed edge is ever on the diagonal
  det <- net_pruning_details(pruned)
  expect_false(any(det$from == det$to))
})

test_that("deprune and reprune are exact inverses without recomputation", {
  net    <- make_net()
  pruned <- net_prune(net, threshold = 0.2)
  undone <- net_deprune(pruned)
  redone <- net_reprune(undone)

  expect_equal(undone$weights, net$weights)
  expect_equal(redone$weights, pruned$weights)
  expect_false(attr(undone, "pruning")$active)
  expect_true(attr(redone, "pruning")$active)
})

test_that("pruning state machine guards against misuse", {
  net    <- make_net()
  pruned <- net_prune(net, threshold = 0.2)

  expect_error(net_prune(pruned), "already been pruned")
  expect_error(net_deprune(net), "must have been pruned")
  expect_error(net_reprune(pruned), "already active")     # still active
  expect_error(net_deprune(net_deprune(pruned)), "already inactive")
})

test_that("net_pruning_details returns a tidy table with metadata", {
  net <- make_net()
  det <- net_pruning_details(net_prune(net, threshold = 0.2))

  expect_true(inherits(det, "net_pruning_details"))
  expect_true(inherits(det, "data.frame"))
  expect_identical(names(det), c("from", "to", "weight"))
  expect_identical(attr(det, "method"), "threshold")
  expect_identical(attr(det, "num_removed"), nrow(det))
  expect_error(net_pruning_details(net), "must have been pruned")
})

test_that("disparity and lowest methods run and stay connected", {
  net <- make_net()
  for (m in c("disparity", "lowest")) {
    pruned <- if (m == "disparity") {
      net_prune(net, method = "disparity", level = 0.5)
    } else {
      net_prune(net, method = "lowest", lowest = 0.3)
    }
    expect_true(inherits(pruned, "netobject"))
    expect_true(all((pruned$weights > 0) <= (net$weights > 0)))
  }
})

test_that("net_prune maps over a netobject_group", {
  net <- make_net()
  grp <- structure(list(`Cluster 1` = net, `Cluster 2` = net),
                   class = "netobject_group")
  pruned <- net_prune(grp, threshold = 0.2)

  expect_true(inherits(pruned, "netobject_group"))
  expect_true(all(vapply(pruned, function(e) !is.null(attr(e, "pruning")),
                         logical(1))))
  det <- net_pruning_details(pruned)
  expect_named(det, c("Cluster 1", "Cluster 2"))
})

test_that("net_prune rejects unsupported input", {
  expect_error(net_prune(1:5), "netobject")
})

test_that("disparity filter keeps the absent edges of a sink or source state at 0, never NA", {
  # D is a sink (never left), S a source (never entered)
  seqs <- data.frame(
    V1 = c("S", "S", "S", "S"), V2 = c("A", "B", "A", "C"),
    V3 = c("B", "C", "C", "A"), V4 = c("C", "A", "B", "B"),
    V5 = c("D", "D", "D", "D"),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "frequency")
  expect_true(all(net$weights["D", ] == 0))
  expect_true(all(net$weights[, "S"] == 0))
  pruned <- net_prune(net, method = "disparity", level = 0.5)
  expect_false(anyNA(pruned$weights))
  expect_true(all(pruned$weights["D", ] == 0))
  expect_true(all(pruned$weights[, "S"] == 0))
  expect_true(all((pruned$weights > 0) <= (net$weights > 0)))
  det <- net_pruning_details(pruned)
  expect_false(anyNA(det$weight))
  expect_identical(attr(det, "num_removed"), nrow(det))
})

test_that("an undirected network is pruned by pairs and stays symmetric", {
  # Co-occurrence counts: a-c (5) and b-d (2) are at most 10. Visited
  # column by column over the pairs, a-c goes first (in both directions):
  # b-d still joins {a, b} to {c, d}. b-d is then the only bridge and stays
  # in both directions. Cell by cell, d->b went and b->d stayed.
  W <- matrix(c(0,  90, 5,  0,
                90, 0,  0,  2,
                5,  0,  0,  80,
                0,  2,  80, 0), 4, byrow = TRUE,
              dimnames = list(letters[1:4], letters[1:4]))
  net <- build_network(W, method = "co_occurrence")
  expect_false(net$directed)
  pr <- net_prune(net, threshold = 10)
  expect_true(isSymmetric(unname(pr$weights)))
  expect_equal(pr$weights["a", "c"], 0)
  expect_equal(pr$weights["c", "a"], 0)
  expect_equal(pr$weights["b", "d"], 2)
  expect_equal(pr$weights["d", "b"], 2)
  det <- net_pruning_details(pr)
  expect_identical(nrow(det), 1L)          # the pair once
  expect_identical(attr(det, "num_removed"), 1L)
  expect_identical(attr(det, "num_retained"), 3L)   # a-b, b-d, c-d
})

test_that("signed networks: every zeroed edge is recorded and counted", {
  # 3 nodes, undirected: A-B 0.6, A-C -0.3, B-C 0.4. Disparity at
  # level 0.5 keeps A-B only: 2 pairs removed (one negative), 1 retained.
  W <- matrix(c(1, 0.6, -0.3,
                0.6, 1, 0.4,
                -0.3, 0.4, 1), 3, byrow = TRUE,
              dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  net <- build_network(W, method = "cor", params = list(n = 100))
  expect_false(net$directed)
  pd <- net_prune(net, method = "disparity", level = 0.5)
  det <- net_pruning_details(pd)
  kept <- sum(pd$weights[upper.tri(pd$weights)] != 0)
  expect_identical(attr(det, "num_removed"), nrow(det))
  expect_identical(attr(det, "num_retained"), kept)
  expect_identical(attr(det, "num_removed"), 2L)
  expect_identical(attr(det, "num_retained"), 1L)
  expect_identical(sum(det$weight < 0), 1L)
})

test_that("bootstrap pruning records non-significant negative edges", {
  cn <- build_network(srl_strategies, method = "cor")
  b <- bootstrap_network(cn, iter = 20L, seed = 1)
  pb <- net_prune(cn, method = "bootstrap", boot = b)
  det <- net_pruning_details(pb)
  up <- upper.tri(cn$weights)
  zeroed <- sum(cn$weights[up] != 0 & pb$weights[up] == 0)
  expect_identical(nrow(det), zeroed)
  expect_identical(sum(det$weight < 0),
                   sum(cn$weights[up] < 0 & pb$weights[up] == 0))
  expect_identical(attr(det, "num_retained"),
                   sum(pb$weights[up] != 0))
})

test_that("a group is pruned by its own bootstraps, matched by name", {
  g <- build_network(group_regulation_long, method = "relative",
                     action = "Action", actor = "Actor", time = "Time",
                     group = "Achiever")
  bg <- bootstrap_network(g, iter = 20L, seed = 1)
  pg <- net_prune(g, method = "bootstrap", boot = bg)
  invisible(lapply(names(g), function(nm) {
    expected <- bg[[nm]]$significant
    diag(expected) <- diag(g[[nm]]$weights)
    expect_equal(unname(pg[[nm]]$weights), unname(expected))
  }))
  expect_error(net_prune(g, method = "bootstrap", boot = bg[[1]]),
               class = "nestimate_bad_input")
  expect_error(net_prune(g[[2]], method = "bootstrap", boot = bg[[1]]),
               class = "nestimate_bad_input")
})

test_that("a depruned network can be pruned again with other settings", {
  net <- make_net()
  undone <- net_deprune(net_prune(net, threshold = 0.2))
  again <- net_prune(undone, threshold = 0.3)
  fresh <- net_prune(net, threshold = 0.3)
  expect_equal(again$weights, fresh$weights)
  expect_equal(net_pruning_details(again), net_pruning_details(fresh))
  expect_error(net_prune(net_prune(net, threshold = 0.2)),
               class = "nestimate_bad_input")
})
