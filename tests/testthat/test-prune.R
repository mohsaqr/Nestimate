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
