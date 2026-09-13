# macro_network() -- the macro layer at mixed resolution.

make_mcml <- function() {
  seqs <- data.frame(
    t1 = c("A", "C", "A", "B", "D", "A"),
    t2 = c("B", "D", "C", "A", "C", "B"),
    t3 = c("C", "A", "D", "C", "A", "D"),
    stringsAsFactors = FALSE
  )
  build_mcml(seqs, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))
}

test_that("with no expand the result is the macro layer", {
  mc <- make_mcml()
  n <- macro_network(mc)
  expect_s3_class(n, "netobject")
  expect_setequal(n$nodes$label, names(mc$cluster_members))
  expect_null(n$expanded)
})

test_that("expanding a cluster replaces its node with its member states", {
  mc <- make_mcml()
  n  <- macro_network(mc, expand = "G2")
  expect_setequal(n$nodes$label, c("G1", "C", "D"))
  expect_false("G2" %in% n$nodes$label)
  expect_equal(n$expanded, "G2")
})

test_that("expanding splits the collapsed column exactly, leaving others alone", {
  mc <- make_mcml()
  a <- macro_network(mc)$weights
  b <- macro_network(mc, expand = "G2")$weights
  # G1's outflow to G2 equals its outflow to C plus D
  expect_equal(unname(a["G1", "G2"]),
               unname(b["G1", "C"] + b["G1", "D"]))
  # and its within-cluster self-loop is untouched
  expect_equal(unname(a["G1", "G1"]), unname(b["G1", "G1"]))
})

test_that("expanded states stay grouped with their parent cluster", {
  mc <- make_mcml()
  g  <- macro_network(mc, expand = "G2")$node_groups
  expect_equal(g$group[g$node == "C"], "G2")
  expect_equal(g$group[g$node == "D"], "G2")
  expect_equal(g$group[g$node == "G1"], "G1")
})

test_that("every cluster can be expanded, giving the node-level alphabet", {
  mc <- make_mcml()
  n  <- macro_network(mc, expand = c("G1", "G2"))
  expect_setequal(n$nodes$label,
                  unlist(mc$cluster_members, use.names = FALSE))
})

test_that("an unknown cluster and a data-less mcml are refused", {
  mc <- make_mcml()
  expect_error(macro_network(mc, expand = "nope"), "Unknown cluster")

  m <- matrix(c(0, 2, 1, 0), 2L, dimnames = list(c("A", "B"), c("A", "B")))
  flat <- build_mcml(m, clusters = list(G1 = "A", G2 = "B"))
  expect_error(macro_network(flat, expand = "G1"),
               class = "nestimate_no_expand_source")
})

test_that("as_tna(expand =) changes only the macro layer", {
  mc <- make_mcml()
  a <- as_tna(mc)
  b <- as_tna(mc, expand = "G2")
  expect_equal(names(a), names(b))
  expect_setequal(b$macro$nodes$label, c("G1", "C", "D"))
  expect_setequal(a$macro$nodes$label, c("G1", "G2"))
  # the per-cluster layers are untouched
  expect_identical(lapply(a[-1L], function(z) z$weights),
                   lapply(b[-1L], function(z) z$weights))
})

test_that("as_tna(expand =) matches macro_network() for the same clusters", {
  mc <- make_mcml()
  expect_equal(as_tna(mc, expand = "G2")$macro$weights,
               macro_network(mc, expand = "G2", method = "relative")$weights)
})

test_that('expand = "all" and TRUE put every state in the macro', {
  seqs <- data.frame(
    t1 = c("A", "C", "A", "B"), t2 = c("B", "D", "C", "A"),
    t3 = c("C", "A", "D", "C"), stringsAsFactors = FALSE
  )
  mc <- build_mcml(seqs, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))

  all_states <- sort(unlist(mc$cluster_members, use.names = FALSE))
  expect_equal(sort(macro_network(mc, expand = "all")$nodes$label), all_states)
  expect_equal(macro_network(mc, expand = TRUE)$weights,
               macro_network(mc, expand = "all")$weights)
  # naming every cluster by hand must give the same network
  expect_equal(macro_network(mc, expand = "all")$weights,
               macro_network(mc, expand = c("G1", "G2"))$weights)
  # FALSE is not a shorthand for "none" -- it is refused, not silently ignored
  expect_error(macro_network(mc, expand = FALSE), "expand")
})
