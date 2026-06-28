# Tests for as_htna(): full node-level network grouped by cluster, rebuilt from
# the original data (the node-level counterpart of build_mcml()).

make_seqs <- function() {
  set.seed(1)
  as.data.frame(
    matrix(sample(LETTERS[1:6], 360, TRUE), nrow = 30),
    stringsAsFactors = FALSE
  )
}
clusters3 <- list(C1 = c("A", "B"), C2 = c("C", "D"), C3 = c("E", "F"))

test_that("as_htna() returns a node-level netobject grouped by cluster", {
  seqs <- make_seqs()
  net <- as_htna(seqs, clusters3)
  expect_s3_class(net, "netobject")
  expect_s3_class(net, "cograph_network")
  # one node per state, with a cluster column
  expect_true("cluster" %in% names(net$nodes))
  expect_setequal(net$nodes$label, LETTERS[1:6])
  expect_identical(attr(net, "cluster_members"), clusters3)
})

test_that("as_htna() is faithful to a direct build (keeps cross-cluster edges)", {
  seqs <- make_seqs()
  net <- as_htna(seqs, clusters3)
  direct <- build_network(seqs, method = "relative")
  expect_equal(net$weights, direct$weights)
  # a genuine between-cluster transition is present and non-zero
  expect_gt(net$weights["A", "C"], 0)
})

test_that("node_groups and cluster column agree with the membership", {
  seqs <- make_seqs()
  net <- as_htna(seqs, clusters3)
  ng <- net$node_groups
  expect_true(is.data.frame(ng))
  expect_setequal(names(ng), c("node", "group"))
  # A,B -> C1 ; C,D -> C2 ; E,F -> C3
  m <- stats::setNames(net$nodes$cluster, net$nodes$label)
  expect_equal(unname(m[c("A", "B")]), c("C1", "C1"))
  expect_equal(unname(m[c("E", "F")]), c("C3", "C3"))
})

test_that("as_htna(mcml) works on its own (stashed source) and with data", {
  # a clearly sequence-detected frame (so build_mcml stashes its source)
  seqs <- data.frame(t1 = c("A", "C", "E", "B"), t2 = c("B", "D", "F", "A"),
                     t3 = c("C", "A", "E", "D"), t4 = c("D", "B", "F", "C"),
                     stringsAsFactors = FALSE)
  m <- build_mcml(seqs, clusters3)
  expect_false(is.null(attr(m, "htna_source")))   # source stashed
  net0 <- as_htna(m)                               # no data needed
  net  <- as_htna(m, data = seqs)
  expect_identical(attr(net0, "cluster_members"), clusters3)
  expect_equal(net0$weights, build_network(seqs, method = "relative")$weights)
  expect_equal(net0$weights, net$weights)
})

test_that("as_htna(mcml) still errors for a matrix-built mcml (no node data)", {
  w <- matrix(c(0, 2, 1, 3, 0, 2, 1, 1, 0), 3, 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  m <- build_mcml(w, list(C1 = c("A", "B"), C2 = "C"))
  expect_error(as_htna(m), "no expandable node-level source")
})

test_that("as_htna() requires clusters when none can be derived", {
  seqs <- make_seqs()
  expect_error(as_htna(seqs), "clusters' is required")
})

test_that("as_htna() rejects a clustering that misses nodes", {
  seqs <- make_seqs()
  expect_error(
    as_htna(seqs, list(C1 = c("A", "B"), C2 = c("C", "D"))),  # E, F unassigned
    regexp = "Unmapped|not assigned|partition"
  )
})

test_that("as_htna() accepts a per-node membership vector", {
  seqs <- make_seqs()
  membership <- stats::setNames(
    c("C1", "C1", "C2", "C2", "C3", "C3"), LETTERS[1:6]
  )
  net <- as_htna(seqs, membership)
  expect_s3_class(net, "netobject")
  expect_true("cluster" %in% names(net$nodes))
})

test_that("as_htna() result is usable for inference", {
  seqs <- make_seqs()
  net <- as_htna(seqs, clusters3)
  expect_no_error(bootstrap_network(net, iter = 20))
})

test_that("as_htna(mcml) robustness across build inputs", {
  cl <- list(C1 = c("A", "B"), C2 = c("C", "D"), C3 = c("E", "F"))

  # wide sequence input -> stashed -> works on its own and matches data path
  seqs <- data.frame(t1 = c("A", "C", "E", "B"), t2 = c("B", "D", "F", "A"),
                     t3 = c("C", "A", "E", "D"), stringsAsFactors = FALSE)
  m_seq <- build_mcml(seqs, cl)
  expect_equal(as_htna(m_seq)$weights, as_htna(m_seq, data = seqs)$weights)

  # long input (widened internally) -> stashed -> works on its own
  long <- data.frame(sess = rep(1:2, each = 4), t = rep(1:4, 2),
                     act = c("A", "B", "C", "D", "E", "F", "A", "B"),
                     stringsAsFactors = FALSE)
  m_long <- build_mcml(long, cl, actor = "sess", action = "act", time = "t")
  expect_s3_class(as_htna(m_long), "netobject")

  # build_mcml(mcml) is idempotent and preserves the stash
  expect_s3_class(as_htna(build_mcml(m_seq)), "netobject")

  # explicit data overrides the stash
  expect_s3_class(as_htna(m_seq, data = seqs), "netobject")

  # edge list and matrix inputs carry no expandable source -> clear error
  el <- data.frame(from = c("A", "B", "C", "E"), to = c("B", "C", "D", "F"),
                   weight = c(2, 1, 3, 1), stringsAsFactors = FALSE)
  expect_error(as_htna(build_mcml(el, cl)), "no expandable node-level source")
  w <- matrix(c(0, 2, 1, 3, 0, 2, 1, 1, 0), 3, 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_error(as_htna(build_mcml(w, list(G1 = c("A", "B"), G2 = "C"))),
               "no expandable node-level source")
})
