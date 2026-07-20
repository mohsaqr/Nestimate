# Downstream contract: htna (CRAN) Imports Nestimate and reaches in via
# `Nestimate::fn()` with no importFrom() and no version pin, so a rename here
# produces no install-time error in htna -- it fails in a user's session.
#
# These tests encode htna's requirements as literal data. They deliberately do
# NOT require htna to be installed, so they guard the contract everywhere,
# including on CRAN. The paired end-to-end check that actually runs htna's
# suite lives in local_testing_and_equivalence/revdep-htna.R.

# ---- Contract A: symbols htna resolves via Nestimate:: ---------------------

test_that("functions htna calls via Nestimate:: are exported", {
  needed <- c(
    "build_network", "bootstrap_network", "sequence_plot", "sequence_compare",
    "state_distribution", "state_frequencies", "plot_state_frequencies",
    "mosaic_plot", "frequencies", "centrality_stability",
    "casedrop_reliability", "network_reliability", "permutation",
    "markov_order_test", "association_rules"
  )
  exported <- getNamespaceExports(asNamespace("Nestimate"))
  expect_true(all(needed %in% exported),
              info = paste("not exported:",
                           paste(setdiff(needed, exported), collapse = ", ")))
})

test_that("datasets htna calls via Nestimate:: are available", {
  needed <- c("ai_long", "human_long")
  lazy <- ls(getNamespaceInfo("Nestimate", "lazydata"))
  expect_true(all(needed %in% lazy),
              info = paste("not lazy-exported:",
                           paste(setdiff(needed, lazy), collapse = ", ")))
})

test_that("argument names htna passes still exist on the generics", {
  # Extracted by parsing htna's R/ for named arguments at Nestimate:: call sites,
  # including those passed through do.call(Nestimate::fn, list(...)).
  # Functions htna calls only positionally are absent here on purpose: a formal
  # rename cannot break a positional call. Their behaviour is covered by the
  # revdep gate (local_testing_and_equivalence/revdep-htna.R), not by this file.
  contract <- list(
    centrality_stability = c("centrality_fn", "certainty", "drop_prop", "iter",
                             "loops", "measures", "method", "seed", "threshold"),
    network_reliability  = c("iter", "scale", "seed", "split"),
    sequence_compare     = c("adjust", "group", "iter", "min_freq", "sub", "test"),
    # passed via do.call(Nestimate::sequence_plot, c(list(x = ., type = .), ...))
    sequence_plot        = c("x", "type")
  )
  invisible(Map(function(fn, args) {
    f <- getExportedValue("Nestimate", fn)
    # Generics delegate to methods; accept the arg on either the generic or,
    # when the generic is `...`-only, on its netobject method.
    have <- names(formals(f))
    if (!all(args %in% have)) {
      m <- tryCatch(getS3method(fn, "netobject", envir = asNamespace("Nestimate")),
                    error = function(e) NULL)
      if (!is.null(m)) have <- union(have, names(formals(m)))
    }
    expect_true(all(args %in% have),
                info = paste0(fn, "() lost: ",
                              paste(setdiff(args, have), collapse = ", ")))
  }, names(contract), contract))
})

# ---- Contract B: Nestimate dispatches on objects of class `htna` -----------

# Minimal object matching what the htna package's build_htna() returns: a
# netobject carrying $node_groups (node, group), with "htna" prepended to the
# class. Constructed here rather than imported, so this file needs no htna.
.fake_htna <- function() {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:4], 40, TRUE), T2 = sample(LETTERS[1:4], 40, TRUE),
    T3 = sample(LETTERS[1:4], 40, TRUE), T4 = sample(LETTERS[1:4], 40, TRUE)
  )
  net <- build_network(seqs, method = "relative")
  net$node_groups <- data.frame(
    node  = net$nodes$label,
    group = rep(c("AI", "Human"), length.out = nrow(net$nodes)),
    stringsAsFactors = FALSE
  )
  net$actor_levels <- c("AI", "Human")
  net$nodes$groups <- factor(
    net$node_groups$group,
    levels = net$actor_levels
  )
  class(net) <- c("htna", class(net))
  net
}

.expect_preserved_htna_group <- function(x, expected_levels) {
  expect_s3_class(x, "htna_group")
  expect_s3_class(x, "netobject_group")
  expect_identical(attr(x, "actor_levels"), expected_levels)
  for (net in x) {
    expect_s3_class(net, "htna")
    expect_identical(net$actor_levels, expected_levels)
    expect_true(all(c("node", "group") %in% names(net$node_groups)))
    expect_identical(as.character(net$node_groups$node),
                     as.character(net$nodes$label))
    expect_false(anyNA(net$node_groups$group))
  }
}

test_that("the htna S3 methods Nestimate registers stay dispatchable", {
  # Read the NAMESPACE registration table, not getS3method(): the latter falls
  # back to any same-named function in the namespace, and devtools::load_all()
  # registers by name convention regardless of NAMESPACE. Both would let a
  # dropped `S3method(gen, htna)` slip through.
  tbl <- getNamespaceInfo(asNamespace("Nestimate"), "S3methods")
  registered <- tbl[tbl[, 2] == "htna", 1]
  for (gen in c("mosaic_plot", "plot_state_frequencies", "state_distribution")) {
    expect_true(gen %in% registered,
                info = paste0("Nestimate no longer registers ", gen, ".htna"))
  }
})

test_that("state_distribution.htna returns a tidy per-group frame", {
  sd_tbl <- state_distribution(.fake_htna())
  expect_s3_class(sd_tbl, "data.frame")
  expect_true(all(c("group", "state", "count", "proportion") %in% names(sd_tbl)))
  expect_setequal(unique(sd_tbl$group), c("AI", "Human"))
  expect_false(any(is.na(sd_tbl$proportion)))
})

test_that("htna objects still carry the netobject fields Nestimate reads", {
  net <- .fake_htna()
  expect_true(all(c("weights", "nodes", "edges", "data", "node_groups") %in% names(net)))
  expect_true(all(c("node", "group") %in% names(net$node_groups)))
  expect_identical(class(net)[1], "htna")
  expect_true(inherits(net, "netobject"))
})

test_that("distance clustering preserves the HTNA actor partition", {
  net <- .fake_htna()
  clustering <- build_clusters(net, k = 2, seed = 1)

  expect_identical(clustering$htna_partition$actor_levels,
                   net$actor_levels)

  grouped <- build_network(clustering)
  .expect_preserved_htna_group(grouped, net$actor_levels)
  expect_s3_class(attr(grouped, "clustering"), "net_clustering")
})

test_that("cluster_network preserves HTNA children and diagnostics", {
  net <- .fake_htna()
  grouped <- cluster_network(net, k = 2, seed = 1)

  .expect_preserved_htna_group(grouped, net$actor_levels)
  expect_s3_class(attr(grouped, "clustering"), "net_clustering")
})

test_that("MMM clustering preserves the HTNA actor partition", {
  net <- .fake_htna()
  grouped <- cluster_mmm(
    net,
    k = 2,
    n_starts = 1,
    max_iter = 20,
    seed = 1
  )

  .expect_preserved_htna_group(grouped, net$actor_levels)
  clustering <- attr(grouped, "clustering")
  expect_s3_class(clustering, "net_mmm_clustering")
  expect_identical(clustering$htna_partition$actor_levels,
                   net$actor_levels)
  expect_equal(nrow(clustering$posterior), nrow(net$data))
})

# ---- Contract B2: the three exports htna's own suite never exercises -------
# htna re-exports these as build-time aliases, e.g.
#   association_rules_htna <- Nestimate::association_rules
# so removing one breaks htna at *install* time (loud), but a behaviour or
# return-shape change passes silently -- and htna's suite never calls them, so
# the revdep gate cannot see it either. These are the only guards they have.

test_that("frequencies() keeps the shape htna's alias re-exports", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:4], 40, TRUE), T2 = sample(LETTERS[1:4], 40, TRUE),
    T3 = sample(LETTERS[1:4], 40, TRUE), T4 = sample(LETTERS[1:4], 40, TRUE)
  )
  f <- frequencies(seqs)
  # Returns a state x state transition-count matrix, not a data.frame.
  expect_true(inherits(f, "nest_transition_counts"))
  expect_true(is.matrix(f))
  expect_identical(dim(f), c(4L, 4L))
})

test_that("state_frequencies() keeps the shape htna's alias re-exports", {
  sf <- state_frequencies(.fake_htna())
  expect_s3_class(sf, "data.frame")
  expect_true(all(c("state", "count", "proportion") %in% names(sf)))
  expect_identical(nrow(sf), 4L)
})

test_that("association_rules() keeps the shape htna's alias re-exports", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:4], 40, TRUE), T2 = sample(LETTERS[1:4], 40, TRUE),
    T3 = sample(LETTERS[1:4], 40, TRUE), T4 = sample(LETTERS[1:4], 40, TRUE)
  )
  ar <- association_rules(seqs, min_support = 0.05, min_confidence = 0.3)
  expect_s3_class(ar, "net_association_rules")
  expect_true(all(c("antecedent", "consequent", "support", "confidence",
                    "lift", "conviction", "count") %in% names(ar$rules)))
})

# ---- Contract C: return-shape stability of centrality_stability ------------
# htna's equivalence suite compares its cs values against Nestimate's. That
# only holds if the *set of measures returned* is a deterministic function of
# (data, requested measures). Commit 0de9168 changed this set and broke htna.

test_that("centrality_stability returns exactly the requested measures when all vary", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:5], 60, TRUE), T2 = sample(LETTERS[1:5], 60, TRUE),
    T3 = sample(LETTERS[1:5], 60, TRUE), T4 = sample(LETTERS[1:5], 60, TRUE)
  )
  net <- build_network(seqs, method = "relative")
  m <- c("InStrength", "OutStrength", "Betweenness")
  st <- centrality_stability(net, measures = m, iter = 5L, seed = 1L)
  expect_identical(names(st$cs), m)
})

test_that("the returned measure set is deterministic across repeated calls", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:5], 60, TRUE), T2 = sample(LETTERS[1:5], 60, TRUE),
    T3 = sample(LETTERS[1:5], 60, TRUE), T4 = sample(LETTERS[1:5], 60, TRUE)
  )
  net <- build_network(seqs, method = "relative")
  m <- c("InStrength", "OutStrength", "Betweenness")
  a <- centrality_stability(net, measures = m, iter = 5L, seed = 1L)
  b <- centrality_stability(net, measures = m, iter = 5L, seed = 1L)
  expect_identical(names(a$cs), names(b$cs))
  expect_equal(a$cs, b$cs)
})

test_that("partially degenerate: zero-variance measures drop, survivors keep order", {
  # Complete digraph over 3 states: every pair is directly connected, so every
  # Betweenness is 0 (sd == 0) while the strengths vary. Exercises the
  # `measures <- measures[keep]` drop branch of centrality_stability().
  set.seed(9)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C"), 30 * 4, TRUE,
                                      prob = c(.6, .3, .1)), nrow = 30))
  names(seqs) <- paste0("T", seq_len(ncol(seqs)))
  net <- build_network(seqs, method = "relative")
  m <- c("InStrength", "OutStrength", "Betweenness")
  st <- suppressWarnings(centrality_stability(net, measures = m, iter = 5L, seed = 1L))
  # Betweenness is degenerate here and must drop; the rest survive, in order.
  expect_identical(names(st$cs), c("InStrength", "OutStrength"))
})

test_that("undefined (NaN) measures do not crash and do not change the default set", {
  # This is the exact condition commit 0de9168 addressed: on a pure 3-cycle
  # Diffusion is NaN, so sd() is NA. The pre-0de9168 filter `sd(x) > 0` put an
  # NA into `keep`, and `if (!any(keep))` then evaluated `if (NA)` -> error.
  #
  # Pinning the DEFAULT measure set (which includes Diffusion) is what makes
  # this test sensitive: htna's equivalence suite compares against Nestimate's
  # returned measure set, so any change to it is a downstream break.
  seqs <- as.data.frame(matrix(rep(c("A", "B", "C"), times = 20 * 2),
                               nrow = 20, byrow = TRUE))
  names(seqs) <- paste0("T", seq_len(ncol(seqs)))
  net <- build_network(seqs, method = "relative")

  # Request Diffusion explicitly: it is not in the default set, and the test is
  # vacuous without it.
  m <- c("InStrength", "Betweenness", "Diffusion")
  st <- suppressWarnings(centrality_stability(net, measures = m, iter = 5L, seed = 1L))
  expect_identical(names(st$cs), m)
  expect_equal(unname(st$cs), rep(0, length(m)))
})

test_that("centrality_stability's DEFAULT measure set is the 0.6.0 trio", {
  # htna's CRAN release calls Nestimate::centrality_stability() with no
  # `measures` and compares the result against its own explicit trio. Any
  # change to this default silently breaks that reverse dependency: commit
  # b0d20cd swapped OutStrength -> Diffusion here and broke 21 htna assertions
  # while Nestimate's own suite stayed green. This is the guard for that.
  expect_identical(eval(formals(centrality_stability)$measures),
                   c("InStrength", "OutStrength", "Betweenness"))
})

test_that("fully degenerate: every requested measure is returned with cs = 0", {
  # Pure 3-cycle A->B->C->A: all three measures have zero variance. This hits
  # the `!any(keep)` branch, which -- unlike the partial branch above -- keeps
  # every requested name. The two branches disagree on shape by design; pin it.
  seqs <- as.data.frame(matrix(rep(c("A", "B", "C"), times = 20 * 2),
                               nrow = 20, byrow = TRUE))
  names(seqs) <- paste0("T", seq_len(ncol(seqs)))
  net <- build_network(seqs, method = "relative")
  m <- c("InStrength", "OutStrength", "Betweenness")
  expect_warning(centrality_stability(net, measures = m, iter = 5L, seed = 1L),
                 "zero variance")
  st <- suppressWarnings(centrality_stability(net, measures = m, iter = 5L, seed = 1L))
  expect_identical(names(st$cs), m)
  expect_equal(unname(st$cs), rep(0, 3))
})

test_that("the restored HTNA partition matches the shape htna itself builds", {
  src <- .fake_htna()
  partition <- .capture_htna_partition(src)
  plain <- build_network(src$data, method = "relative")
  restored <- .restore_htna_partition(plain, partition)

  # htna's build_htna() stores $node_groups$group as character and mirrors it
  # as a factor on $nodes$groups. A factor in both places is not that shape.
  expect_type(restored$node_groups$group, "character")
  expect_s3_class(restored$nodes$groups, "factor")
  expect_identical(levels(restored$nodes$groups), src$actor_levels)
  expect_identical(attr(restored$node_groups, "actor_levels"),
                   src$actor_levels)
  expect_identical(restored$node_groups$group, src$node_groups$group)
  expect_false(anyNA(restored$node_groups$group))

  # Capturing the restored object again must reproduce the same partition.
  expect_identical(.capture_htna_partition(restored), partition)
})

test_that("actor levels that omit an observed actor error instead of NA-ing", {
  src <- .fake_htna()
  # A stale level set: factor() would silently coerce "Human" to NA.
  src$actor_levels <- "AI"
  partition <- .capture_htna_partition(src)
  plain <- build_network(src$data, method = "relative")

  expect_error(
    .restore_htna_partition(plain, partition),
    "absent from 'actor_levels'"
  )
})

test_that("restoring a group keeps the clustering attributes intact", {
  src <- .fake_htna()
  nets <- cluster_network(src, k = 2)

  expect_s3_class(nets, "htna_group")
  expect_false(is.null(attr(nets, "clustering")))
  expect_identical(attr(nets, "group_col"), "cluster")
  expect_identical(attr(nets, "actor_levels"), src$actor_levels)
})
