# build_mcml(), cluster_summary(), summary.mcml(): multi-cluster multi-level networks of a partition of the states.
# The macro network is the transition network of the sequences recoded to their clusters; each within-cluster
# network is the conditional chain of the cluster's states (consecutive pairs with both states in the cluster).
# References: tna::tna() on the recoded and on the NA-masked sequences, and the transition counts written from the
# definition. Re-expressed from local_testing_and_equivalence/test-equiv-mcml.R, test-mcml-within-equivalence.R,
# test-equiv-mcml-repartition.R and test-mcml-equiv-extracted.R. The layers are read from the result's fields
# `macro` and `clusters` (API gap: no accessor returns one layer).

set.seed(1)
states <- c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "c3")
clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2"), C = c("c1", "c2", "c3"))
wide <- as.data.frame(matrix(sample(states, 80L * 25L, replace = TRUE), 80L, 25L), stringsAsFactors = FALSE)
names(wide) <- paste0("T", seq_len(25L))
# Sequences of different lengths: trailing positions missing.
ragged <- wide
lengths_kept <- sample(5:25, 80L, replace = TRUE)
ragged[col(as.matrix(ragged)) > lengths_kept] <- NA_character_

test_that("the macro network and its initial probabilities equal tna::tna() on the sequences recoded to clusters", {
  skip_without_reference("tna")
  invisible(lapply(list(full = wide, ragged = ragged), \(data) {
    found <- build_mcml(data, clusters = clusters, type = "tna")
    reference <- tna::tna(sq_recode(data, clusters))
    label <- if (anyNA(data)) "ragged" else "full"
    expect_reference(found$macro$weights, reference$weights, reference = "tna", covers = "build_mcml()",
                     quantity = sprintf("macro transition probability, %s", label), tolerance = 1e-12)
    expect_reference(found$macro$inits, reference$inits, reference = "tna", covers = "build_mcml()",
                     quantity = sprintf("macro initial probability, %s", label), tolerance = 1e-12)
  }))
})

test_that("each within-cluster network equals tna::tna() on the sequences masked to that cluster", {
  skip_without_reference("tna")
  found <- build_mcml(wide, clusters = clusters, type = "tna")
  invisible(lapply(names(clusters), \(cluster) {
    reference <- tna::tna(sq_mask(wide, clusters[[cluster]]))
    expect_reference(found$clusters[[cluster]]$weights, reference$weights[clusters[[cluster]], clusters[[cluster]]],
                     reference = "tna", covers = "build_mcml()",
                     quantity = sprintf("within-cluster transition probability, cluster %s", cluster), tolerance = 1e-12)
  }))
})

test_that("within-cluster networks equal the conditional chain counted from the definition", {
  found <- build_mcml(ragged, clusters = clusters, type = "tna")
  invisible(lapply(names(clusters), \(cluster) {
    expect_reference(found$clusters[[cluster]]$weights, sq_transitions(ragged, clusters[[cluster]]),
                     reference = "independent: consecutive pairs within the cluster, counted and row-normalised",
                     covers = "build_mcml()", quantity = sprintf("within-cluster transition probability, ragged, cluster %s", cluster),
                     tolerance = 1e-12)
  }))
})

test_that("the cluster summary table follows from tna's macro network", {
  skip_without_reference("tna")
  macro <- tna::tna(sq_recode(wide, clusters))$weights
  expected <- data.frame(cluster = names(clusters), size = lengths(clusters),
                         between_out = rowSums(macro) - diag(macro), between_in = colSums(macro) - diag(macro))
  expect_reference(summary(build_mcml(wide, clusters = clusters, type = "tna")), expected, reference = "tna",
                   covers = c("build_mcml()", "summary.mcml()"), quantity = "between-cluster out and in weight",
                   tolerance = 1e-12, key = "cluster")
})

test_that("cluster_summary(method = 'sum') of a frequency network equals block sums of tna::ftna() counts", {
  skip_without_reference("tna")
  counts <- tna::ftna(wide)$weights
  found <- cluster_summary(build_network(wide, method = "frequency"), clusters = clusters, method = "sum")
  blocks <- outer(names(clusters), names(clusters), Vectorize(\(a, b) sum(counts[clusters[[a]], clusters[[b]]])))
  dimnames(blocks) <- list(names(clusters), names(clusters))
  expect_reference(found$macro$weights, blocks, reference = "tna", covers = "cluster_summary()",
                   quantity = "summed transitions between clusters", tolerance = 0,
                   note = "block sums written in the test on tna::ftna()")
  expect_reference(found$macro$inits, colSums(blocks) / sum(blocks), reference = "tna", covers = "cluster_summary()",
                   quantity = "share of transitions into each cluster", tolerance = 1e-12)
  expect_reference(found$clusters$A$weights, counts[clusters$A, clusters$A], reference = "tna", covers = "cluster_summary()",
                   quantity = "within-cluster counts, cluster A", tolerance = 0)
})
