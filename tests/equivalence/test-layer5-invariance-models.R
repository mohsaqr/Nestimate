# Layer 5, invariance under transformation for the models fitted on sequences: build_mcml(), build_mmm(),
# sequence_patterns(), markov_order_test() and build_hypa() are refitted after the rows are permuted, the states
# renamed and every sequence duplicated, and every output must equal the original carried through the same
# transformation. Counts are integers and are compared exactly. A mixture refitted from a fixed seed reaches the
# same optimum up to a relabelling of its clusters, which is undone before comparing; it stops when the
# log-likelihood changes by less than `tol` (1e-6), so its parameters agree to 1e-8, not bit for bit.

wide <- li_two_regimes(n_per = 30L, length = 12L, seed = 961L)
set.seed(962L)
row_order <- sample(nrow(wide))
undo_order <- order(row_order)
state_map <- c(A = "zeta", B = "alpha", C = "mu", D = "beta")
back <- li_invert(state_map)
renamed <- li_rename_states(wide, state_map)
doubled <- rbind(wide, wide)

test_that("build_mcml() is invariant to row order, renamed states and renamed or reordered clusters", {
  clusters <- list(g1 = c("A", "B"), g2 = c("C", "D"))
  original <- build_mcml(wide, clusters = clusters)
  reference <- li_ref("build_mcml() on the original input")
  check <- function(fit, quantity, map = NULL, cluster_map = NULL) {
    expect_reference(li_relabel(fit$macro$weights, cluster_map %||% character()), original$macro$weights,
                     reference = reference, covers = "build_mcml()", quantity = paste(quantity, "macro weights"),
                     tolerance = 0, authority = "internal")
    Map(\(name) {
      layer <- fit$clusters[[if (is.null(cluster_map)) name else li_invert(cluster_map)[[name]]]]
      expect_reference(li_relabel(layer$weights, map %||% character()), original$clusters[[name]]$weights,
                       reference = reference, covers = "build_mcml()",
                       quantity = sprintf("%s cluster %s weights", quantity, name), tolerance = 0, authority = "internal")
    }, names(clusters))
  }
  check(build_mcml(wide[row_order, ], clusters = clusters), "rows permuted:")
  check(build_mcml(renamed, clusters = list(g1 = c("zeta", "alpha"), g2 = c("mu", "beta"))), "states renamed:", map = back)
  check(build_mcml(wide, clusters = list(second = c("C", "D"), first = c("A", "B"))), "clusters renamed and reordered:",
        cluster_map = c(second = "g2", first = "g1"))
  check(build_mcml(doubled, clusters = clusters), "every sequence twice:")
  counts <- build_mcml(wide, clusters = clusters, type = "frequency")
  expect_reference(build_mcml(doubled, clusters = clusters, type = "frequency")$macro$weights, 2 * counts$macro$weights,
                   reference = li_ref("twice build_mcml() on the sequences once"), covers = "build_mcml()",
                   quantity = "frequency macro weights, every sequence twice", tolerance = 0, authority = "internal")
})

test_that("build_mmm() reaches the same mixture after the rows are permuted, the states renamed or the data duplicated", {
  original <- build_mmm(wide, k = 2L, seed = 3L, n_starts = 10L, smooth = 0)
  note <- "EM stops when the log-likelihood changes by less than 1e-6, so a refit agrees to 1e-8."
  compare_fit <- function(refit, assignments, quantity, map = NULL, scale = 1) {
    matched <- li_match_clusters(original$assignments, assignments)
    reference <- li_ref("build_mmm() on the original input, clusters matched")
    expect_reference(assignments, matched[original$assignments],
                     reference = reference, covers = "build_mmm()", quantity = paste(quantity, "assignments"),
                     tolerance = 0, authority = "internal")
    expect_reference(refit$log_likelihood / scale, original$log_likelihood, reference = reference, covers = "build_mmm()",
                     quantity = paste(quantity, "log-likelihood"), tolerance = 1e-8, measure = "relative",
                     authority = "internal", note = note)
    expect_reference(refit$mixing[matched], original$mixing, reference = reference, covers = "build_mmm()",
                     quantity = paste(quantity, "mixing proportions"), tolerance = 1e-8, authority = "internal", note = note)
    Map(\(cluster) {
      model <- refit$models[[matched[[cluster]]]]
      expect_reference(li_relabel(model$weights, map %||% character()), original$models[[cluster]]$weights,
                       reference = reference, covers = "build_mmm()",
                       quantity = sprintf("%s cluster %d transition probabilities", quantity, cluster),
                       tolerance = 1e-8, authority = "internal", note = note)
      expect_reference(li_relabel(model$initial, map %||% character()), original$models[[cluster]]$initial,
                       reference = reference, covers = "build_mmm()",
                       quantity = sprintf("%s cluster %d initial probabilities", quantity, cluster),
                       tolerance = 1e-8, authority = "internal", note = note)
    }, seq_len(original$k))
  }
  permuted <- build_mmm(wide[row_order, ], k = 2L, seed = 3L, n_starts = 10L, smooth = 0)
  # each refit's assignments are taken in the original row order
  compare_fit(permuted, permuted$assignments[undo_order], "rows permuted:")
  renamed_fit <- build_mmm(renamed, k = 2L, seed = 3L, n_starts = 10L, smooth = 0)
  compare_fit(renamed_fit, renamed_fit$assignments, "states renamed:", map = back)
  doubled_fit <- build_mmm(doubled, k = 2L, seed = 3L, n_starts = 10L, smooth = 0)
  compare_fit(doubled_fit, doubled_fit$assignments[seq_len(nrow(wide))], "every sequence twice:", scale = 2)
})

test_that("sequence_patterns() is invariant to row order, follows renamed states and scales with duplication", {
  original <- sequence_patterns(wide, len = 2:3)
  reference <- li_ref("sequence_patterns() on the original input")
  columns <- c("pattern", "length", "occurrences", "proportion", "n_sequences", "support", "lift")
  expect_reference(subset(sequence_patterns(wide[row_order, ], len = 2:3), select = columns), subset(original, select = columns),
                   reference = reference, covers = "sequence_patterns()", quantity = "patterns, rows permuted",
                   key = "pattern", tolerance = 0, authority = "internal")
  renamed_patterns <- transform(sequence_patterns(renamed, len = 2:3), pattern = li_map_patterns(pattern, back))
  expect_reference(subset(renamed_patterns, select = columns), subset(original, select = columns), reference = reference,
                   covers = "sequence_patterns()", quantity = "patterns, states renamed", key = "pattern", tolerance = 0,
                   authority = "internal")
  # min_freq counts occurrences, which duplication doubles: with min_freq = 1 both fits keep the same patterns
  once <- sequence_patterns(wide, len = 2:3, min_freq = 1)
  twice <- sequence_patterns(doubled, len = 2:3, min_freq = 1)
  expect_reference(subset(twice, select = c("pattern", "proportion", "support", "lift")),
                   subset(once, select = c("pattern", "proportion", "support", "lift")), reference = reference,
                   covers = "sequence_patterns()", quantity = "shares, every sequence twice", key = "pattern", tolerance = 0,
                   authority = "internal")
  expect_reference(subset(twice, select = c("pattern", "occurrences", "n_sequences")),
                   transform(subset(once, select = c("pattern", "occurrences", "n_sequences")),
                             occurrences = 2L * occurrences, n_sequences = 2L * n_sequences),
                   reference = li_ref("twice sequence_patterns() on the sequences once"), covers = "sequence_patterns()",
                   quantity = "counts, every sequence twice", key = "pattern", tolerance = 0, authority = "internal")
})
