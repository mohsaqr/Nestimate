# Part 2 of test-layer5-invariance-models.R, split so that one Prova case holds at most 7 references;
# the shared set-up is repeated.

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

test_that("markov_order_test() is invariant to row order and renamed states, and its likelihoods double with the data", {
  original <- markov_order_test(wide, max_order = 2L, n_perm = 50L, seed = 1L)$test_table
  deterministic <- c("order", "loglik", "AIC", "BIC", "df", "g2", "p_asymptotic")
  note <- "A log-likelihood summed over the windows in another order; 1e-10 allows for it."
  permuted <- markov_order_test(wide[row_order, ], max_order = 2L, n_perm = 50L, seed = 1L)$test_table
  expect_reference(subset(permuted, select = deterministic), subset(original, select = deterministic),
                   reference = li_ref("markov_order_test() on the original rows"), covers = "markov_order_test()",
                   quantity = "deterministic columns, rows permuted", key = "order", tolerance = 1e-10, measure = "relative",
                   authority = "internal", note = note)
  renamed_test <- markov_order_test(renamed, max_order = 2L, n_perm = 50L, seed = 1L)$test_table
  expect_reference(renamed_test, original, reference = li_ref("markov_order_test() on the original state names"),
                   covers = "markov_order_test()", quantity = "every column, states renamed", key = "order",
                   tolerance = 1e-10, measure = "relative", authority = "internal", note = note)
  twice <- markov_order_test(doubled, max_order = 2L, n_perm = 50L, seed = 1L)$test_table
  expect_reference(subset(twice, select = c("order", "loglik", "g2", "df")),
                   transform(subset(original, select = c("order", "loglik", "g2", "df")), loglik = 2 * loglik, g2 = 2 * g2),
                   reference = li_ref("twice the log-likelihoods of the sequences once"), covers = "markov_order_test()",
                   quantity = "log-likelihood, G2 and df, every sequence twice", key = "order", tolerance = 1e-10,
                   measure = "relative", authority = "internal", note = note)
})

test_that("build_hypa() is invariant to row order, follows renamed states and doubles its counts with the data", {
  original <- build_hypa(wide, order = 2L)$scores
  columns <- c("path", "observed", "expected", "ratio", "p_value", "p_under", "p_over", "anomaly")
  expect_reference(subset(build_hypa(wide[row_order, ], order = 2L)$scores, select = columns), subset(original, select = columns),
                   reference = li_ref("build_hypa() on the original rows"), covers = "build_hypa()",
                   quantity = "path scores, rows permuted", key = "path", tolerance = 0, authority = "internal")
  renamed_scores <- transform(build_hypa(renamed, order = 2L)$scores, path = li_map_patterns(path, back, sep = " -> "))
  expect_reference(subset(renamed_scores, select = columns), subset(original, select = columns),
                   reference = li_ref("build_hypa() on the original state names"), covers = "build_hypa()",
                   quantity = "path scores, states renamed", key = "path", tolerance = 1e-12, measure = "relative",
                   authority = "internal", note = "Renaming reorders the nodes of the null model's sums; 1e-12 allows for it.")
  twice <- build_hypa(doubled, order = 2L)$scores
  expect_reference(subset(twice, select = c("path", "observed", "expected", "ratio")),
                   transform(subset(original, select = c("path", "observed", "expected", "ratio")),
                             observed = 2L * observed, expected = 2 * expected),
                   reference = li_ref("twice build_hypa() on the sequences once"), covers = "build_hypa()",
                   quantity = "observed and expected counts, every sequence twice", key = "path", tolerance = 1e-12,
                   measure = "relative", authority = "internal")
})
