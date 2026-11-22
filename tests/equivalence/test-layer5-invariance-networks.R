# Layer 5, invariance under transformation (metamorphic tests): a network is refitted after a transformation of
# its input whose effect on every output is known in advance, and the refit must equal the original carried
# through the same transformation. Rows permuted, actors' blocks of a long log shuffled, states, actors and groups
# renamed, factor levels reversed: nothing changes. Every sequence duplicated: relative weights and initial
# probabilities unchanged, counts doubled. Association data with every variable multiplied by a positive constant:
# correlations, partial correlations and the graphical lasso unchanged.
# Transition counts are integers, so those checks are exact (tolerance 0). A correlation of permuted or rescaled
# rows sums the same terms in another order; those checks allow 1e-12.

wide <- tn_simulate(n_actors = 50L, n_states = 5L, seq_length = 10L, seed = 951L, ragged = TRUE)
set.seed(952L)
row_order <- sample(nrow(wide))
state_map <- c(A = "zeta", B = "alpha", C = "mu", D = "beta", E = "omega")
renamed <- li_rename_states(wide, state_map)
back <- li_invert(state_map)
doubled <- rbind(wide, wide)
sequence_methods <- c("relative", "frequency", "co_occurrence", "attention", "ngram", "gap", "reverse")
count_methods <- setdiff(sequence_methods, "relative")

test_that("sequence networks are invariant to the order of the rows", {
  lapply(sequence_methods, \(method) {
    li_expect_same_network(build_network(wide[row_order, ], method = method), build_network(wide, method = method),
                           route = "fit on the rows in their original order", covers = "build_network()",
                           quantity = sprintf("%s, rows permuted", method))
  })
  li_expect_same_network(build_network(wide[row_order, ], method = "relative", start = TRUE, end = TRUE, threshold = 0.1),
                         build_network(wide, method = "relative", start = TRUE, end = TRUE, threshold = 0.1),
                         route = "fit on the rows in their original order", covers = "build_network()",
                         quantity = "relative with start, end and threshold 0.1, rows permuted")
})

test_that("sequence networks are invariant to renaming the states", {
  lapply(sequence_methods, \(method) {
    original <- build_network(wide, method = method)
    refit <- build_network(renamed, method = method)
    reference <- li_ref("fit on the original state names")
    expect_reference(li_relabel(refit$weights, back), original$weights, reference = reference, covers = "build_network()",
                     quantity = sprintf("%s weights, states renamed", method), tolerance = 0, authority = "internal")
    expect_reference(li_relabel(li_initial(refit), back), li_initial(original), reference = reference,
                     covers = "build_network()", quantity = sprintf("%s initial probabilities, states renamed", method),
                     tolerance = 0, authority = "internal")
  })
})

test_that("relative networks are unchanged and count networks doubled when every sequence is duplicated", {
  li_expect_same_network(build_network(doubled, method = "relative"), build_network(wide, method = "relative"),
                         route = "fit on the sequences once", covers = "build_network()",
                         quantity = "relative, every sequence twice", parts = c("weights", "edges", "initial"))
  lapply(count_methods, \(method) {
    original <- build_network(wide, method = method)
    expect_reference(build_network(doubled, method = method)$weights, 2 * original$weights,
                     reference = li_ref("twice the fit on the sequences once"), covers = "build_network()",
                     quantity = sprintf("%s weights, every sequence twice", method), tolerance = 0, authority = "internal")
  })
})

test_that("a long event log is invariant to shuffling the actors' blocks and renaming the actors", {
  long <- transform(tn_long(wide), time = time + 1000 * actor)
  blocks <- split(long, long$actor)
  set.seed(953L)
  shuffled <- do.call(rbind, blocks[sample(length(blocks))])
  renamed_actors <- transform(long, actor = sprintf("id_%03d", 1000L - actor))
  fit <- function(data, method) build_network(data, method = method, actor = "actor", action = "action", order = "time")
  lapply(c("relative", "frequency", "attention"), \(method) {
    original <- fit(long, method)
    li_expect_same_network(fit(shuffled, method), original, route = "fit on the log in its original order",
                           covers = "build_network()", quantity = sprintf("%s, actors' blocks shuffled", method))
    li_expect_same_network(fit(renamed_actors, method), original, route = "fit with the original actor names",
                           covers = "build_network()", quantity = sprintf("%s, actors renamed", method))
  })
})

test_that("a factor coding of the states with reversed levels gives the same network", {
  reversed <- wide
  reversed[] <- lapply(wide, \(column) factor(column, levels = rev(LETTERS[1:5])))
  li_expect_same_network(build_network(reversed, method = "relative"), build_network(wide, method = "relative"),
                         route = "fit on character states", covers = "build_network()",
                         quantity = "relative, factor states with reversed levels")
})
