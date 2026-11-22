# Part 3 of test-layer5-invariance-networks.R, split so that one Prova case holds at most 7 references;
# the shared set-up is repeated.

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

test_that("grouped networks follow renamed groups, and the order of the groups does not depend on the rows", {
  grouped <- transform(wide, g = rep(c("x", "y", "z"), length.out = nrow(wide)))
  original <- build_network(grouped, method = "relative", group = "g")
  group_map <- c(x = "Zed", y = "Alf", z = "Mid")
  renamed_groups <- build_network(transform(grouped, g = unname(group_map[g])), method = "relative", group = "g")
  Map(\(level, new) {
    li_expect_same_network(renamed_groups[[new]], original[[level]], route = "fit with the original group names",
                           covers = "build_network()", quantity = sprintf("group %s renamed %s", level, new),
                           drop_params = c("cols", "alphabet"))
  }, names(group_map), unname(group_map))
  expect_reference(names(renamed_groups), unname(group_map[names(original)]), reference = li_ref("original group order, renamed"),
                   covers = "build_network()", quantity = "group order after renaming", tolerance = 0, authority = "internal")
  order_defect <- "build_network(group = ) orders the groups by first appearance in the rows (unique()), so permuting rows reorders them and factor levels are ignored"
  reversed_rows <- build_network(grouped[rev(seq_len(nrow(grouped))), ], method = "relative", group = "g")
  Map(\(level) {
    li_expect_same_network(reversed_rows[[level]], original[[level]], route = "fit on the rows in their original order",
                           covers = "build_network()", quantity = sprintf("group %s, rows reversed", level),
                           drop_params = c("cols", "alphabet"))
  }, c("x", "y", "z"))
  expect_reference(names(reversed_rows), names(original), reference = li_ref("fit on the rows in their original order"),
                   covers = "build_network()", quantity = "group order, rows reversed", tolerance = 0,
                   authority = "internal", defect = order_defect)
  levelled <- build_network(transform(grouped, g = factor(g, levels = c("z", "x", "y"))), method = "relative", group = "g")
  expect_reference(names(levelled), c("z", "x", "y"), reference = li_ref("the declared factor levels"),
                   covers = "build_network()", quantity = "group order of a factor with levels z, x, y", tolerance = 0,
                   authority = "internal", defect = order_defect)
  # downstream, the pair direction (and so the sign of every difference) follows that order
  expect_reference(unique(edge_differences(compare_networks(reversed_rows), pair = c("x", "y"))$network_a),
                   unique(edge_differences(compare_networks(original), pair = c("x", "y"))$network_a),
                   reference = li_ref("comparison of the fit on the rows in their original order"),
                   covers = c("compare_networks()", "build_network()"), quantity = "first network of the pair x, y, rows reversed",
                   tolerance = 0, authority = "internal", defect = order_defect)
})
