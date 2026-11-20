# net_prune(), net_pruning_details(), net_deprune(), net_reprune(): pruning of a transition network against
# tna::prune() and tna::pruning_details() (threshold, lowest share, disparity filter). Nestimate never prunes a
# self-loop; tna prunes the diagonal like any other edge, so the off-diagonal edges are compared as a reference and
# the diagonal as a known difference. Re-expressed from local_testing_and_equivalence/test-equiv-prune.R (40 random
# configs, maximum of differences).

inputs <- list(rs_random_sequences(150L, 4L, 6L, 1L), rs_random_sequences(150L, 5L, 6L, 2L),
               rs_random_sequences(150L, 6L, 6L, 3L))
settings <- list(list(method = "threshold", threshold = 0.15), list(method = "threshold", threshold = 0.25),
                 list(method = "lowest", lowest = 0.20), list(method = "lowest", lowest = 0.40),
                 list(method = "disparity", level = 0.5))
setting_label <- function(setting) paste(names(setting), unlist(setting), sep = " = ", collapse = ", ")

test_that("pruned off-diagonal weights equal tna::prune() for three methods on three sequence sets", {
  skip_without_reference("tna")
  invisible(lapply(seq_along(inputs), \(i) lapply(settings, \(setting) {
    found <- do.call(net_prune, c(list(build_network(inputs[[i]], method = "relative")), setting))
    expected <- do.call(tna::prune, c(list(tna::tna(inputs[[i]])), setting))
    expect_reference(rs_edge_table(weight = extract_transition_matrix(found), drop_diagonal = TRUE),
                     rs_edge_table(weight = expected$weights, drop_diagonal = TRUE), reference = "tna",
                     covers = c("net_prune()", "net_prune.netobject()"),
                     quantity = sprintf("random set %d, %s, off-diagonal", i, setting_label(setting)),
                     tolerance = 1e-10, key = c("from", "to"))
  })))
})

test_that("self-loops are kept by net_prune() and pruned by tna::prune()", {
  skip_without_reference("tna")
  found <- net_prune(build_network(inputs[[1L]], method = "relative"), method = "threshold", threshold = 0.25)
  expected <- tna::prune(tna::tna(inputs[[1L]]), method = "threshold", threshold = 0.25)
  expect_known_difference(diag(extract_transition_matrix(found)), diag(expected$weights), reference = "tna",
                          covers = "net_prune()", quantity = "random set 1, threshold = 0.25, self-loops",
                          tolerance = 1e-10,
                          reason = "Nestimate never removes a self-loop when it prunes; tna::prune() prunes the diagonal like any other edge.")
})

test_that("removed off-diagonal edges equal tna::pruning_details()", {
  skip_without_reference("tna")
  invisible(lapply(settings[c(1L, 3L, 5L)], \(setting) {
    found <- do.call(net_prune, c(list(build_network(inputs[[2L]], method = "relative")), setting))
    removed <- as.data.frame(tna::pruning_details(do.call(tna::prune, c(list(tna::tna(inputs[[2L]])), setting))))
    expect_reference(net_pruning_details(found), subset(removed, from != to), reference = "tna",
                     covers = c("net_pruning_details()", "net_pruning_details.netobject()"),
                     quantity = sprintf("random set 2, %s, removed edges", setting_label(setting)),
                     tolerance = 1e-10, key = c("from", "to"))
  }))
})

test_that("deprune restores the weights of tna::deprune(); reprune those of tna::reprune()", {
  skip_without_reference("tna")
  pruned <- net_prune(build_network(inputs[[3L]], method = "relative"), method = "threshold", threshold = 0.2)
  pruned_tna <- tna::prune(tna::tna(inputs[[3L]]), method = "threshold", threshold = 0.2)
  restored <- net_deprune(pruned)
  expect_reference(extract_transition_matrix(restored), tna::deprune(pruned_tna)$weights, reference = "tna",
                   covers = c("net_deprune()", "net_deprune.netobject()"), quantity = "random set 3, depruned weights",
                   tolerance = 1e-10)
  expect_reference(rs_edge_table(weight = extract_transition_matrix(net_reprune(restored)), drop_diagonal = TRUE),
                   rs_edge_table(weight = tna::reprune(tna::deprune(pruned_tna))$weights, drop_diagonal = TRUE),
                   reference = "tna", covers = c("net_reprune()", "net_reprune.netobject()"),
                   quantity = "random set 3, repruned weights, off-diagonal", tolerance = 1e-10, key = c("from", "to"))
})

test_that("pruning a grouped network equals tna::prune() of each group", {
  skip_without_reference("tna")
  sequences <- inputs[[1L]]
  group <- rep(c("x", "y"), length.out = nrow(sequences))
  pruned <- net_prune(build_network(cbind(sequences, grp = group), method = "relative", group = "grp"),
                      method = "threshold", threshold = 0.25)
  details <- net_pruning_details(pruned)
  restored <- net_deprune(pruned)
  # A group of networks has no tidy accessor: each group is read by name (API gap).
  invisible(lapply(c("x", "y"), \(g) {
    expected <- tna::prune(tna::tna(sequences[group == g, ]), method = "threshold", threshold = 0.25)
    expect_reference(rs_edge_table(weight = extract_transition_matrix(pruned[[g]]), drop_diagonal = TRUE),
                     rs_edge_table(weight = expected$weights, drop_diagonal = TRUE), reference = "tna",
                     covers = "net_prune.netobject_group()", quantity = sprintf("group %s, off-diagonal", g),
                     tolerance = 1e-10, key = c("from", "to"))
    expect_reference(details[[g]], subset(as.data.frame(tna::pruning_details(expected)), from != to),
                     reference = "tna", covers = "net_pruning_details.netobject_group()",
                     quantity = sprintf("group %s, removed edges", g), tolerance = 1e-10, key = c("from", "to"))
    expect_reference(extract_transition_matrix(restored[[g]]), tna::deprune(expected)$weights, reference = "tna",
                     covers = "net_deprune.netobject_group()", quantity = sprintf("group %s, depruned weights", g),
                     tolerance = 1e-10)
    expect_reference(rs_edge_table(weight = extract_transition_matrix(net_reprune(restored)[[g]]), drop_diagonal = TRUE),
                     rs_edge_table(weight = tna::reprune(tna::deprune(expected))$weights, drop_diagonal = TRUE),
                     reference = "tna", covers = "net_reprune.netobject_group()",
                     quantity = sprintf("group %s, repruned weights, off-diagonal", g), tolerance = 1e-10,
                     key = c("from", "to"))
  }))
})
