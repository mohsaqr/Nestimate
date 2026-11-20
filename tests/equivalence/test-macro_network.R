# macro_network(): the cluster-level network of an mcml, with clusters collapsed or expanded, recounted from the
# sequences. The reference recodes each state to its cluster (or keeps it, for an expanded cluster) and
# row-normalises the table of consecutive pairs (helper-transitions.R): the definition the documentation states.

wide <- tn_simulate(n_actors = 40L, n_states = 6L, seq_length = 12L, seed = 501L, ragged = TRUE)
clusters <- list(g1 = c("A", "B"), g2 = c("C", "D"), g3 = c("E", "F"))

# Sequences with every state of the collapsed clusters replaced by its cluster's name.
tn_recode <- function(sequences, clusters, expand = character()) {
  collapsed <- clusters[setdiff(names(clusters), expand)]
  lookup <- stats::setNames(rep(names(collapsed), lengths(collapsed)), unlist(collapsed, use.names = FALSE))
  lapply(sequences, \(s) ifelse(s %in% names(lookup), lookup[s], s))
}

test_that("the collapsed and the partly expanded macro networks equal the recoded pair tables", {
  fit <- build_mcml(wide, clusters = clusters)
  sequences <- tn_sequences(wide)
  expect_reference(extract_transition_matrix(macro_network(fit)),
                   tn_row_normalise(tn_pair_counts(tn_recode(sequences, clusters))),
                   reference = "independent: table of consecutive pairs of the recoded sequences, base R",
                   covers = "macro_network()", quantity = "all clusters collapsed", tolerance = 1e-10)
  expect_reference(extract_transition_matrix(macro_network(fit, expand = "g2")),
                   tn_row_normalise(tn_pair_counts(tn_recode(sequences, clusters, expand = "g2"))),
                   reference = "independent: table of consecutive pairs of the recoded sequences, base R",
                   covers = "macro_network()", quantity = "cluster g2 expanded", tolerance = 1e-10)
  expect_reference(extract_transition_matrix(macro_network(fit, expand = "all", method = "frequency")),
                   tn_pair_counts(sequences),
                   reference = "independent: table of consecutive pairs, base R",
                   covers = "macro_network()", quantity = "every cluster expanded, frequency", tolerance = 0)
})
