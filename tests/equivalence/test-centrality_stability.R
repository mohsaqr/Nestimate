# centrality_stability(): the case-dropping stability of centralities (CS coefficient, Epskamp, Borsboom & Fried
# 2018) against tna::estimate_cs(). Both drop sequences from the same random stream (FINDINGS.md), so under one
# seed the correlation of every draw is the same number, and with it the CS coefficients and the summary per drop
# proportion (the mean, SD and share of draws at or above the threshold, written here from tna's correlations).
# Re-expressed from the centrality_stability parts of local_testing_and_equivalence/test-equiv-tna-functions.R.
# The comparison with bootnet::corStability() in test-equiv-centrality-stability.R accepted differences of 0.3 on
# 80 % of the configurations (different estimator, centralities and random stream): not a declarable comparison,
# it stays a standing study. The actor-level comparison needs tna >= 1.5.4 (estimate_cs(actor =)). The windowed
# comparison of test-equiv-tna-functions.R is not carried over: tna 1.3.1 ignores params = list(window_size = )
# for a relative model (its weights equal the lag-1 model's), so its reference was not a windowed network.

# summary.net_stability()'s table written from tna's matrices of correlations (draws x drop proportions).
rs_stability_summary <- function(estimate, threshold = 0.7) {
  drop_prop <- attr(estimate, "drop_prop")
  do.call(rbind, lapply(names(estimate), \(m) {
    correlations <- estimate[[m]]$correlations
    data.frame(measure = m, drop_prop = drop_prop,
               mean_cor = colMeans(correlations, na.rm = TRUE),
               sd_cor = apply(correlations, 2L, stats::sd, na.rm = TRUE),
               prop_above = colMeans(!is.na(correlations) & correlations >= threshold))
  }))
}

inputs <- list(rs_random_sequences(20L, 4L, 15L, 11L), rs_random_sequences(30L, 3L, 12L, 12L),
               rs_random_sequences(15L, 5L, 25L, 13L))

test_that("CS coefficients and per-drop correlations equal tna::estimate_cs() on three sequence sets", {
  skip_without_reference("tna")
  invisible(lapply(seq_along(inputs), \(i) {
    set.seed(1)
    reference <- tna::estimate_cs(tna::tna(inputs[[i]]), iter = 30L)
    found <- centrality_stability(build_network(inputs[[i]], method = "relative"), iter = 30L, seed = 1)
    # The CS coefficients are a field of the result (API gap, FINDINGS.md: `stability$cs`).
    expect_reference(found$cs, vapply(reference, \(m) m$cs_coefficient, numeric(1)), reference = "tna",
                     covers = "centrality_stability()", quantity = sprintf("random set %d, CS coefficients", i),
                     tolerance = 1e-10)
    expect_reference(summary(found), rs_stability_summary(reference), reference = "tna",
                     covers = c("centrality_stability()", "summary.net_stability()"),
                     quantity = sprintf("random set %d, per drop proportion", i), tolerance = 1e-10,
                     key = c("measure", "drop_prop"))
  }))
})

test_that("a grouped network gives tna::estimate_cs() of each group under the same seed", {
  skip_without_reference("tna")
  sequences <- rs_random_sequences(40L, 4L, 12L, 14L)
  group <- rep(c("x", "y"), length.out = nrow(sequences))
  # As for the bootstrap, each group is drawn from the same seed.
  expected <- do.call(rbind, lapply(c("x", "y"), \(g) {
    set.seed(1)
    cbind(group = g, rs_stability_summary(tna::estimate_cs(tna::tna(sequences[group == g, ]), iter = 10L)))
  }))
  found <- centrality_stability(build_network(cbind(sequences, grp = group), method = "relative", group = "grp"),
                                iter = 10L, seed = 1)
  expect_reference(summary(found), expected, reference = "tna",
                   covers = c("centrality_stability()", "summary.net_stability_group()"),
                   quantity = "two groups, per drop proportion", tolerance = 1e-10,
                   key = c("group", "measure", "drop_prop"))
})
