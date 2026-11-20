# bootstrap_network(): the stability bootstrap of a transition network against tna::bootstrap(). Both resample
# whole sequences with replacement from the same random stream (FINDINGS.md: "Bootstrap and case-dropping
# stability do share tna's random stream"), so under one seed every bootstrap mean, SD, percentile interval,
# consistency range, p-value and significance flag is the same number.
# Re-expressed from local_testing_and_equivalence/test-equiv-bootstrap.R (20 random configs, one expect_equal per
# cell) and test-bootstrap_network-equiv-extracted.R (group_regulation, 1000 draws). The coverage and convergence
# checks of test-equiv-bootstrap.R compare Nestimate with itself and stay there as a standing study.

# tna's summary with the bootstrap means and SDs of its matrices added, one row per edge of its summary.
rs_tna_bootstrap <- function(boot) {
  moments <- rs_edge_table(mean = boot$weights_mean, sd = boot$weights_sd)
  merge(as.data.frame(boot$summary), moments, by = c("from", "to"))
}

test_that("bootstrap_network() equals tna::bootstrap() on group_regulation, 1000 draws", {
  skip_without_reference("tna")
  sequences <- tna::group_regulation
  set.seed(42)
  reference <- tna::bootstrap(tna::tna(sequences), iter = 1000L)
  found <- bootstrap_network(build_network(sequences, method = "relative"), iter = 1000L, seed = 42)
  # Same draws, same arithmetic: only the last bits of a mean or a quantile can differ.
  expect_reference(summary(found), rs_tna_bootstrap(reference), reference = "tna",
                   covers = c("bootstrap_network()", "summary.net_bootstrap()"),
                   quantity = "group_regulation, 1000 draws", tolerance = 1e-10, key = c("from", "to"),
                   note = "tna::bootstrap(tna::tna(group_regulation), iter = 1000) after set.seed(42)")
})

test_that("bootstrap_network() equals tna::bootstrap() on three random sequence sets", {
  skip_without_reference("tna")
  inputs <- list(rs_random_sequences(20L, 4L, 15L, 101L), rs_random_sequences(30L, 3L, 12L, 102L),
                 rs_random_sequences(25L, 5L, 20L, 103L))
  invisible(lapply(seq_along(inputs), \(i) {
    set.seed(1)
    reference <- tna::bootstrap(tna::tna(inputs[[i]]), iter = 200L)
    found <- bootstrap_network(build_network(inputs[[i]], method = "relative"), iter = 200L, seed = 1)
    expect_reference(summary(found), rs_tna_bootstrap(reference), reference = "tna",
                     covers = c("bootstrap_network()", "summary.net_bootstrap()"),
                     quantity = sprintf("random set %d, 200 draws", i), tolerance = 1e-10, key = c("from", "to"))
  }))
})

test_that("bootstrap_network() on a grouped network equals tna::bootstrap() of each group", {
  skip_without_reference("tna")
  sequences <- tna::group_regulation
  group <- rep(c("A", "B"), length.out = nrow(sequences))
  # Nestimate draws every group from the same seed; tna::bootstrap() of a group_tna() continues one random stream
  # from group to group. The reference is therefore tna::bootstrap() of each group's model after set.seed(3).
  expected <- do.call(rbind, lapply(c("A", "B"), \(g) {
    set.seed(3)
    cbind(group = g, rs_tna_bootstrap(tna::bootstrap(tna::tna(sequences[group == g, ]), iter = 100L)))
  }))
  found <- bootstrap_network(build_network(cbind(sequences, grp = group), method = "relative", group = "grp"),
                             iter = 100L, seed = 3)
  expect_reference(summary(found), expected, reference = "tna",
                   covers = c("bootstrap_network()", "summary.net_bootstrap_group()"),
                   quantity = "two groups of group_regulation, 100 draws", tolerance = 1e-10,
                   key = c("group", "from", "to"))
})
