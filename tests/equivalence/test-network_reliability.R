# network_reliability(): split-half reliability of a transition network against tna::reliability(). Both split
# the sequences in halves from the same random stream, so under one seed the four metrics Nestimate reports (mean,
# median and maximum absolute difference of the two halves' weights, and their Pearson correlation) have the same
# mean and SD over the draws. Re-expressed from the reliability parts of
# local_testing_and_equivalence/test-equiv-tna-functions.R (tolerance 1e-8, kept). Its windowed comparison is not
# carried over: tna 1.3.1 ignores params = list(window_size = ) for a relative model, so the reference was not
# windowed.

# tna's summary rows of the four metrics, under Nestimate's metric names.
rs_reliability_summary <- function(reliability, model) {
  names_tna <- c(mean_dev = "Mean Abs. Diff.", median_dev = "Median Abs. Diff.", cor = "Pearson",
                 max_dev = "Max Abs. Diff.")
  table <- as.data.frame(reliability$summary)
  rows <- match(names_tna, table$metric)
  data.frame(model = model, metric = names(names_tna), mean = table$mean[rows], sd = table$sd[rows])
}

inputs <- list(rs_random_sequences(20L, 4L, 15L, 21L), rs_random_sequences(30L, 3L, 12L, 22L),
               rs_random_sequences(25L, 5L, 20L, 23L))

test_that("split-half metrics equal tna::reliability() on three random sequence sets", {
  skip_without_reference("tna")
  invisible(lapply(seq_along(inputs), \(i) {
    set.seed(1)
    reference <- tna::reliability(tna::tna(inputs[[i]]), iter = 30L)
    found <- network_reliability(build_network(inputs[[i]], method = "relative"), iter = 30L, seed = 1)
    expect_reference(summary(found), rs_reliability_summary(reference, "relative"), reference = "tna",
                     covers = c("network_reliability()", "summary.net_reliability()"),
                     quantity = sprintf("random set %d", i), tolerance = 1e-8, key = c("model", "metric"))
  }))
})

