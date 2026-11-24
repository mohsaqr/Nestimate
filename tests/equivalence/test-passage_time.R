# passage_time(), markov_stability() and their summaries: mean first passage times, return times and the
# stationary distribution of an ergodic Markov chain (Kemeny and Snell 1960), against markovchain 1.2's
# meanFirstPassageTime() and steadyStates(). The diagonal of passage_time() holds the mean return time 1 / pi
# (Kac's lemma), where markovchain puts 0: the reference matrix carries 1 / pi there.
# Re-expressed from local_testing_and_equivalence/test-equiv-markov.R (100 random chains by their largest
# difference), test-equiv-moved3-markov.R and test-markov-equiv-extracted.R.
# The unrounded matrix and stationary distribution are read from the result's fields (API gap: summary() rounds,
# and no accessor returns the matrix).

inputs <- list(ho_stochastic(3L, 11L), ho_stochastic(5L, 12L), ho_stochastic(7L, 13L))

# markovchain's mean first passage times with the mean return times 1 / pi on the diagonal.
ho_passage <- function(matrix_p) {
  passage <- markovchain::meanFirstPassageTime(ho_chain(matrix_p))
  passage <- passage[rownames(matrix_p), colnames(matrix_p)]
  diag(passage) <- 1 / ho_stationary(matrix_p)[rownames(matrix_p)]
  passage
}

test_that("mean first passage times and the stationary distribution equal markovchain", {
  skip_without_reference("markovchain")
  invisible(lapply(seq_along(inputs), \(i) {
    found <- passage_time(inputs[[i]])
    expect_reference(found$matrix, ho_passage(inputs[[i]]), reference = "markovchain", covers = "passage_time()",
                     quantity = sprintf("mean first passage time, chain %d", i), tolerance = 1e-8, measure = "relative",
                     note = "markovchain::meanFirstPassageTime(), diagonal 1 / steadyStates()")
    expect_reference(found$stationary, ho_stationary(inputs[[i]]), reference = "markovchain", covers = "passage_time()",
                     quantity = sprintf("stationary distribution, chain %d", i), tolerance = 1e-10)
  }))
})

test_that("the rounded summaries of passage_time() and markov_stability() equal markovchain's values", {
  skip_without_reference("markovchain")
  invisible(lapply(seq_along(inputs), \(i) {
    matrix_p <- inputs[[i]]
    passage <- ho_passage(matrix_p)
    off <- passage
    diag(off) <- NA
    stationary <- ho_stationary(matrix_p)[rownames(matrix_p)]
    expected <- data.frame(state = rownames(matrix_p), persistence = diag(matrix_p), stationary_prob = stationary,
                           return_time = 1 / stationary, sojourn_time = 1 / (1 - diag(matrix_p)),
                           avg_time_to_others = rowMeans(off, na.rm = TRUE),
                           avg_time_from_others = colMeans(off, na.rm = TRUE))
    printed <- utils::capture.output(stability <- summary(markov_stability(matrix_p)))
    # markov_stability() rounds persistence and the stationary probability to 4 decimals and the times to 2:
    # the tolerance is half a unit in the last kept decimal.
    expect_reference(subset(stability, select = c(state, persistence, stationary_prob)),
                     subset(expected, select = c(state, persistence, stationary_prob)), reference = "markovchain",
                     covers = c("markov_stability()", "summary.net_markov_stability()"),
                     quantity = sprintf("persistence and stationary probability, chain %d", i), tolerance = 5e-5,
                     key = "state")
    expect_reference(subset(stability, select = -c(persistence, stationary_prob)),
                     subset(expected, select = -c(persistence, stationary_prob)), reference = "markovchain",
                     covers = c("markov_stability()", "summary.net_markov_stability()"),
                     quantity = sprintf("return, sojourn and passage times, chain %d", i), tolerance = 5e-3,
                     key = "state")
    table <- summary(passage_time(matrix_p))
    expect_reference(subset(table, select = c(state, stationary)),
                     data.frame(state = rownames(matrix_p), stationary = stationary), reference = "markovchain",
                     covers = c("passage_time()", "summary.net_mpt()"),
                     quantity = sprintf("summary stationary, chain %d", i), tolerance = 5e-5, key = "state")
    expect_reference(subset(table, select = -stationary),
                     data.frame(state = rownames(matrix_p), return_time = 1 / stationary,
                                mean_out = rowMeans(off, na.rm = TRUE), mean_in = colMeans(off, na.rm = TRUE)),
                     reference = "markovchain", covers = c("passage_time()", "summary.net_mpt()"),
                     quantity = sprintf("summary times, chain %d", i), tolerance = 5e-3, key = "state")
  }))
})
