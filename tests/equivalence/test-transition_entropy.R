# transition_entropy(): per-state (row) entropies of a transition matrix, the entropy rate h = sum pi_i H(P_i.),
# the entropy of the stationary distribution and the redundancy H(pi) - h (Cover and Thomas 2006, ch. 4).
# References: entropy::entropy.plugin() for each row and the stationary distribution from
# markovchain::steadyStates(). Re-expressed from local_testing_and_equivalence/test-equiv-transition_entropy.R
# (100 random chains by their largest difference). The result's fields are read with `$` where summary() has no
# unrounded table for them (API gap).

inputs <- list(list(p = ho_stochastic(3L, 21L), base = 2), list(p = ho_stochastic(5L, 22L), base = exp(1)),
               list(p = ho_stochastic(7L, 23L), base = 10))
units <- c("2" = "log2", "2.71828182845905" = "log", "10" = "log10")

ho_entropy <- function(matrix_p, base) {
  unit <- units[[as.character(base)]]
  rows <- vapply(rownames(matrix_p), \(state) entropy::entropy.plugin(matrix_p[state, ], unit = unit), numeric(1))
  stationary <- ho_stationary(matrix_p)[rownames(matrix_p)]
  rate <- sum(stationary * rows)
  stationary_entropy <- entropy::entropy.plugin(stationary, unit = unit)
  list(rows = rows, stationary = stationary,
       chain = c(entropy_rate = rate, stationary_entropy = stationary_entropy, redundancy = stationary_entropy - rate))
}

test_that("row entropies, entropy rate and stationary entropy equal entropy::entropy.plugin() with markovchain", {
  skip_without_reference("entropy")
  skip_without_reference("markovchain")
  invisible(lapply(inputs, \(input) {
    found <- transition_entropy(input$p, base = input$base)
    expected <- ho_entropy(input$p, input$base)
    label <- sprintf("%d states, base %.4g", nrow(input$p), input$base)
    expect_reference(found$row_entropy, expected$rows, reference = "entropy", covers = "transition_entropy()",
                     quantity = sprintf("row entropy, %s", label), tolerance = 1e-10)
    expect_reference(found$stationary, expected$stationary, reference = "markovchain", covers = "transition_entropy()",
                     quantity = sprintf("stationary distribution, %s", label), tolerance = 1e-8)
    expect_reference(c(entropy_rate = found$entropy_rate, stationary_entropy = found$stationary_entropy,
                       redundancy = found$redundancy), expected$chain, reference = "entropy",
                     covers = "transition_entropy()", quantity = sprintf("chain entropies, %s", label), tolerance = 1e-10,
                     note = "stationary distribution from markovchain::steadyStates()")
    # summary() returns a list holding a per-state table (API gap: read with `$`).
    expect_reference(summary(found)$table,
                     data.frame(state = rownames(input$p), stationary = expected$stationary, row_entropy = expected$rows,
                                contribution = expected$stationary * expected$rows),
                     reference = "entropy", covers = c("transition_entropy()", "summary.net_transition_entropy()"),
                     quantity = sprintf("summary per state, %s", label), tolerance = 1e-8, key = "state")
  }))
})
