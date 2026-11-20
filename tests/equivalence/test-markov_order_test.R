# markov_order_test(): likelihood-ratio (G2) tests of Markov order k against k - 1 (Anderson and Goodman 1957).
# References: markovchain::verifyMarkovProperty() (the order-2 test of one sequence, G2 method) and the G2 of
# conditional independence written from its definition: the sum over the k - 1 intermediate states of the G2 of the
# (first state x next state) table. The permutation p-values are drawn inside the function (seed only) and are not
# compared. Re-expressed from local_testing_and_equivalence/test-equiv-markov_order.R (30 + 20 random inputs by
# their largest difference).

# G2 and degrees of freedom of order k against k - 1, summed over the intermediate contexts.
ho_order_g2 <- function(sequences, k) {
  tuples <- do.call(rbind, lapply(sequences, \(trajectory) {
    positions <- seq_len(length(trajectory) - k)
    data.frame(first = trajectory[positions],
               context = vapply(positions, \(p) paste(trajectory[p + seq_len(k - 1L)], collapse = "|"), character(1)),
               last = trajectory[positions + k])
  }))
  parts <- lapply(split(tuples, tuples$context), \(rows) {
    counts <- table(rows$first, rows$last)
    if (nrow(counts) < 2L || ncol(counts) < 2L) return(c(g2 = 0, df = 0))
    expected <- outer(rowSums(counts), colSums(counts)) / sum(counts)
    observed <- counts > 0
    c(g2 = 2 * sum(counts[observed] * log(counts[observed] / expected[observed])),
      df = (nrow(counts) - 1) * (ncol(counts) - 1))
  })
  Reduce(`+`, parts)
}

test_that("the order-2 G2, degrees of freedom and asymptotic p-value equal markovchain::verifyMarkovProperty()", {
  skip_without_reference("markovchain")
  invisible(lapply(c(1L, 2L, 3L, 4L), \(seed) {
    sequence <- unlist(ho_sequences(1L, 600L, ho_stochastic(4L, seed), seed + 100L))
    found <- subset(summary(markov_order_test(list(sequence), max_order = 2L, n_perm = 1L, seed = seed)), order == 2L,
                    select = c(order, g2, df, p_asymptotic))
    reference <- markovchain::verifyMarkovProperty(sequence, method = "G", verbose = FALSE)
    expect_reference(found, data.frame(order = 2L, g2 = unname(reference$statistic), df = unname(reference$parameter),
                                       p_asymptotic = reference$p.value),
                     reference = "markovchain", covers = c("markov_order_test()", "summary.net_markov_order()"),
                     quantity = sprintf("order 2 test, seed %d", seed), tolerance = 1e-10, measure = "relative",
                     key = "order")
  }))
})

test_that("G2 and degrees of freedom of orders 1 to 3 equal the stratified conditional-independence G2", {
  invisible(lapply(c(5L, 6L), \(seed) {
    wide <- ho_sequences(20L, 80L, ho_stochastic(4L, seed), seed + 100L)
    sequences <- lapply(seq_len(nrow(wide)), \(i) unlist(wide[i, ], use.names = FALSE))
    statistics <- lapply(1:3, \(k) ho_order_g2(sequences, k))
    expected <- data.frame(order = 1:3, g2 = vapply(statistics, \(s) s[["g2"]], numeric(1)),
                           df = vapply(statistics, \(s) s[["df"]], numeric(1)))
    found <- subset(summary(markov_order_test(sequences, max_order = 3L, n_perm = 1L, seed = seed)), order >= 1L,
                    select = c(order, g2, df))
    expect_reference(found, expected, reference = "independent: stratified G2 of conditional independence, base R",
                     covers = c("markov_order_test()", "summary.net_markov_order()"),
                     quantity = sprintf("orders 1 to 3, seed %d", seed), tolerance = 1e-10, measure = "relative",
                     key = "order")
  }))
})
