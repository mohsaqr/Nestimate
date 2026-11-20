# sequence_patterns(), sequence_compare(), summary.net_sequence_comparison(), pattern_outcome(): contiguous
# patterns (n-grams) of sequences, their comparison between groups, and their association with an outcome.
# References: the n-grams counted from the definition in base R; stats::chisq.test() on each pattern's count
# against the other n-grams of the group; stats::glm() of the outcome on the presence of the selected patterns.
# New: local_testing_and_equivalence/test-equiv-sequence-patterns.R compared with codyna, which is abandoned and is
# not carried over.

regulation <- tna::group_regulation[seq_len(200L), ]

test_that("n-gram occurrences, sequences, support and share equal a count of the n-grams", {
  invisible(lapply(2:3, \(len) {
    expected <- sq_pattern_counts(sq_ngrams(regulation, len), nrow(regulation))
    expect_reference(sequence_patterns(regulation, type = "ngram", len = len, min_freq = 1, min_support = 0), expected,
                     reference = "independent: contiguous n-grams counted in base R", covers = "sequence_patterns()",
                     quantity = sprintf("%d-grams", len), tolerance = 1e-12, key = "pattern")
  }))
})

test_that("group comparison of bigrams equals stats::chisq.test() of each pattern against the other bigrams", {
  data <- tna::group_regulation[seq_len(300L), ]
  groups <- rep(c("A", "B", "C"), 100L)
  found <- sequence_compare(build_network(cbind(data, grp = groups), method = "relative", group = "grp"), sub = 2,
                            test = "chisq", adjust = "none", min_freq = 5L)
  grams <- lapply(c("A", "B", "C"), \(g) sq_ngrams(data[groups == g, ], 2L))
  totals <- vapply(grams, nrow, numeric(1))
  counts <- lapply(grams, \(x) table(x$pattern))
  patterns <- summary(found)$pattern
  rows <- lapply(patterns, \(p) {
    freq <- vapply(counts, \(x) if (p %in% names(x)) as.numeric(x[[p]]) else 0, numeric(1))
    test <- stats::chisq.test(rbind(freq, totals - freq), correct = FALSE)
    data.frame(pattern = p, freq_A = freq[[1L]], freq_B = freq[[2L]], freq_C = freq[[3L]],
               prop_A = freq[[1L]] / totals[[1L]], prop_B = freq[[2L]] / totals[[2L]], prop_C = freq[[3L]] / totals[[3L]],
               resid_A = test$stdres[1L, 1L], resid_B = test$stdres[1L, 2L], resid_C = test$stdres[1L, 3L],
               statistic = unname(test$statistic), p_value = test$p.value)
  })
  expect_reference(summary(found), do.call(rbind, rows), reference = "stats",
                   covers = c("sequence_compare()", "summary.net_sequence_comparison()"),
                   quantity = "bigram counts and chi-square per pattern, three groups", tolerance = 1e-10,
                   measure = "relative", key = "pattern",
                   note = "bigrams counted in base R; stats::chisq.test(correct = FALSE) of the pattern against the rest")
})

test_that("pattern effects equal one stats::glm() of the outcome on the selected patterns' presence", {
  grams <- sq_ngrams(regulation, 2L)
  set.seed(1)
  outcome <- ifelse(stats::runif(nrow(regulation)) < 0.5, "hi", "lo")
  found <- pattern_outcome(regulation, outcome = outcome, len = 2, n = 3)
  presence <- as.data.frame(lapply(stats::setNames(found$pattern, paste0("p", seq_along(found$pattern))),
                                   \(p) as.numeric(seq_len(nrow(regulation)) %in% grams$sequence[grams$pattern == p])))
  fit <- summary(stats::glm(factor(outcome) ~ ., data = cbind(outcome = outcome, presence), family = stats::binomial))
  table_glm <- fit$coefficients[-1L, , drop = FALSE]
  supports <- lapply(c("hi", "lo"), \(level) colMeans(presence[outcome == level, , drop = FALSE]))
  expected <- data.frame(pattern = found$pattern, estimate = table_glm[, "Estimate"], std_error = table_glm[, "Std. Error"],
                         p_value = table_glm[, "Pr(>|z|)"], odds_ratio = exp(table_glm[, "Estimate"]),
                         support_hi = supports[[1L]], support_lo = supports[[2L]])
  # The same glm on both sides; relative 1e-8 for IRLS run to the same convergence.
  expect_reference(found, expected, reference = "stats", covers = "pattern_outcome()",
                   quantity = "logistic regression on three bigrams", tolerance = 1e-8, measure = "relative", key = "pattern",
                   note = "the event is the second outcome level, as glm() takes a factor")
})
