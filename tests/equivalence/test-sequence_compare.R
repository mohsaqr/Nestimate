# sequence_compare(test = "chisq"): for each pattern, the 2 x G table of its count against the rest of each group's
# patterns, its standardised (adjusted) residuals, chi-square statistic and p-value. Reference: stats::chisq.test()
# on the same 2 x G tables (stdres, statistic, p.value). Moved from tests/testthat/test-houserules-D.R (1.2.0).

test_that("pattern residuals, statistics and p-values equal stats::chisq.test() on each 2 x G table", {
  halves <- rep(c("x", "y"), length.out = nrow(trajectories))
  res <- summary(sequence_compare(trajectories, group = halves, sub = 2, min_freq = 1, test = "chisq",
                                  adjust = "none"))
  counts <- cbind(res$freq_x, res$freq_y)
  totals <- colSums(counts)
  reference <- do.call(rbind, lapply(seq_len(nrow(counts)), \(i) {
    fit <- suppressWarnings(stats::chisq.test(rbind(counts[i, ], totals - counts[i, ])))
    data.frame(pattern = res$pattern[i], resid_x = fit$stdres[1L, 1L], resid_y = fit$stdres[1L, 2L],
               statistic = unname(fit$statistic), p_value = fit$p.value)
  }))
  expect_reference(res[c("pattern", "resid_x", "resid_y", "statistic", "p_value")], reference, reference = "stats",
                   covers = c("sequence_compare()", "summary.net_sequence_comparison()"),
                   quantity = "chi-square test of a pattern, trajectories halves", tolerance = 1e-12, key = "pattern")
})
