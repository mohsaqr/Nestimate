# Layer 4, recovery of known parameters: mixtures of first-order Markov chains. Sequences are drawn from three
# chains with known transition matrices, initial distributions and mixing proportions; every sequence's component is
# known. build_mmm() must recover the parameters (after relabelling its components by the permutation that brings
# them closest to the generating ones) and classify the sequences about as well as the Bayes rule under the true
# parameters; compare_mmm()'s BIC must pick the generating number of components.
# Simulators and the criteria's arithmetic: helper-layer-recovery.R.

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
# Bias: the mean of each matched estimate is within MC_Z Monte Carlo SEs of its generating value (Morris, White
# and Crowther 2019). The mixture MLE is consistent but not unbiased; with 300 sequences of 20 states and
# components this far apart its bias is expected to be well inside the Monte Carlo error of R = 200 replicates.
MC_Z <- 3.5
# Classification: the generating memberships are known, so the yardstick is the Bayes (oracle) classification
# under the generating parameters, the best any estimate can do on the same data. The fitted model's mean accuracy
# may fall short of the oracle's by at most ACCURACY_GAP and its mean adjusted Rand index by at most ARI_GAP.
ACCURACY_GAP <- 0.02
ARI_GAP <- 0.05
# Selection: among k = 2, 3, 4, the lowest BIC is at the generating k = 3 in at least SELECT_MIN of the replicates.
SELECT_MIN <- 0.9
R_RECOVERY <- 200L
R_SELECTION <- 100L
N_STARTS <- 10L

# ---- Generating process -----------------------------------------------------------------------------------------
states <- c("A", "B", "C")
named <- \(values) matrix(values, 3L, byrow = TRUE, dimnames = list(states, states))
transitions <- list(named(c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.3, 0.2, 0.5)),
                    named(c(0.1, 0.8, 0.1, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1)),
                    named(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8)))
initials <- list(c(A = 0.6, B = 0.2, C = 0.2), c(A = 0.2, B = 0.6, C = 0.2), c(A = 0.2, B = 0.2, C = 0.6))
mixing <- c(0.5, 0.3, 0.2)
n_seq <- 300L
len <- 20L
reference <- "independent: generating parameters and memberships of the simulated mixture of Markov chains"

# One replicate: the data, its true components and the oracle's classification.
rc_mixture_draw <- function(seed) {
  rc_with_seed(seed, {
    component <- sample.int(3L, n_seq, replace = TRUE, prob = mixing)
    wide <- rc_markov_sequences(n_seq, len, transitions, initials, component)
    list(wide = wide, component = component, oracle = rc_oracle_classes(wide, transitions, initials, mixing))
  })
}

test_that("build_mmm() recovers the generating transition matrices, initial distributions and mixing proportions", {
  # Serial: build_mmm() runs its random starts in parallel itself.
  replicates <- rc_replicate(seq_len(R_RECOVERY), \(r) {
    draw <- rc_mixture_draw(40000L + r)
    fit <- build_mmm(draw$wide, k = 3L, n_starts = N_STARTS, seed = r)
    # models, mixing, assignments: fields of the result (API gap: no accessor returns the component parameters).
    estimated <- lapply(fit$models, \(m) m$weights[states, states])
    match_to <- rc_match_components(estimated, transitions)
    relabel <- match(fit$assignments, match_to)
    list(transitions = estimated[match_to],
         initials = lapply(fit$models[match_to], \(m) as.numeric(m$initial[states])),
         mixing = as.numeric(fit$mixing[match_to]),
         accuracy = mean(relabel == draw$component), oracle_accuracy = mean(draw$oracle == draw$component),
         ari = rc_adjusted_rand(fit$assignments, draw$component),
         oracle_ari = rc_adjusted_rand(draw$oracle, draw$component))
  }, cores = 1L)
  invisible(lapply(seq_along(transitions), \(j) {
    bias <- rc_standardised_bias(lapply(replicates, \(e) e$transitions[[j]]), transitions[[j]])
    expect_reference(rc_cells(z = bias$z), rc_cells(z = 0 * transitions[[j]]), reference = reference,
                     covers = "build_mmm()",
                     quantity = sprintf("component %d transitions: (mean - generating) / MC SE, R = %d", j, R_RECOVERY),
                     tolerance = MC_Z, key = c("from", "to"),
                     note = sprintf("largest |mean - P| %.4f", max(abs(bias$mean - transitions[[j]]))))
    first <- rc_standardised_bias(lapply(replicates, \(e) e$initials[[j]]), unname(initials[[j]]))
    expect_reference(stats::setNames(first$z, states), stats::setNames(c(0, 0, 0), states), reference = reference,
                     covers = "build_mmm()",
                     quantity = sprintf("component %d initial: (mean - generating) / MC SE, R = %d", j, R_RECOVERY),
                     tolerance = MC_Z, note = sprintf("largest |mean - initial| %.4f", max(abs(first$mean - initials[[j]]))))
  }))
  shares <- rc_standardised_bias(lapply(replicates, \(e) e$mixing), mixing)
  expect_reference(stats::setNames(shares$z, paste("component", 1:3)), stats::setNames(c(0, 0, 0), paste("component", 1:3)),
                   reference = reference, covers = "build_mmm()",
                   quantity = sprintf("mixing proportions: (mean - generating) / MC SE, R = %d", R_RECOVERY),
                   tolerance = MC_Z, note = sprintf("mean estimates %s", toString(round(shares$mean, 4))))
  mean_of <- \(field) mean(vapply(replicates, \(e) e[[field]], numeric(1L)))
  expect_reference(mean_of("accuracy"), mean_of("oracle_accuracy"),
                   reference = "independent: Bayes classification under the generating mixture",
                   covers = "build_mmm()", quantity = sprintf("mean classification accuracy, R = %d", R_RECOVERY),
                   tolerance = ACCURACY_GAP)
  expect_reference(mean_of("ari"), mean_of("oracle_ari"),
                   reference = "independent: Bayes classification under the generating mixture",
                   covers = "build_mmm()", quantity = sprintf("mean adjusted Rand index, R = %d", R_RECOVERY),
                   tolerance = ARI_GAP)
})

test_that("compare_mmm() BIC selects the generating number of components", {
  chosen <- rc_replicate(seq_len(R_SELECTION), \(r) {
    draw <- rc_mixture_draw(50000L + r)
    table <- compare_mmm(draw$wide, k = 2:4, n_starts = N_STARTS, seed = r)
    subset(summary(table), best == "BIC")$k
  }, cores = 1L)
  expect_reference(mean(unlist(chosen) == 3L), 1, reference = reference, covers = c("compare_mmm()", "summary.mmm_compare()"),
                   quantity = sprintf("share of replicates whose lowest BIC is at k = 3 (of 2, 3, 4), R = %d", R_SELECTION),
                   tolerance = 1 - SELECT_MIN,
                   note = sprintf("chosen k: %s", toString(sprintf("%s: %d", names(table(unlist(chosen))),
                                                                  as.vector(table(unlist(chosen)))))))
})
