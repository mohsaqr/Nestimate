# Layer 4, recovery of known parameters: tests of two groups of sequences. Both groups are drawn from one
# first-order chain (the null is true), or the second from a chain that differs in one row (it is false).
# permutation() and sequence_compare() permute whole sequences (or their per-sequence counts) between the groups, so
# under the null every edge's, pattern's and global p-value must reject at its nominal level; under the
# alternative the global test must reject with high power.
# Simulators and the criteria's arithmetic: helper-layer-recovery.R.

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
# Rejection rates under the null lie within the central 99.9% binomial range of ALPHA for R independent replicates,
# for each edge, each pattern and each global statistic separately (the unadjusted p-values; no multiplicity
# correction is involved). The harness judges |rate - ALPHA| against the larger half-width of that range.
BAND_LEVEL <- 0.999
ALPHA <- 0.05
# 199 permutations: p = (b + 1) / 200, so a valid p-value has P(p <= 0.05) = 0.05 exactly; a test "rejects" when
# p <= ALPHA.
ITER <- 199L
# Power: the global test rejects in at least POWER_MIN of the replicates when one row of the chain differs by 0.2 in
# two cells (about five standard errors of the difference at these sizes).
POWER_MIN <- 0.9
R_NULL <- 500L
R_POWER <- 100L

# ---- Generating process -----------------------------------------------------------------------------------------
states <- c("A", "B", "C")
transition <- matrix(c(0.5, 0.3, 0.2,
                       0.2, 0.5, 0.3,
                       0.3, 0.2, 0.5), 3L, byrow = TRUE, dimnames = list(states, states))
shifted <- transition
shifted["A", ] <- c(0.3, 0.5, 0.2)
initial <- c(A = 0.6, B = 0.3, C = 0.1)
n_group <- 50L
len <- 20L
reference <- "independent: two groups drawn from one generating chain (null) or from two (alternative)"

# Two groups of sequences, the second from `second`.
rc_two_groups <- function(seed, second = transition) {
  rc_with_seed(seed, list(x = rc_markov_sequences(n_group, len, transition, initial),
                          y = rc_markov_sequences(n_group, len, second, initial)))
}

test_that("permutation() edge and global p-values reject at the nominal level when the groups share a chain", {
  band <- rc_binomial_band(R_NULL, ALPHA, BAND_LEVEL)
  p_values <- rc_replicate(seq_len(R_NULL), \(r) {
    groups <- rc_two_groups(100000L + r)
    found <- permutation(build_network(groups$x, method = "relative"), build_network(groups$y, method = "relative"),
                         iter = ITER, seed = r)
    edges <- summary(found)
    # global: a field of the result (API gap: no accessor returns the global test table).
    c(stats::setNames(edges$p_value, paste(edges$from, "->", edges$to)),
      stats::setNames(found$global$p_value, paste("global", found$global$statistic)))
  })
  rates <- rowMeans(simplify2array(p_values) <= ALPHA)
  note <- sprintf("99.9%% binomial range [%.3f, %.3f]; %d + %d sequences of %d states, %d permutations",
                  band[["lower"]], band[["upper"]], n_group, n_group, len, ITER)
  expect_reference(rates, stats::setNames(rep(ALPHA, length(rates)), names(rates)), reference = reference,
                   covers = c("permutation()", "summary.net_permutation()"),
                   quantity = sprintf("rejection rate under the null, alpha = %.2f, R = %d", ALPHA, R_NULL),
                   tolerance = band[["tolerance"]], note = note)
})

test_that("permutation() global test detects a shifted row of the chain", {
  p_values <- rc_replicate(seq_len(R_POWER), \(r) {
    groups <- rc_two_groups(110000L + r, second = shifted)
    found <- permutation(build_network(groups$x, method = "relative"), build_network(groups$y, method = "relative"),
                         iter = ITER, seed = r)
    stats::setNames(found$global$p_value, paste("global", found$global$statistic))
  })
  power <- rowMeans(simplify2array(p_values) <= ALPHA)
  expect_reference(power[["global M"]], 1, reference = reference, covers = "permutation()",
                   quantity = sprintf("global M: rejection rate under the alternative (power), R = %d", R_POWER),
                   tolerance = 1 - POWER_MIN,
                   note = sprintf("row A: (0.5, 0.3, 0.2) in group 1, (0.3, 0.5, 0.2) in group 2; global S power %.2f",
                                  power[["global S"]]))
})

test_that("sequence_compare() pattern p-values reject at the nominal level when the groups share a chain", {
  band <- rc_binomial_band(R_NULL, ALPHA, BAND_LEVEL)
  p_values <- rc_replicate(seq_len(R_NULL), \(r) {
    groups <- rc_two_groups(120000L + r)
    # sequence_compare() has no seed argument: its permutations are drawn from the session stream, seeded here.
    patterns <- summary(rc_with_seed(r, sequence_compare(rbind(groups$x, groups$y),
                                                         group = rep(c("x", "y"), each = n_group), sub = 2L,
                                                         iter = ITER, adjust = "none")))
    stats::setNames(patterns$p_value, patterns$pattern)[sort(patterns$pattern)]
  })
  rates <- rowMeans(simplify2array(p_values) <= ALPHA)
  expect_reference(rates, stats::setNames(rep(ALPHA, length(rates)), names(rates)), reference = reference,
                   covers = c("sequence_compare()", "summary.net_sequence_comparison()"),
                   quantity = sprintf("rejection rate under the null, bigrams, alpha = %.2f, R = %d", ALPHA, R_NULL),
                   tolerance = band[["tolerance"]],
                   note = sprintf("99.9%% binomial range [%.3f, %.3f]; unadjusted p-values, %d permutations",
                                  band[["lower"]], band[["upper"]], ITER))
})
