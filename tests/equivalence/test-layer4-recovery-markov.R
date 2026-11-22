# Layer 4, recovery of known parameters: first-order Markov chains. Sequences are drawn from a chain whose
# transition matrix and initial distribution are set below, and the estimates are judged against them:
# build_network(method = "relative") and extract_initial_probs() must be unbiased for the transition matrix and the
# initial distribution, build_network(method = "frequency") for the expected transition counts; markov_order_test()
# must reject at its nominal level when the tested order is not needed and with high power when it is (order 1 of a
# first-order chain, order 2 of a chain with a planted second-order rule); bootstrap_network()'s percentile
# intervals must cover the generating probabilities at their nominal rate.
# Simulators and the criteria's arithmetic: helper-layer-recovery.R.

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
# Bias: the mean estimate over R replicates lies within MC_Z Monte Carlo standard errors (sd / sqrt(R)) of the
# generating value, per parameter (Morris, White and Crowther 2019). 3.5 MC SE keeps the chance that an unbiased
# estimator fails one of the 9 cells of a 3-state matrix below 0.5%.
MC_Z <- 3.5
# Rates: a rejection rate under the null, or a coverage, lies within the central 99.9% binomial range of its
# nominal value for R independent replicates. The harness judges |observed - nominal| against the larger of the
# range's two half-widths (rc_binomial_band()); the exact bounds are in each comparison's note.
BAND_LEVEL <- 0.999
ALPHA <- 0.05
# 199 permutations: p = (b + 1) / 200, so a valid permutation p-value has P(p <= 0.05) = 0.05 exactly. A test
# "rejects" when its p-value is at most ALPHA (the definition of a valid p-value); the functions' own
# `significant` flags use p < alpha, whose exact level with 199 permutations is 9 / 200 = 0.045.
N_PERM <- 199L
# Power: where the alternative is far from the null (an order the data plainly need), the test rejects in at least
# POWER_MIN of the replicates.
POWER_MIN <- 0.9
# Replicates.
R_ESTIMATION <- 300L
R_SIZE <- 300L
R_POWER <- 100L
R_COVERAGE <- 200L
BOOT_ITER <- 200L
CI_LEVEL <- 0.05

# ---- Generating process -----------------------------------------------------------------------------------------
states <- c("A", "B", "C")
transition <- matrix(c(0.5, 0.3, 0.2,
                       0.2, 0.5, 0.3,
                       0.3, 0.2, 0.5), 3L, byrow = TRUE, dimnames = list(states, states))
initial <- c(A = 0.6, B = 0.3, C = 0.1)
independent <- matrix(initial, 3L, 3L, byrow = TRUE, dimnames = list(states, states))
n_seq <- 100L
len <- 20L
reference <- "independent: generating parameters of the simulated Markov chain"

test_that("relative transition probabilities, initial probabilities and transition counts are unbiased", {
  estimates <- rc_replicate(seq_len(R_ESTIMATION), \(r) {
    wide <- rc_with_seed(1000L + r, rc_markov_sequences(n_seq, len, transition, initial))
    list(relative = extract_transition_matrix(build_network(wide, method = "relative"))[states, states],
         frequency = extract_transition_matrix(build_network(wide, method = "frequency"))[states, states],
         initial = extract_initial_probs(build_network(wide, method = "relative"))[states])
  })
  relative <- rc_standardised_bias(lapply(estimates, \(e) e$relative), transition)
  frequency <- rc_standardised_bias(lapply(estimates, \(e) e$frequency),
                                    rc_expected_counts(n_seq, len, transition, initial))
  first <- rc_standardised_bias(lapply(estimates, \(e) as.numeric(e$initial)), unname(initial))
  zero_cells <- rc_cells(z = 0 * transition)
  expect_reference(rc_cells(z = relative$z), zero_cells, reference = reference,
                   covers = c("build_network()", "extract_transition_matrix()"),
                   quantity = sprintf("relative: (mean estimate - generating P) / MC SE, R = %d", R_ESTIMATION),
                   tolerance = MC_Z, key = c("from", "to"),
                   note = sprintf("largest |mean - P| %.4f; largest MC SE %.4f", max(abs(relative$mean - transition)),
                                  max(relative$mcse)))
  expect_reference(rc_cells(z = frequency$z), zero_cells, reference = reference,
                   covers = c("build_network()", "extract_transition_matrix()"),
                   quantity = sprintf("frequency: (mean count - expected count) / MC SE, R = %d", R_ESTIMATION),
                   tolerance = MC_Z, key = c("from", "to"),
                   note = "expected count: n_seq * sum over t of P(X_t = i) P_ij, from the initial distribution")
  expect_reference(stats::setNames(first$z, states), stats::setNames(c(0, 0, 0), states), reference = reference,
                   covers = "extract_initial_probs()",
                   quantity = sprintf("initial: (mean estimate - generating probability) / MC SE, R = %d", R_ESTIMATION),
                   tolerance = MC_Z,
                   note = sprintf("largest |mean - initial| %.4f", max(abs(first$mean - initial))))
})

test_that("markov_order_test() rejects an unneeded order at its nominal level", {
  band <- rc_binomial_band(R_SIZE, ALPHA, BAND_LEVEL)
  # Order 1 of independent draws (order 0 is the truth) and order 2 of a first-order chain.
  p_values <- rc_replicate(seq_len(R_SIZE), \(r) {
    iid <- rc_with_seed(2000L + r, rc_markov_sequences(50L, 30L, independent, initial))
    chain <- rc_with_seed(3000L + r, rc_markov_sequences(50L, 30L, transition, initial))
    c(order_1_of_iid = subset(summary(markov_order_test(iid, max_order = 1L, n_perm = N_PERM, seed = r)),
                              order == 1L)$p_permutation,
      order_2_of_chain = subset(summary(markov_order_test(chain, max_order = 2L, n_perm = N_PERM, seed = r)),
                                order == 2L)$p_permutation)
  })
  rates <- rowMeans(simplify2array(p_values) <= ALPHA)
  expect_reference(rates, c(order_1_of_iid = ALPHA, order_2_of_chain = ALPHA), reference = reference,
                   covers = "markov_order_test()",
                   quantity = sprintf("rejection rate under the null, alpha = %.2f, R = %d", ALPHA, R_SIZE),
                   tolerance = band[["tolerance"]],
                   note = sprintf("99.9%% binomial range [%.3f, %.3f]; 50 sequences of 30 states, %d permutations",
                                  band[["lower"]], band[["upper"]], N_PERM))
})

test_that("markov_order_test() detects a needed first order and a planted second-order rule", {
  p_values <- rc_replicate(seq_len(R_POWER), \(r) {
    chain <- rc_with_seed(4000L + r, rc_markov_sequences(50L, 30L, transition, initial))
    planted <- rc_with_seed(5000L + r, rc_planted_sequences(100L, 30L, transition, initial,
                                                             context = c("A", "B"), target = "C", strength = 0.9))
    c(order_1_of_chain = subset(summary(markov_order_test(chain, max_order = 1L, n_perm = N_PERM, seed = r)),
                                order == 1L)$p_permutation,
      order_2_of_planted = subset(summary(markov_order_test(planted, max_order = 2L, n_perm = N_PERM, seed = r)),
                                  order == 2L)$p_permutation)
  })
  power <- rowMeans(simplify2array(p_values) <= ALPHA)
  expect_reference(power, c(order_1_of_chain = 1, order_2_of_planted = 1), reference = reference,
                   covers = "markov_order_test()",
                   quantity = sprintf("rejection rate under the alternative (power), alpha = %.2f, R = %d", ALPHA, R_POWER),
                   tolerance = 1 - POWER_MIN,
                   note = "planted rule: after A, B the next state is C with probability 0.9 (first-order P(C | B) = 0.3)")
})

test_that("bootstrap_network() percentile intervals cover the generating probabilities at the nominal rate", {
  band <- rc_binomial_band(R_COVERAGE, 1 - CI_LEVEL, BAND_LEVEL)
  covered <- rc_replicate(seq_len(R_COVERAGE), \(r) {
    wide <- rc_with_seed(6000L + r, rc_markov_sequences(n_seq, len, transition, initial))
    intervals <- summary(bootstrap_network(build_network(wide, method = "relative"), iter = BOOT_ITER,
                                           ci_level = CI_LEVEL, seed = r))
    truth <- merge(intervals, rc_cells(truth = transition), by = c("from", "to"))
    stats::setNames(truth$ci_lower <= truth$truth & truth$truth <= truth$ci_upper, paste(truth$from, truth$to))
  })
  coverage <- rowMeans(simplify2array(covered))
  expect_reference(coverage, stats::setNames(rep(1 - CI_LEVEL, length(coverage)), names(coverage)),
                   reference = reference, covers = c("bootstrap_network()", "summary.net_bootstrap()"),
                   quantity = sprintf("coverage of the 95%% percentile interval, %d draws, R = %d", BOOT_ITER, R_COVERAGE),
                   tolerance = band[["tolerance"]],
                   note = sprintf("99.9%% binomial range [%.3f, %.3f]; %d sequences of %d states", band[["lower"]],
                                  band[["upper"]], n_seq, len))
})
