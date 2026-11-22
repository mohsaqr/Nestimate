# Layer 4, recovery of known parameters: the Network Comparison Test (van Borkulo et al. 2022). Two groups are drawn
# from multivariate normal distributions with known precision matrices: the same ring of partial correlations in
# both (the null is true), or a ring with three of its six edges removed in the second (it is false). nct()
# permutes rows between the groups, so under the null its structure (M), global strength (S) and edge p-values
# must reject at their nominal level; under the alternative M and S must reject with high power.
# Simulators and the criteria's arithmetic: helper-layer-recovery.R.

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
# Rejection rates under the null lie within the central 99.9% binomial range of ALPHA for R independent replicates,
# for M, S and each edge separately (unadjusted p-values). The harness judges |rate - ALPHA| against the larger
# half-width of that range. Edges absent from the ring are estimated zero in most permutations, so their tests are
# expected to be conservative (rates below ALPHA, inside the range's lower part).
BAND_LEVEL <- 0.999
ALPHA <- 0.05
# 39 permutations: p = (b + 1) / 40, so a valid p-value has P(p <= 0.05) = 2 / 40 = 0.05 exactly. Each permutation
# fits two EBIC graphical lassos (about 0.07 s each on the CI machine), which sets the number.
ITER <- 39L
# Power: with three ring edges of 0.3 removed in group 2 and n = 500 per group, M and S reject in at least
# POWER_MIN of the replicates.
POWER_MIN <- 0.9
R_NULL <- 100L
R_POWER <- 40L

# ---- Generating processes ---------------------------------------------------------------------------------------
full <- rc_ring_precision(6L, rep(0.3, 6L))
sparse <- rc_ring_precision(6L, rep(0.3, 6L), drop = c(1L, 3L, 5L))
reference <- "independent: two groups drawn from one generating precision matrix (null) or from two (alternative)"

# nct() has no seed argument: its permutations are drawn from the session stream, seeded here.
rc_nct_p_values <- function(seed, n, second) {
  rc_with_seed(seed, {
    x <- rc_normal(n, full)
    y <- rc_normal(n, second)
    found <- nct(x, y, iter = ITER)
    global <- summary(found)
    edges <- as.data.frame(found)
    c(stats::setNames(global$p_value, global$statistic), stats::setNames(edges$p_value, paste(edges$from, "--", edges$to)))
  })
}

test_that("nct() structure, strength and edge p-values reject at the nominal level when the groups share a network", {
  band <- rc_binomial_band(R_NULL, ALPHA, BAND_LEVEL)
  p_values <- rc_replicate(seq_len(R_NULL), \(r) rc_nct_p_values(130000L + r, 200L, full))
  rates <- rowMeans(simplify2array(p_values) <= ALPHA)
  expect_reference(rates, stats::setNames(rep(ALPHA, length(rates)), names(rates)), reference = reference,
                   covers = c("nct()", "summary.net_nct()", "as.data.frame.net_nct()"),
                   quantity = sprintf("rejection rate under the null, alpha = %.2f, R = %d", ALPHA, R_NULL),
                   tolerance = band[["tolerance"]],
                   note = sprintf("99.9%% binomial range [%.3f, %.3f]; n = 200 per group, p = 6 ring of 0.3, %d permutations",
                                  band[["lower"]], band[["upper"]], ITER))
})

test_that("nct() detects three removed edges", {
  p_values <- rc_replicate(seq_len(R_POWER), \(r) rc_nct_p_values(140000L + r, 500L, sparse))
  power <- rowMeans(simplify2array(p_values) <= ALPHA)
  expect_reference(power[c("M", "S")], c(M = 1, S = 1), reference = reference,
                   covers = c("nct()", "summary.net_nct()"),
                   quantity = sprintf("M and S: rejection rate under the alternative (power), R = %d", R_POWER),
                   tolerance = 1 - POWER_MIN,
                   note = "group 2 ring without edges V1--V2, V3--V4, V5--V6; n = 500 per group")
})
