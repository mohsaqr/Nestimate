# Layer 4, recovery of known parameters: psychometric networks. Continuous data are drawn from a multivariate
# normal distribution whose precision matrix is a ring (every variable conditionally related to its two neighbours,
# independent of the rest given them), binary data exactly from an Ising model on a ring. The unregularised partial
# correlations (build_network(method = "pcor")) must be unbiased for the generating ones; the EBIC graphical lasso
# (method = "glasso") and the eLasso Ising fit (method = "ising") must recover which edges exist; boot_glasso()'s
# percentile intervals are judged by how often they cover the generating partial correlation.
# Simulators and the criteria's arithmetic: helper-layer-recovery.R.

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
# Bias of the partial correlations: within MC_Z Monte Carlo SEs of the generating value, every pair (Morris, White
# and Crowther 2019). The sample partial correlation is biased by about -rho (1 - rho^2) / (2 n), 3e-4 here, a
# tenth of the Monte Carlo SE of R = 300 replicates.
MC_Z <- 3.5
# Structure: mean sensitivity (share of the generating edges estimated non-zero) at least SENSITIVITY_MIN and mean
# specificity (share of the absent edges estimated zero) at least SPECIFICITY_MIN, over the replicates. With n = 500
# rows and partial correlations of 0.25-0.3 (glasso), or interactions of 1 on the logit scale and n = 1000 (Ising),
# every generating edge is several standard errors from zero, so an estimator selecting by EBIC (Foygel and Drton
# 2010; van Borkulo et al. 2014) should find nearly all of them; the specificity bound leaves room for the
# occasional weak false edge these selectors admit.
SENSITIVITY_MIN <- 0.9
SPECIFICITY_MIN <- 0.8
# Coverage of boot_glasso()'s 95% percentile intervals: the shortfall below 95%, max(0, 0.95 - coverage), is at
# most the distance from 95% to the lower end of the central 99.9% binomial range for R replicates. One-sided: an
# absent edge's interval is pinned at zero by the lasso, so its coverage of the true zero is expected near 1, which
# is not a failure.
BAND_LEVEL <- 0.999
CI_ALPHA <- 0.05
R_PCOR <- 300L
R_STRUCTURE <- 200L
R_ISING <- 100L
R_BOOT <- 100L
BOOT_ITER <- 100L

# ---- Generating processes ---------------------------------------------------------------------------------------
precision <- rc_ring_precision(8L, c(0.3, -0.25, 0.3, 0.3, -0.25, 0.3, 0.3, 0.3))
partial <- rc_partial(precision)
precision_six <- rc_ring_precision(6L, rep(0.3, 6L))
partial_six <- rc_partial(precision_six)
interactions <- 1 * (partial_six != 0)
thresholds <- rep(-1, 6L)
reference <- "independent: generating precision matrix (partial correlations) or Ising interactions"

# The upper-triangle pairs of a symmetric matrix, keyed from < to (names V1..V9).
rc_upper <- function(...) subset(rc_cells(...), from < to)

test_that("build_network(method = 'pcor') partial correlations are unbiased for the generating ones", {
  estimates <- rc_replicate(seq_len(R_PCOR), \(r) {
    data <- rc_with_seed(60000L + r, rc_normal(500L, precision))
    extract_transition_matrix(build_network(data, method = "pcor"))[rownames(partial), colnames(partial)]
  })
  bias <- rc_standardised_bias(estimates, partial)
  expect_reference(rc_upper(z = bias$z), rc_upper(z = 0 * partial), reference = reference,
                   covers = c("build_network()", "build_pcor()"),
                   quantity = sprintf("pcor: (mean estimate - generating partial correlation) / MC SE, R = %d", R_PCOR),
                   tolerance = MC_Z, key = c("from", "to"),
                   note = sprintf("n = 500, p = 8 ring; largest |mean - truth| %.4f", max(abs(bias$mean - partial))))
})

test_that("build_network(method = 'glasso') and method = 'ising' recover the generating edge set", {
  glasso <- rc_replicate(seq_len(R_STRUCTURE), \(r) {
    data <- rc_with_seed(70000L + r, rc_normal(500L, precision))
    rc_edge_recovery(extract_transition_matrix(build_network(data, method = "glasso"))[rownames(partial),
                                                                                     colnames(partial)], partial)
  })
  ising <- rc_replicate(seq_len(R_ISING), \(r) {
    data <- rc_with_seed(80000L + r, rc_ising(1000L, interactions, thresholds))
    rc_edge_recovery(extract_transition_matrix(build_network(data, method = "ising"))[rownames(interactions),
                                                                                    colnames(interactions)],
                     interactions)
  })
  recovery <- list(glasso = rowMeans(simplify2array(glasso)), ising = rowMeans(simplify2array(ising)))
  covers <- list(glasso = c("build_network()", "build_glasso()"), ising = c("build_network()", "build_ising()"))
  replicates <- c(glasso = R_STRUCTURE, ising = R_ISING)
  invisible(lapply(names(recovery), \(method) {
    expect_reference(recovery[[method]][["sensitivity"]], 1, reference = reference, covers = covers[[method]],
                     quantity = sprintf("%s: mean sensitivity, R = %d", method, replicates[[method]]),
                     tolerance = 1 - SENSITIVITY_MIN)
    expect_reference(recovery[[method]][["specificity"]], 1, reference = reference, covers = covers[[method]],
                     quantity = sprintf("%s: mean specificity, R = %d", method, replicates[[method]]),
                     tolerance = 1 - SPECIFICITY_MIN)
  }))
})

test_that("boot_glasso() percentile intervals cover the generating partial correlations", {
  band <- rc_binomial_band(R_BOOT, 1 - CI_ALPHA, BAND_LEVEL)
  truth <- rc_upper(truth = partial_six)
  truth$edge <- paste(truth$from, "--", truth$to)
  replicates <- rc_replicate(seq_len(R_BOOT), \(r) {
    data <- rc_with_seed(90000L + r, rc_normal(500L, precision_six))
    intervals <- merge(summary(boot_glasso(data, iter = BOOT_ITER, cs_iter = 2L, alpha = CI_ALPHA, nlambda = 50L,
                                           centrality = "strength", seed = r)), truth, by = "edge")
    list(covered = stats::setNames(intervals$ci_lower <= intervals$truth & intervals$truth <= intervals$ci_upper,
                                   intervals$edge),
         weight = stats::setNames(intervals$weight, intervals$edge),
         contains_estimate = stats::setNames(intervals$ci_lower <= intervals$weight &
                                               intervals$weight <= intervals$ci_upper, intervals$edge))
  })
  summarise <- \(field) rowMeans(simplify2array(lapply(replicates, \(e) e[[field]])))
  coverage <- summarise("covered")
  present <- truth$edge[truth$truth != 0]
  absent <- setdiff(names(coverage), present)
  shortfall <- pmax((1 - CI_ALPHA) - coverage, 0)
  diagnosis <- sprintf(paste0("coverage %s; lower end of the 99.9%% binomial range %.3f; mean EBIC-glasso estimate of",
                              " the present edges (truth 0.3) %.3f; share of intervals containing their own estimate",
                              " %.3f; n = 500, p = 6 ring of 0.3, %d draws"),
                       toString(sprintf("%s %.2f", names(coverage), coverage)), band[["lower"]],
                       mean(summarise("weight")[present]), mean(summarise("contains_estimate")), BOOT_ITER)
  # A property of the method, not of the code (diagnosed on the first run of this file): the lasso shrinks every
  # present edge toward zero, and the percentile bootstrap re-applies the same shrinkage to each resample, so the
  # interval sits around the shrunken estimate, not around the generating value. Bootstrap intervals of
  # lasso-regularised estimates are not confidence intervals for the true parameter (Chatterjee and Lahiri 2011,
  # JASA 106:608-625); bootnet's authors say the same of these edge intervals (Epskamp, Borsboom and Fried 2018,
  # Behavior Research Methods 50:195-212). The criterion is unchanged; the result is recorded as known.
  expect_known_difference(shortfall[present], stats::setNames(rep(0, length(present)), present), reference = reference,
                          covers = c("boot_glasso()", "summary.boot_glasso()"),
                          quantity = sprintf("present edges: shortfall of coverage below 95%%, R = %d", R_BOOT),
                          tolerance = (1 - CI_ALPHA) - band[["lower"]], note = diagnosis,
                          reason = paste("percentile bootstrap of an L1-penalised estimator: the intervals follow the",
                                         "shrunken estimate, not the generating partial correlation",
                                         "(Chatterjee and Lahiri 2011; Epskamp, Borsboom and Fried 2018)"))
  expect_reference(shortfall[absent], stats::setNames(rep(0, length(absent)), absent), reference = reference,
                   covers = c("boot_glasso()", "summary.boot_glasso()"),
                   quantity = sprintf("absent edges: shortfall of coverage below 95%%, R = %d", R_BOOT),
                   tolerance = (1 - CI_ALPHA) - band[["lower"]], note = diagnosis)
})
