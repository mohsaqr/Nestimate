# Layer 4, recovery of known parameters: HYPA path anomalies (LaRock et al. 2020) on data whose second-order
# structure is known. Under a first-order Markov chain no order-2 path is anomalous, so the share of paths
# build_hypa() flags is its false-positive rate; with one planted second-order rule (after A, B the next state is C
# far more often than the first-order chain says) the planted path must be flagged over-represented.
# CRAN 0.9.24 flagged 174 paths of first-order data as anomalous (its null model was wrong); this file keeps the
# corrected behaviour under watch. Simulators and the criteria's arithmetic: helper-layer-recovery.R.

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
ALPHA <- 0.05
# False positives: each order-2 path is tested in two directions (over and under), each at level ALPHA. A valid
# test flags at most ALPHA of the paths per direction, up to Monte Carlo error. The paths of one replicate share
# their data, so the Monte Carlo SE of the mean flagged share is taken from the replicate-to-replicate spread
# (sd of the per-replicate shares / sqrt(R)), not from a binomial count of paths. The criterion is one-sided:
# the excess of the mean share over ALPHA, max(0, share - ALPHA), must be within MC_Z Monte Carlo SEs of 0. A
# share below ALPHA (a conservative null) passes; HYPA's hypergeometric null is wider than the conditional
# multinomial variation of a first-order chain, so a share below ALPHA is expected.
MC_Z <- 3.5
# Power: the planted path is flagged over-represented in at least POWER_MIN of the replicates.
POWER_MIN <- 0.9
R_NULL <- 500L
R_POWER <- 100L

# ---- Generating process -----------------------------------------------------------------------------------------
states <- c("A", "B", "C")
transition <- matrix(c(0.5, 0.3, 0.2,
                       0.2, 0.5, 0.3,
                       0.3, 0.2, 0.5), 3L, byrow = TRUE, dimnames = list(states, states))
initial <- c(A = 0.6, B = 0.3, C = 0.1)
four <- c("A", "B", "C", "D")
transition_four <- matrix(c(0.4, 0.3, 0.2, 0.1,
                            0.1, 0.4, 0.3, 0.2,
                            0.2, 0.1, 0.4, 0.3,
                            0.3, 0.2, 0.1, 0.4), 4L, byrow = TRUE, dimnames = list(four, four))
initial_four <- c(A = 0.25, B = 0.25, C = 0.25, D = 0.25)
reference <- "independent: generating first-order chain (no order-2 path is anomalous) and its planted rule"

# The flagged shares of one design, and the one-sided excess over ALPHA with its Monte Carlo SE.
rc_hypa_excess <- function(shares) {
  stack <- simplify2array(shares)
  share <- rowMeans(stack)
  mcse <- apply(stack, 1L, stats::sd) / sqrt(ncol(stack))
  list(share = share, mcse = mcse, excess = pmax(share - ALPHA, 0))
}

test_that("build_hypa() flags at most alpha of the order-2 paths of a first-order chain, per direction", {
  designs <- list(
    large_counts = list(n_seq = 200L, len = 30L, transition = transition, initial = initial, seed = 10000L),
    small_counts = list(n_seq = 40L, len = 12L, transition = transition_four, initial = initial_four, seed = 20000L)
  )
  invisible(lapply(names(designs), \(name) {
    design <- designs[[name]]
    shares <- rc_replicate(seq_len(R_NULL), \(r) {
      wide <- rc_with_seed(design$seed + r, rc_markov_sequences(design$n_seq, design$len, design$transition,
                                                                design$initial))
      fit <- build_hypa(wide, order = 2L, alpha = ALPHA)
      # n_over, n_under, n_edges: fields of the result (API gap: no accessor returns the counts as a table).
      c(over = fit$n_over / fit$n_edges, under = fit$n_under / fit$n_edges)
    })
    found <- rc_hypa_excess(shares)
    note <- sprintf("mean flagged share: over %.4f (MC SE %.4f), under %.4f (MC SE %.4f); %d sequences of %d states",
                    found$share[["over"]], found$mcse[["over"]], found$share[["under"]], found$mcse[["under"]],
                    design$n_seq, design$len)
    expect_reference(found$excess[["over"]], 0, reference = reference, covers = "build_hypa()",
                     quantity = sprintf("%s: excess of the over-represented share over alpha, R = %d", name, R_NULL),
                     tolerance = MC_Z * found$mcse[["over"]], note = note)
    expect_reference(found$excess[["under"]], 0, reference = reference, covers = "build_hypa()",
                     quantity = sprintf("%s: excess of the under-represented share over alpha, R = %d", name, R_NULL),
                     tolerance = MC_Z * found$mcse[["under"]], note = note)
  }))
})

test_that("build_hypa() flags a planted second-order path as over-represented", {
  flagged <- rc_replicate(seq_len(R_POWER), \(r) {
    wide <- rc_with_seed(30000L + r, rc_planted_sequences(200L, 30L, transition, initial, context = c("A", "B"),
                                                          target = "C", strength = 0.9))
    scores <- build_hypa(wide, order = 2L, alpha = ALPHA)$scores
    identical(subset(scores, path == "A -> B -> C")$anomaly, "over")
  })
  expect_reference(mean(unlist(flagged)), 1, reference = reference, covers = "build_hypa()",
                   quantity = sprintf("share of replicates flagging A -> B -> C over-represented, R = %d", R_POWER),
                   tolerance = 1 - POWER_MIN,
                   note = "after A, B the next state is C with probability 0.9 (first-order P(C | B) = 0.3)")
})
