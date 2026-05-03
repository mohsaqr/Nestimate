# ---- Equivalence: centrality_stability() vs bootnet::corStability() ----
#
# Validates Nestimate::centrality_stability() against bootnet::corStability()
# run on bootnet::bootnet(..., type = "case"). Both implement the same
# case-dropping CS-coefficient (Epskamp et al. 2018, Eq. 4): the maximum
# drop proportion at which >95% of bootstrap subsets correlate with the
# full-sample centrality at >= 0.7.
#
# Configuration is matched as tightly as possible:
#   * same continuous data (multivariate normal, AR(1)-like Sigma)
#   * same estimator: method = "glasso" in Nestimate (EBICglasso in bootnet)
#   * same drop grid: seq(0.1, 0.9, by = 0.1) vs caseMin/Max/N = 0.1/0.9/9
#   * same nBoots (Nestimate iter)
#   * same threshold (0.7) and certainty (0.95)
#
# What we DO NOT match:
#   * RNG stream. Nestimate uses sample.int(n_seq, n_keep, replace = FALSE)
#     per iteration; bootnet uses a different internal sampling order.
#   * Centrality internals. Nestimate uses Floyd-Warshall-based
#     closeness/betweenness on abs(weights) with inverse-weight distances;
#     bootnet uses qgraph/igraph's centrality_auto which relies on
#     igraph::distances + betweenness. These differ slightly on ties, on
#     disconnected components, and on edge-sign handling (bootnet signs are
#     squared via centralityTable; Nestimate takes |w|). Both give the same
#     centrality ORDERING on the full sample but may diverge on small
#     subsets, which adds a small amount of noise to the Spearman/Pearson
#     correlations inside the CS calculation.
#
# Because of (1) + (2), we cannot expect TOL = 1e-8. Empirical probing on
# 20+ configs shows CS-coefficients agree to within ~0.1-0.3. We therefore
# use TOL = 0.3 (roughly 3 drop-grid steps) and require that at least
# 80% of (config, measure) pairs come in under that tolerance. We also log
# exact deltas to CSV so regressions show up quantitatively.
#
# Note: build_network() stores $data as a raw numeric matrix for
# association methods (glasso/pcor/cor). .prepare_association_input was
# previously hard-failing on non-square matrices, which broke
# centrality_stability's resample-and-re-estimate loop. Fixed 2026-04-21
# by teaching .prepare_association_input to coerce non-square matrices
# to data.frame and reuse the data-frame branch (see R/estimators.R).

TOL_CS <- 0.3           # drop-grid tolerance; documented rationale above
TOL_PASS_RATE <- 0.80   # require >= 80% of checks within TOL_CS

test_that("centrality_stability CS-coefficients match bootnet::corStability (20 configs)", {
  skip_equiv_tests()
  skip_if_not_installed("bootnet")
  skip_if_not_installed("qgraph")

  rep <- equiv_report()

  drop_grid <- seq(0.1, 0.9, by = 0.1)
  measures <- c("InStrength", "Closeness", "Betweenness")
  boot_names <- c(InStrength = "strength",
                  Closeness  = "closeness",
                  Betweenness = "betweenness")

  # 20 random configs spanning n, p, rho, seed. Expensive: each config
  # runs two case-drop bootstraps at nBoots = 60, so ~10-20s per config.
  set.seed(20260420L)
  configs <- lapply(seq_len(20L), function(i) {
    list(
      n    = sample(c(150L, 200L, 250L, 300L), 1L),
      p    = sample(c(5L, 6L, 7L), 1L),
      rho  = sample(c(0.4, 0.5, 0.6, 0.7), 1L),
      seed = 1000L + i
    )
  })

  run_config <- function(cfg) {
    # Generate AR(1)-structured continuous data
    set.seed(cfg$seed)
    Sigma <- matrix(0, cfg$p, cfg$p)
    for (i in seq_len(cfg$p)) {
      for (j in seq_len(cfg$p)) {
        Sigma[i, j] <- cfg$rho ^ abs(i - j)
      }
    }
    L <- chol(Sigma)
    mat <- matrix(rnorm(cfg$n * cfg$p), cfg$n, cfg$p) %*% L
    df  <- as.data.frame(mat)
    colnames(df) <- paste0("V", seq_len(cfg$p))

    # Nestimate side
    net <- build_network(df, method = "glasso")
cs_nest <- suppressWarnings(suppressMessages(
      centrality_stability(
        net,
        measures = measures,
        iter = 60L,
        drop_prop = drop_grid,
        threshold = 0.7,
        certainty = 0.95,
        seed = cfg$seed
      )
    ))

    # bootnet reference
    net_b <- suppressWarnings(suppressMessages(
      bootnet::estimateNetwork(df, default = "EBICglasso", verbose = FALSE)
    ))
    cs_boot <- suppressWarnings(suppressMessages({
      b <- bootnet::bootnet(
        net_b,
        nBoots = 60L,
        type = "case",
        caseMin = 0.1, caseMax = 0.9, caseN = 9L,
        statistics = c("strength", "closeness", "betweenness"),
        verbose = FALSE, nCores = 1L
      )
      bootnet::corStability(b, verbose = FALSE)
    }))

    # Collate CS per measure (default 0 when missing, matching bootnet's
    # max0() when corStability returns nothing for a statistic)
    get_cs <- function(x, nm) {
      v <- tryCatch(x[[nm]], error = function(e) NA_real_)
      if (is.null(v) || !is.finite(v)) 0 else as.numeric(v)
    }

    data.frame(
      config      = sprintf("n=%d p=%d rho=%.1f seed=%d",
                            cfg$n, cfg$p, cfg$rho, cfg$seed),
      measure     = measures,
      cs_nest     = vapply(measures, function(m) get_cs(cs_nest$cs, m), numeric(1)),
      cs_boot     = vapply(measures, function(m) get_cs(cs_boot, boot_names[m]),
                           numeric(1)),
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, lapply(configs, run_config))
  results$delta <- abs(results$cs_nest - results$cs_boot)

  # Per-(config, measure) log
  invisible(lapply(seq_len(nrow(results)), function(i) {
    r <- results[i, ]
    rep$log(
      func           = "centrality_stability",
      config         = sprintf("%s | %s", r$config, r$measure),
      n_checked      = 1L,
      n_failed       = as.integer(r$delta > TOL_CS),
      max_abs_err    = r$delta,
      mean_abs_err   = r$delta,
      median_abs_err = r$delta,
      p95_abs_err    = r$delta,
      reference      = "bootnet::corStability (type='case')",
      notes          = sprintf("nest=%.2f boot=%.2f", r$cs_nest, r$cs_boot)
    )
  }))

  rep$write_csv("centrality_stability")
  message(rep$summary())

  # Structural invariants -- must ALL hold regardless of Monte Carlo noise.
  # CS is the MAX drop proportion in drop_prop (0.1..0.9) where stability
  # holds, so the valid range is [0, max(drop_grid)] = [0, 0.9].
  expect_true(all(results$cs_nest >= 0 & results$cs_nest <= 0.9),
              info = "CS-coefficients must lie in [0, 0.9] (max drop grid value)")
  expect_true(all(results$cs_boot >= 0 & results$cs_boot <= 0.9),
              info = "bootnet CS-coefficients must lie in [0, 0.9]")

  # Numerical agreement: require the majority to agree within TOL_CS.
  # We do NOT require all because two independent Monte Carlo runs at
  # nBoots = 60 on a discrete 0.1-step grid have an unavoidable +/- one
  # grid step of jitter, plus centrality-algorithm divergence on small
  # subsets. See header comment for full rationale.
  n_pass <- sum(results$delta <= TOL_CS)
  pass_rate <- n_pass / nrow(results)
  expect_gte(
    pass_rate, TOL_PASS_RATE,
    label = sprintf(
      "CS agreement rate within TOL=%.2f: %d/%d = %.1f%% (max delta=%.2f, median=%.2f)",
      TOL_CS, n_pass, nrow(results), 100 * pass_rate,
      max(results$delta), stats::median(results$delta)
    )
  )
})


# ---- Structural / invariant checks (cheap, no bootnet) -----------------
#
# These verify properties of Nestimate::centrality_stability() that must
# hold regardless of whether bootnet is available. Guards against
# regressions in the CS-coefficient formula and the correlation grid.

test_that("centrality_stability mean-correlation curve is near-monotone decreasing", {
  skip_equiv_tests()

  set.seed(7L)
  p <- 6L; n <- 250L
  Sigma <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) Sigma[i, j] <- 0.6 ^ abs(i - j)
  L <- chol(Sigma)
  df <- as.data.frame(matrix(rnorm(n * p), n, p) %*% L)
  colnames(df) <- paste0("V", seq_len(p))

  net <- build_network(df, method = "glasso")

  cs <- suppressWarnings(suppressMessages(
    centrality_stability(
      net,
      measures  = c("InStrength", "Closeness"),
      iter      = 80L,
      drop_prop = seq(0.1, 0.9, by = 0.1),
      seed      = 7L
    )
  ))

  summ <- summary(cs)
  # Check that mean correlation at drop=0.1 exceeds drop=0.8 for every measure.
  # (Not strict monotonicity -- Monte Carlo can create tiny local bumps.)
  endpoints <- do.call(rbind, lapply(split(summ, summ$measure), function(s) {
    data.frame(measure = s$measure[1],
               low_drop  = s$mean_cor[s$drop_prop == 0.1],
               high_drop = s$mean_cor[s$drop_prop == 0.8],
               stringsAsFactors = FALSE)
  }))
  expect_true(all(endpoints$low_drop > endpoints$high_drop),
              info = sprintf("Endpoints: %s",
                             paste(sprintf("%s: %.3f->%.3f",
                                           endpoints$measure,
                                           endpoints$low_drop,
                                           endpoints$high_drop),
                                   collapse = "; ")))
})


test_that("centrality_stability CS-coefficient matches manual .calculate_cs on stored correlations", {
  skip_equiv_tests()

  set.seed(99L)
  p <- 5L; n <- 200L
  Sigma <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) Sigma[i, j] <- 0.6 ^ abs(i - j)
  L <- chol(Sigma)
  df <- as.data.frame(matrix(rnorm(n * p), n, p) %*% L)
  colnames(df) <- paste0("V", seq_len(p))

  net <- build_network(df, method = "glasso")

  cs <- suppressWarnings(suppressMessages(
    centrality_stability(
      net,
      measures  = c("InStrength", "Closeness"),
      iter      = 60L,
      drop_prop = seq(0.1, 0.9, by = 0.1),
      threshold = 0.7,
      certainty = 0.95,
      seed      = 42L
    )
  ))

  # Recompute CS directly from stored correlation matrices using the
  # public definition: max drop_prop where >= 95% of iters exceed 0.7.
  recomputed <- vapply(names(cs$cs), function(m) {
    corr_mat   <- cs$correlations[[m]]
    prop_above <- colMeans(corr_mat >= 0.7, na.rm = TRUE)
    valid      <- which(prop_above >= 0.95)
    if (length(valid) == 0L) 0 else cs$drop_prop[max(valid)]
  }, numeric(1))

  expect_equal(unname(cs$cs), unname(recomputed), tolerance = 1e-12,
               info = "CS-coefficient must match direct recomputation from $correlations")
})
