# ---- Equivalence Tests: Bootstrap & Permutation ----
#
# Tests bootstrap/permutation against EXTERNAL references:
#   1-2. Statistical property tests (CI coverage, mean convergence)
#   3-4. Permutation property tests (null uniformity, alternative detection)
#   5.   Bootstrap vs tna::bootstrap() — means, SDs, CIs (20 configs)
#   6.   Permutation vs tna::permutation_test() — true diffs (10 configs)
#
# Requires: NESTIMATE_EQUIV_TESTS=true

# ---- Config generation ----
set.seed(7777)
N_BOOT <- 20L
N_PERM <- 10L
boot_configs <- lapply(seq_len(N_BOOT), function(i) list(seed = sample.int(100000, 1)))
perm_configs <- lapply(seq_len(N_PERM), function(i) list(seed = sample.int(100000, 1)))


# ---- 1. Bootstrap CI coverage ----

test_that("bootstrap CI coverage >= 0.80 across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  coverage_rates <- vapply(seq_len(N_BOOT), function(i) {
    cfg <- boot_configs[[i]]
    data <- simulate_sequences(n_actors = 30, n_states = 4,
                               seq_length = 20, seed = cfg$seed)
    net <- build_network(data, method = "relative")
    boot <- bootstrap_network(net, iter = 200, seed = 1)

    s <- summary(boot)
    true_mat <- net$weights
    node_labels <- rownames(true_mat)

    # Extract true weight for each edge in summary
    covered <- vapply(seq_len(nrow(s)), function(j) {
      tw <- true_mat[s$from[j], s$to[j]]
      tw >= s$ci_lower[j] & tw <= s$ci_upper[j]
    }, logical(1))

    n_edges <- length(covered)
    rate <- mean(covered)
    ci_widths <- s$ci_upper - s$ci_lower

    report$log(
      func = "bootstrap_ci_coverage",
      config = sprintf("seed=%d", cfg$seed),
      n_checked = n_edges,
      n_failed = sum(!covered),
      max_abs_err = max(ci_widths),
      mean_abs_err = mean(ci_widths),
      median_abs_err = stats::median(ci_widths),
      p95_abs_err = stats::quantile(ci_widths, 0.95),
      reference = "self",
      notes = sprintf("coverage=%.3f", rate)
    )

    rate
  }, numeric(1))

  report$write_csv("bootstrap")

  # Each config should achieve >= 0.80 coverage (conservative for iter=200)
  lapply(seq_len(N_BOOT), function(i) {
    expect_true(
      coverage_rates[i] >= 0.80,
      info = sprintf("Config %d (seed=%d): coverage=%.3f < 0.80",
                     i, boot_configs[[i]]$seed, coverage_rates[i])
    )
  })
})


# ---- 2. Bootstrap mean convergence ----

test_that("bootstrap means converge to full-data weights across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  max_devs <- vapply(seq_len(N_BOOT), function(i) {
    cfg <- boot_configs[[i]]
    data <- simulate_sequences(n_actors = 30, n_states = 4,
                               seq_length = 20, seed = cfg$seed)
    net <- build_network(data, method = "relative")
    boot <- bootstrap_network(net, iter = 500, seed = 1)

    # Compare bootstrap mean matrix to original weights
    true_w <- net$weights
    boot_mean <- boot$mean

    abs_dev <- abs(boot_mean - true_w)
    max_dev <- max(abs_dev)
    mean_dev <- mean(abs_dev)

    report$log(
      func = "bootstrap_mean_convergence",
      config = sprintf("seed=%d", cfg$seed),
      n_checked = length(abs_dev),
      n_failed = sum(abs_dev > 0.15),
      max_abs_err = max_dev,
      mean_abs_err = mean_dev,
      median_abs_err = stats::median(abs_dev),
      p95_abs_err = stats::quantile(abs_dev, 0.95),
      reference = "self",
      notes = sprintf("max_dev=%.4f, mean_dev=%.4f", max_dev, mean_dev)
    )

    max_dev
  }, numeric(1))

  report$write_csv("bootstrap")

  # Max absolute deviation between bootstrap mean and original should be < 0.15
  lapply(seq_len(N_BOOT), function(i) {
    expect_true(
      max_devs[i] < 0.15,
      info = sprintf("Config %d (seed=%d): max_dev=%.4f >= 0.15",
                     i, boot_configs[[i]]$seed, max_devs[i])
    )
  })
})


# ---- 3. Permutation null distribution ----

test_that("permutation p-values are not clustered near 0 under the null (10 configs)", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  null_props <- vapply(seq_len(N_PERM), function(i) {
    cfg <- perm_configs[[i]]

    # Two datasets from the SAME distribution (same params, different seeds)
    seed_a <- cfg$seed
    seed_b <- cfg$seed + 50000L
    data1 <- simulate_sequences(n_actors = 30, n_states = 4,
                                seq_length = 20, seed = seed_a)
    data2 <- simulate_sequences(n_actors = 30, n_states = 4,
                                seq_length = 20, seed = seed_b)

    net1 <- build_network(data1, method = "relative")
    net2 <- build_network(data2, method = "relative")
    perm <- permutation_test(net1, net2, iter = 200, seed = 1)

    pvals <- as.vector(perm$p_values)
    # Exclude diagonal (self-loops always 0)
    n_nodes <- nrow(perm$p_values)
    diag_idx <- seq(1, n_nodes * n_nodes, by = n_nodes + 1)
    pvals_off <- pvals[-diag_idx]

    prop_nonsig <- mean(pvals_off > 0.05)
    min_p <- min(pvals_off)
    max_p <- max(pvals_off)
    med_p <- stats::median(pvals_off)

    report$log(
      func = "permutation_null",
      config = sprintf("seed=%d", cfg$seed),
      n_checked = length(pvals_off),
      n_failed = sum(pvals_off <= 0.05),
      max_abs_err = 1 - prop_nonsig,
      mean_abs_err = mean(pvals_off),
      median_abs_err = med_p,
      p95_abs_err = stats::quantile(pvals_off, 0.05),
      reference = "self",
      notes = sprintf("prop_nonsig=%.3f, min_p=%.4f, max_p=%.4f, med_p=%.4f",
                       prop_nonsig, min_p, max_p, med_p)
    )

    prop_nonsig
  }, numeric(1))

  report$write_csv("bootstrap")

  # At least 70% of edges should have p > 0.05 under the null
  lapply(seq_len(N_PERM), function(i) {
    expect_true(
      null_props[i] >= 0.70,
      info = sprintf("Config %d (seed=%d): prop(p>0.05)=%.3f < 0.70",
                     i, perm_configs[[i]]$seed, null_props[i])
    )
  })
})


# ---- 4. Permutation alternative detection ----

test_that("permutation test detects genuine differences (10 configs)", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  sig_counts <- vapply(seq_len(N_PERM), function(i) {
    cfg <- perm_configs[[i]]

    # Group 1: uniform distribution across states
    set.seed(cfg$seed)
    states <- LETTERS[seq_len(4)]
    data1 <- as.data.frame(
      matrix(sample(states, 30 * 20, replace = TRUE), nrow = 30, ncol = 20),
      stringsAsFactors = FALSE
    )
    colnames(data1) <- paste0("T", seq_len(20))

    # Group 2: heavily biased distribution (first state dominates)
    set.seed(cfg$seed + 50000L)
    biased_probs <- c(0.70, 0.10, 0.10, 0.10)
    data2 <- as.data.frame(
      matrix(sample(states, 30 * 20, replace = TRUE, prob = biased_probs),
             nrow = 30, ncol = 20),
      stringsAsFactors = FALSE
    )
    colnames(data2) <- paste0("T", seq_len(20))

    net1 <- build_network(data1, method = "relative")
    net2 <- build_network(data2, method = "relative")
    perm <- permutation_test(net1, net2, iter = 200, seed = 1)

    pvals <- as.vector(perm$p_values)
    n_nodes <- nrow(perm$p_values)
    diag_idx <- seq(1, n_nodes * n_nodes, by = n_nodes + 1)
    pvals_off <- pvals[-diag_idx]

    n_sig <- sum(pvals_off < 0.05)
    min_p <- min(pvals_off)

    report$log(
      func = "permutation_alternative",
      config = sprintf("seed=%d", cfg$seed),
      n_checked = length(pvals_off),
      n_failed = length(pvals_off) - n_sig,
      max_abs_err = min_p,
      mean_abs_err = mean(pvals_off),
      median_abs_err = stats::median(pvals_off),
      p95_abs_err = stats::quantile(pvals_off, 0.05),
      reference = "self",
      notes = sprintf("n_sig=%d, min_p=%.6f", n_sig, min_p)
    )

    n_sig
  }, numeric(1))

  report$write_csv("bootstrap")

  # At least some edges should show significant differences
  lapply(seq_len(N_PERM), function(i) {
    expect_true(
      sig_counts[i] >= 1L,
      info = sprintf("Config %d (seed=%d): n_sig=%d, expected >= 1",
                     i, perm_configs[[i]]$seed, sig_counts[i])
    )
  })
})


# ---- 5. Bootstrap vs tna::bootstrap() (20 configs) ----

# ---- 5. Bootstrap vs tna::bootstrap() — per-edge (20 configs) ----

test_that("bootstrap: every edge mean/SD/CI matches tna::bootstrap() across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  report <- equiv_report()
  TOL <- 1e-10

  invisible(lapply(seq_len(N_BOOT), function(i) {
    cfg  <- boot_configs[[i]]
    data <- simulate_sequences(n_actors = 20, n_states = 4,
                               seq_length = 15, seed = cfg$seed)
    iter <- 200L

    tna_model <- tna::tna(data)
    set.seed(1)
    tna_boot <- tna::bootstrap(tna_model, iter = iter)

    nest_net  <- build_network(data, method = "relative")
    nest_boot <- bootstrap_network(nest_net, iter = iter, seed = 1)

    states <- rownames(nest_boot$mean)
    n <- length(states)
    lbl <- sprintf("cfg%d(seed=%d)", i, cfg$seed)

    # Per-edge: means
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        our_val <- nest_boot$mean[ii, jj]
        tna_val <- tna_boot$weights_mean[ii, jj]
        if (is.na(our_val)) our_val <- 0
        if (is.na(tna_val)) tna_val <- 0
        expect_equal(our_val, tna_val, tolerance = TOL,
                     label = sprintf("boot mean %s [%s -> %s]",
                                     lbl, states[ii], states[jj]))
      })
    }))

    # Per-edge: SDs
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        our_val <- nest_boot$sd[ii, jj]
        tna_val <- tna_boot$weights_sd[ii, jj]
        if (is.na(our_val)) our_val <- 0
        if (is.na(tna_val)) tna_val <- 0
        expect_equal(our_val, tna_val, tolerance = TOL,
                     label = sprintf("boot SD %s [%s -> %s]",
                                     lbl, states[ii], states[jj]))
      })
    }))

    # Per-edge: CI lower
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        our_val <- nest_boot$ci_lower[ii, jj]
        tna_val <- tna_boot$ci_lower[ii, jj]
        if (is.na(our_val)) our_val <- 0
        if (is.na(tna_val)) tna_val <- 0
        expect_equal(our_val, tna_val, tolerance = TOL,
                     label = sprintf("boot CI_lo %s [%s -> %s]",
                                     lbl, states[ii], states[jj]))
      })
    }))

    # Per-edge: CI upper
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        our_val <- nest_boot$ci_upper[ii, jj]
        tna_val <- tna_boot$ci_upper[ii, jj]
        if (is.na(our_val)) our_val <- 0
        if (is.na(tna_val)) tna_val <- 0
        expect_equal(our_val, tna_val, tolerance = TOL,
                     label = sprintf("boot CI_hi %s [%s -> %s]",
                                     lbl, states[ii], states[jj]))
      })
    }))

    # Per-edge: significant edges
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        our_val <- nest_boot$significant[ii, jj]
        tna_val <- tna_boot$weights_sig[ii, jj]
        if (is.na(our_val)) our_val <- 0
        if (is.na(tna_val)) tna_val <- 0
        expect_equal(our_val, tna_val, tolerance = TOL,
                     label = sprintf("boot sig %s [%s -> %s]",
                                     lbl, states[ii], states[jj]))
      })
    }))

    n_checked <- n * n * 5L  # mean + sd + ci_lo + ci_hi + sig
    delta_all <- c(
      abs(replace(nest_boot$mean, is.na(nest_boot$mean), 0) -
          replace(tna_boot$weights_mean, is.na(tna_boot$weights_mean), 0)),
      abs(replace(nest_boot$sd, is.na(nest_boot$sd), 0) -
          replace(tna_boot$weights_sd, is.na(tna_boot$weights_sd), 0))
    )
    report$log("bootstrap_vs_tna", lbl, n_checked, sum(delta_all >= TOL),
               max(delta_all), mean(delta_all), median(delta_all),
               quantile(delta_all, 0.95, names = FALSE), "tna::bootstrap()")
  }))

  report$write_csv("bootstrap")
})


# ---- 6. Permutation vs tna::permutation_test() — per-edge (10 configs) ----

test_that("permutation: every edge diff matches tna across 10 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  report <- equiv_report()
  TOL <- 1e-12

  invisible(lapply(seq_len(N_PERM), function(i) {
    cfg <- perm_configs[[i]]
    data1 <- simulate_sequences(n_actors = 20, n_states = 4,
                                seq_length = 15, seed = cfg$seed)
    data2 <- simulate_sequences(n_actors = 20, n_states = 4,
                                seq_length = 15, seed = cfg$seed + 50000L)
    iter <- 200L

    m1_tna <- tna::tna(data1)
    m2_tna <- tna::tna(data2)
    set.seed(1)
    tna_perm <- tna::permutation_test(m1_tna, m2_tna, iter = iter)

    n1 <- build_network(data1, method = "relative")
    n2 <- build_network(data2, method = "relative")
    nest_perm <- permutation_test(n1, n2, iter = iter, seed = 1)

    tna_diff <- tna_perm$edges$diffs_true
    states <- rownames(nest_perm$diff)
    n <- length(states)
    lbl <- sprintf("cfg%d(seed=%d)", i, cfg$seed)

    # Per-edge: true differences
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        our_val <- nest_perm$diff[ii, jj]
        tna_val <- tna_diff[ii, jj]
        if (is.na(our_val)) our_val <- 0
        if (is.na(tna_val)) tna_val <- 0
        expect_equal(our_val, tna_val, tolerance = TOL,
                     label = sprintf("perm diff %s [%s -> %s]",
                                     lbl, states[ii], states[jj]))
      })
    }))

    delta <- abs(replace(nest_perm$diff, is.na(nest_perm$diff), 0) -
                 replace(tna_diff, is.na(tna_diff), 0))
    report$log("permutation_vs_tna", lbl, n * n, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "tna::permutation_test()")
  }))

  report$write_csv("bootstrap")
})
