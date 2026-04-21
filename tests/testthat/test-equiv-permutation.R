# ---- Equivalence: permutation() against independent references ----
#
# Validates Nestimate::permutation() numerics against two independent
# references:
#   (1) p.adjust correspondence: the adjusted p-value matrix returned by
#       permutation(..., adjust = method) must equal
#       stats::p.adjust(raw_p, method) applied to the raw p-values from
#       a twin call with adjust = "none". Verified across 8 adjust methods
#       x 12 random configurations.
#   (2) Hand-coded permutation reference: with a fixed RNG seed, we
#       reimplement the unpaired permutation loop in pure base R and
#       confirm per-edge p-values, observed diffs, and perm_sd match.
#       Uses method = "relative" (simplest estimator).
# TOL = 1e-8 for floats; exact integer match for counts.

TOL <- 1e-8

# -- Reference 1: p.adjust correspondence (8 methods x 12 configs) --

test_that("permutation() adjusted p-values match stats::p.adjust (8 methods x 12 configs)", {
  skip_equiv_tests()

  rep <- equiv_report()

  adjust_methods <- c("holm", "hochberg", "hommel", "bonferroni",
                      "BH", "BY", "fdr", "none")
  seeds <- seq_len(12L)

  configs <- expand.grid(seed = seeds, adjust = adjust_methods,
                         stringsAsFactors = FALSE)

  res <- lapply(seq_len(nrow(configs)), function(i) {
    seed_i <- configs$seed[i]
    adjust_i <- configs$adjust[i]

    seq1 <- simulate_sequences(n_actors = 8L, n_states = 4L,
                               seq_length = 15L, seed = seed_i * 7L)
    seq2 <- simulate_sequences(n_actors = 8L, n_states = 4L,
                               seq_length = 15L, seed = seed_i * 7L + 3L)

    net1 <- build_network(seq1, method = "relative")
    net2 <- build_network(seq2, method = "relative")

    # Twin calls with same seed → same raw permutation null
    perm_raw <- permutation(net1, net2, iter = 50L,
                            adjust = "none", seed = 123L + seed_i)
    perm_adj <- permutation(net1, net2, iter = 50L,
                            adjust = adjust_i, seed = 123L + seed_i)

    raw_p <- as.vector(perm_raw$p_values)
    ref_adj <- stats::p.adjust(raw_p, method = adjust_i)
    our_adj <- as.vector(perm_adj$p_values)

    diffs <- abs(our_adj - ref_adj)
    n_fail <- sum(diffs > TOL)

    rep$log(
      func           = "permutation",
      config         = sprintf("seed=%d adjust=%s", seed_i, adjust_i),
      n_checked      = length(diffs),
      n_failed       = n_fail,
      max_abs_err    = max(diffs),
      mean_abs_err   = mean(diffs),
      median_abs_err = stats::median(diffs),
      p95_abs_err    = stats::quantile(diffs, 0.95),
      reference      = "stats::p.adjust",
      notes          = "adjusted p_values vs p.adjust(raw_p, method)"
    )

    expect_true(
      n_fail == 0L,
      info = sprintf("seed=%d adjust=%s: max delta=%.2e",
                     seed_i, adjust_i, max(diffs))
    )

    n_fail
  })

  rep$write_csv("permutation_padjust")
  message(rep$summary())
})


# -- Reference 2: Hand-coded permutation reference --

# Helper: compute relative transition matrix from a wide data.frame
# of sequences using the same byrow/colSums logic the package uses.
.relative_from_wide <- function(df, states) {
  n <- length(states)
  nbins <- n * n
  # Per-row pair counts (t, t+1) encoded as (from-1)*n + to
  idx_list <- lapply(seq_len(nrow(df)), function(i) {
    row <- match(as.character(unlist(df[i, ])), states)
    from_ <- row[-length(row)]
    to_ <- row[-1L]
    ok <- !is.na(from_) & !is.na(to_)
    (from_[ok] - 1L) * n + to_[ok]
  })
  all_idx <- unlist(idx_list)
  counts <- tabulate(all_idx, nbins = nbins)
  # byrow layout: rows = from, cols = to
  mat <- matrix(counts, nrow = n, ncol = n, byrow = TRUE)
  rs <- rowSums(mat)
  nz <- rs > 0
  mat[nz, ] <- mat[nz, ] / rs[nz]
  dimnames(mat) <- list(states, states)
  mat
}

# Per-row count vector (length n*n) in byrow order used by pooled approach
.per_row_counts <- function(df, states) {
  n <- length(states)
  nbins <- n * n
  t(vapply(seq_len(nrow(df)), function(i) {
    row <- match(as.character(unlist(df[i, ])), states)
    from_ <- row[-length(row)]
    to_ <- row[-1L]
    ok <- !is.na(from_) & !is.na(to_)
    idx <- (from_[ok] - 1L) * n + to_[ok]
    tabulate(idx, nbins = nbins)
  }, numeric(nbins)))
}

test_that("permutation() matches hand-coded permutation reference (10 configs)", {
  skip_equiv_tests()

  rep <- equiv_report()

  seeds <- seq_len(10L)

  lapply(seeds, function(seed_i) {
    # small, deterministic data
    seq1 <- simulate_sequences(n_actors = 6L, n_states = 4L,
                               seq_length = 10L, seed = seed_i * 11L)
    seq2 <- simulate_sequences(n_actors = 6L, n_states = 4L,
                               seq_length = 10L, seed = seed_i * 11L + 5L)

    net1 <- build_network(seq1, method = "relative")
    net2 <- build_network(seq2, method = "relative")

    states <- net1$nodes$label
    n_nodes <- length(states)

    # Package result with known seed
    perm_obj <- permutation(net1, net2, iter = 100L,
                            adjust = "none", paired = FALSE,
                            seed = 2026L + seed_i)

    # Hand-coded reference using the SAME per-iteration RNG state.
    # Nestimate calls set.seed(seed) BEFORE the precompute step, then
    # the precompute itself does NOT consume RNG (it's deterministic
    # tabulate + indexing), then the permutation loop consumes RNG via
    # sample.int. We therefore call set.seed with the same value and
    # advance the RNG in lockstep.
    pooled_counts <- rbind(
      .per_row_counts(net1$data, states),
      .per_row_counts(net2$data, states)
    )
    n_x <- nrow(net1$data)
    n_y <- nrow(net2$data)
    n_total <- n_x + n_y
    nbins <- n_nodes * n_nodes
    iter <- 100L

    obs_x <- .relative_from_wide(net1$data, states)
    obs_y <- .relative_from_wide(net2$data, states)
    obs_diff <- obs_x - obs_y
    obs_flat <- as.vector(obs_diff)

    set.seed(2026L + seed_i)
    exceed_counts <- integer(nbins)
    sum_diffs <- numeric(nbins)
    sum_diffs_sq <- numeric(nbins)

    for (i in seq_len(iter)) {
      idx_x <- sample.int(n_total, n_x)
      cnt_x <- colSums(pooled_counts[idx_x, , drop = FALSE])
      cnt_y <- colSums(pooled_counts[-idx_x, , drop = FALSE])
      mx <- matrix(cnt_x, n_nodes, n_nodes, byrow = TRUE)
      my <- matrix(cnt_y, n_nodes, n_nodes, byrow = TRUE)
      rs_x <- rowSums(mx); nzx <- rs_x > 0
      rs_y <- rowSums(my); nzy <- rs_y > 0
      mx[nzx, ] <- mx[nzx, ] / rs_x[nzx]
      my[nzy, ] <- my[nzy, ] / rs_y[nzy]
      perm_diff <- as.vector(mx) - as.vector(my)
      exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
      sum_diffs <- sum_diffs + perm_diff
      sum_diffs_sq <- sum_diffs_sq + perm_diff^2
    }

    ref_p_flat <- (exceed_counts + 1L) / (iter + 1L)
    perm_mean <- sum_diffs / iter
    ref_sd <- sqrt(pmax(sum_diffs_sq / iter - perm_mean^2, 0))

    # Compare p-values
    our_p <- as.vector(perm_obj$p_values)
    diffs_p <- abs(our_p - ref_p_flat)

    # Compare observed diff
    our_diff <- as.vector(perm_obj$diff)
    diffs_d <- abs(our_diff - obs_flat)

    # Compare effect size implicitly via perm_sd:
    # es = obs_flat / perm_sd; if ref_sd==0, both set to 0
    our_es <- as.vector(perm_obj$effect_size)
    ref_es <- ifelse(ref_sd == 0, 0, obs_flat / ref_sd)
    diffs_es <- abs(our_es - ref_es)

    all_diffs <- c(diffs_p, diffs_d, diffs_es)
    n_checked <- length(all_diffs)
    n_fail <- sum(all_diffs > TOL)

    rep$log(
      func           = "permutation",
      config         = sprintf("seed=%d hand-coded", seed_i),
      n_checked      = n_checked,
      n_failed       = n_fail,
      max_abs_err    = max(all_diffs),
      mean_abs_err   = mean(all_diffs),
      median_abs_err = stats::median(all_diffs),
      p95_abs_err    = stats::quantile(all_diffs, 0.95),
      reference      = "hand-coded base-R permutation",
      notes          = "p_values, diff, effect_size combined"
    )

    expect_true(
      n_fail == 0L,
      info = sprintf("seed=%d: max delta=%.2e (p=%.2e d=%.2e es=%.2e)",
                     seed_i, max(all_diffs),
                     max(diffs_p), max(diffs_d), max(diffs_es))
    )
  })

  rep$write_csv("permutation_handcoded")
  message(rep$summary())
})
