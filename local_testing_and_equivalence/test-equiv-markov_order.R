# Equivalence tests for markov_order_test()
#
# Layer 1 (deterministic, <1e-10): G^2 matches from-scratch manual formula.
# Layer 2 (deterministic, <1e-10): Pearson chi-square variant matches
#   markovchain::verifyMarkovProperty() byte-for-byte (statistic, dof, p.value).
# Layer 3 (Monte Carlo): permutation p-value is uniform under true null
#   (Kolmogorov-Smirnov test vs U(0,1) over 150 replicates).


# ---------------------------------------------------------------------------
# Manual G^2 (from-scratch)
# ---------------------------------------------------------------------------

.mo_manual_g2 <- function(seqs, k) {
  tuples <- do.call(rbind, lapply(seqs, function(traj) {
    n <- length(traj)
    if (n < k + 1L) return(NULL)
    pos <- seq_len(n - k)
    data.frame(
      x = traj[pos],
      w = if (k == 1L) rep("", length(pos)) else
          vapply(pos, function(p) paste(traj[(p + 1L):(p + k - 1L)],
                                        collapse = "\x01"),
                 character(1L)),
      s = traj[pos + k],
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(tuples) || nrow(tuples) == 0L) return(c(stat = 0, df = 0))
  total_g2 <- 0; total_df <- 0L
  for (wv in unique(tuples$w)) {
    sub <- tuples[tuples$w == wv, , drop = FALSE]
    tab <- table(sub$x, sub$s)
    if (nrow(tab) < 2L || ncol(tab) < 2L) next
    rs <- rowSums(tab); cs <- colSums(tab); tot <- sum(tab)
    expt <- outer(rs, cs) / tot
    nz <- tab > 0 & expt > 0
    total_g2 <- total_g2 + 2 * sum(tab[nz] * log(tab[nz] / expt[nz]))
    total_df <- total_df + (nrow(tab) - 1L) * (ncol(tab) - 1L)
  }
  c(stat = total_g2, df = total_df)
}


# ---------------------------------------------------------------------------
# verifyMarkovProperty() byte-for-byte replica
#
# Replicates markovchain::verifyMarkovProperty() exactly:
#   - transMatrix = first-order MLE on whole sequence
#   - subSample trimmed to length divisible by 3
#   - statistic = sum_{i,j,k} (N_ijk - N_ij*P(k|j))^2 / (N_ij*P(k|j))
#   - dof = #unique triples - #unique doubles + #unique states - 1
# ---------------------------------------------------------------------------

.mo_verify_replica <- function(sequence) {
  # MLE first-order transition matrix (matches markovchainFit$estimate)
  states <- sort(unique(sequence))
  n <- length(sequence)
  src <- sequence[-n]; dst <- sequence[-1]
  cnt <- table(factor(src, levels = states), factor(dst, levels = states))
  rs <- rowSums(cnt)
  P <- cnt / ifelse(rs == 0, 1, rs)

  # subSample trimmed to length %% 3 == 0
  sub <- sequence[1:(n - (n %% 3L))]
  # Triples (stride 1, NOT 3 — matches the package's sliding window)
  m <- length(sub)
  i <- sub[1:(m - 2L)]
  j <- sub[2:(m - 1L)]
  kk <- sub[3:m]

  # Nijk via aggregate-equivalent
  keys_ijk <- paste(i, j, kk, sep = "\x02")
  tab_ijk <- table(keys_ijk)
  split_ijk <- strsplit(names(tab_ijk), "\x02", fixed = TRUE)
  Nijk <- data.frame(
    V1 = vapply(split_ijk, `[`, character(1L), 1L),
    V2 = vapply(split_ijk, `[`, character(1L), 2L),
    V3 = vapply(split_ijk, `[`, character(1L), 3L),
    n  = as.integer(tab_ijk),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Nij
  keys_ij <- paste(i, j, sep = "\x02")
  tab_ij <- table(keys_ij)
  split_ij <- strsplit(names(tab_ij), "\x02", fixed = TRUE)
  Nij <- data.frame(
    V1 = vapply(split_ij, `[`, character(1L), 1L),
    V2 = vapply(split_ij, `[`, character(1L), 2L),
    n  = as.integer(tab_ij),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Per-triple contribution: (Nijk - Nij * P(k|j))^2 / (Nij * P(k|j))
  term <- vapply(seq_len(nrow(Nijk)), function(idx) {
    i_ <- Nijk$V1[idx]; j_ <- Nijk$V2[idx]; k_ <- Nijk$V3[idx]
    nij <- Nij$n[Nij$V1 == i_ & Nij$V2 == j_]
    pjk <- P[j_, k_]
    expected <- nij * pjk
    if (expected <= 0) 0 else (Nijk$n[idx] - expected)^2 / expected
  }, numeric(1L))
  statistic <- sum(term)

  # DOF (verifyMarkovProperty's formula — operates on the ORIGINAL sequence)
  doubles <- paste(sequence[-n], sequence[-1L], sep = "")
  triples <- paste(sequence[1:(n - 2L)], sequence[2:(n - 1L)],
                   sequence[3:n], sep = "")
  dof <- length(unique(triples)) - length(unique(doubles)) +
         length(unique(sequence)) - 1L

  list(statistic = statistic, dof = as.numeric(dof),
       p.value   = 1 - stats::pchisq(statistic, df = dof))
}


# ---------------------------------------------------------------------------
# TEST 1: G^2 = manual formula, byte-for-byte
# ---------------------------------------------------------------------------

test_that("G^2 statistic and df match manual formula to machine precision", {
  skip_if_not_installed("Nestimate")
  skip_equiv_tests()

  states <- letters[1:5]
  configs <- expand.grid(n_seqs = c(10L, 30L),
                         len    = c(50L, 150L),
                         seed   = 1:5)

  results <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
    cfg <- configs[i, ]
    set.seed(cfg$seed)
    tm <- matrix(runif(5 * 5), 5, 5, dimnames = list(states, states))
    tm <- tm / rowSums(tm)
    seqs <- lapply(seq_len(cfg$n_seqs), function(.) {
      s <- character(cfg$len); s[1] <- sample(states, 1L)
      for (ii in 2:cfg$len) s[ii] <- sample(states, 1L, prob = tm[s[ii - 1L], ])
      s
    })

    res <- markov_order_test(seqs, max_order = 3L, n_perm = 1L, seed = cfg$seed)

    do.call(rbind, lapply(1:3, function(k) {
      ours <- c(stat = res$test_table$g2[k + 1L],
                df   = res$test_table$df[k + 1L])
      theirs <- .mo_manual_g2(seqs, k)
      data.frame(seed = cfg$seed, n_seqs = cfg$n_seqs, len = cfg$len, k = k,
                 delta_stat = abs(ours["stat"] - theirs["stat"]),
                 delta_df   = abs(ours["df"] - theirs["df"]))
    }))
  }))

  expect_lt(max(results$delta_stat), 1e-10)
  expect_equal(max(results$delta_df), 0)
})


# ---------------------------------------------------------------------------
# TEST 2: Pearson chi-square replica matches markovchain byte-for-byte
# ---------------------------------------------------------------------------

test_that("Pearson chi-square replica matches markovchain::verifyMarkovProperty byte-for-byte", {
  skip_if_not_installed("Nestimate")
  skip_if_not_installed("markovchain")
  skip_equiv_tests()

  states <- letters[1:4]

  results <- do.call(rbind, lapply(1:30, function(seed) {
    set.seed(seed)
    tm <- matrix(runif(16), 4, 4, dimnames = list(states, states))
    tm <- tm / rowSums(tm)
    n <- 600L
    s <- character(n); s[1] <- sample(states, 1L)
    for (i in 2:n) s[i] <- sample(states, 1L, prob = tm[s[i - 1L], ])

    ours  <- .mo_verify_replica(s)
    theirs <- suppressWarnings(markovchain::verifyMarkovProperty(s, verbose = FALSE))

    data.frame(
      seed          = seed,
      delta_stat    = abs(ours$statistic - theirs$statistic),
      delta_dof     = abs(ours$dof - theirs$dof),
      delta_pvalue  = abs(ours$p.value - theirs$p.value)
    )
  }))

  expect_lt(max(results$delta_stat),   1e-10)
  expect_equal(max(results$delta_dof), 0)
  expect_lt(max(results$delta_pvalue), 1e-10)
})


# ---------------------------------------------------------------------------
# TEST 3: permutation p-value is uniform under true null
# ---------------------------------------------------------------------------

test_that("permutation p-value distribution matches U(0,1) under true H0 (KS test)", {
  skip_if_not_installed("Nestimate")
  skip_equiv_tests()

  states <- letters[1:4]
  set.seed(2026)
  tm <- matrix(runif(16), 4, 4, dimnames = list(states, states))
  tm <- tm / rowSums(tm)

  # 150 independent pure order-1 datasets; test order-2-vs-1
  pvals <- vapply(1:150, function(seed) {
    set.seed(seed)
    seqs <- lapply(1:8, function(.) {
      s <- character(200); s[1] <- sample(states, 1L)
      for (i in 2:200) s[i] <- sample(states, 1L, prob = tm[s[i - 1L], ])
      s
    })
    res <- markov_order_test(seqs, max_order = 2L, n_perm = 200L, seed = seed)
    res$test_table$p_permutation[3L]
  }, numeric(1L))

  # Kolmogorov-Smirnov: p-values should come from U(0,1)
  ks <- suppressWarnings(stats::ks.test(pvals, "punif"))
  expect_gt(ks$p.value, 0.05)

  # Mean close to 0.5
  expect_gt(mean(pvals), 0.40)
  expect_lt(mean(pvals), 0.60)

  # Type-I rate at alpha=0.05 close to nominal
  rej_rate <- mean(pvals < 0.05)
  expect_gte(rej_rate, 0.01)
  expect_lte(rej_rate, 0.10)
})
