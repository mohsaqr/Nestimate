# Numerical equivalence: build_mmm() internal identities + seqHMM parity + speed
#
# Nestimate's build_mmm fits a mixture of Markov chains via EM. This file does
# three things:
#
#   (A) Internal identities — formula-level checks needing no external dep.
#       (1) posterior rows sum to 1
#       (2) assignments == argmax(posterior)
#       (3) mixing ~= colMeans(posterior)
#       (4) BIC and AIC match textbook formulas
#       (5) recomputed log-likelihood from {mixing, P_k, init_k} == stored ll
#
#   (B) External parity vs seqHMM — gated on `skip_if_not_installed("seqHMM")`.
#       Per package policy, seqHMM is NOT in Suggests; install locally.
#       The core difficulty: EM is non-convex, label-switching is unidentified,
#       and seqHMM's modern API requires explicit per-cluster transition_probs
#       and initial_probs. Strategy:
#         - run both packages with multi-start (N_STARTS), keep the best
#           log-likelihood on each side
#         - if the two best log-likelihoods agree (absolute or relative), align
#           components via Hungarian-equivalent permutation enumeration
#           (closed form for k <= 4) on mixing proportions
#         - then compare aligned mixing, transition matrices, initial probs
#       Tolerances are loose because EM converges to numerically distinct but
#       structurally identical optima.
#
#   (C) Speed benchmark — single-start wall-clock of Nestimate vs seqHMM at
#       the same (n, p, T, k). Reports per-config ratio + aggregate stats.
#       Logged via equiv_report so runs land in the validation dashboard.

set.seed(20260426)
N_CONFIGS    <- 100L
N_STARTS     <- 5L
MAX_ITER     <- 200L
TOL          <- 1e-8
TOL_LL_ABS   <- 1.0       # absolute log-likelihood tolerance for "same optimum"
TOL_LL_REL   <- 5e-3      # relative log-likelihood tolerance
TOL_PARAM    <- 5e-2      # max abs delta on aligned transitions/initials
PARITY_PASS_RATE <- 0.90  # fraction of configs that must hit the same optimum


.gen_seqs <- function(n_actors, n_states, seq_length, seed) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- replicate(n_actors, sample(states, seq_length, replace = TRUE),
                    simplify = FALSE)
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_length))
  df
}

# Generate sequences from a true K-component MMM with separated transition
# matrices. Used by the parity test so that EM has an identifiable optimum
# (uniform-random noise gives a flat likelihood surface where many parameter
# settings tie at the same logLik, defeating point comparison).
.gen_seqs_mixture <- function(n_actors, n_states, seq_length, k, seed,
                              concentration = 4) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]

  # Generate k well-separated transition matrices via Dirichlet(concentration)
  # — higher concentration on diagonal-shifted alpha gives distinct dynamics.
  trans_list <- lapply(seq_len(k), function(c_) {
    M <- matrix(0, n_states, n_states)
    for (i in seq_len(n_states)) {
      alpha <- rep(1, n_states)
      # Shift mass to a cluster-specific "preferred next state"
      pref <- ((i + c_ - 2L) %% n_states) + 1L
      alpha[pref] <- concentration
      g <- stats::rgamma(n_states, alpha)
      M[i, ] <- g / sum(g)
    }
    M
  })
  init_list <- lapply(seq_len(k), function(c_) {
    g <- stats::rgamma(n_states, 1)
    g / sum(g)
  })
  mixing <- {
    g <- stats::rgamma(k, 4)
    g / sum(g)
  }

  z <- sample.int(k, n_actors, replace = TRUE, prob = mixing)
  rows <- lapply(seq_len(n_actors), function(i) {
    cluster <- z[i]
    seq_idx <- integer(seq_length)
    seq_idx[1L] <- sample.int(n_states, 1L, prob = init_list[[cluster]])
    if (seq_length >= 2L) {
      for (t in 2:seq_length) {
        seq_idx[t] <- sample.int(n_states, 1L,
                                 prob = trans_list[[cluster]][seq_idx[t - 1L], ])
      }
    }
    states[seq_idx]
  })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_length))
  df
}

.recompute_ll <- function(fit, data) {
  states    <- fit$states
  n_states  <- length(states)
  mixing    <- fit$mixing
  k         <- fit$k

  inits <- vapply(fit$models, function(m) m$initial[states], numeric(n_states))
  Ps    <- lapply(fit$models, function(m) m$weights[states, states])

  seq_mat <- as.matrix(data)
  seq_lls <- apply(seq_mat, 1L, function(seq_row) {
    seq_chr <- as.character(seq_row)
    seq_chr <- seq_chr[!is.na(seq_chr) & seq_chr != ""]
    if (length(seq_chr) < 1L) return(NA_real_)
    idx <- match(seq_chr, states)

    per_comp <- vapply(seq_len(k), function(c_) {
      lp <- log(inits[idx[1L], c_])
      if (length(idx) >= 2L) {
        from_i <- idx[-length(idx)]
        to_i   <- idx[-1L]
        lp <- lp + sum(log(Ps[[c_]][cbind(from_i, to_i)]))
      }
      log(mixing[c_]) + lp
    }, numeric(1))

    m_max <- max(per_comp)
    m_max + log(sum(exp(per_comp - m_max)))
  })
  sum(seq_lls)
}

# Brute-force best permutation by sum of squared mixing distance.
# k = 2..4 covered; for k > 4 fall back to identity (rare in our config grid).
.all_perms <- function(k) {
  if (k <= 1L) return(list(seq_len(k)))
  if (k == 2L) return(list(1:2, 2:1))
  base <- .all_perms(k - 1L)
  out <- list()
  for (p in base) {
    for (pos in seq_len(k)) {
      out[[length(out) + 1L]] <- append(p, k, pos - 1L)
    }
  }
  out
}

.best_permutation <- function(mix_a, mix_b) {
  k <- length(mix_a)
  if (k > 4L) return(seq_len(k))
  perms <- .all_perms(k)
  errs <- vapply(perms, function(p) sum((mix_a - mix_b[p])^2), numeric(1))
  perms[[which.min(errs)]]
}

# Pick permutation minimizing total parameter distance (mixing + transitions
# + initials). Robust to near-symmetric mixing where mixing-only matching
# would pick the wrong cluster pairing.
.best_full_permutation <- function(nx, sx) {
  k <- length(nx$mixing)
  if (k > 4L) return(seq_len(k))
  perms <- .all_perms(k)
  errs <- vapply(perms, function(p) {
    e_mix   <- sum(abs(nx$mixing - sx$mixing[p]))
    e_trans <- sum(vapply(seq_len(k), function(c_) {
      sum(abs(nx$transitions[[c_]] - sx$transitions[[p[c_]]]))
    }, numeric(1)))
    e_init  <- sum(vapply(seq_len(k), function(c_) {
      sum(abs(nx$initials[[c_]] - sx$initials[[p[c_]]]))
    }, numeric(1)))
    e_mix + e_trans + e_init
  }, numeric(1))
  perms[[which.min(errs)]]
}

# Run seqHMM::build_mmm + fit_model with `n_starts` random restarts; keep best.
# Returns NULL on total failure.
.fit_seqhmm_best <- function(seqs_df, alphabet, k, n_starts, max_iter, tol, seed) {
  n_states <- length(alphabet)
  seq_obj <- TraMineR::seqdef(seqs_df, alphabet = alphabet, void = "%",
                              nr = "*", silent = TRUE)
  set.seed(seed)
  start_seeds <- sample.int(.Machine$integer.max, n_starts)

  fits <- lapply(start_seeds, function(s) {
    set.seed(s)
    init_list <- replicate(k, {
      v <- runif(n_states); v / sum(v)
    }, simplify = FALSE)
    trans_list <- replicate(k, {
      M <- matrix(runif(n_states^2), n_states, n_states)
      M / rowSums(M)
    }, simplify = FALSE)
    init_list  <- lapply(init_list,  function(v) { names(v) <- alphabet; v })
    trans_list <- lapply(trans_list, function(M) { dimnames(M) <- list(alphabet, alphabet); M })

    mmm <- tryCatch(
      seqHMM::build_mmm(seq_obj, n_clusters = k,
                        transition_probs = trans_list,
                        initial_probs = init_list),
      error = function(e) NULL
    )
    if (is.null(mmm)) return(NULL)
    tryCatch(
      suppressMessages(seqHMM::fit_model(
        mmm, em_step = TRUE, global_step = FALSE, local_step = FALSE,
        control_em = list(maxeval = max_iter, reltol = tol),
        threads = 1L
      )),
      error = function(e) NULL
    )
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0L) return(NULL)
  lls <- vapply(fits, function(f) f$logLik, numeric(1))
  fits[[which.max(lls)]]
}

# Extract aligned (Nestimate-style) mixing / transitions / initials from a
# fitted seqHMM model.
.seqhmm_extract <- function(seqhmm_fit, alphabet) {
  mod <- seqhmm_fit$model
  coef <- mod$coefficients               # (n_covariates+1) x n_clusters
  log_pi <- coef[1L, ]                   # row 1 = intercept
  pi <- exp(log_pi) / sum(exp(log_pi))
  names(pi) <- mod$cluster_names

  trans <- lapply(mod$transition_probs, function(M) {
    M2 <- M[alphabet, alphabet]
    M2
  })
  names(trans) <- mod$cluster_names

  inits <- lapply(mod$initial_probs, function(v) v[alphabet])
  names(inits) <- mod$cluster_names

  list(mixing = unname(pi), transitions = unname(trans),
       initials = unname(inits), logLik = as.numeric(seqhmm_fit$logLik))
}

# Same shape from Nestimate fit
.nestimate_extract <- function(nfit, alphabet) {
  mixing <- unname(nfit$mixing)
  trans <- lapply(nfit$models, function(m) m$weights[alphabet, alphabet])
  inits <- lapply(nfit$models, function(m) unname(m$initial[alphabet]))
  list(mixing = mixing, transitions = trans, initials = inits,
       logLik = nfit$log_likelihood)
}

configs <- replicate(N_CONFIGS, list(
  n_actors   = sample(c(40L, 60L, 100L), 1L),
  n_states   = sample(3:4, 1L),
  seq_length = sample(c(8L, 14L, 20L), 1L),
  k          = sample(2:3, 1L),
  seed       = sample.int(1e5, 1L)
), simplify = FALSE)


# ---- (1) Posterior row sums ---------------------------------------------

test_that("build_mmm posterior rows sum to 1 (100 configs)", {
  skip_on_cran()
  skip_equiv_tests()
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = MAX_ITER,
                      seed = cfg$seed)
    max(abs(rowSums(fit$posterior) - 1))
  }, numeric(1))
  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e", max(deltas)))
})


# ---- (2) Assignments = argmax posterior ---------------------------------

test_that("build_mmm assignments equal argmax(posterior)", {
  skip_on_cran()
  skip_equiv_tests()
  mismatches <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = MAX_ITER,
                      seed = cfg$seed)
    sum(fit$assignments != apply(fit$posterior, 1L, which.max))
  }, integer(1))
  expect_true(all(mismatches == 0L),
              info = sprintf("total mismatches = %d", sum(mismatches)))
})


# ---- (3) Mixing ~= colMeans(posterior) ----------------------------------

test_that("build_mmm mixing equals colMeans(posterior) (loose)", {
  skip_on_cran()
  skip_equiv_tests()
  TOL_LOOSE <- 1e-3
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = MAX_ITER,
                      seed = cfg$seed)
    max(abs(fit$mixing - colMeans(fit$posterior)))
  }, numeric(1))
  expect_true(all(deltas < TOL_LOOSE),
              info = sprintf("max delta = %.3e", max(deltas)))
})


# ---- (4) BIC / AIC identities -------------------------------------------

test_that("build_mmm BIC and AIC match their textbook formulas", {
  skip_on_cran()
  skip_equiv_tests()
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = MAX_ITER,
                      seed = cfg$seed)
    N <- nrow(seqs)
    bic_ref <- -2 * fit$log_likelihood + fit$n_params * log(N)
    aic_ref <- -2 * fit$log_likelihood + 2 * fit$n_params
    max(abs(fit$BIC - bic_ref), abs(fit$AIC - aic_ref))
  }, numeric(1))
  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e", max(deltas)))
})


# ---- (5) Recomputed log-likelihood matches stored ll --------------------

test_that("build_mmm stored log_likelihood equals recomputation from parameters", {
  skip_on_cran()
  skip_equiv_tests()
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = MAX_ITER,
                      seed = cfg$seed)
    ll_ref <- .recompute_ll(fit, seqs)
    abs(fit$log_likelihood - ll_ref)
  }, numeric(1))
  expect_true(all(deltas < TOL),
              info = sprintf("max ll delta = %.3e", max(deltas)))
})


# ---- (6) seqHMM parity --------------------------------------------------
#
# For each of N_CONFIGS, fit Nestimate and seqHMM with N_STARTS random
# restarts, take best logLik on each side, then:
#   (a) check whether logLik values agree (absolute OR relative)
#   (b) if yes, align components via permutation enumeration on mixing,
#       and compare aligned transitions/initials.
# Pass criterion: at least PARITY_PASS_RATE of configs converge to the
# same logLik basin.

test_that("build_mmm parity vs seqHMM (multi-start, Hungarian alignment)", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("TraMineR")

  agree <- logical(N_CONFIGS)
  ll_deltas <- numeric(N_CONFIGS)
  param_deltas <- rep(NA_real_, N_CONFIGS)

  for (i in seq_along(configs)) {
    cfg <- configs[[i]]
    seqs <- .gen_seqs_mixture(cfg$n_actors, cfg$n_states, cfg$seq_length,
                              cfg$k, cfg$seed)
    alphabet <- LETTERS[seq_len(cfg$n_states)]

    nfit <- tryCatch(
      build_mmm(seqs, k = cfg$k, n_starts = N_STARTS, max_iter = MAX_ITER,
                tol = 1e-08, smooth = 0, seed = cfg$seed),
      error = function(e) NULL
    )
    sfit <- .fit_seqhmm_best(seqs, alphabet, cfg$k, N_STARTS, MAX_ITER,
                              tol = 1e-08, seed = cfg$seed + 1L)

    if (is.null(nfit) || is.null(sfit)) {
      agree[i] <- FALSE
      ll_deltas[i] <- NA_real_
      next
    }

    nx <- .nestimate_extract(nfit, alphabet)
    sx <- .seqhmm_extract(sfit, alphabet)

    ll_diff <- abs(nx$logLik - sx$logLik)
    ll_rel  <- ll_diff / max(abs(nx$logLik), abs(sx$logLik))
    ll_deltas[i] <- ll_diff
    same_optimum <- (ll_diff < TOL_LL_ABS) || (ll_rel < TOL_LL_REL)
    agree[i] <- same_optimum

    if (same_optimum) {
      perm <- .best_full_permutation(nx, sx)
      mix_d <- max(abs(nx$mixing - sx$mixing[perm]))
      trans_d <- max(vapply(seq_along(perm), function(c_) {
        max(abs(nx$transitions[[c_]] - sx$transitions[[perm[c_]]]))
      }, numeric(1)))
      init_d <- max(vapply(seq_along(perm), function(c_) {
        max(abs(nx$initials[[c_]] - sx$initials[[perm[c_]]]))
      }, numeric(1)))
      param_deltas[i] <- max(mix_d, trans_d, init_d)
    }
  }

  pass_rate <- mean(agree)
  agreed_idx <- which(agree)

  message(sprintf(
    "seqHMM parity: %d/%d configs (%.1f%%) hit same optimum; ll-delta median %.2e p95 %.2e; param-delta median %.2e p95 %.2e",
    sum(agree), N_CONFIGS, 100 * pass_rate,
    stats::median(ll_deltas[agreed_idx]),
    stats::quantile(ll_deltas[agreed_idx], 0.95, names = FALSE),
    stats::median(param_deltas[agreed_idx], na.rm = TRUE),
    stats::quantile(param_deltas[agreed_idx], 0.95, names = FALSE, na.rm = TRUE)
  ))

  expect_gte(pass_rate, PARITY_PASS_RATE)
  if (length(agreed_idx) > 0L) {
    expect_lt(stats::quantile(param_deltas[agreed_idx], 0.95,
                              names = FALSE, na.rm = TRUE),
              TOL_PARAM)
  }
})


# ---- (7) Speed benchmark vs seqHMM --------------------------------------
#
# Fair single-start comparison: n_starts = 1 on Nestimate (so its parallel
# mclapply does nothing), 1 init on seqHMM. Warm both packages once on a
# tiny problem before timing.

test_that("build_mmm speed vs seqHMM (single-start, wall-clock)", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("TraMineR")

  # Warm up
  warm <- .gen_seqs(20L, 3L, 8L, 1L)
  invisible(build_mmm(warm, k = 2L, n_starts = 1L, max_iter = 50L, seed = 1L))
  invisible(.fit_seqhmm_best(warm, LETTERS[1:3], 2L, 1L, 50L, 1e-08, 1L))

  bench_idx <- sample.int(N_CONFIGS, min(40L, N_CONFIGS))
  speedups <- numeric(length(bench_idx))
  t_n <- numeric(length(bench_idx))
  t_s <- numeric(length(bench_idx))

  for (j in seq_along(bench_idx)) {
    cfg <- configs[[bench_idx[j]]]
    seqs <- .gen_seqs_mixture(cfg$n_actors, cfg$n_states, cfg$seq_length,
                              cfg$k, cfg$seed)
    alphabet <- LETTERS[seq_len(cfg$n_states)]

    t0 <- Sys.time()
    invisible(build_mmm(seqs, k = cfg$k, n_starts = 1L, max_iter = MAX_ITER,
                        tol = 1e-08, smooth = 0, seed = cfg$seed))
    t_n[j] <- as.numeric(Sys.time() - t0, units = "secs")

    t0 <- Sys.time()
    invisible(.fit_seqhmm_best(seqs, alphabet, cfg$k, 1L, MAX_ITER, 1e-08,
                               cfg$seed + 1L))
    t_s[j] <- as.numeric(Sys.time() - t0, units = "secs")

    speedups[j] <- t_s[j] / t_n[j]
  }

  median_speedup <- stats::median(speedups)
  message(sprintf(
    "Speed benchmark (%d configs, single start each): Nestimate median %.3fs, seqHMM median %.3fs, speedup median %.2fx geomean %.2fx",
    length(bench_idx),
    stats::median(t_n), stats::median(t_s),
    median_speedup,
    exp(mean(log(speedups)))
  ))

  # Soft requirement: Nestimate is at least as fast on the median.
  expect_gt(median_speedup, 1.0)
})


# ---- (8) Delta report + validation dashboard emit ----------------------

test_that("build_mmm equivalence report (CSV + CVS JSON)", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()
  has_seqhmm <- requireNamespace("seqHMM", quietly = TRUE) &&
                requireNamespace("TraMineR", quietly = TRUE)

  invisible(lapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = MAX_ITER,
                      seed = cfg$seed)

    ll_err  <- abs(fit$log_likelihood - .recompute_ll(fit, seqs))
    bic_err <- abs(fit$BIC - (-2 * fit$log_likelihood +
                              fit$n_params * log(nrow(seqs))))
    post_err <- max(abs(rowSums(fit$posterior) - 1))
    errs <- c(ll_err = ll_err, bic_err = bic_err, posterior_err = post_err)

    report$log(
      func = "build_mmm (internal identities)",
      config = sprintf("n_actors=%d n_states=%d seq_len=%d k=%d",
                       cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$k),
      n_checked = length(errs),
      n_failed  = sum(errs >= TOL),
      max_abs_err    = max(errs),
      mean_abs_err   = mean(errs),
      median_abs_err = stats::median(errs),
      p95_abs_err    = stats::quantile(errs, 0.95, names = FALSE),
      reference = "recomputed ll from {mixing, P_k, init_k}; BIC formula; posterior row-sum",
      notes = "internal-identities sweep; 100 configs"
    )

    if (has_seqhmm) {
      alphabet <- LETTERS[seq_len(cfg$n_states)]
      seqs_mix <- .gen_seqs_mixture(cfg$n_actors, cfg$n_states,
                                    cfg$seq_length, cfg$k, cfg$seed)
      sfit <- .fit_seqhmm_best(seqs_mix, alphabet, cfg$k, N_STARTS, MAX_ITER,
                                1e-08, cfg$seed + 1L)
      nfit_ms <- build_mmm(seqs_mix, k = cfg$k, n_starts = N_STARTS,
                           max_iter = MAX_ITER, tol = 1e-08, smooth = 0,
                           seed = cfg$seed)
      if (!is.null(sfit)) {
        nx <- .nestimate_extract(nfit_ms, alphabet)
        sx <- .seqhmm_extract(sfit, alphabet)
        ll_diff <- abs(nx$logLik - sx$logLik)
        ll_rel  <- ll_diff / max(abs(nx$logLik), abs(sx$logLik))
        same_optimum <- (ll_diff < TOL_LL_ABS) || (ll_rel < TOL_LL_REL)
        if (same_optimum) {
          perm <- .best_full_permutation(nx, sx)
          mix_d <- max(abs(nx$mixing - sx$mixing[perm]))
          trans_d <- max(vapply(seq_along(perm), function(c_) {
            max(abs(nx$transitions[[c_]] - sx$transitions[[perm[c_]]]))
          }, numeric(1)))
          init_d <- max(vapply(seq_along(perm), function(c_) {
            max(abs(nx$initials[[c_]] - sx$initials[[perm[c_]]]))
          }, numeric(1)))
          errs2 <- c(ll = ll_diff, mixing = mix_d, transitions = trans_d,
                     initials = init_d)
          notes2 <- "same EM optimum; aligned via Hungarian permutation"
          n_failed2 <- sum(errs2 >= TOL_PARAM)
        } else {
          errs2 <- c(ll = ll_diff)
          notes2 <- sprintf("different optima (ll_rel=%.2e); not aligned",
                            ll_rel)
          n_failed2 <- 1L
        }

        report$log(
          func = "build_mmm vs seqHMM (parity)",
          config = sprintf("n_actors=%d n_states=%d seq_len=%d k=%d n_starts=%d",
                           cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$k,
                           N_STARTS),
          n_checked = length(errs2),
          n_failed  = n_failed2,
          max_abs_err    = max(errs2),
          mean_abs_err   = mean(errs2),
          median_abs_err = stats::median(errs2),
          p95_abs_err    = stats::quantile(errs2, 0.95, names = FALSE),
          reference = "seqHMM::build_mmm + fit_model (em_step), N_STARTS restarts",
          notes = notes2
        )
      }
    }
  }))

  report$write_csv("build_mmm")
  report$write_cvs("build_mmm")
  expect_true(length(report$rows) > 0L)
})
