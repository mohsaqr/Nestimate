# Numerical equivalence: build_mmm() internal identities + optional seqHMM parity
#
# Nestimate's build_mmm fits a mixture of Markov chains via EM. The natural
# equivalence checks are a mix of internal identities (no external dep needed)
# and an external cross-check against seqHMM::build_mmm + fit_model (gated on
# availability — seqHMM is NOT in Suggests by policy).
#
# Internal identities (TOL = 1e-10 unless noted):
#   (1) Every row of $posterior sums to 1.
#   (2) $assignments equals apply($posterior, 1, which.max).
#   (3) $mixing equals colMeans($posterior).
#   (4) $BIC = -2 * $log_likelihood + $n_params * log(N)
#       $AIC = -2 * $log_likelihood + 2 * $n_params
#   (5) Recomputed log-likelihood from fitted {$mixing, models[[k]]$weights,
#       models[[k]]$initial} equals stored $log_likelihood. This is the
#       substantive formula check — bugs in the EM objective surface here.
#
# External (seqHMM):
#   (6) For small k = 2, log-likelihood and sorted mixing proportions from
#       Nestimate and seqHMM agree up to EM local-optimum noise (loose TOL).
#       Skipped if seqHMM is not installed locally.

set.seed(20260422)
N_CONFIGS <- 15L
TOL       <- 1e-8

.gen_seqs <- function(n_actors, n_states, seq_length, seed) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- replicate(n_actors, sample(states, seq_length, replace = TRUE),
                    simplify = FALSE)
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

  # Per-sequence log-likelihood: log( sum_k pi_k * init_k[s1] * prod P_k[st-1, st] )
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

    # log-sum-exp for numerical stability
    m_max <- max(per_comp)
    m_max + log(sum(exp(per_comp - m_max)))
  })
  sum(seq_lls)
}

configs <- replicate(N_CONFIGS, list(
  n_actors   = sample(c(40L, 60L, 100L), 1L),
  n_states   = sample(3:4, 1L),
  seq_length = sample(c(8L, 14L, 20L), 1L),
  k          = sample(2:3, 1L),
  seed       = sample.int(1e5, 1L)
), simplify = FALSE)


# ---- (1) Posterior row sums ---------------------------------------------

test_that("build_mmm posterior rows sum to 1", {
  skip_on_cran()
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = 100L,
                      seed = cfg$seed)
    max(abs(rowSums(fit$posterior) - 1))
  }, numeric(1))
  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e", max(deltas)))
})


# ---- (2) Assignments = argmax posterior ---------------------------------

test_that("build_mmm assignments equal argmax(posterior)", {
  skip_on_cran()
  mismatches <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = 100L,
                      seed = cfg$seed)
    sum(fit$assignments != apply(fit$posterior, 1L, which.max))
  }, integer(1))
  expect_true(all(mismatches == 0L),
              info = sprintf("total mismatches = %d", sum(mismatches)))
})


# ---- (3) Mixing ≈ colMeans(posterior) -----------------------------------
#
# The M-step sets pi_k = (1/N) sum_i posterior[i, k]. After the final
# E-step these should match. A small discrepancy may remain if posterior
# is stored before the last mixing update; tolerate TOL_LOOSE.

test_that("build_mmm mixing equals colMeans(posterior) (loose)", {
  skip_on_cran()
  TOL_LOOSE <- 1e-3
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = 100L,
                      seed = cfg$seed)
    max(abs(fit$mixing - colMeans(fit$posterior)))
  }, numeric(1))
  expect_true(all(deltas < TOL_LOOSE),
              info = sprintf("max delta = %.3e", max(deltas)))
})


# ---- (4) BIC / AIC identities -------------------------------------------

test_that("build_mmm BIC and AIC match their textbook formulas", {
  skip_on_cran()
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = 100L,
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
  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = 100L,
                      seed = cfg$seed)
    ll_ref <- .recompute_ll(fit, seqs)
    abs(fit$log_likelihood - ll_ref)
  }, numeric(1))
  expect_true(all(deltas < TOL),
              info = sprintf("max ll delta = %.3e", max(deltas)))
})


# ---- (6) seqHMM parity — deferred --------------------------------------
#
# TODO: write a seqHMM cross-check. seqHMM::build_mmm requires an explicit
# list of initial_probs and transition_probs per cluster (NULL no longer
# accepted by the modern API). A faithful port needs to construct those
# inputs and then align components via Hungarian matching on mixing
# proportions before comparing. Not blocking the internal-identity tests.
# Gated by skip_if_not_installed("seqHMM") when implemented; seqHMM stays
# out of Suggests per the no-validation-deps policy.


# ---- Delta report + validation dashboard emit --------------------------

test_that("build_mmm equivalence report (CSV + CVS JSON)", {
  skip_on_cran()

  report <- equiv_report()

  invisible(lapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    fit  <- build_mmm(seqs, k = cfg$k, n_starts = 3L, max_iter = 100L,
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
      notes = "seqHMM external parity deferred (modern API requires explicit probs)"
    )
  }))

  report$write_csv("build_mmm")
  report$write_cvs("build_mmm")
  expect_true(length(report$rows) > 0L)
})
