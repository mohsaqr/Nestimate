# ---- Equivalence: build_mlvar() vs mlVAR::mlVAR(estimator="lmer", temporal="fixed", contemporaneous="fixed") ----
#
# Validates bit-for-bit equivalence across 20+ random VAR-1 panel configurations on three
# networks returned by build_mlvar(): temporal B (lagged fixed-effect coefficients),
# contemporaneous partial correlations of lmer residuals, and between-person partial
# correlations derived via corpcor::pseudoinverse + cor2pcor. Four byte-level subtleties
# documented in CLAUDE.md are implicitly exercised here: (1) random intercept ONLY
# (`(1 | id)`), not random slopes, for temporal="fixed"; (2) cor2pcor applied directly
# to the residual correlation for contemporaneous (no EBIC-glasso); (3) mlVAR's
# forcePositive scalar-recycling quirk (a 0.001 bump injected into every off-diagonal)
# is replicated byte-for-byte; (4) `.mlvar_aveLag()` uses `rep(NA, lag)` (logical)
# rather than `rep(NA_real_, lag)` so integer columns stay integer and base R's
# `mean()` skips its two-pass pass for integer input — matching mlVAR's code path.
# Reference: mlVAR::mlVAR(..., estimator = "lmer", temporal = "fixed",
# contemporaneous = "fixed", scale = FALSE). Prior validation in tmp/ scripts reached
# ~1e-15 max delta; TOL here is 1e-8 with effective machine precision expected.

TOL <- 1e-8

strip_dimnames <- function(m) {
  dimnames(m) <- NULL
  m
}

# ---- 20-config equivalence sweep ----

test_that("build_mlvar() temporal/contemporaneous/between match mlVAR::mlVAR() (20 configs)", {
  skip_equiv_tests()
  skip_if_not_installed("mlVAR")
  skip_if_not_installed("lme4")
  skip_if_not_installed("corpcor")

  rep <- equiv_report()

  configs <- expand.grid(
    n_subjects = c(15L, 20L, 25L),
    d          = c(3L, 4L, 5L),
    n_obs      = c(30L, 40L),
    stringsAsFactors = FALSE
  )
  # Take 20 configs via seeds
  seeds <- seq(301L, 320L)
  configs <- configs[rep_len(seq_len(nrow(configs)), length(seeds)), , drop = FALSE]
  configs$seed <- seeds

  results <- lapply(seq_len(nrow(configs)), function(i) {
    cfg    <- configs[i, , drop = FALSE]
    seed   <- cfg$seed
    ns     <- cfg$n_subjects
    d      <- cfg$d
    no     <- cfg$n_obs
    cfg_id <- sprintf("seed=%d n_subjects=%d d=%d n_obs=%d", seed, ns, d, no)

    d_df <- simulate_data("mlvar",
                          seed       = seed,
                          n_subjects = ns,
                          d          = d,
                          n_obs      = no)
    vars <- attr(d_df, "vars")

    ours <- suppressMessages(suppressWarnings(
      build_mlvar(
        d_df,
        vars        = vars,
        id          = "id",
        day         = "day",
        beep        = "beep",
        standardize = FALSE
      )
    ))

    ref <- suppressWarnings(mlVAR::mlVAR(
      d_df,
      vars            = vars,
      idvar           = "id",
      dayvar          = "day",
      beepvar         = "beep",
      lags            = 1,
      estimator       = "lmer",
      temporal        = "fixed",
      contemporaneous = "fixed",
      scale           = FALSE,
      verbose         = FALSE
    ))

    ref_B       <- ref$results$Beta$mean[, , 1]
    ref_theta   <- ref$results$Theta$pcor$mean;    diag(ref_theta)   <- 0
    ref_between <- ref$results$Omega_mu$pcor$mean; diag(ref_between) <- 0

    our_B       <- strip_dimnames(ours$temporal$weights)
    our_theta   <- strip_dimnames(ours$contemporaneous$weights)
    our_between <- strip_dimnames(ours$between$weights)

    nets <- list(
      temporal        = list(o = our_B,       r = strip_dimnames(ref_B)),
      contemporaneous = list(o = our_theta,   r = strip_dimnames(ref_theta)),
      between         = list(o = our_between, r = strip_dimnames(ref_between))
    )

    per_net <- lapply(names(nets), function(nm) {
      diffs  <- abs(nets[[nm]]$o - nets[[nm]]$r)
      n_fail <- sum(diffs > TOL)
      rep$log(
        func           = paste0("build_mlvar:", nm),
        config         = cfg_id,
        n_checked      = length(diffs),
        n_failed       = n_fail,
        max_abs_err    = max(diffs),
        mean_abs_err   = mean(diffs),
        median_abs_err = stats::median(diffs),
        p95_abs_err    = stats::quantile(diffs, 0.95),
        reference      = "mlVAR::mlVAR(estimator='lmer', temporal='fixed', contemporaneous='fixed', scale=FALSE)",
        notes          = switch(
          nm,
          temporal        = "Beta$mean[,,1]",
          contemporaneous = "Theta$pcor$mean (diag zeroed)",
          between         = "Omega_mu$pcor$mean (diag zeroed)"
        )
      )
      list(network = nm, n_fail = n_fail, max_delta = max(diffs))
    })

    total_fail <- sum(vapply(per_net, `[[`, numeric(1), "n_fail"))
    max_delta  <- max(vapply(per_net, `[[`, numeric(1), "max_delta"))
    expect_true(
      total_fail == 0L,
      info = sprintf("%s: total_fail=%d max_delta=%.2e", cfg_id, total_fail, max_delta)
    )
    list(cfg_id = cfg_id, total_fail = total_fail, max_delta = max_delta)
  })

  rep$write_csv("mlvar")
  message(rep$summary())

  all_ok      <- sum(vapply(results, `[[`, numeric(1), "total_fail")) == 0L
  worst_delta <- max(vapply(results, `[[`, numeric(1), "max_delta"))
  expect_true(
    all_ok,
    info = sprintf("Some configs failed; worst max_delta across all nets = %.2e",
                   worst_delta)
  )
})
