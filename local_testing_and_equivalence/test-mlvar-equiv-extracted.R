# Equivalence test_that() blocks extracted from
# tests/testthat/test-mlvar.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ---- build_mlvar() -------------------------------------------------------
# Tests focus on:
#   1. S3 interface (class = c("net_mlvar","netobject_group"), aliases,
#      print/summary dispatch, coefs() accessor)
#   2. Constituent netobjects (fit$temporal, fit$contemporaneous, fit$between
#      are each c("netobject","cograph_network"))
#   3. Input validation
#   4. Numerical equivalence against `mlVAR::mlVAR()` on simulated data
#
# Equivalence target: bit-identical temporal B, contemporaneous <2e-15,
# between <1e-10 across 20 simulated configurations (seeds 201-220).
# Full real-ESM validation (25 datasets) lives in tmp/mlvar_equivalence_real20.R
# and is not run under R CMD check because it depends on external data.
#
# Plotting is cograph's responsibility. Nestimate never calls cograph and
# does not define a plot method for net_mlvar; each constituent is a
# standard cograph_network that cograph::splot dispatches on directly.

test_that("mlVAR equivalence - temporal/contemporaneous/between", {
  skip_if_not_installed("mlVAR")
  skip_if_not_installed("lme4")
  skip_if_not_installed("corpcor")

  # Run 5 simulated configurations — each fit takes ~0.2s so the whole
  # block is well under a second. The full 20-seed validation runs in
  # tmp/mlvar_equivalence_20seeds.R.
  seeds <- c(201, 205, 210, 215, 220)
  for (s in seeds) {
    d <- simulate_data("mlvar", seed = s)
    vars <- attr(d, "vars")

    ours <- build_mlvar(d, vars = vars,
                        id = "id", day = "day", beep = "beep")
    ref <- suppressWarnings(mlVAR::mlVAR(
      d, vars = vars, idvar = "id", dayvar = "day", beepvar = "beep",
      lags = 1, estimator = "lmer",
      temporal = "fixed", contemporaneous = "fixed",
      scale = FALSE, verbose = FALSE
    ))

    ref_B     <- ref$results$Beta$mean[, , 1]
    ref_theta <- ref$results$Theta$pcor$mean
    diag(ref_theta) <- 0
    ref_betw  <- ref$results$Omega_mu$pcor$mean
    diag(ref_betw)  <- 0

    strip <- function(m) { dimnames(m) <- NULL; m }
    expect_equal(max(abs(strip(ours$temporal$weights)        - strip(ref_B))),
                 0, tolerance = 1e-12,
                 label = sprintf("seed %d temporal", s))
    expect_equal(max(abs(strip(ours$contemporaneous$weights) - strip(ref_theta))),
                 0, tolerance = 1e-12,
                 label = sprintf("seed %d contemporaneous", s))
    expect_equal(max(abs(strip(ours$between$weights)         - strip(ref_betw))),
                 0, tolerance = 1e-10,
                 label = sprintf("seed %d between", s))
  }
})


# ---- Convergence diagnostics ----

test_that("build_mlvar numerical equivalence preserved with convergence checks", {
  # Verify that the convergence checking code does not change numerical output
  skip_if_not_installed("mlVAR")
  skip_if_not_installed("lme4")
  skip_if_not_installed("corpcor")

  # Use seed 210 — known to work from the 20-seed validation suite
  d <- simulate_data("mlvar", seed = 210)
  vars <- attr(d, "vars")

  ours <- suppressWarnings(build_mlvar(
    d, vars = vars, id = "id", day = "day", beep = "beep"
  ))
  ref  <- suppressWarnings(mlVAR::mlVAR(
    d, vars = vars, idvar = "id", dayvar = "day", beepvar = "beep",
    lags = 1, estimator = "lmer",
    temporal = "fixed", contemporaneous = "fixed",
    scale = FALSE, verbose = FALSE
  ))

  ref_B     <- ref$results$Beta$mean[, , 1]
  ref_theta <- ref$results$Theta$pcor$mean
  ref_betw  <- ref$results$Omega_mu$pcor$mean
  diag(ref_theta) <- 0; diag(ref_betw) <- 0

  strip <- function(m) { dimnames(m) <- NULL; m }
  expect_equal(max(abs(strip(ours$temporal$weights) - strip(ref_B))),
               0, tolerance = 1e-12)
  expect_equal(max(abs(strip(ours$contemporaneous$weights) - strip(ref_theta))),
               0, tolerance = 1e-12)
  expect_equal(max(abs(strip(ours$between$weights) - strip(ref_betw))),
               0, tolerance = 1e-10)
})