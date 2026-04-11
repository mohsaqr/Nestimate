# ---- build_mlvar() -------------------------------------------------------
# Tests focus on:
#   1. S3 interface (class, aliases, print/summary/plot dispatch)
#   2. Input validation
#   3. Numerical equivalence against `mlVAR::mlVAR()` on simulated data
#
# Equivalence target: bit-identical temporal B, contemporaneous <2e-15,
# between <1e-10 across 20 simulated configurations (seeds 201-220).
# Full real-ESM validation (25 datasets) lives in tmp/mlvar_equivalence_real20.R
# and is not run under R CMD check because it depends on external data.

test_that("build_mlvar returns a net_mlvar / cograph_network object", {
  d <- simulate_data("mlvar", seed = 1)
  vars <- attr(d, "vars")

  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  expect_true(inherits(fit, "net_mlvar"))
  expect_true(inherits(fit, "cograph_network"))
  expect_true(all(c("temporal", "contemporaneous", "between",
                    "coefs", "n_obs", "n_subjects", "lag",
                    "weights", "nodes", "edges", "directed") %in% names(fit)))
  expect_equal(dim(fit$temporal),        c(length(vars), length(vars)))
  expect_equal(dim(fit$contemporaneous), c(length(vars), length(vars)))
  expect_equal(dim(fit$between),         c(length(vars), length(vars)))
  # Default cograph view is the temporal network
  expect_identical(fit$weights, fit$temporal)
  expect_true(fit$directed)
})

test_that("mlvar() is an alias for build_mlvar()", {
  d <- simulate_data("mlvar", seed = 2)
  vars <- attr(d, "vars")

  fit1 <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")
  fit2 <- mlvar(d,       vars = vars, id = "id", day = "day", beep = "beep")

  expect_identical(fit1, fit2)
  expect_identical(body(build_mlvar), body(mlvar))
})

test_that("build_mlvar validates required arguments", {
  d <- simulate_data("mlvar", seed = 3)
  vars <- attr(d, "vars")

  expect_error(build_mlvar(list(), vars = vars, id = "id"),
               "is.data.frame")
  expect_error(build_mlvar(d, vars = "only_one", id = "id"),
               "length")
  expect_error(build_mlvar(d, vars = vars, id = c("a", "b")),
               "length")
  expect_error(build_mlvar(d, vars = vars, id = "nonexistent_column"),
               "not found in data")
  expect_error(build_mlvar(d, vars = vars, id = "id", lag = 0),
               "lag")
})

test_that("print.net_mlvar prints matrix dimensions honestly", {
  d <- simulate_data("mlvar", seed = 4)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  joined <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(joined,
               sprintf("%d x %d directed", length(vars), length(vars)))
  expect_match(joined, "Contemporaneous network")
  expect_match(joined, "Between network")
  # Significant-edge report uses the tidy coefs
  expect_match(joined, "edges significant at p<0.05")
})

test_that("summary.net_mlvar reports B, Theta, and between blocks", {
  d <- simulate_data("mlvar", seed = 5)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  out <- capture.output(summary(fit))
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "Temporal Network")
  expect_match(joined, "Contemporaneous Network")
  expect_match(joined, "Between-Subjects Network")
})

test_that("coefs is a tidy data.frame with one row per (outcome, predictor)", {
  d <- simulate_data("mlvar", seed = 6)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  expect_s3_class(fit$coefs, "data.frame")
  expected_cols <- c("outcome", "predictor", "beta", "se", "t", "p",
                     "ci_lower", "ci_upper", "significant")
  expect_true(all(expected_cols %in% names(fit$coefs)))
  # d outcomes x d predictors = d^2 rows
  expect_equal(nrow(fit$coefs), length(vars)^2)
  # Every (outcome, predictor) pair appears exactly once
  expect_equal(nrow(unique(fit$coefs[, c("outcome", "predictor")])),
               length(vars)^2)
  # `significant` is a proper logical
  expect_type(fit$coefs$significant, "logical")
})

# ---- Numerical equivalence against mlVAR --------------------------------

test_that("mlVAR equivalence — temporal/contemporaneous/between", {
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
    expect_equal(max(abs(strip(ours$temporal)        - strip(ref_B))),
                 0, tolerance = 1e-12,
                 label = sprintf("seed %d temporal", s))
    expect_equal(max(abs(strip(ours$contemporaneous) - strip(ref_theta))),
                 0, tolerance = 1e-12,
                 label = sprintf("seed %d contemporaneous", s))
    expect_equal(max(abs(strip(ours$between)         - strip(ref_betw))),
                 0, tolerance = 1e-10,
                 label = sprintf("seed %d between", s))
  }
})

# Plotting is cograph's responsibility. Nestimate never calls cograph and
# does not define a plot method for `net_mlvar`. The cograph-side splot
# dispatch and per-type styling spec lives in cograph's CLAUDE.md
# ("Nestimate net_mlvar" section) — equivalence and output structure are
# tested here; rendering is validated on the cograph side.
