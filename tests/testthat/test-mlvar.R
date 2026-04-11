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

test_that("build_mlvar returns a net_mlvar / netobject_group", {
  d <- simulate_data("mlvar", seed = 1)
  vars <- attr(d, "vars")

  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  expect_true(inherits(fit, "net_mlvar"))
  expect_true(inherits(fit, "netobject_group"))
  expect_false(inherits(fit, "cograph_network"))  # the GROUP is not itself a network
  expect_equal(names(fit), c("temporal", "contemporaneous", "between"))
  expect_length(fit, 3L)
})

test_that("each constituent is a full netobject/cograph_network", {
  d <- simulate_data("mlvar", seed = 1)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  for (nm in c("temporal", "contemporaneous", "between")) {
    net <- fit[[nm]]
    expect_true(inherits(net, "netobject"),
                label = paste0("fit$", nm, " should be a netobject"))
    expect_true(inherits(net, "cograph_network"),
                label = paste0("fit$", nm, " should be a cograph_network"))
    expect_true(is.matrix(net$weights),
                label = paste0("fit$", nm, "$weights should be a matrix"))
    expect_equal(dim(net$weights), c(length(vars), length(vars)),
                 label = paste0("fit$", nm, "$weights shape"))
    expect_s3_class(net$nodes, "data.frame")
    expect_s3_class(net$edges, "data.frame")
  }

  expect_true(fit$temporal$directed)
  expect_false(fit$contemporaneous$directed)
  expect_false(fit$between$directed)

  expect_equal(fit$temporal$method,        "mlvar_temporal")
  expect_equal(fit$contemporaneous$method, "mlvar_contemporaneous")
  expect_equal(fit$between$method,         "mlvar_between")
})

test_that("model metadata is stored as attributes", {
  d <- simulate_data("mlvar", seed = 1)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  expect_type(attr(fit, "n_obs"), "integer")
  expect_type(attr(fit, "n_subjects"), "integer")
  expect_type(attr(fit, "lag"), "integer")
  expect_type(attr(fit, "standardize"), "logical")
  expect_s3_class(attr(fit, "coefs"), "data.frame")
})

test_that("coefs(fit) returns a tidy data.frame", {
  d <- simulate_data("mlvar", seed = 6)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  cf <- coefs(fit)
  expect_s3_class(cf, "data.frame")
  expected_cols <- c("outcome", "predictor", "beta", "se", "t", "p",
                     "ci_lower", "ci_upper", "significant")
  expect_true(all(expected_cols %in% names(cf)))
  # d outcomes x d predictors = d^2 rows
  expect_equal(nrow(cf), length(vars)^2)
  # Every (outcome, predictor) pair appears exactly once
  expect_equal(nrow(unique(cf[, c("outcome", "predictor")])),
               length(vars)^2)
  expect_type(cf$significant, "logical")
})

test_that("coefs() errors cleanly on unsupported classes", {
  expect_error(coefs(1L),           "No coefs")
  expect_error(coefs(data.frame()), "No coefs")
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

test_that("print.net_mlvar shows the three-network structure", {
  d <- simulate_data("mlvar", seed = 4)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  joined <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(joined, "\\$temporal")
  expect_match(joined, "\\$contemporaneous")
  expect_match(joined, "\\$between")
  expect_match(joined,
               sprintf("%d x %d directed", length(vars), length(vars)))
  expect_match(joined, "edges significant at p<0.05")
})

test_that("summary.net_mlvar reports B, Theta, and between blocks", {
  d <- simulate_data("mlvar", seed = 5)
  vars <- attr(d, "vars")
  fit <- build_mlvar(d, vars = vars, id = "id", day = "day", beep = "beep")

  joined <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(joined, "Temporal Network")
  expect_match(joined, "Contemporaneous Network")
  expect_match(joined, "Between-Subjects Network")
})

# ---- Numerical equivalence against mlVAR --------------------------------

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
