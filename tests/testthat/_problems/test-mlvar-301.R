# Extracted from test-mlvar.R:301

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
cors <- vapply(1:5, function(seed) {
    d <- simulate_data("mlvar", seed = seed)
    fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
                 day = "day", beep = "beep")
    true_B <- attr(d, "true_temporal")
    cor(as.vector(fit$temporal), as.vector(true_B))
  }, numeric(1L))
