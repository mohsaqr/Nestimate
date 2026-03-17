# Extracted from test-mlvar.R:337

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- simulate_data("mlvar", seed = 3)
fit <- mlvar(d, vars = attr(d, "vars"), id = "id",
               day = "day", beep = "beep")
expect_true(all(diag(fit$temporal) > 0))
