# Extracted from test-mlvar.R:509

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- simulate_data("mlvar", seed = 1)
B <- attr(d, "true_temporal")
vars <- attr(d, "vars")
expect_equal(dim(B), c(length(vars), length(vars)))
expect_true(all(diag(B) > 0))
expect_equal(rownames(B), vars)
