# Extracted from test-mlvar.R:498

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- simulate_data("mlvar", seed = 1)
expect_equal(attr(d, "type"), "mlvar")
expect_true(!is.null(attr(d, "vars")))
expect_true(!is.null(attr(d, "true_temporal")))
expect_true(!is.null(attr(d, "true_contemporaneous")))
