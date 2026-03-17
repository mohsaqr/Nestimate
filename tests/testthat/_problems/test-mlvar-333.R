# Extracted from test-mlvar.R:333

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- simulate_data("mlvar", seed = 3)
