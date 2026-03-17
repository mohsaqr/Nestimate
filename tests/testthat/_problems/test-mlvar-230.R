# Extracted from test-mlvar.R:230

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(42)
d <- simulate_data("mlvar", seed = 42)
