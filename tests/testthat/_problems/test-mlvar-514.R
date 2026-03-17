# Extracted from test-mlvar.R:514

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d1 <- simulate_data("mlvar", seed = 42)
