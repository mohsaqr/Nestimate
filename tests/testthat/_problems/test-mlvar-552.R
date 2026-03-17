# Extracted from test-mlvar.R:552

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- simulate_data("mlvar", seed = 100, n_subjects = 40, d = 4,
                     n_days = 5, beeps_per_day = 10)
