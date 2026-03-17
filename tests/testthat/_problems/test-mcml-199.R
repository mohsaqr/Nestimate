# Extracted from test-mcml.R:199

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.make_test_matrix <- function(seed = 42) {
  sim <- simulate_mtna(n_nodes = 4, n_types = 3, seed = seed)
  sim
}
.make_wide_sequences <- function(seed = 1) {
  set.seed(seed)
  states <- c("plan", "monitor", "adapt", "discuss", "evaluate", "reflect")
  as.data.frame(matrix(
    sample(states, 500, replace = TRUE), nrow = 100, ncol = 5,
    dimnames = list(NULL, paste0("T", 1:5))
  ))
}
.test_clusters_6 <- list(
  Regulation = c("plan", "adapt"),
  Cognition = c("monitor", "evaluate"),
  Social = c("discuss", "reflect")
)

# test -------------------------------------------------------------------------
sim <- .make_test_matrix()
