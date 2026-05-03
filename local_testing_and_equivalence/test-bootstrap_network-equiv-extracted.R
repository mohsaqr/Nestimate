# Equivalence test_that() blocks extracted from
# tests/testthat/test-bootstrap_network.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ---- bootstrap_network() Tests ----

# Helper: generate wide sequence data
.make_boot_wide <- function(n = 50, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# Helper: generate frequency-like data for association methods
.make_boot_assoc <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}


# ---- Basic functionality: transition methods ----

test_that("bootstrap_network dispatches over netobject_group", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  df$grp <- rep(c("A", "B"), length.out = nrow(df))
  group_net <- build_network(df, method = "relative", group = "grp")
  expect_s3_class(group_net, "netobject_group")

  results <- bootstrap_network(group_net, iter = 20L, seed = 1)
  expect_true(is.list(results))
  expect_equal(length(results), 2L)
  expect_s3_class(results[[1]], "net_bootstrap")
  expect_s3_class(results[[2]], "net_bootstrap")
})


# ---- cograph_network input (L96) ----

test_that("bootstrap_network matches tna::bootstrap numerically", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  iter <- 1000L
  seed <- 42

  # tna approach
  tna_model <- tna::tna(df)
  set.seed(seed)
  tna_boot <- tna::bootstrap(tna_model, iter = iter)

  # Nestimate approach
  nest_net <- build_network(df, method = "relative")
  nest_boot <- bootstrap_network(nest_net, iter = iter, seed = seed)

  # 1. Original weights must be identical
  expect_equal(nest_boot$original$weights, tna_model$weights,
               tolerance = 1e-12)

  # 2. Bootstrap means must match to machine precision
  expect_equal(nest_boot$mean, tna_boot$weights_mean, tolerance = 1e-10)

  # 3. Bootstrap SDs must match
  expect_equal(nest_boot$sd, tna_boot$weights_sd, tolerance = 1e-10)

  # 4. CI bounds must match
  expect_equal(nest_boot$ci_lower, tna_boot$ci_lower, tolerance = 1e-10)
  expect_equal(nest_boot$ci_upper, tna_boot$ci_upper, tolerance = 1e-10)

  # 5. CR bounds must match
  expect_equal(nest_boot$cr_lower, tna_boot$cr_lower, tolerance = 1e-10)
  expect_equal(nest_boot$cr_upper, tna_boot$cr_upper, tolerance = 1e-10)

  # 6. P-values: RNG consumption patterns may differ slightly between
  #    tna and Nestimate bootstrap implementations. With 1000 iterations,
  #    p-values should correlate > 0.9999 and max diff < 0.01.
  #    Zero-weight edges excluded (tna sets p=1, Nestimate computes actual).
  non_zero <- nest_boot$original$weights != 0
  expect_gt(cor(as.vector(nest_boot$p_values[non_zero]),
                as.vector(tna_boot$p_values[non_zero])), 0.9999)
  expect_lt(max(abs(nest_boot$p_values[non_zero] -
                    tna_boot$p_values[non_zero])), 0.01)

  # 7. Significant edges must match
  expect_equal(nest_boot$significant, tna_boot$weights_sig, tolerance = 1e-10)
})


# ---- Group and mixed dispatches ----

