# Tests for magnitude_difference(): structure, metrics, scalings, S3 methods.
# Compares the frequency (FTNA) and probability (TNA) views of a transition
# network and quantifies the per-edge discrepancy.

test_that("magnitude_difference returns the documented structure", {
  data(group_regulation_long, package = "Nestimate")
  fit <- magnitude_difference(group_regulation_long,
                              actor = "Actor", action = "Action",
                              time = "Time")

  expect_true(inherits(fit, "magnitude_difference"),
              info = "class is magnitude_difference")
  expect_named(fit, c("edges", "metric", "scale", "weights_ftna",
                      "weights_tna", "states"))
  expect_equal(fit$metric, "abs_diff")
  expect_equal(fit$scale, "tna_range")

  # n_states^2 edges (full dense expand.grid over the state set).
  n <- length(fit$states)
  expect_equal(nrow(fit$edges), n * n)
  expect_named(fit$edges, c("from", "to", "ftna", "tna", "signed", "value"))

  # signed is exactly tna - ftna; value (abs_diff) is |tna - ftna| == |signed|.
  expect_equal(fit$edges$signed, fit$edges$tna - fit$edges$ftna)
  expect_equal(fit$edges$value, abs(fit$edges$signed))
  expect_false(anyNA(fit$edges$value))
})

test_that("tna_range leaves TNA untouched and rescales FTNA into its range", {
  data(group_regulation_long, package = "Nestimate")
  net_t <- build_network(group_regulation_long, actor = "Actor",
                         action = "Action", method = "relative")
  fit <- magnitude_difference(group_regulation_long, actor = "Actor",
                              action = "Action", scale = "tna_range")

  # TNA matrix is passed through unchanged.
  expect_equal(unname(fit$weights_tna), unname(net_t$weights))
  # Rescaled FTNA shares TNA's [min, max] range.
  expect_equal(range(fit$weights_ftna), range(net_t$weights))
})

test_that("all five metrics compute without error and are non-negative", {
  data(group_regulation_long, package = "Nestimate")
  metrics <- c("abs_diff", "chord_dist", "atanh_diff",
               "geom_norm_diff", "cv_inflation")
  ok <- vapply(metrics, function(m) {
    f <- magnitude_difference(group_regulation_long, actor = "Actor",
                              action = "Action", metric = m)
    f$metric == m && all(f$edges$value >= 0) && !anyNA(f$edges$value)
  }, logical(1))
  expect_true(all(ok))
})

test_that("all four scalings compute without error", {
  data(group_regulation_long, package = "Nestimate")
  scales_arg <- c("tna_range", "rank_minmax", "minmax", "none")
  ok <- vapply(scales_arg, function(s) {
    f <- magnitude_difference(group_regulation_long, actor = "Actor",
                              action = "Action", scale = s)
    inherits(f, "magnitude_difference") && f$scale == s
  }, logical(1))
  expect_true(all(ok))
})

test_that("invalid arguments error cleanly", {
  data(group_regulation_long, package = "Nestimate")
  expect_error(magnitude_difference(group_regulation_long, metric = "nope"))
  expect_error(magnitude_difference(group_regulation_long, scale = "nope"))
  expect_error(magnitude_difference(list(1, 2, 3)))
})

test_that("print returns its argument invisibly", {
  data(group_regulation_long, package = "Nestimate")
  fit <- magnitude_difference(group_regulation_long, actor = "Actor",
                              action = "Action")
  expect_output(print(fit), "magnitude_difference object")
  expect_identical(withVisible(print(fit))$visible, FALSE)
})

test_that("both plot types return a warning-free ggplot", {
  skip_if_not_installed("ggplot2")
  data(group_regulation_long, package = "Nestimate")
  fit <- magnitude_difference(group_regulation_long, actor = "Actor",
                              action = "Action")
  expect_warning(p_stacked <- plot(fit), NA)
  p_circular <- plot(fit, type = "circular")
  expect_s3_class(p_stacked, "ggplot")
  expect_s3_class(p_circular, "ggplot")
})
