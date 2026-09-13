# outcome_model() -- unit-level outcome regression with honest inference.

make_binary <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(
    hint     = rbinom(n, 1L, 0.4),
    think    = rbinom(n, 1L, 0.3),
    noise    = rbinom(n, 1L, 0.5),
    constant = 1L,
    n_events = rpois(n, 20L),
    actor    = rep(letters[1:10], length.out = n)
  )
  d$success <- rbinom(n, 1L, stats::plogis(-0.5 + 1.2 * d$hint))
  d
}

test_that("a two-valued outcome is modelled as binomial and gives odds ratios", {
  fit <- outcome_model(make_binary(), outcome = "success",
                       predictors = c("hint", "think"))
  expect_s3_class(fit, "net_outcome_model")
  expect_equal(fit$family, "binomial")
  e <- effects_table(fit)
  expect_true(all(c("odds_ratio", "or_lower", "or_upper") %in% names(e)))
  # the planted effect is positive and its interval excludes 1
  expect_gt(e$odds_ratio[e$term == "hint"], 1)
  expect_gt(e$or_lower[e$term == "hint"], 1)
})

test_that("a numeric outcome is modelled as gaussian, with no odds ratios", {
  d <- make_binary()
  d$score <- 2 * d$hint + rnorm(nrow(d))
  fit <- outcome_model(d, outcome = "score", predictors = c("hint", "think"))
  expect_equal(fit$family, "gaussian")
  expect_false("odds_ratio" %in% names(effects_table(fit)))
  expect_equal(effects_table(fit)$estimate[1], 2, tolerance = 0.2)
})

test_that("p-values carry the multiplicity correction and intercept is excluded", {
  fit <- outcome_model(make_binary(), outcome = "success",
                       predictors = c("hint", "think", "noise"))
  e <- effects_table(fit, intercept = TRUE)
  expect_true(is.na(e$p_adj[e$term == "(Intercept)"]))
  term_rows <- e[e$term != "(Intercept)", ]
  expect_true(all(term_rows$p_adj >= term_rows$p_value - 1e-12))
  expect_equal(fit$correction, "BH")
})

test_that("constant predictors are dropped rather than breaking the fit", {
  fit <- outcome_model(make_binary(), outcome = "success",
                       predictors = c("hint", "constant"))
  expect_equal(fit$dropped, "constant")
  expect_false("constant" %in% effects_table(fit)$term)
})

test_that("select = 'split' fits on rows the selection did not see", {
  d <- make_binary(n = 600L)
  fit <- outcome_model(d, outcome = "success",
                       predictors = c("hint", "think", "noise"),
                       select = "split", n_select = 2L, seed = 42L)
  expect_equal(length(fit$selected), 2L)
  expect_lt(fit$n, nrow(d))              # only the held-out half is fitted
  expect_equal(fit$n, nrow(d) - floor(nrow(d) / 2))
  # same seed, same split
  again <- outcome_model(d, outcome = "success",
                         predictors = c("hint", "think", "noise"),
                         select = "split", n_select = 2L, seed = 42L)
  expect_equal(fit$selected, again$selected)
})

test_that("effects_table filters by significance and rounds all but p-values", {
  fit <- outcome_model(make_binary(), outcome = "success",
                       predictors = c("hint", "think", "noise"))
  sig <- effects_table(fit, significant = TRUE, alpha = 0.05)
  expect_true(all(sig$p_adj < 0.05))
  e <- effects_table(fit, digits = 2)
  expect_equal(e$estimate, round(e$estimate, 2))
  expect_false(identical(e$p_value, round(e$p_value, 2)))
})

test_that("summary returns the tidy table and plot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- outcome_model(make_binary(), outcome = "success",
                       predictors = c("hint", "think"))
  expect_s3_class(summary(fit), "data.frame")
  expect_s3_class(plot(fit), "ggplot")
  expect_silent(ggplot2::ggplot_build(plot(fit)))
})

test_that("degenerate and malformed input raise classed errors", {
  d <- make_binary()
  expect_error(outcome_model(d, "success", "nope"), "not found")
  expect_error(outcome_model(d, "success", "constant"),
               class = "nestimate_no_variance")
  empty <- d[0, ]
  expect_error(outcome_model(empty, "success", "hint"),
               class = "nestimate_empty_model_frame")
  expect_error(outcome_model(d, "success", "hint", ci_level = 2), "`ci_level`")
})
