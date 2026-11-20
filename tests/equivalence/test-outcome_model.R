# outcome_model(), effects_table(), summary.net_outcome_model(): a regression of an outcome on predictors (logistic
# for a binary outcome, linear otherwise; with `group`, a random intercept per group). References: stats::glm()
# with Wald intervals (stats::confint.default()) and stats::p.adjust(method = "BH"), stats::lm(), and
# lme4::glmer(). New: these functions had no comparison with a reference.

set.seed(2)
data <- data.frame(y = stats::rbinom(150L, 1L, 0.4), x1 = stats::rnorm(150L), x2 = stats::rnorm(150L),
                   g = sample(c("a", "b", "c"), 150L, replace = TRUE))
data$z <- data$x1 + stats::rnorm(150L)

# Coefficient table of a fitted model without its intercept, in effects_table()'s columns.
sq_effects <- function(coefficients, intervals, log_scale) {
  table <- coefficients[-1L, , drop = FALSE]
  bounds <- intervals[-1L, , drop = FALSE]
  out <- data.frame(term = rownames(table), estimate = table[, 1L], std_error = table[, 2L], statistic = table[, 3L],
                    ci_lower = bounds[, 1L], ci_upper = bounds[, 2L], p_value = table[, 4L],
                    p_adj = stats::p.adjust(table[, 4L], method = "BH"))
  if (log_scale) {
    out$odds_ratio <- exp(out$estimate)
    out$or_lower <- exp(out$ci_lower)
    out$or_upper <- exp(out$ci_upper)
  }
  out
}

test_that("a logistic outcome model equals stats::glm() with Wald intervals and BH-adjusted p-values", {
  fit <- stats::glm(y ~ x1 + x2, family = stats::binomial, data = data)
  expected <- sq_effects(summary(fit)$coefficients, stats::confint.default(fit), log_scale = TRUE)
  found <- outcome_model(data, outcome = "y", predictors = c("x1", "x2"))
  # The same IRLS on both sides.
  expect_reference(effects_table(found, digits = 15), expected, reference = "stats",
                   covers = c("outcome_model()", "effects_table()"), quantity = "logistic regression", tolerance = 1e-10,
                   measure = "relative", key = "term")
  # summary() rounds every column but the p-values to three decimals: half a unit in the third decimal.
  expect_reference(summary(found), expected, reference = "stats", covers = c("outcome_model()", "summary.net_outcome_model()"),
                   quantity = "logistic regression, summary()", tolerance = 5e-4, key = "term")
})

test_that("a linear outcome model equals stats::lm()", {
  fit <- stats::lm(z ~ x2, data = data)
  expected <- sq_effects(summary(fit)$coefficients, stats::confint.default(fit), log_scale = FALSE)
  expect_reference(effects_table(outcome_model(data, outcome = "z", predictors = "x2"), digits = 15),
                   subset(expected, select = c(term, estimate, std_error, statistic, p_value)), reference = "stats",
                   covers = c("outcome_model()", "effects_table()"), quantity = "linear regression", tolerance = 1e-10,
                   measure = "relative", key = "term")
})

test_that("a grouped logistic outcome model equals lme4::glmer() with a random intercept", {
  skip_without_reference("lme4")
  fit <- lme4::glmer(y ~ x1 + x2 + (1 | g), family = stats::binomial, data = data)
  coefficients <- summary(fit)$coefficients
  expected <- data.frame(term = rownames(coefficients)[-1L], estimate = coefficients[-1L, 1L],
                         std_error = coefficients[-1L, 2L], p_value = coefficients[-1L, 4L])
  # The same lme4 fit on both sides (Laplace approximation, default optimiser).
  expect_reference(effects_table(outcome_model(data, outcome = "y", predictors = c("x1", "x2"), group = "g"), digits = 15),
                   expected, reference = "lme4", covers = c("outcome_model()", "effects_table()"),
                   quantity = "mixed logistic regression", tolerance = 1e-8, measure = "relative", key = "term")
})
