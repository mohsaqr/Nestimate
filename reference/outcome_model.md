# Model Unit-Level Outcomes from Sequence or Network Predictors

`outcome_model()` fits a regression of a unit-level outcome on predictor
columns (pattern indicators, topological features from
[`simplicial_features()`](https://pak.dynasite.org/Nestimate/reference/simplicial_features.md),
or any other numeric covariates) and returns a tidy effect table with
confidence intervals and multiplicity-corrected p-values. A two-valued
outcome gets a logistic regression and a numeric one a linear
regression; a `group` column adds a random intercept per group, and
`select = "split"` chooses the predictors on one half of the data and
fits them on the other.

Each effect is conditional on the other terms in the model; for a
two-valued outcome, the odds ratio is the factor by which the odds of
the outcome value coded 1 change per unit of the predictor. The p-values
are corrected across every term except the intercept, the `adjust`
covariates included. A predictor that counts events also measures how
long the unit was observed; enter that exposure through `adjust`.

## Usage

``` r
outcome_model(
  data,
  outcome,
  predictors,
  group = NULL,
  adjust = NULL,
  family = c("auto", "binomial", "gaussian"),
  select = c("none", "split"),
  n_select = 10L,
  correction = "BH",
  ci_level = 0.95,
  seed = NULL
)

# S3 method for class 'net_outcome_model'
print(x, ...)

# S3 method for class 'net_outcome_model'
summary(object, ...)

# S3 method for class 'net_outcome_model'
plot(x, ...)
```

## Arguments

- data:

  A data frame with one row per unit of analysis. Rows with a missing
  value in the outcome, a predictor, an `adjust` covariate or the
  `group` column are left out; when none is left, an error of class
  `nestimate_empty_model_frame` is raised. A column that is not in
  `data` raises an error.

- outcome:

  The outcome column. With `family = "auto"`, an outcome with two
  distinct values is modelled with the binomial family (the second value
  in sort order, level order for a factor, or `TRUE`, coded 1) and any
  other with the gaussian family, which needs a numeric outcome.

- predictors:

  The predictor columns. A predictor that is not numeric (a logical
  column included; convert it with
  [`as.integer()`](https://rdrr.io/r/base/integer.html)) or has no
  variance in the complete rows is dropped and listed in the result;
  when every predictor is dropped, an error of class
  `nestimate_no_variance` is raised.

- group:

  A grouping column, such as the actor the units are nested in. When
  lme4 is installed, a random intercept per group is added, the right
  treatment for units nested in actors. lme4 is a suggested package:
  without it the random intercept is dropped, a warning of class
  `nestimate_no_lme4` is raised, and a plain
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) is fitted. Default
  `NULL`.

- adjust:

  Covariates entered before the predictors. Use it for exposure: a unit
  observed longer contains more of every pattern, so an unadjusted
  effect can be volume in disguise. Default `NULL`.

- family:

  The model family: `"auto"` (default) reads it from the outcome,
  `"binomial"` and `"gaussian"` force one. `"binomial"` codes the
  second-smallest outcome value (or `TRUE`) as 1 and every other value
  as 0, so give it a two-valued outcome.

- select:

  How the predictors are chosen: `"none"` (default) fits all of them;
  `"split"` ranks them on a random half of the rows by the absolute z or
  t statistic of a one-predictor
  [`glm()`](https://rdrr.io/r/stats/glm.html) (without `adjust` or
  `group`), keeps the best `n_select`, and fits the reported model on
  the other half.

- n_select:

  The number of predictors kept when `select = "split"`. Default `10`.
  Ignored otherwise.

- correction:

  The multiplicity correction, one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html),
  passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"`.

- ci_level:

  The confidence level of the intervals, between 0 and 1. Default
  `0.95`.

- seed:

  An integer that makes the split reproducible; used only with
  `select = "split"`. A random-number state that existed before the call
  is restored afterwards. Default `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  `outcome_model()` result.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: not
  used.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): an
  `outcome_model()` result.

## Value

A `net_outcome_model` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the outcome and
family, the number of units (and groups), the number of selected terms,
the adjustments, the dropped predictors and the correction, followed by
the effect table, and returns its input invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`effects_table()`](https://pak.dynasite.org/Nestimate/reference/effects_table.md)
return the effect table without the intercept, rounded to three decimals
except the p-values;
[`effects_table()`](https://pak.dynasite.org/Nestimate/reference/effects_table.md)
also keeps the intercept, or only the significant terms, on request.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` forest plot of the effects (odds ratios for a binomial fit,
against a dashed line at no effect), with filled markers for terms whose
corrected p-value is below 0.05. The object's fields:

- effects:

  The effect table: a data frame with one row per model term, the
  intercept and the `adjust` covariates included, and columns `term`,
  `estimate`, `std_error`, `statistic`, `ci_lower`, `ci_upper`,
  `p_value` and `p_adj` (`NA` on the intercept row, which is left out of
  the correction), plus `odds_ratio`, `or_lower` and `or_upper` for a
  binomial fit.

- model:

  The fitted model: a `glm`, or an lme4 fit when a random intercept was
  added.

- family:

  `"binomial"` or `"gaussian"`.

- n:

  The number of rows the reported model was fitted on.

- n_groups:

  The number of groups; `NA` unless a random intercept was fitted.

- selected:

  The predictors in the reported model.

- dropped:

  The predictors dropped as non-numeric or constant.

- adjust, select, correction, ci_level, outcome:

  The settings used.

- mixed:

  Whether a random intercept was fitted.

## Details

Ordinary gaussian fits use Student t confidence limits with the model's
residual degrees of freedom; binomial and mixed fits use normal (Wald)
limits. A linear mixed model
([`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)) reports no
p-values, so they are computed from the normal distribution of its t
statistics.

## Honest inference

Choosing predictors by their association with the outcome and then
testing them on the same rows invalidates the p-values. With
`select = "split"` the data is halved: predictors are ranked on one half
and the reported model is fitted on the other, so the returned inference
is valid for the selected set. `select = "none"` (default) fits every
supplied predictor and needs no split.

## See also

[`simplicial_features`](https://pak.dynasite.org/Nestimate/reference/simplicial_features.md),
[`effects_table`](https://pak.dynasite.org/Nestimate/reference/effects_table.md)

## Examples

``` r
set.seed(1)
d <- data.frame(
  hint    = rbinom(300, 1, 0.4),
  think   = rbinom(300, 1, 0.3),
  n_events = rpois(300, 20),
  actor   = rep(letters[1:10], each = 30)
)
d$success <- rbinom(300, 1, plogis(-0.5 + 0.8 * d$hint))
fit <- outcome_model(d, outcome = "success",
                     predictors = c("hint", "think"),
                     adjust = "n_events")
effects_table(fit)
#>       term estimate std_error statistic ci_lower ci_upper     p_value
#> 1 n_events    0.064     0.026     2.434    0.012    0.115 0.014937387
#> 2     hint    0.635     0.245     2.591    0.155    1.116 0.009562848
#> 3    think   -0.048     0.258    -0.188   -0.555    0.458 0.851074699
#>        p_adj odds_ratio or_lower or_upper
#> 1 0.02240608      1.066    1.012    1.122
#> 2 0.02240608      1.888    1.167    3.053
#> 3 0.85107470      0.953    0.574    1.580
```
