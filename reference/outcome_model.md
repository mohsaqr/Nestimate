# Model Unit-Level Outcomes from Sequence or Network Predictors

Fits a regression of an outcome on predictor columns – pattern
indicators, topological features from
[`simplicial_features`](https://saqr.me/Nestimate/reference/simplicial_features.md),
or any other numeric covariates – and returns a tidy effect table with
confidence intervals and multiplicity-corrected p-values.

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

  A `data.frame` with one row per unit of analysis.

- outcome:

  Name of the outcome column. A two-valued outcome is modelled with a
  binomial family, a numeric one with gaussian; override with `family`.

- predictors:

  Character vector of predictor column names.

- group:

  Optional column name giving a grouping factor. When supplied and lme4
  is installed, a random intercept per group is added – the right
  treatment for units nested in actors. lme4 is a suggested package:
  when it is not installed the random intercept is dropped, a
  `"nestimate_no_lme4"` warning is raised, and a plain
  [`glm`](https://rdrr.io/r/stats/glm.html) is fitted instead.

- adjust:

  Optional character vector of covariates entered before the predictors.
  Use it for exposure: a unit observed longer contains more of every
  pattern, so an unadjusted effect can be volume in disguise.

- family:

  `"auto"` (default), `"binomial"` or `"gaussian"`.

- select:

  `"none"` (default) fits all predictors; `"split"` ranks them on half
  the data and fits on the other half.

- n_select:

  Number of predictors kept when `select = "split"`. Default `10`.

- correction:

  Multiplicity correction passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`.

- ci_level:

  Confidence level for the intervals. Default `0.95`.

- seed:

  Optional integer seed for the split, so the result is reproducible.
  The RNG state is restored on exit.

- x:

  A `net_outcome_model`, for the `print` and `plot` methods.

- ...:

  Unused.

- object:

  A `net_outcome_model`.

## Value

An object of class `net_outcome_model`: a list whose `$effects` element
is the tidy `data.frame`, one row per model term (the intercept
included), with columns `term`, `estimate`, `std_error`, `statistic`,
`ci_lower`, `ci_upper`, `p_value` and `p_adj` (`NA` on the intercept
row, which is excluded from the correction), plus `odds_ratio`,
`or_lower` and `or_upper` for a binomial fit. The remaining elements are
the fitted `$model` (a `glm`, or an lme4 fit when a random intercept was
added), `$family`, `$n` (rows the reported model was fitted on),
`$n_groups` (`NA` unless mixed), `$selected`, `$dropped` (zero-variance
predictors), `$adjust`, `$select`, `$correction`, `$ci_level`, `$mixed`
and `$outcome`. Retrieve the table with
[`effects_table`](https://saqr.me/Nestimate/reference/effects_table.md).

`print` returns its input invisibly.

`summary` returns the tidy effect table.

`plot` returns a `ggplot` forest of the effects.

## Honest inference

Choosing predictors by their association with the outcome and then
testing them on the same rows invalidates the p-values. With
`select = "split"` the data is halved: predictors are ranked on one half
and the reported model is fitted on the other, so the returned inference
is valid for the selected set. `select = "none"` (default) fits every
supplied predictor and needs no split.

## See also

[`simplicial_features`](https://saqr.me/Nestimate/reference/simplicial_features.md),
[`effects_table`](https://saqr.me/Nestimate/reference/effects_table.md)

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
