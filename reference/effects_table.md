# Effect Table of a Fitted Outcome Model

The tidy one-row-per-term table of estimates, confidence intervals and
corrected p-values.

## Usage

``` r
effects_table(
  x,
  intercept = FALSE,
  significant = FALSE,
  alpha = 0.05,
  digits = 3
)
```

## Arguments

- x:

  A `net_outcome_model` from
  [`outcome_model`](https://saqr.me/Nestimate/reference/outcome_model.md).

- intercept:

  Keep the intercept row? Default `FALSE`.

- significant:

  Keep only terms whose corrected p-value is below `alpha`? Default
  `FALSE`.

- alpha:

  Threshold used by `significant`. Default `0.05`.

- digits:

  Rounding for the numeric columns. Default `3`; `p_value` and `p_adj`
  are never rounded.

## Value

A `data.frame` with the same columns as the model's effect table, one
row per retained term. The intercept row is dropped unless
`intercept = TRUE`, and every term is kept unless `significant = TRUE`
restricts them to `p_adj < alpha`.

## Examples

``` r
set.seed(1)
d <- data.frame(hint = rbinom(200, 1, 0.5))
d$success <- rbinom(200, 1, plogis(-0.3 + 0.9 * d$hint))
effects_table(outcome_model(d, outcome = "success", predictors = "hint"))
#>   term estimate std_error statistic ci_lower ci_upper      p_value        p_adj
#> 1 hint    1.522     0.305     4.994    0.925     2.12 5.910034e-07 5.910034e-07
#>   odds_ratio or_lower or_upper
#> 1      4.583    2.522     8.33
```
