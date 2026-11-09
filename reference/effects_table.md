# Effect Table of a Fitted Outcome Model

`effects_table()` returns the effect table of an
[`outcome_model()`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
fit: one row per term, with its estimate, confidence interval and
corrected p-value, rounded for reading. Arguments keep the intercept or
only the significant terms.

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

  An
  [`outcome_model()`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  result; anything else raises an error.

- intercept:

  Whether to keep the intercept row. Default `FALSE`.

- significant:

  Whether to keep only the terms whose corrected p-value (`p_adj`) is
  below `alpha`. Default `FALSE`.

- alpha:

  The threshold for `significant`, between 0 and 1. Default `0.05`.

- digits:

  The decimal places of the numeric columns. Default `3`; `p_value` and
  `p_adj` are never rounded.

## Value

A data frame with the columns of the model's effect table (`term`,
`estimate`, `std_error`, `statistic`, `ci_lower`, `ci_upper`, `p_value`,
`p_adj`, and for a binomial fit `odds_ratio`, `or_lower` and
`or_upper`), one row per retained term. The intercept row is dropped
unless `intercept = TRUE`, and every term is kept unless
`significant = TRUE` restricts them to `p_adj < alpha`.

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
