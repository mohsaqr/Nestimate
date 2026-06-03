# Print method for tidy covariate output

Prints a one-line header naming the estimator, then the data.frame. The
full human-readable view (per-cluster stats, profiles, OR/test tables)
was already printed by
[`summary()`](https://rdrr.io/r/base/summary.html) when this object was
produced, so this method intentionally stays minimal to avoid
duplication. Auto-prints when the user types the variable at the REPL.

## Usage

``` r
# S3 method for class 'tidy_covariates'
print(x, ...)
```

## Arguments

- x:

  A `tidy_covariates` data.frame.

- ...:

  Ignored.

## Value

The input invisibly.
