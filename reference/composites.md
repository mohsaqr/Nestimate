# Cluster Scores From a Psychometric MCML Fit

The per-observation cluster scores that a re-estimated
[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
macro network was fitted on: each respondent's weighted, sign-corrected
score on every cluster. These are the scores to carry into a profile
analysis, a regression, or any downstream model that needs one number
per cluster per respondent.

## Usage

``` r
composites(x, ...)

# S3 method for class 'mcml_pc'
composites(x, ...)
```

## Arguments

- x:

  An object carrying cluster scores.

- ...:

  Ignored.

## Value

A data frame with one row per row of the input data, in input order and
with the input's row names, and one numeric column per cluster, named by
the cluster. A row whose members of a cluster are all missing is `NA` in
that column (items missing only in part are averaged over the observed
ones). The macro network is estimated on the complete rows, so
`build_network(composites(fit), method = ...)` reproduces it.

The `mcml_pc` method errors with class `"nestimate_no_composites"` for
the descriptive aggregations (`"average"`, `"escoufier"`, `"cancor"`),
which relate clusters without ever forming a score.

## See also

[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
to create the fit,
[`item_loadings`](https://saqr.me/Nestimate/reference/item_loadings.md)
for the item weights behind these scores.

## Examples

``` r
set.seed(1)
f <- stats::rnorm(200)
g <- stats::rnorm(200)
df <- data.frame(a1 = f + stats::rnorm(200), a2 = f + stats::rnorm(200),
                 a3 = f + stats::rnorm(200), b1 = g + stats::rnorm(200),
                 b2 = g + stats::rnorm(200), b3 = g + stats::rnorm(200))
clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
fit <- build_mcml_pc(df, clusters, aggregation = "loadings",
                     method = "cor")
head(composites(fit))
#>             A          B
#> 1 -0.51616641  0.5188386
#> 2  0.46126437  1.9244471
#> 3 -0.34848366  0.3766971
#> 4  1.12989912 -0.8649034
#> 5 -0.03460362 -1.5703036
#> 6 -0.93775141  1.3283922
```
