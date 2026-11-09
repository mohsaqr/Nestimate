# Cluster Scores From a Psychometric MCML Fit

`composites()` returns the per-observation cluster scores that the macro
network of a re-estimated
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
fit was estimated from: each respondent's weighted, sign-corrected score
on every cluster. These are the scores to carry into a profile analysis,
a regression or any other model that needs one number per cluster per
respondent.

A score is only as good as the items behind it: check the item weights
and misfit flags with
[`item_loadings()`](https://pak.dynasite.org/Nestimate/reference/item_loadings.md)
and the stability of the weights with
[`loading_stability()`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
before using the scores elsewhere.

## Usage

``` r
composites(x, ...)

# S3 method for class 'mcml_pc'
composites(x, ...)
```

## Arguments

- x:

  An `mcml_pc` fit with a score aggregation (`"scaled"`, `"composite"`,
  `"mean"`, `"median"` or `"loadings"`). The descriptive aggregations
  (`"average"`, `"escoufier"`, `"cancor"`) relate the clusters without
  ever forming a score and raise an error of class
  `nestimate_no_composites`.

- ...:

  Ignored.

## Value

A data frame with one row per row of the input data, in input order and
with the input's row names, and one numeric column per cluster, named by
the cluster (as a syntactic name, so `"Clu A"` becomes `Clu.A`). A row
whose members of a cluster are all missing is `NA` in that column; items
missing only in part are averaged over the observed ones. The macro
network is estimated on the complete rows, so
`build_network(composites(fit), method = ...)` with the fit's estimator
reproduces it.

## See also

[`build_mcml_pc`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
to create the fit,
[`item_loadings`](https://pak.dynasite.org/Nestimate/reference/item_loadings.md)
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
