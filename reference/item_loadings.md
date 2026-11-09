# Item Diagnostics From a Psychometric MCML Fit

`item_loadings()` returns the item table behind a
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
fit: one row per item, with how strongly the item connects to its own
cluster, the weight it carries into that cluster's score, and whether it
connects more strongly to some other cluster. The loadings are network
loadings, read from the node-level network (an item's mean absolute
connection to the members of a cluster), not factor loadings.

An item flagged as `misfit` connects more strongly to another cluster
than to its own, which is evidence that it is assigned to the wrong
cluster; every aggregation of the fit inherits that membership. A `sign`
of -1 marks an item treated as reverse-keyed. How stable the weights are
under resampling is checked with
[`loading_stability()`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md).

## Usage

``` r
item_loadings(x, ...)

# S3 method for class 'mcml_pc'
item_loadings(x, misfit = NULL, ...)
```

## Arguments

- x:

  An `mcml_pc` object from
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md).
  A fit that carries no item table raises an error of class
  `nestimate_no_loadings`.

- ...:

  Ignored.

- misfit:

  Which items to return: `NULL` (default) returns every item, `TRUE`
  only those whose strongest cross-cluster connection exceeds their
  own-cluster loading, and `FALSE` only those that fit where they were
  assigned.

## Value

A data frame with one row per item (only the misfitting or the fitting
items when `misfit` is set) and columns `node`, `cluster`, `loading`
(the mean absolute connection to the other members of its own cluster,
signed by the leading eigenvector of the cluster's block; under a
sign-carrying weighting the `sign` column holds the sign actually used,
which can differ), `weight` (its weight in the cluster score, or its
strength weight for a descriptive aggregation; the weights of a cluster
sum to 1), `sign` (+1, or -1 for a reverse-keyed item), `max_cross` (its
strongest mean absolute connection to any other cluster),
`cross_cluster` (which cluster that is) and `misfit` (logical). An item
alone in its cluster has loading 0, weight 1 and is never flagged.

## See also

[`build_mcml_pc`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
to create the fit,
[`loading_stability`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
for the weights' sampling uncertainty,
[`composites`](https://pak.dynasite.org/Nestimate/reference/composites.md)
for the scores these weights produce.

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
item_loadings(fit)
#>   node cluster   loading    weight sign  max_cross cross_cluster misfit
#> 1   a1       A 0.4535516 0.3288962    1 0.02471784             B  FALSE
#> 2   a2       A 0.4572043 0.3315450    1 0.04824045             B  FALSE
#> 3   a3       A 0.4682554 0.3395588    1 0.05652080             B  FALSE
#> 4   b1       B 0.5363571 0.3493775    1 0.05405353             A  FALSE
#> 5   b2       B 0.4977645 0.3242387    1 0.01233061             A  FALSE
#> 6   b3       B 0.5010576 0.3263838    1 0.06309494             A  FALSE
item_loadings(fit, misfit = TRUE)
#> [1] node          cluster       loading       weight        sign         
#> [6] max_cross     cross_cluster misfit       
#> <0 rows> (or 0-length row.names)
```
