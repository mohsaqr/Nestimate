# Item Diagnostics From a Psychometric MCML Fit

The item table behind
[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md):
one row per node, reporting how strongly the item connects to its own
cluster, the weight it carries into that cluster's composite, and
whether it connects more strongly to some other cluster.

Reading the table is `item_loadings(fit)` - never a reach into the fit's
internals. The name says *item*: these are network loadings of items on
their own cluster, not factor loadings.

## Usage

``` r
item_loadings(x, ...)

# S3 method for class 'mcml_pc'
item_loadings(x, misfit = NULL, ...)
```

## Arguments

- x:

  An `mcml_pc` object from
  [`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md).

- ...:

  Ignored.

- misfit:

  Logical or NULL. `NULL` (default) returns every item; `TRUE` returns
  only the items whose strongest cross-cluster connection exceeds their
  own-cluster loading; `FALSE` returns only the items that fit where
  they were assigned.

## Value

A data frame with one row per node (one row per misfitting or fitting
node when `misfit` is set) and columns `node`, `cluster`, `loading`
(signed mean connection to its own cluster), `weight` (its composite
weight), `sign` (+1, or -1 for a reverse-keyed item), `max_cross`
(strongest connection to any other cluster), `cross_cluster` (which
cluster that is), and `misfit` (logical).

## See also

[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
to create the fit,
[`loading_stability`](https://saqr.me/Nestimate/reference/loading_stability.md)
for the weights' sampling uncertainty,
[`composites`](https://saqr.me/Nestimate/reference/composites.md) for
the scores these weights produce.

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
