# Promote a psychometric MCML result to a network group

`as_networks()` gathers the cluster-level (macro) network and the
within-cluster networks of a
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
fit into one `netobject_group`, so they go through the same verbs as any
other group of networks: [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md),
and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
for drawing. It is the psychometric-network counterpart of
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md).

[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
converts *transition* networks (directed, row-normalised, with initial
probabilities) and re-wraps raw matrices. The macro and within-cluster
networks of an `mcml_pc` object are already full netobjects of
*psychometric* networks (undirected; correlation, partial correlation or
glasso) carrying their estimator and directedness, and their data when
they were re-estimated, so `as_networks()` only assembles them into a
group.

## Usage

``` r
as_networks(x)

# S3 method for class 'mcml_pc'
as_networks(x)

# Default S3 method
as_networks(x)
```

## Arguments

- x:

  The object to convert: an `mcml_pc` fit from
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md),
  or a `netobject_group`, which is returned unchanged. Any other object
  raises an error.

## Value

A `netobject_group`, read by
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
and drawn by
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html):
a named list whose first element, `macro`, is the cluster-level network,
followed by one netobject per cluster with more than one member. A
singleton cluster has no within-cluster network and is left out, with a
warning naming it.

## See also

[`build_mcml_pc`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
to create the input,
[`as_tna`](https://pak.dynasite.org/Nestimate/reference/as_tna.md) for
the transition-network counterpart.

## Examples

``` r
set.seed(1)
f <- stats::rnorm(200)
g <- stats::rnorm(200)
df <- data.frame(a1 = f + stats::rnorm(200), a2 = f + stats::rnorm(200),
                 a3 = f + stats::rnorm(200), b1 = g + stats::rnorm(200),
                 b2 = g + stats::rnorm(200), b3 = g + stats::rnorm(200))
clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
fit <- build_mcml_pc(df, clusters, aggregation = "composite", method = "cor")
nets <- as_networks(fit)
nets
#> Group Networks (3 groups)
#> 
#>   Group  Nodes  Edges  Weights
#>   macro  2      1      [0.007, 0.007]
#>   A      3      3      [0.443, 0.472]
#>   B      3      3      [0.462, 0.540]
summary(nets)
#> Network metrics by group:
#>                       metric   macro       A       B
#>                   Node Count       2       3       3
#>                   Edge Count       1       3       3
#>              Network Density       1       1       1
#>                Mean Distance 0.00676  0.4597  0.5117
#>            Mean Out-Strength 0.00676  0.9193   1.023
#>              SD Out-Strength       0 0.01531 0.04279
#>             Mean In-Strength 0.00676  0.9193   1.023
#>               SD In-Strength       0 0.01531 0.04279
#>              Mean Out-Degree       1       2       2
#>                SD Out-Degree       0       0       0
#>  Centralization (Out-Degree)       0       0       0
#>   Centralization (In-Degree)       0       0       0
#>                  Reciprocity       1       1       1
```
