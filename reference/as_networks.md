# Promote a psychometric MCML result to a network group

`as_networks()` is the psychometric-network counterpart of
[`as_tna`](https://pak.dynasite.org/Nestimate/reference/as_tna.md). It
promotes the cluster-level (macro) and within-cluster networks produced
by
[`build_mcml_pc`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
into a single `netobject_group`, so the result flows into the same
downstream verbs as any other group of networks
([`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`net_centrality`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)).

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

  An object to convert. The `mcml_pc` method (from
  [`build_mcml_pc`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md))
  is the primary path.

## Value

A `netobject_group`: a named list whose first element is `macro` (the
cluster-level network), followed by one netobject per non-singleton
cluster.

The `mcml_pc` method returns a `netobject_group`; singleton clusters (no
within-network) are dropped with a
[`warning()`](https://rdrr.io/r/base/warning.html).

The default method returns the input unchanged if it is already a
`netobject_group`, otherwise it errors.

## Details

Where
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
promotes *transition* networks (directed, row-normalised, with initial
probabilities) and re-wraps raw matrices, `as_networks()` promotes
*psychometric* networks (undirected; correlation / partial-correlation /
glasso). The macro and within-cluster components of an `mcml_pc` object
are already full netobjects carrying their estimator, directedness and
data, so this function assembles them into a group rather than
re-wrapping matrices.

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
#>                   Edge Count       2       6       6
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
