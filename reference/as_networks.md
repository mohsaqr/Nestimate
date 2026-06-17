# Promote a psychometric MCML result to a network group

`as_networks()` is the psychometric-network counterpart of
[`as_tna`](https://saqr.me/Nestimate/reference/as_tna.md). It promotes
the cluster-level (macro) and within-cluster networks produced by
[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
into a single `netobject_group`, so the result flows into the same
downstream verbs as any other group of networks
([`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`net_centrality`](https://saqr.me/Nestimate/reference/net_centrality.md)).

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
  [`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md))
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

Where [`as_tna()`](https://saqr.me/Nestimate/reference/as_tna.md)
promotes *transition* networks (directed, row-normalised, with initial
probabilities) and re-wraps raw matrices, `as_networks()` promotes
*psychometric* networks (undirected; correlation / partial-correlation /
glasso). The macro and within-cluster components of an `mcml_pc` object
are already full netobjects carrying their estimator, directedness and
data, so this function assembles them into a group rather than
re-wrapping matrices.

## See also

[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
to create the input,
[`as_tna`](https://saqr.me/Nestimate/reference/as_tna.md) for the
transition-network counterpart.

## Examples

``` r
set.seed(1)
df <- as.data.frame(matrix(stats::rnorm(200 * 6), 200, 6))
names(df) <- c("a1", "a2", "a3", "b1", "b2", "b3")
clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
fit <- build_mcml_pc(df, clusters, aggregation = "composite", method = "cor")
#> Warning: Item(s) more strongly connected to another cluster than their own (possible misassignment): a1, a2, b1. See $loadings (misfit, cross_cluster).
#> Warning: Reverse-keyed item(s) flipped in composites: a2, b3. See $loadings (sign).
nets <- as_networks(fit)
nets
#> Group Networks (3 groups)
#> 
#>   Group  Nodes  Edges  Weights
#>   macro  2      1      [0.028, 0.028]
#>   A      3      3      [-0.026, 0.068]
#>   B      3      3      [-0.070, 0.038]
nets$macro$weights
#>            A          B
#> A 0.00000000 0.02847388
#> B 0.02847388 0.00000000
```
