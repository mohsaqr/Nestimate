# Build a Partial Correlation Network

Convenience wrapper for `build_network(method = "pcor")`. Computes
partial correlations from numeric data.

## Usage

``` r
build_pcor(data, ...)
```

## Arguments

- data:

  Data frame (sequences or per-observation frequencies) or a square
  symmetric matrix (correlation or covariance). A fitted
  `net_clustering` or `net_mmm` object is also accepted: the per-cluster
  networks are (re)built and a `netobject_group` is returned.

- ...:

  Additional arguments passed to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).

## Value

A `netobject` (see
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)).

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)

## Examples

``` r
data(srl_strategies)
net <- build_pcor(srl_strategies)
```
