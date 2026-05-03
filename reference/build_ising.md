# Build an Ising Network

Convenience wrapper for `build_network(method = "ising")`. Computes
L1-regularized logistic regression network for binary data.

## Usage

``` r
build_ising(data, ...)
```

## Arguments

- data:

  Data frame (sequences or per-observation frequencies) or a square
  symmetric matrix (correlation or covariance).

- ...:

  Additional arguments passed to
  [`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).

## Value

A `netobject` (see
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)).

## See also

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)

## Examples

``` r
# \donttest{
bin_data <- data.frame(matrix(rbinom(200, 1, 0.5), ncol = 5))
net <- build_ising(bin_data)
# }
```
