# Build an Attention-Weighted Transition Network (ATNA)

Convenience wrapper for `build_network(method = "attention")`. Computes
decay-weighted transitions from sequence data.

## Usage

``` r
build_atna(data, ...)
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
seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
net <- build_atna(seqs)
```
