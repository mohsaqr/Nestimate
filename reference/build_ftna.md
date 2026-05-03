# Build a Frequency Transition Network (FTNA)

Convenience wrapper for `build_network(method = "frequency")`. Computes
raw transition counts from sequence data.

## Usage

``` r
build_ftna(data, ...)
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
net <- build_ftna(seqs)
```
