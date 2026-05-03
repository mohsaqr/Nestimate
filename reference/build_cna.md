# Build a Co-occurrence Network (CNA)

Convenience wrapper for `build_network(method = "co_occurrence")`.
Computes co-occurrence counts from binary or sequence data.

## Usage

``` r
build_cna(data, ...)
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

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
[`cooccurrence`](https://mohsaqr.github.io/Nestimate/reference/cooccurrence.md)
for delimited-field, bipartite, and other non-sequence co-occurrence
formats.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
net <- build_cna(seqs)
```
