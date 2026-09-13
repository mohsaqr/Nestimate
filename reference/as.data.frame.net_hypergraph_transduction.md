# Coerce a net_hypergraph_transduction to a data.frame

Coerce a net_hypergraph_transduction to a data.frame

## Usage

``` r
# S3 method for class 'net_hypergraph_transduction'
as.data.frame(x, what = c("predictions", "scores"), ...)
```

## Arguments

- x:

  A `net_hypergraph_transduction` object.

- what:

  Character. `"predictions"` (default) for the one-row-per-node table,
  `"scores"` for the tidy long score table (one row per node x class:
  `node`, `class`, `score`).

- ...:

  Additional arguments (ignored).

## Value

A data.frame as selected by `what`.
