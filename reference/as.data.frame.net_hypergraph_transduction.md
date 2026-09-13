# Coerce a net_hypergraph_transduction to a data.frame

Coerce a net_hypergraph_transduction to a data.frame

## Usage

``` r
# S3 method for class 'net_hypergraph_transduction'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("predictions", "scores"),
  ...
)
```

## Arguments

- x:

  A `net_hypergraph_transduction` object.

- row.names:

  `NULL` (default) or a character vector of row names for the returned
  data frame.

- optional:

  Ignored; present so the method matches the signature of the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  generic.

- what:

  Character. `"predictions"` (default) for the one-row-per-node table,
  `"scores"` for the tidy long score table (one row per node x class:
  `node`, `class`, `score`).

- ...:

  Additional arguments (ignored).

## Value

A data.frame selected by `what`: for `"predictions"`, one row per node
with columns `node`, `label` (the given label, `NA` if unlabeled),
`predicted`, `score` and `margin`; for `"scores"`, one row per node x
class with columns `node`, `class` and `score`.
