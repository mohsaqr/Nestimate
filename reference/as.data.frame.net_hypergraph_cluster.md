# Coerce a net_hypergraph_cluster to a data.frame

Coerce a net_hypergraph_cluster to a data.frame

## Usage

``` r
# S3 method for class 'net_hypergraph_cluster'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `net_hypergraph_cluster` object.

- ...:

  Additional arguments (ignored).

## Value

The tidy assignment table: one row per node, columns `node`, `cluster`,
`pi` (stationary probability of the node under the Laplacian's random
walk) and the spectral-embedding coordinates `dim1..dimk`.
