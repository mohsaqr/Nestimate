# Simplicial Degree

`simplicial_degree()` counts, for each node of a simplicial complex, how
many simplices of each dimension contain it: edges (dimension 1),
triangles (dimension 2) and so on. It is the higher-order counterpart of
a node's degree, and says which nodes take part in many group-level
structures rather than only in pairs.

The counts describe the complex as built, so they move with the
threshold or pathways
[`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
was given. A node in a clique of m nodes is in `choose(m - 1, k)` of its
k-simplices, so one large clique raises the node's counts in every
dimension up to the clique's own; compare nodes within a dimension.

## Usage

``` r
simplicial_degree(sc, normalized = FALSE)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md).

- normalized:

  Whether to divide each count of dimension k by the largest possible
  one, `choose(n_nodes - 1, k)`, the number of k-simplices a node can be
  in. Default `FALSE`.

## Value

A data frame with one row per node, sorted by decreasing `total`:
`node`, then `d0`, `d1`, ... up to the dimension of the complex (the
number of simplices of that dimension containing the node; `d0` is
always 1), and `total`, the sum of `d1` onwards.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
simplicial_degree(sc)
#>   node d0 d1 d2 total
#> 1    A  1  2  1     3
#> 2    B  1  2  1     3
#> 3    C  1  2  1     3
```
