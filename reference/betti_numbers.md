# Betti Numbers

`betti_numbers()` counts the holes of a simplicial complex in each
dimension: \\\beta_0\\ is the number of connected components,
\\\beta_1\\ the number of independent loops (cycles not filled in by
triangles), \\\beta_2\\ the number of enclosed voids, and so on up to
the dimension of the complex.

Betti numbers describe the complex as it was built; a clique complex
built at a different threshold can give different numbers.
[`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
follows each hole across all thresholds and shows how long it lasts.

## Usage

``` r
betti_numbers(sc)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md).
  When it was built with a `max_dim` below its largest clique, the top
  Betti number counts cycles whose filling simplices were not built, so
  it is an upper bound, and a warning of class
  `nestimate_truncated_complex` says so.

## Value

A named integer vector `c(b0 = , b1 = , ...)`, one entry per dimension
from 0 to the dimension of the complex.

## Details

\\\beta_k\\ is the dimension of the k-th homology group, the number of
k-simplices minus the rank of the k-th boundary matrix, minus the rank
of the (k + 1)-th. The ranks are computed with
[`qr()`](https://rdrr.io/r/base/qr.html) on the signed boundary
matrices.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
betti_numbers(sc)
#> b0 b1 b2 
#>  1  0  0 
```
