# Euler Characteristic

`euler_characteristic()` returns the Euler characteristic of a
simplicial complex, \\\chi = \sum\_{k=0}^{d} (-1)^k f_k\\, where \\f_k\\
is the number of k-simplices and \\d\\ the dimension of the complex. By
the Euler-Poincare theorem it also equals \\\sum\_{k} (-1)^k \beta_k\\,
the alternating sum of the
[`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md).

## Usage

``` r
euler_characteristic(sc)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md).

## Value

A single integer.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
euler_characteristic(sc)
#> [1] 1
```
