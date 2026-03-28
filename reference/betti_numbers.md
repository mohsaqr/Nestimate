# Betti Numbers

Computes Betti numbers: \\\beta_0\\ (components), \\\beta_1\\ (loops),
\\\beta_2\\ (voids), etc.

## Usage

``` r
betti_numbers(sc)
```

## Arguments

- sc:

  A `simplicial_complex` object.

## Value

Named integer vector `c(b0 = ..., b1 = ..., ...)`.

## Examples

``` r
# \donttest{
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
betti_numbers(sc)
#> b0 b1 b2 
#>  1  0  0 
# }
```
