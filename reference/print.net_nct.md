# Print Method for net_nct

Print Method for net_nct

## Usage

``` r
# S3 method for class 'net_nct'
print(x, ...)
```

## Arguments

- x:

  A `net_nct` object.

- ...:

  Ignored.

## Value

The input object, invisibly.

## Examples

``` r
set.seed(1)
x1 <- matrix(rnorm(200 * 5), 200, 5)
x2 <- matrix(rnorm(200 * 5), 200, 5)
colnames(x1) <- colnames(x2) <- paste0("V", 1:5)
# iter = 20 keeps the example fast; a real analysis uses 1000 or more.
res <- nct(x1, x2, iter = 20)
res
#> Network Comparison Test  [20 permutations | unpaired]
#>   Global strength (M):  observed = 0.0000   p = 0.5714
#>   Network structure (S): observed = 0.0000   p = 0.6667
#>   Edge tests (E):       10 edges, 0 significant at p < 0.05
summary(res)
#>    from to diff_observed   p_value significant
#> 1    V1 V2  0.000000e+00 1.0000000       FALSE
#> 2    V1 V3  0.000000e+00 1.0000000       FALSE
#> 3    V2 V3  0.000000e+00 1.0000000       FALSE
#> 4    V1 V4  1.387779e-17 0.1428571       FALSE
#> 5    V2 V4  0.000000e+00 1.0000000       FALSE
#> 6    V3 V4  0.000000e+00 1.0000000       FALSE
#> 7    V1 V5  0.000000e+00 1.0000000       FALSE
#> 8    V2 V5  0.000000e+00 1.0000000       FALSE
#> 9    V3 V5  0.000000e+00 1.0000000       FALSE
#> 10   V4 V5  0.000000e+00 1.0000000       FALSE
```
