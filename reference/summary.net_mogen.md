# Summary Method for net_mogen

Summary Method for net_mogen

## Usage

``` r
# S3 method for class 'net_mogen'
summary(object, ...)
```

## Arguments

- object:

  A `net_mogen` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
mog <- build_mogen(seqs, max_order = 2L)
summary(mog)
#> Multi-Order Generative Model (MOGen) Summary
#> 
#>   States: A, B, C
#>   Paths: 5 | Observations: 15
#> 
#>  order layer_dof cum_dof loglik   aic   bic selected
#>      0         2       2 -16.48 36.96 38.37         
#>      1         0       2  -5.49 14.99 16.40      <--
#>      2         0       2  -5.49 14.99 16.40         
#> 
#>   Optimal order: 1 (by aic)
# }
```
