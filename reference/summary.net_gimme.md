# Summary Method for net_gimme

Summary Method for net_gimme

## Usage

``` r
# S3 method for class 'net_gimme'
summary(object, ...)
```

## Arguments

- object:

  A `net_gimme` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
set.seed(1)
panel <- data.frame(
  id = rep(1:5, each = 20),
  t  = rep(seq_len(20), 5),
  A  = rnorm(100), B = rnorm(100), C = rnorm(100)
)
gm <- build_gimme(panel, vars = c("A","B","C"), id = "id", time = "t")
summary(gm)
#> GIMME Network Analysis -- Summary
#> ======================================== 
#> 
#> FIT INDICES (per subject)
#> ------------------------------ 
#>  metric   mean     sd     min   max
#>   RMSEA  0.107 0.1219  0.0000 0.279
#>    SRMR  0.128 0.0345  0.0966 0.180
#>    NNFI -0.598 1.5306 -3.3036 0.326
#>     CFI  0.562 0.4528  0.0000 1.000
#> 
#> 
#> GROUP-LEVEL PATHS ( 0 )
#> ------------------------------ 
#>   (none)
#> 
#> AVERAGE TEMPORAL COEFFICIENTS
#> ------------------------------ 
#>        A      B      C
#> A -0.008  0.000  0.000
#> B  0.000 -0.225  0.000
#> C  0.000  0.000 -0.175
#> 
#> AVERAGE CONTEMPORANEOUS COEFFICIENTS
#> ------------------------------ 
#>   A B C
#> A 0 0 0
#> B 0 0 0
#> C 0 0 0
# }
```
