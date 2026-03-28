# Summary Method for wtna_boot_mixed

Summary Method for wtna_boot_mixed

## Usage

``` r
# S3 method for class 'wtna_boot_mixed'
summary(object, ...)
```

## Arguments

- object:

  A `wtna_boot_mixed` object.

- ...:

  Additional arguments (ignored).

## Value

A list with `$transition` and `$cooccurrence` summary data frames.

## Examples

``` r
# \donttest{
set.seed(1)
oh <- data.frame(
  A = c(1,0,1,0,1,0,1,0),
  B = c(0,1,0,1,0,1,0,1),
  C = c(1,1,0,0,1,1,0,0)
)
mixed <- wtna(oh, method = "both")
boot  <- bootstrap_network(mixed, iter = 20)
summary(boot)
#> $transition
#>   from to weight mean        sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  B      4 1.65 0.6708204 0.9523810 FALSE    0.475    2.525     3.00
#> 2    A  C      2 1.65 1.2680279 0.7619048 FALSE    0.000    4.000     1.50
#> 3    B  A      3 1.50 0.6882472 0.9523810 FALSE    0.475    2.525     2.25
#> 4    B  C      1 1.95 1.2763022 0.6666667 FALSE    0.475    4.575     0.75
#> 5    C  A      2 1.65 1.4608937 0.8571429 FALSE    0.000    4.525     1.50
#> 6    C  B      2 2.20 1.3218806 0.5714286 FALSE    0.000    4.525     1.50
#> 7    C  C      2 2.00 2.0000000 0.8095238 FALSE    0.000    6.525     1.50
#>   cr_upper
#> 1     5.00
#> 2     2.50
#> 3     3.75
#> 4     1.25
#> 5     2.50
#> 6     2.50
#> 7     2.50
#> 
#> $cooccurrence
#>   from to weight mean       sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  A      4 4.15 1.531253 0.3333333 FALSE    1.475    6.525      3.0
#> 2    A  C      2 2.40 1.142481 0.6666667 FALSE    1.000    4.000      1.5
#> 3    B  B      4 3.85 1.531253 0.3333333 FALSE    1.475    6.525      3.0
#> 4    B  C      2 1.90 1.333772 0.8571429 FALSE    0.000    4.525      1.5
#> 5    C  C      4 4.30 1.380313 0.3809524 FALSE    2.000    6.000      3.0
#>   cr_upper
#> 1      5.0
#> 2      2.5
#> 3      5.0
#> 4      2.5
#> 5      5.0
#> 
# }
```
