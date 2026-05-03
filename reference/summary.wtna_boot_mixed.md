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
oh <- data.frame(A = c(1,0,1,0), B = c(0,1,0,1), C = c(1,1,0,0))
mixed <- wtna(oh, method = "both")
boot  <- bootstrap_network(mixed, iter = 10)
summary(boot)
#> $transition
#>   from to weight mean        sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  B      2  1.2 0.6324555 0.7272727 FALSE    0.225    2.000     1.50
#> 2    A  C      1  1.2 0.7888106 0.3636364 FALSE    0.225    2.775     0.75
#> 3    B  A      1  0.9 0.7378648 0.5454545 FALSE    0.000    2.000     0.75
#> 4    C  A      1  0.6 0.6992059 0.6363636 FALSE    0.000    1.775     0.75
#> 5    C  B      1  1.1 0.5676462 0.3636364 FALSE    0.225    2.000     0.75
#> 6    C  C      1  0.8 0.7888106 0.6363636 FALSE    0.000    2.000     0.75
#>   cr_upper
#> 1     2.50
#> 2     1.25
#> 3     1.25
#> 4     1.25
#> 5     1.25
#> 6     1.25
#> 
#> $cooccurrence
#>   from to weight mean        sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  A      2  2.4 0.9660918 0.7272727 FALSE    1.000    3.775     1.50
#> 2    A  C      1  0.8 0.4216370 0.2727273 FALSE    0.000    1.000     0.75
#> 3    B  B      2  1.6 0.9660918 0.7272727 FALSE    0.225    3.000     1.50
#> 4    B  C      1  0.7 0.8232726 0.7272727 FALSE    0.000    2.000     0.75
#> 5    C  C      2  1.5 0.9718253 0.8181818 FALSE    0.225    3.000     1.50
#>   cr_upper
#> 1     2.50
#> 2     1.25
#> 3     2.50
#> 4     1.25
#> 5     2.50
#> 
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
#> 1    A  B      4 1.65 0.6708204 1.0000000 FALSE    0.475    2.525     3.00
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
#> 1    A  A      4 4.15 1.531253 0.8571429 FALSE    1.475    6.525      3.0
#> 2    A  C      2 2.40 1.142481 0.6666667 FALSE    1.000    4.000      1.5
#> 3    B  B      4 3.85 1.531253 0.8571429 FALSE    1.475    6.525      3.0
#> 4    B  C      2 1.90 1.333772 0.8571429 FALSE    0.000    4.525      1.5
#> 5    C  C      4 4.30 1.380313 0.8571429 FALSE    2.000    6.000      3.0
#>   cr_upper
#> 1      5.0
#> 2      2.5
#> 3      5.0
#> 4      2.5
#> 5      5.0
#> 
# }
```
