# Summary Method for net_bootstrap

Summary Method for net_bootstrap

## Usage

``` r
# S3 method for class 'net_bootstrap'
summary(object, ...)
```

## Arguments

- object:

  A `net_bootstrap` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with edge-level bootstrap statistics.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
  V3 = c("C","A","C","B","A")
)
net  <- build_network(seqs, method = "relative")
boot <- bootstrap_network(net, iter = 20)
summary(boot)
#>   from to weight mean        sd    p_value   sig ci_lower ci_upper cr_lower
#> 1    A  B      1 0.95 0.2236068 0.09523810 FALSE    0.475        1     0.75
#> 2    B  C      1 1.00 0.0000000 0.04761905  TRUE    1.000        1     0.75
#> 3    C  A      1 1.00 0.0000000 0.04761905  TRUE    1.000        1     0.75
#>   cr_upper
#> 1     1.25
#> 2     1.25
#> 3     1.25
# }
```
