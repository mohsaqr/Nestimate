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
net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
  method = "relative")
boot <- bootstrap_network(net, iter = 10)
summary(boot)
#>   from to weight mean        sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  B      1  0.7 0.4830459 0.3636364 FALSE    0.000        1     0.75
#> 2    B  C      1  0.9 0.3162278 0.1818182 FALSE    0.225        1     0.75
#> 3    C  A      1  0.6 0.5163978 0.4545455 FALSE    0.000        1     0.75
#>   cr_upper
#> 1     1.25
#> 2     1.25
#> 3     1.25
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
