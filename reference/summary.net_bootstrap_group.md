# Summary Method for net_bootstrap_group

Summary Method for net_bootstrap_group

## Usage

``` r
# S3 method for class 'net_bootstrap_group'
summary(object, ...)
```

## Arguments

- object:

  A `net_bootstrap_group` object.

- ...:

  Ignored.

## Value

A data frame with group, edge, and bootstrap statistics columns.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = c("A","B","A","C","B","A"),
  V2 = c("B","C","B","A","C","B"),
  V3 = c("C","A","C","B","A","C"),
  grp = c("X","X","X","Y","Y","Y")
)
nets <- build_network(seqs, method = "relative", group = "grp")
boot <- bootstrap_network(nets, iter = 20)
summary(boot)
#>   group from to weight mean        sd    p_value   sig ci_lower ci_upper
#> 1     X    A  B      1 0.80 0.4103913 0.23809524 FALSE    0.000        1
#> 2     X    B  C      1 1.00 0.0000000 0.04761905  TRUE    1.000        1
#> 3     X    C  A      1 0.80 0.4103913 0.23809524 FALSE    0.000        1
#> 4     Y    A  B      1 0.95 0.2236068 0.09523810 FALSE    0.475        1
#> 5     Y    B  C      1 0.95 0.2236068 0.09523810 FALSE    0.475        1
#> 6     Y    C  A      1 0.90 0.3077935 0.14285714 FALSE    0.000        1
#>   cr_lower cr_upper
#> 1     0.75     1.25
#> 2     0.75     1.25
#> 3     0.75     1.25
#> 4     0.75     1.25
#> 5     0.75     1.25
#> 6     0.75     1.25
# }
```
