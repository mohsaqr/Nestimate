# Print Method for net_bootstrap_group

Print Method for net_bootstrap_group

## Usage

``` r
# S3 method for class 'net_bootstrap_group'
print(x, ...)
```

## Arguments

- x:

  A `net_bootstrap_group` object.

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","C","A"),
  V3 = c("C","A","B","B"), grp = c("X","X","Y","Y"))
nets <- build_network(seqs, method = "relative", group = "grp")
boot <- bootstrap_network(nets, iter = 10)
print(boot)
#> Grouped Bootstrap  [2 groups | 10 iterations | 95% CI]
#>   X                     0 sig / 3 total
#>   Y                     0 sig / 4 total
#>   Shared (all groups)   0 edges
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
print(boot)
#> Grouped Bootstrap  [2 groups | 20 iterations | 95% CI]
#>   X                     1 sig / 3 total
#>   Y                     0 sig / 3 total
#>   Shared (all groups)   0 edges
# }
```
