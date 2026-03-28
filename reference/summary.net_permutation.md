# Summary Method for net_permutation

Summary Method for net_permutation

## Usage

``` r
# S3 method for class 'net_permutation'
summary(object, ...)
```

## Arguments

- object:

  A `net_permutation` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with edge-level permutation test results.

## Examples

``` r
# \donttest{
set.seed(1)
d1 <- data.frame(V1 = c("A","B","A"), V2 = c("B","C","B"),
                 V3 = c("C","A","C"))
d2 <- data.frame(V1 = c("C","A","C"), V2 = c("A","B","A"),
                 V3 = c("B","C","B"))
net1 <- build_network(d1, method = "relative")
net2 <- build_network(d2, method = "relative")
perm <- permutation_test(net1, net2, iter = 20, seed = 1)
summary(perm)
#>   from to weight_x weight_y diff effect_size p_value   sig
#> 1    A  B        1        1    0           0       1 FALSE
#> 2    B  C        1        1    0           0       1 FALSE
#> 3    C  A        1        1    0           0       1 FALSE
# }
```
