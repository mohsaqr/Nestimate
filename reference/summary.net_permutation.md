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
s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
n1 <- build_network(s1, method = "relative")
n2 <- build_network(s2, method = "relative")
perm <- permutation(n1, n2, iter = 10)
summary(perm)
#>   from to weight_x weight_y diff effect_size   p_value   sig
#> 1    A  B        1        0    1    1.333333 0.5454545 FALSE
#> 2    A  C        0        1   -1   -1.267449 0.5454545 FALSE
#> 3    B  A        0        1   -1   -1.360828 0.4545455 FALSE
#> 4    B  C        1        0    1    1.360828 0.4545455 FALSE
#> 5    C  A        1        0    1    1.237969 0.7272727 FALSE
#> 6    C  B        0        1   -1   -1.184698 0.7272727 FALSE
# \donttest{
set.seed(1)
d1 <- data.frame(V1 = c("A","B","A"), V2 = c("B","C","B"),
                 V3 = c("C","A","C"))
d2 <- data.frame(V1 = c("C","A","C"), V2 = c("A","B","A"),
                 V3 = c("B","C","B"))
net1 <- build_network(d1, method = "relative")
net2 <- build_network(d2, method = "relative")
perm <- permutation(net1, net2, iter = 20, seed = 1)
summary(perm)
#>   from to weight_x weight_y diff effect_size p_value   sig
#> 1    A  B        1        1    0           0       1 FALSE
#> 2    B  C        1        1    0           0       1 FALSE
#> 3    C  A        1        1    0           0       1 FALSE
# }
```
