# Summary Method for net_permutation_group

Returns a combined summary data frame across all groups.

## Usage

``` r
# S3 method for class 'net_permutation_group'
summary(object, ...)
```

## Arguments

- object:

  A `net_permutation_group` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with group, edge, p_value, and sig columns.

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
  V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
nets1 <- build_network(s1, method = "relative", group = "grp")
nets2 <- build_network(s2, method = "relative", group = "grp")
perm  <- permutation(nets1, nets2, iter = 10)
summary(perm)
#>   group from to weight_x weight_y diff effect_size p_value   sig
#> 1     X    A  B        1        1    0           0       1 FALSE
#> 2     X    B  C        1        1    0           0       1 FALSE
#> 3     X    C  A        1        1    0           0       1 FALSE
#> 4     Y    A  B        1        1    0           0       1 FALSE
#> 5     Y    B  C        1        1    0           0       1 FALSE
#> 6     Y    C  A        1        1    0           0       1 FALSE
# \donttest{
set.seed(1)
s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
                 V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
                 V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
nets1 <- build_network(s1, method = "relative", group = "grp")
nets2 <- build_network(s2, method = "relative", group = "grp")
perm  <- permutation(nets1, nets2, iter = 20, seed = 1)
summary(perm)
#>   group from to weight_x weight_y diff effect_size p_value   sig
#> 1     X    A  B        1        1    0           0       1 FALSE
#> 2     X    B  C        1        1    0           0       1 FALSE
#> 3     X    C  A        1        1    0           0       1 FALSE
#> 4     Y    A  B        1        1    0           0       1 FALSE
#> 5     Y    B  C        1        1    0           0       1 FALSE
#> 6     Y    C  A        1        1    0           0       1 FALSE
# }
```
