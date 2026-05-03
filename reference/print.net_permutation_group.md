# Print Method for net_permutation_group

Print Method for net_permutation_group

## Usage

``` r
# S3 method for class 'net_permutation_group'
print(x, ...)
```

## Arguments

- x:

  A `net_permutation_group` object.

- ...:

  Additional arguments (ignored).

## Value

`x` invisibly.

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
  V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
nets1 <- build_network(s1, method = "relative", group = "grp")
nets2 <- build_network(s2, method = "relative", group = "grp")
perm  <- permutation(nets1, nets2, iter = 10)
print(perm)
#> Grouped Permutation Test
#> Groups: X, Y 
#> Use summary() on each element for edge-level results.
# \donttest{
set.seed(1)
s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
                 V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
                 V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
nets1 <- build_network(s1, method = "relative", group = "grp")
nets2 <- build_network(s2, method = "relative", group = "grp")
perm  <- permutation(nets1, nets2, iter = 20, seed = 1)
print(perm)
#> Grouped Permutation Test
#> Groups: X, Y 
#> Use summary() on each element for edge-level results.
# }
```
