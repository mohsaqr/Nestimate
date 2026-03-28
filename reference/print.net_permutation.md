# Print Method for net_permutation

Print Method for net_permutation

## Usage

``` r
# S3 method for class 'net_permutation'
print(x, ...)
```

## Arguments

- x:

  A `net_permutation` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

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
print(perm)
#> Permutation Test:Transition Network (relative probabilities) [directed]
#>   Iterations: 20  |  Alpha: 0.05
#>   Nodes: 3  |  Edges tested: 3  |  Significant: 0
# }
```
