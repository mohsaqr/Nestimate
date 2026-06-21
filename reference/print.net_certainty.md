# Print Method for net_certainty

Print Method for net_certainty

## Usage

``` r
# S3 method for class 'net_certainty'
print(x, ...)
```

## Arguments

- x:

  A `net_certainty` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
print(certainty(build_network(seqs, method = "relative")))
#> Certainty (Dirichlet)  [Transition Network (relative) | directed]
#>   Prior      : Dirichlet(0.50)  |  Nodes : 3
#>   Edges      : 0 certain / 3 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
```
