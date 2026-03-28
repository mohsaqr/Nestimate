# Print Method for net_bootstrap

Print Method for net_bootstrap

## Usage

``` r
# S3 method for class 'net_bootstrap'
print(x, ...)
```

## Arguments

- x:

  A `net_bootstrap` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

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
print(boot)
#>   Edge                   Mean     95% CI          p
#>   -----------------------------------------------
#>   B → C                1.000  [1.000, 1.000]  *  
#>   C → A                1.000  [1.000, 1.000]  *  
#> 
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 20  |  Nodes : 3
#>   Edges      : 2 significant / 3 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
# }
```
