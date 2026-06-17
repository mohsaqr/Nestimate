# Print Method for Multilevel Network Object

Print Method for Multilevel Network Object

## Usage

``` r
# S3 method for class 'netobject_ml'
print(x, ...)
```

## Arguments

- x:

  A `netobject_ml`.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
set.seed(1)
obs <- data.frame(id = rep(1:3, each = 5),
                  A = rnorm(15), B = rnorm(15), C = rnorm(15))
net_ml <- build_network(obs, method = "cor",
                         params = list(id = "id"), level = "both")
print(net_ml)
#> Multilevel Network (method: cor)
#> -- Between-person --
#>   Nodes: 3  |  Edges: 3
#>   Sample size: 3 (unique persons)
#> -- Within-person --
#>   Nodes: 3  |  Edges: 3
#>   Sample size: 15 (observations)
# \donttest{
set.seed(1)
obs <- data.frame(
  id  = rep(1:5, each = 8),
  A   = rnorm(40), B = rnorm(40),
  C   = rnorm(40), D = rnorm(40)
)
net_ml <- build_network(obs, method = "cor",
                         params = list(id = "id"), level = "both")
print(net_ml)
#> Multilevel Network (method: cor)
#> -- Between-person --
#>   Nodes: 4  |  Edges: 6
#>   Sample size: 5 (unique persons)
#> -- Within-person --
#>   Nodes: 4  |  Edges: 6
#>   Sample size: 40 (observations)
# }
```
