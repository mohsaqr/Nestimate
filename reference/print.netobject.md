# Print Method for Network Object

Print Method for Network Object

## Usage

``` r
# S3 method for class 'netobject'
print(x, ...)
```

## Arguments

- x:

  A `netobject`.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
net <- build_network(seqs, method = "relative")
print(net)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [1.000, 1.000]  |  mean: 1.000
#> 
#>   Weight matrix:
#>     A B C
#>   A 0 1 0
#>   B 0 0 1
#>   C 1 0 0 
#> 
#>   Initial probabilities:
#>   A             0.500  ████████████████████████████████████████
#>   B             0.250  ████████████████████
#>   C             0.250  ████████████████████
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B")
)
net <- build_network(seqs, method = "relative")
print(net)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [1.000, 1.000]  |  mean: 1.000
#> 
#>   Weight matrix:
#>     A B C
#>   A 0 1 0
#>   B 0 0 1
#>   C 1 0 0 
#> 
#>   Initial probabilities:
#>   A             0.500  ████████████████████████████████████████
#>   B             0.250  ████████████████████
#>   C             0.250  ████████████████████
# }
```
