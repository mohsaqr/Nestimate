# Print Method for net_honem

Print Method for net_honem

## Usage

``` r
# S3 method for class 'net_honem'
print(x, ...)
```

## Arguments

- x:

  A `net_honem` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
hon   <- build_hon(seqs, max_order = 2L)
honem <- build_honem(hon, dim = 2L)
print(honem)
#> HONEM: Higher-Order Network Embedding
#>   Nodes:      3
#>   Dimensions: 2
#>   Max power:  10
#>   Variance explained: 82.6%
# }
```
