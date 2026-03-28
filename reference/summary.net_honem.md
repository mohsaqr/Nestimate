# Summary Method for net_honem

Summary Method for net_honem

## Usage

``` r
# S3 method for class 'net_honem'
summary(object, ...)
```

## Arguments

- object:

  A `net_honem` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon <- build_hon(seqs, max_order = 3)
he <- build_honem(hon, dim = 2)
summary(he)
#> HONEM Summary
#> 
#>   Nodes: 4 | Dimensions: 2
#>   Variance explained: 78.2%
#>   Top singular values: 1.023, 0.6
#>   Embedding range: [-0.573, 0.430]
# }
```
