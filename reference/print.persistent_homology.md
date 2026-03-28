# Print persistent homology results

Print persistent homology results

## Usage

``` r
# S3 method for class 'persistent_homology'
print(x, ...)
```

## Arguments

- x:

  A `persistent_homology` object.

- ...:

  Additional arguments (unused).

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
net <- build_network(seqs, method = "relative")
ph  <- persistent_homology(net)
print(ph)
#> Persistent Homology
#>   20 filtration steps [1.0000 → 0.0100]
#>   Features: b0: 3 (1 persistent) 
#>   Longest-lived:
#>     b0: 1.0000 → 0.0000 (life: 1.0000)
#>     b0: 1.0000 → 0.9479 (life: 0.0521)
#>     b0: 1.0000 → 0.9479 (life: 0.0521)
# }
```
