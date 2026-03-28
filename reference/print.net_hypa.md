# Print Method for net_hypa

Print Method for net_hypa

## Usage

``` r
# S3 method for class 'net_hypa'
print(x, ...)
```

## Arguments

- x:

  A `net_hypa` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B","C","A","B","C","A"),
  V2 = c("B","C","A","B","C","A","B","C","A","B"),
  V3 = c("C","A","B","C","A","B","C","A","B","C"),
  V4 = c("A","B","C","A","B","C","A","B","C","A")
)
hypa <- build_hypa(seqs, k = 2L)
print(hypa)
#> HYPA: Path Anomaly Detection
#>   Order k:      2
#>   Edges:        3
#>   Anomalous:    0 (alpha=0.05)
#>     Over-repr:  0
#>     Under-repr: 0
# }
```
