# Summary Method for mcml

Summary Method for mcml

## Usage

``` r
# S3 method for class 'mcml'
summary(object, ...)
```

## Arguments

- object:

  An `mcml` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  T1 = c("A","B","A"), T2 = c("B","C","B"),
  T3 = c("C","A","C"), T4 = c("A","B","A")
)
clusters <- c("Alpha", "Beta", "Alpha")
cs <- build_mcml(seqs, clusters, type = "raw")
summary(cs)
#> MCML Network
#> ============
#> Type: raw  | Method: sum 
#> Nodes: 3  | Clusters: 2 
#> Transitions: 9 
#>   Macro: 6  | Per-cluster: 3 
#> 
#> Clusters:
#>   Alpha (2): A, C
#>   Beta (1): B
#> 
#> Macro (cluster-level) weights:
#>       Alpha Beta
#> Alpha     3    3
#> Beta      3    0
# }
```
