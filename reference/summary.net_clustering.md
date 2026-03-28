# Summary Method for net_clustering

Summary Method for net_clustering

## Usage

``` r
# S3 method for class 'net_clustering'
summary(object, ...)
```

## Arguments

- object:

  A `net_clustering` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 20, TRUE),
  V2 = sample(c("A","B","C"), 20, TRUE),
  V3 = sample(c("A","B","C"), 20, TRUE)
)
cl <- cluster_data(seqs, k = 2)
summary(cl)
#> Sequence Clustering Summary
#>   Method:        pam 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.3407 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1   11         1.600000
#>        2    9         1.444444
# }
```
