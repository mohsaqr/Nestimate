# Print Method for net_clustering

Print Method for net_clustering

## Usage

``` r
# S3 method for class 'net_clustering'
print(x, ...)
```

## Arguments

- x:

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
print(cl)
#> Sequence Clustering
#>   Method:        pam 
#>   Dissimilarity: hamming  
#>   Clusters:      2 
#>   Silhouette:    0.3407 
#>   Cluster sizes: 11, 9 
#>   Medoids:       14, 18 
# }
```
