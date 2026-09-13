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

  Unsupported. Supplying unused arguments raises an error.

## Value

A data frame of per-cluster statistics, one row per cluster, with
columns `cluster`, `size` and `mean_within_dist`, returned *visibly*.
When the clustering was fitted with `covariates`, a
`tidy_covariates`/`data.frame` (the tidied covariate table, with cluster
sizes, fit statistics and profiles attached as attributes) is returned
*invisibly* instead. In both cases the printed summary is a side effect.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
cl <- build_clusters(seqs, k = 2)
summary(cl)
#> Sequence Clustering Summary
#>   Method:        pam 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.6 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1    2                0
#>        2    3                2
#>   cluster size mean_within_dist
#> 1       1    2                0
#> 2       2    3                2
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 20, TRUE),
  V2 = sample(c("A","B","C"), 20, TRUE),
  V3 = sample(c("A","B","C"), 20, TRUE)
)
cl <- build_clusters(seqs, k = 2)
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
#>   cluster size mean_within_dist
#> 1       1   11         1.600000
#> 2       2    9         1.444444
# }
```
