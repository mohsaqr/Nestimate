# Cluster Choice – sweep k, dissimilarity and method

One-call sweep across any combination of k, dissimilarity metric, and
clustering algorithm for distance-based sequence clustering. Mirrors
[`compare_mmm`](https://mohsaqr.github.io/Nestimate/reference/compare_mmm.md)
for model-based clustering: returns a data frame with one row per swept
configuration, a `best` marker on the silhouette-max row in the print
method, and a [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
that adapts to the swept axes.

## Usage

``` r
cluster_choice(
  data,
  k = 2:5,
  dissimilarity = "hamming",
  method = "ward.D2",
  ...
)
```

## Arguments

- data:

  Sequence data (data frame or matrix) – forwarded to
  [`build_clusters`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md).

- k:

  Integer vector of cluster counts to sweep. Default `2:5`. Each value
  must be \>= 2 and \<= n - 1.

- dissimilarity:

  Character vector of dissimilarity metrics. Use `"all"` to expand to
  every supported metric:
  `c("hamming", "osa", "lv", "dl", "lcs", "qgram", "cosine",`
  `"jaccard", "jw")`. Default `"hamming"`.

- method:

  Character vector of clustering algorithms. Use `"all"` to expand to
  every supported method:
  `c("pam", "ward.D2", "ward.D", "complete", "average", "single",`
  `"mcquitty", "median", "centroid")`. Default `"ward.D2"`.

- ...:

  Other arguments forwarded to
  [`build_clusters`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md)
  (`weighted`, `lambda`, `q`, `p`, `seed`, `na_syms`, `covariates`).
  Note: `weighted = TRUE` only works with `dissimilarity = "hamming"`
  and is rejected up-front when sweeping mixed dissimilarities.

## Value

A `cluster_choice` object (a data.frame subclass) with one row per (k,
dissimilarity, method) combination and columns:

- k, dissimilarity, method:

  The configuration for that row.

- silhouette:

  Overall average silhouette width (from
  [`cluster::silhouette`](https://rdrr.io/pkg/cluster/man/silhouette.html),
  computed inside
  [`build_clusters`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md)).

- mean_within_dist:

  Size-weighted mean of within-cluster distances, in the units of the
  row's dissimilarity.

- min_size, max_size, size_ratio:

  Cluster-size balance bounds and their ratio (`max / min`).

## See also

[`build_clusters`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md),
[`compare_mmm`](https://mohsaqr.github.io/Nestimate/reference/compare_mmm.md)
for the model-based equivalent,
[`cluster_diagnostics`](https://mohsaqr.github.io/Nestimate/reference/cluster_diagnostics.md)
for the post-fit diagnostic surface on a single clustering.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 40, TRUE),
                   V2 = sample(c("A","B","C"), 40, TRUE))
cluster_choice(seqs, k = 2:4)
#> Cluster Choice (sweep: k)
#> 
#>  k silhouette within_dist sizes    ratio best    
#>  2 0.421      0.965       [18, 22] 1.222         
#>  3 0.561      0.697       [8, 18]  2.250 <-- best
#>  4 0.556      0.532       [7, 14]  2.000         
# \donttest{
# Sweep dissimilarities at fixed k
cluster_choice(seqs, k = 3, dissimilarity = c("hamming", "lcs", "jaccard"))
#> Cluster Choice (sweep: dissimilarity)
#> 
#>  dissimilarity silhouette within_dist sizes   ratio best    
#>  hamming       0.561      0.697       [8, 18] 2.250 <-- best
#>  lcs           0.434      1.254       [7, 17] 2.429         
#>  jaccard       0.413      0.587       [6, 27] 4.500         

# Full grid of k x dissimilarity
cluster_choice(seqs, k = 2:4, dissimilarity = c("hamming", "lcs"))
#> Cluster Choice (sweep: k x dissimilarity)
#> 
#>  k dissimilarity silhouette within_dist sizes    ratio best    
#>  2 hamming       0.421      0.965       [18, 22] 1.222         
#>  3 hamming       0.561      0.697       [8, 18]  2.250 <-- best
#>  4 hamming       0.556      0.532       [7, 14]  2.000         
#>  2 lcs           0.371      1.712       [17, 23] 1.353         
#>  3 lcs           0.434      1.254       [7, 17]  2.429         
#>  4 lcs           0.553      0.960       [6, 17]  2.833         

# "all" sentinel
cluster_choice(seqs, k = 3, dissimilarity = "all")
#> Cluster Choice (sweep: dissimilarity)
#> 
#>  dissimilarity silhouette within_dist sizes   ratio best    
#>  hamming       0.561      0.697       [8, 18] 2.250         
#>  osa           0.500      0.697       [8, 18] 2.250         
#>  lv            0.561      0.697       [8, 18] 2.250         
#>  dl            0.500      0.697       [8, 18] 2.250         
#>  lcs           0.434      1.254       [7, 17] 2.429         
#>  qgram         0.413      1.173       [6, 27] 4.500         
#>  cosine        0.413      0.587       [6, 27] 4.500         
#>  jaccard       0.413      0.587       [6, 27] 4.500         
#>  jw            0.711      0.209       [8, 18] 2.250 <-- best
# }
```
