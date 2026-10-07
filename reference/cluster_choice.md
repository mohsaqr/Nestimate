# Cluster Choice – sweep k, dissimilarity and method

One-call sweep across any combination of k, dissimilarity metric, and
clustering algorithm for distance-based sequence clustering. Mirrors
[`compare_mmm`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
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

# S3 method for class 'cluster_choice'
print(x, digits = 3L, ...)

# S3 method for class 'cluster_choice'
summary(object, ...)

# S3 method for class 'cluster_choice'
plot(
  x,
  type = c("auto", "lines", "bars", "heatmap", "tradeoff", "facet"),
  abbrev = FALSE,
  combined = TRUE,
  ...
)
```

## Arguments

- data:

  Sequence data (data frame or matrix) – forwarded to
  [`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md).

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
  [`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  (`weighted`, `lambda`, `q`, `p`, `seed`, `na_syms`, `covariates`,
  `estimator`). Note: `weighted = TRUE` only works with
  `dissimilarity = "hamming"` and is rejected up-front when sweeping
  mixed dissimilarities. In `plot.cluster_choice()`,
  `print.cluster_choice()` and `summary.cluster_choice()`: Unsupported.
  Supplying unused arguments raises an error.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `cluster_choice`.

- digits:

  Integer. Decimal places for floating-point columns. Default `3L`.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `cluster_choice`.

- type:

  Character. One of `"auto"` (default), `"lines"`, `"bars"`,
  `"heatmap"`, `"tradeoff"`, `"facet"`.

- abbrev:

  Logical. If `TRUE`, dissimilarity and method names shown on tick
  labels and point labels are shortened (e.g. `"hamming"` -\> `"ham"`,
  `"ward.D2"` -\> `"wD2"`). The legend shows the full canonical name.
  Default `FALSE`.

- combined:

  Only meaningful for `type = "facet"`. When `TRUE` (default), all
  methods are shown in one ggplot via `facet_wrap(~ method)`. When
  `FALSE`, returns a named list of single-panel ggplots, one per method.

## Value

A `cluster_choice` object (a data.frame subclass) with one row per (k,
dissimilarity, method) combination and columns:

- k, dissimilarity, method:

  The configuration for that row.

- silhouette:

  Overall average silhouette width (from
  [`cluster::silhouette`](https://rdrr.io/pkg/cluster/man/silhouette.html),
  computed inside
  [`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)).

- mean_within_dist:

  Size-weighted mean of within-cluster distances, in the units of the
  row's dissimilarity.

- min_size, max_size, size_ratio:

  Cluster-size balance bounds and their ratio (`max / min`).

In `print.cluster_choice()`: The input object, invisibly.

In `summary.cluster_choice()`: A data frame with the swept
configurations, all metrics, and a `best` character column flagging the
silhouette-max row.

In `plot.cluster_choice()`: A `ggplot` object, invisibly; for
`type = "facet"` with `combined = FALSE`, a named list of ggplots.

## Methods

- `plot.cluster_choice()`: Six explicit chart types plus a smart
  `"auto"` default. The user picks the shape; the function does not
  editorialise (no "best" annotation, no interpretive subtitles, no
  inferred recommendation).

## Plot types

Type cheat-sheet:

- `"auto"`:

  Default. Picks one of the others based on which axes were swept.
  k-only -\> `"lines"`; one categorical axis swept -\> `"bars"`; k plus
  one categorical -\> `"lines"`; k plus two categoricals -\> `"facet"`;
  both categoricals without k -\> `"heatmap"`.

- `"lines"`:

  Silhouette across k (and `mean_within_dist` when `k` is the only swept
  axis), one line per non-k axis when present.

- `"bars"`:

  Horizontal bar chart of silhouette per axis level. Bars sorted by
  silhouette.

- `"heatmap"`:

  Tiled silhouette across two categorical axes. Requires both
  `dissimilarity` and `method` swept.

- `"tradeoff"`:

  Scatter: silhouette (y) vs `size_ratio` (x). Works for any sweep;
  labels each point.

- `"facet"`:

  Lines vs k, colour by one categorical axis, facet by another. Requires
  `k` plus two categoricals.

Asking for a type the data can't support raises an error pointing at the
alternatives.

## See also

[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
[`compare_mmm`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
for the model-based equivalent,
[`cluster_diagnostics`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
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
