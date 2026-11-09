# Compare Clusterings Across k, Dissimilarity and Method

`cluster_choice()` clusters the same sequences under every combination
of the cluster counts, dissimilarities and algorithms given, each with
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
and returns one row per combination with its average silhouette width,
mean within-cluster distance and cluster-size balance. It is the
distance-based counterpart of
[`compare_mmm()`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md),
which compares mixtures of Markov chains across k. `"all"` expands
`dissimilarity` or `method` to every supported option.

The silhouette of a sequence (Rousseeuw, 1987) compares its mean
distance to its own cluster with that to the nearest other cluster: near
1 it sits well inside its cluster, near 0 between two, and below 0 it is
closer to another cluster. Rows that share a dissimilarity compare k and
algorithms on one distance scale; across dissimilarities the silhouettes
and within-cluster distances rest on different distances and compare
less directly. The row with the highest silhouette is marked, but that
is not a decision: also weigh the cluster sizes (`min_size`,
`size_ratio`), as a high silhouette can come with a cluster of very few
sequences; then fit the chosen configuration with
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
and read it with its [`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md).

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

  The sequences, passed to
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  unchanged: any input it takes, such as a wide data frame of sequences
  or a network built from sequences.

- k:

  The cluster counts to try. Default `2:5`. Each must be a whole finite
  number, or an error is raised;
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  also requires each to be at least 2 and at most the number of
  sequences minus 1.

- dissimilarity:

  The dissimilarity measures to try. Default `"hamming"`. `"all"`
  expands to every supported measure: `"hamming"`, `"osa"`, `"lv"`,
  `"dl"`, `"lcs"`, `"qgram"`, `"cosine"`, `"jaccard"` and `"jw"`, and
  may stand beside named ones (each is fitted once).

- method:

  The clustering algorithms to try. Default `"ward.D2"` (where
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  defaults to `"pam"`). `"all"` expands to every supported algorithm:
  `"pam"`, `"ward.D2"`, `"ward.D"`, `"complete"`, `"average"`,
  `"single"`, `"mcquitty"`, `"median"` and `"centroid"`.

- ...:

  Further arguments of
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
  used for every combination: `weighted`, `lambda`, `q`, `p`, `seed`,
  `na_syms`, `covariates` and `estimator`; any other name raises an
  error. `weighted = TRUE` works only with `dissimilarity = "hamming"`
  and is refused before any fit when the sweep holds another measure.
  The [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods take
  no further arguments and raise an error on any.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `cluster_choice`.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: the
  decimal places of the silhouette, within-cluster distance and size
  ratio. Default `3`. Must be a single non-negative whole number.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `cluster_choice`.

- type:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the kind of plot, one of `"auto"` (default), `"lines"`,
  `"bars"`, `"heatmap"`, `"tradeoff"` and `"facet"` (see *Plot types*).
  A type the sweep cannot support raises an error that names the
  alternatives, and so does a table with a single configuration, which
  has nothing to compare.

- abbrev:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: whether to shorten the dissimilarity and method names (e.g.
  `"hamming"` to `"ham"`, `"ward.D2"` to `"wD2"`) on the axes, point
  labels and legend. Default `FALSE`.

- combined:

  For `plot(type = "facet")`: `TRUE` (default) draws one panel per
  method in a single plot; `FALSE` returns a named list of one-method
  plots instead. Ignored by the other types, but must still be a single
  logical value.

## Value

A `cluster_choice` (a data frame subclass) with one row per combination,
`k` varying fastest, then `dissimilarity`, then `method`.
[`print()`](https://rdrr.io/r/base/print.html) shows the table, with
only the settings that vary as columns, and flags the row with the
highest silhouette in a `best` column;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the full
table as a plain data frame with a character `best` column that holds
`"silhouette"` on that row;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
silhouette across the swept settings (see *Plot types*) and returns the
`ggplot` invisibly (with `type = "facet"` and `combined = FALSE`, a
named list of `ggplot`s, returned invisibly and not drawn). The
[`print()`](https://rdrr.io/r/base/print.html) method returns its input
invisibly. The columns:

- k, dissimilarity, method:

  The configuration of the row.

- silhouette:

  The average silhouette width over all sequences, computed by
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  (from [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) for
  `"pam"`, with
  [`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html)
  otherwise).

- mean_within_dist:

  The mean distance between members of the same cluster, averaged over
  clusters weighted by their size (a one-member cluster counts as 0), in
  the units of the row's dissimilarity.

- min_size, max_size, size_ratio:

  The sizes of the smallest and the largest cluster, and their ratio
  `max_size / min_size`.

## Plot types

The plots show the values only: they carry no "best" marker and no
recommendation. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
needs at least one swept setting.

- `"auto"`:

  The default. Picks a type from the swept settings: `"lines"` when k is
  swept with at most one other setting, `"bars"` when only
  `dissimilarity` or only `method` is swept, `"heatmap"` when both are
  swept at a single k, and `"facet"` when all three are swept.

- `"lines"`:

  The silhouette across k, one line per value of the other swept
  setting. When only k is swept, `mean_within_dist` is drawn as a second
  line on the same axis. With all three swept, the lines follow
  `dissimilarity` and the methods share them; use `"facet"` instead.
  Needs k swept.

- `"bars"`:

  Horizontal bars of the silhouette per value of the swept setting, or
  per combination when several are swept, sorted by silhouette, with the
  values printed.

- `"heatmap"`:

  Tiles of the silhouette over `dissimilarity` and `method`, with the
  values printed. Needs both swept, and is meant for a single k.

- `"tradeoff"`:

  The silhouette against `size_ratio`, each point labelled with its
  swept settings.

- `"facet"`:

  The silhouette across k, one line per dissimilarity and one panel per
  method. Needs all three swept.

## References

Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the
interpretation and validation of cluster analysis. *Journal of
Computational and Applied Mathematics*, 20, 53–65.
[doi:10.1016/0377-0427(87)90125-7](https://doi.org/10.1016/0377-0427%2887%2990125-7)

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
#>  2 0.403      0.995       [17, 23] 1.353         
#>  3 0.563      0.711       [11, 17] 1.545 <-- best
#>  4 0.541      0.553       [7, 12]  1.714         
# \donttest{
# Sweep dissimilarities at fixed k
cluster_choice(seqs, k = 3, dissimilarity = c("hamming", "lcs", "jaccard"))
#> Cluster Choice (sweep: dissimilarity)
#> 
#>  dissimilarity silhouette within_dist sizes    ratio best    
#>  hamming       0.563      0.711       [11, 17] 1.545 <-- best
#>  lcs           0.403      1.363       [9, 19]  2.111         
#>  jaccard       0.381      0.619       [5, 28]  5.600         

# Full grid of k x dissimilarity
cluster_choice(seqs, k = 2:4, dissimilarity = c("hamming", "lcs"))
#> Cluster Choice (sweep: k x dissimilarity)
#> 
#>  k dissimilarity silhouette within_dist sizes    ratio best    
#>  2 hamming       0.403      0.995       [17, 23] 1.353         
#>  3 hamming       0.563      0.711       [11, 17] 1.545 <-- best
#>  4 hamming       0.541      0.553       [7, 12]  1.714         
#>  2 lcs           0.368      1.709       [19, 21] 1.105         
#>  3 lcs           0.403      1.363       [9, 19]  2.111         
#>  4 lcs           0.519      1.044       [5, 19]  3.800         

# "all" sentinel
cluster_choice(seqs, k = 3, dissimilarity = "all")
#> Cluster Choice (sweep: dissimilarity)
#> 
#>  dissimilarity silhouette within_dist sizes    ratio best    
#>  hamming       0.563      0.711       [11, 17] 1.545         
#>  osa           0.463      0.714       [10, 17] 1.700         
#>  lv            0.563      0.711       [11, 17] 1.545         
#>  dl            0.463      0.714       [10, 17] 1.700         
#>  lcs           0.403      1.363       [9, 19]  2.111         
#>  qgram         0.381      1.237       [5, 28]  5.600         
#>  cosine        0.381      0.619       [5, 28]  5.600         
#>  jaccard       0.381      0.619       [5, 28]  5.600         
#>  jw            0.718      0.209       [12, 14] 1.167 <-- best
# }
```
