# Cluster Diagnostics

`cluster_diagnostics()` collects the quality measures of a clustering in
one `net_cluster_diagnostics` object, laid out the same way whichever
route produced the clusters: a distance clustering from
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
a mixed Markov model from
[`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
or
[`cluster_mmm()`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md),
the `netobject_group` that
[`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
or
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
builds from one of them, or the `net_mmm_clustering` attribute of such a
group. It packages what the fit already computed and adds only the
per-cluster summaries the fit does not store (silhouette means and mean
within-cluster distances for a distance clustering, the classification
error within each cluster for a mixture); it does not refit.

For a distance clustering it reports the average silhouette width,
overall and per cluster, and each cluster's mean within-cluster
distance. For a mixture it reports each cluster's average posterior
probability (AvePP) and mixing proportion, the overall AvePP, the
normalised entropy of the posterior probabilities (0 when every sequence
is assigned with certainty, 1 when all assignments are equally likely;
the reverse of the relative entropy often reported for latent class
models), the classification error (the expected share of sequences
misclassified when each is assigned to its most likely cluster: one
minus the largest posterior, averaged over all sequences and over each
cluster's members, so it equals one minus the AvePP; Vermunt & Magidson,
2005), and the log-likelihood with BIC, AIC and ICL.

These measures describe how cleanly one solution separates the
sequences; none of them picks the number of clusters. Compare solutions
with
[`cluster_choice()`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
for distance clustering or
[`compare_mmm()`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
for mixtures, and build the clusters' networks with
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
on the clustering to see how they differ, with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
for how stable their edges are.

## Usage

``` r
cluster_diagnostics(x, ...)

# S3 method for class 'net_cluster_diagnostics'
print(x, digits = 3L, ...)

# S3 method for class 'net_cluster_diagnostics'
plot(x, type = NULL, ...)

# S3 method for class 'net_cluster_diagnostics'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  The clustering: a `net_clustering`, a `net_mmm`, a `netobject_group`
  carrying a clustering (built by
  [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  or by
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  from a `net_clustering` or `net_mmm`), or a `net_mmm_clustering`. A
  `netobject_group` without a clustering, and any other object, raise an
  error. For the [`print()`](https://rdrr.io/r/base/print.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods: a `net_cluster_diagnostics`.

- ...:

  Not used: an extra argument raises an error, here and in the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods. For [`plot()`](https://rdrr.io/r/graphics/plot.default.html):
  further arguments for the clustering's own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method, such
  as `combined`.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html): decimal places of
  the silhouette widths, distances, AvePP, entropy and information
  criteria, a non-negative whole number. Default `3`. Percentages are
  always shown to one decimal.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  view, passed to the clustering's own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.
  Default `NULL`, which draws `"silhouette"` for a distance clustering
  and `"posterior"` for a mixture. A distance clustering also takes
  `"mds"`, `"heatmap"` and `"predictors"` (see
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md));
  a mixture takes `"covariates"`, with `"predictors"` accepted for it
  when the source is a `net_mmm_clustering` (see
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)).

- row.names, optional:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
  accepted for the generic and ignored.

## Value

A `net_cluster_diagnostics` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the numbers of
sequences and clusters (with the method and dissimilarity of a distance
clustering, or the number of states of a mixture), the overall quality
measures, the information criteria of a mixture, one row per cluster,
and the covariates of the fit, if any;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
clustering's own views. The object's fields:

- family:

  `"distance"` or `"mmm"`.

- k, n, sizes:

  Number of clusters, number of sequences, and the cluster sizes.

- per_cluster:

  Data frame with one row per cluster; its columns depend on the family.
  Distance: `cluster`, `size`, `pct`, `mean_within_dist`, `sil_mean`.
  Mixture: `cluster`, `size`, `pct`, `mix_pct`, `avepp`,
  `class_err_pct`. `pct`, `mix_pct` and `class_err_pct` are in percent.

- overall:

  List of the family's summary measures: `silhouette` for distance;
  `avepp_overall`, `entropy` and `classification_error` (a proportion)
  for a mixture.

- ics:

  For a mixture, a list with `BIC`, `AIC`, `ICL` and `log_likelihood`;
  NULL for distance.

- metadata:

  The settings of the fit. Distance: `method`, `dissimilarity`,
  `weighted`, `lambda`, `medoids`, `covariates`. Mixture: `states`,
  `n_states`, `iterations`, `converged`, `covariates`.

- source:

  The clustering object itself, which
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws without
  recomputing anything.

For [`print()`](https://rdrr.io/r/base/print.html): the input,
invisibly. For [`plot()`](https://rdrr.io/r/graphics/plot.default.html):
what the clustering's own
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method returns,
a `ggplot` object (invisibly) or, for the covariate views with
`combined = FALSE`, a list of `ggplot` objects named by cluster
(invisibly). For
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): the
per-cluster table described under `per_cluster`.

## References

Vermunt, J. K., & Magidson, J. (2005). *Latent GOLD 4.0 User's Guide*.
Statistical Innovations.

## See also

[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
[`build_mmm`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md),
[`cluster_choice`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md),
[`compare_mmm`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
for k-sweep model selection (MMM only).

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
cl <- build_clusters(seqs, k = 2, method = "ward.D2")
cluster_diagnostics(cl)
#> Cluster Diagnostics (distance) [ward.D2 / hamming]
#>   Sequences: 30  |  Clusters: 2
#>   Quality: silhouette = 0.463
#> 
#>   Cluster  N           Mean within-dist  Silhouette
#>   1        14 (46.7%)  1.132             0.283
#>   2        16 (53.3%)  0.608             0.620
# \donttest{
fit <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 20, seed = 1)
cluster_diagnostics(fit)
#> Cluster Diagnostics (mmm) [k = 2]
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   Quality: AvePP = 0.946  |  Entropy = 0.287  |  Class.Err = 5.4%
#>   ICs: LL = -56.481  |  BIC = 170.782  |  AIC = 146.962  |  ICL = 174.190
#> 
#>   Cluster  N           Mix%   AvePP  Class.Err%
#>   1        21 (70.0%)  67.5%  0.943   5.7%
#>   2        9 (30.0%)   32.5%  0.952   4.8%
as.data.frame(cluster_diagnostics(fit))
#>   cluster size pct  mix_pct     avepp class_err_pct
#> 1       1   21  70 67.48869 0.9428695      5.713045
#> 2       2    9  30 32.51131 0.9519734      4.802660
# }
```
