# Plot Method for net_cluster_diagnostics

Delegates to the original clustering object's plot method
([`plot.net_clustering`](https://mohsaqr.github.io/Nestimate/reference/plot.net_clustering.md)
for distance-based diagnostics,
[`plot.net_mmm_clustering`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mmm_clustering.md)
or
[`plot.net_mmm`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mmm.md)
for model-based). The diagnostics object itself stores no plot geometry
– it just keeps a reference to the source so the existing visual layer
is reused.

## Usage

``` r
# S3 method for class 'net_cluster_diagnostics'
plot(x, type = NULL, ...)
```

## Arguments

- x:

  A `net_cluster_diagnostics` object.

- type:

  Character. Forwarded to the underlying plot method. Valid values for
  distance: `"silhouette"` (default), `"mds"`, `"heatmap"`,
  `"predictors"`. Valid values for mmm: `"posterior"` (default),
  `"covariates"` / `"predictors"`.

- ...:

  Forwarded to the underlying plot method.

## Value

A `ggplot` object, invisibly.
