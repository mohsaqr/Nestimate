# Plot Method for MMM Clustering Attribute

Plot routines for the MMM clustering metadata attached to the
`netobject_group` that
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
materializes from a
[`cluster_mmm`](https://saqr.me/Nestimate/reference/cluster_mmm.md) fit
(or that
[`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md)
returns directly with `cluster_by = "mmm"`). Mirrors the type-driven
surface of
[`plot.net_clustering`](https://saqr.me/Nestimate/reference/plot.net_clustering.md)
but covers only the metrics the EM fit produces – there is no distance
matrix on an MMM clustering, so `"silhouette"` / `"mds"` / `"heatmap"`
aren't defined here and the dispatcher raises a clear error if you ask
for one of those on an MMM result.

## Usage

``` r
# S3 method for class 'net_mmm_clustering'
plot(
  x,
  type = c("posterior", "covariates", "predictors"),
  combined = TRUE,
  ...
)
```

## Arguments

- x:

  A `net_mmm_clustering` object.

- type:

  Character. One of `"posterior"` (default; histogram of max posterior
  probability per sequence, coloured by cluster), `"covariates"` or its
  alias `"predictors"` (covariate forest plot when
  [`cluster_mmm()`](https://saqr.me/Nestimate/reference/cluster_mmm.md)
  was run with `covariates`).

- combined:

  Logical. For `type` in `"covariates"` or `"predictors"` only: when
  `TRUE` (default), forest panels are combined into a single faceted
  plot; when `FALSE`, a list of separate ggplots is returned.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A `ggplot` object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(V1 = sample(c("A","B","C"), 40, TRUE),
                   V2 = sample(c("A","B","C"), 40, TRUE))
fit <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 20, seed = 1)
grp <- build_network(fit)
plot(attr(grp, "clustering"), type = "posterior")

# }
```
