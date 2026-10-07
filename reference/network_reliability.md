# Split-Half Reliability for Network Estimates

Assesses the stability of network estimates by repeatedly splitting
sequences into two halves, building networks from each half, and
comparing them. Supports single-model reliability assessment and
multi-model comparison with optional scaling for cross-method
comparability.

For transition methods (`"relative"`, `"frequency"`, `"co_occurrence"`),
uses pre-computed per-sequence count matrices for fast resampling (same
infrastructure as
[`bootstrap_network`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)).

## Usage

``` r
network_reliability(
  ...,
  iter = 1000L,
  split = 0.5,
  scale = "none",
  seed = NULL
)

# S3 method for class 'net_reliability'
print(x, ...)

# S3 method for class 'net_reliability'
summary(object, ...)

# S3 method for class 'net_reliability'
plot(x, bins = 60L, combined = TRUE, ...)
```

## Arguments

- ...:

  One or more `netobject`s (from
  [`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)).
  If unnamed, each model is auto-named from its `$method`; duplicate
  names are made unique with
  [`make.unique()`](https://rdrr.io/r/base/make.unique.html). A
  `netobject_group` is flattened into its constituent models (named by
  group), and an `mcml` or `cograph_network` is converted first. In
  `plot.net_reliability()` and `print.net_reliability()`: Additional
  arguments (ignored). In `summary.net_reliability()`: Ignored.

- iter:

  Integer. Number of split-half iterations (default: 1000).

- split:

  Numeric. Fraction of sequences assigned to the first half (default:
  0.5).

- scale:

  Character. Scaling applied to both split-half matrices before
  computing metrics. One of `"none"` (default), `"minmax"`,
  `"standardize"`, or `"proportion"`. Use scaling when comparing models
  on different scales (e.g. frequency vs relative).

- seed:

  Integer or NULL. RNG seed for reproducibility.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `net_reliability`.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_reliability`.

- bins:

  Integer. Number of histogram bins per panel (default 60).

- combined:

  When `TRUE` (default), all four metrics are shown in one ggplot via
  `facet_wrap(~ metric)`. When `FALSE`, returns a named list of four
  single-panel ggplots, one per metric.

## Value

An object of class `"net_reliability"` containing:

- iterations:

  Data frame with columns `model`, `mean_dev`, `median_dev`, `cor`,
  `max_dev` (one row per iteration per model).

- summary:

  Data frame with columns `model`, `metric`, `mean`, `sd`.

- models:

  Named list of the original `netobject`s.

- iter:

  Number of iterations.

- split:

  Split fraction.

- scale:

  Scaling method used.

In `print.net_reliability()`: The input object, invisibly.

In `summary.net_reliability()`: A tidy data frame with columns `model`,
`metric`, `mean`, `sd` summarising the split-half iterations.

In `plot.net_reliability()`: A `ggplot` object (invisibly), or a named
list of four ggplots when `combined = FALSE`.

## Methods

- `plot.net_reliability()`: Density plots of split-half metrics faceted
  by metric type. Multi-model comparisons show overlaid densities
  colored by model.

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`bootstrap_network`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
rel <- network_reliability(net, iter = 10)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
rel <- network_reliability(net, iter = 100, seed = 42)
print(rel)
#> Split-Half Reliability (100 iterations, split = 50%)
#>   Mean Abs. Diff.     mean = 0.1537  sd = 0.0300
#>   Median Abs. Diff.   mean = 0.1329  sd = 0.0322
#>   Pearson             mean = 0.2006  sd = 0.1882
#>   Max Abs. Diff.      mean = 0.3760  sd = 0.0972
# }
```
