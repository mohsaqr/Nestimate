# Split-Half Reliability for Network Estimates

Assesses the stability of network estimates by repeatedly splitting
sequences into two halves, building networks from each half, and
comparing them. Supports single-model reliability assessment and
multi-model comparison with optional scaling for cross-method
comparability.

For transition methods (`"relative"`, `"frequency"`, `"co_occurrence"`),
uses pre-computed per-sequence count matrices for fast resampling (same
infrastructure as
[`bootstrap_network`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)).

## Usage

``` r
network_reliability(
  ...,
  iter = 1000L,
  split = 0.5,
  scale = "none",
  seed = NULL
)
```

## Arguments

- ...:

  One or more `netobject`s (from
  [`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)).
  If unnamed, each model is auto-named from its `$method`. A
  `netobject_group` is flattened into its constituent models.

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

## See also

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
[`bootstrap_network`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
rel <- network_reliability(net, iter = 10)
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
rel <- network_reliability(net, iter = 100, seed = 42)
print(rel)
#> Split-Half Reliability (100 iterations, split = 50%)
#>   Mean Abs. Diff.     mean = 0.1454  sd = 0.0297
#>   Median Abs. Diff.   mean = 0.1297  sd = 0.0378
#>   Pearson             mean = -0.0740  sd = 0.1979
#>   Max Abs. Diff.      mean = 0.3669  sd = 0.0929
# }
```
