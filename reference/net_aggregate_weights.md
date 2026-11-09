# Aggregate Edge Weights

`net_aggregate_weights()` reduces a vector of edge weights to one
number: their sum, mean, median, maximum, minimum, product, density or
geometric mean. Missing and zero weights are removed first, so every
method works on the edges that are present.
[`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
and
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
use it to combine the edges between two clusters into one cluster-level
weight.

## Usage

``` r
net_aggregate_weights(w, method = "sum", n_possible = NULL)
```

## Arguments

- w:

  The edge weights, a numeric vector. `NA` is allowed and removed with
  the zeros; an infinite value raises an error.

- method:

  How to combine the weights: `"sum"` (default), `"mean"`, `"median"`,
  `"max"`, `"min"`, `"prod"`, `"density"` or `"geomean"`. A single
  value; any other raises an error. `"density"` divides the sum by
  `n_possible`, and `"mean"` by the number of non-zero weights, so the
  two agree when every possible edge is present, and `"density"` is the
  smaller when some are zero or missing. `"geomean"` is the geometric
  mean of the positive weights; negative weights are left out, and the
  result is `0` when none is positive.

- n_possible:

  For `method = "density"`: the number of possible edges, such as the
  block size `n_i * n_j` of two clusters, a single finite number.
  Default `NULL`, under which (as with a value of 0 or less) `"density"`
  returns the mean of the non-zero weights. Other methods ignore it.

## Value

A single number: the chosen aggregate of the non-zero, non-missing
weights, or `0` when none remain.

## Details

`"sum"`, `"prod"`, `"min"`, `"max"`, `"mean"` and `"median"` mean what
they mean in igraph's `edge.attr.comb` argument.

## Examples

``` r
w <- c(0.5, 0.8, 0.3, 0.9)
net_aggregate_weights(w, "sum")   # 2.5
#> [1] 2.5
net_aggregate_weights(w, "mean")  # 0.625
#> [1] 0.625
net_aggregate_weights(w, "max")   # 0.9
#> [1] 0.9
net_aggregate_weights(w, "density", n_possible = 9)  # 2.5 / 9
#> [1] 0.2777778
```
