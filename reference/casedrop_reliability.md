# Edge-weight Case-dropping Stability

`casedrop_reliability()` measures how stable a network's edge weights
are when cases are left out. It re-estimates the network on random
subsets that drop each proportion in `drop_prop` of the cases, compares
each re-estimated edge vector with the original, and reports a
CS-coefficient for the edge weights: the largest proportion of cases
that can be dropped while the correlation between the original and the
re-estimated edge weights stays at or above `threshold` in at least
`certainty` of the iterations (the case-dropping CS-coefficient of
Epskamp, Borsboom & Fried, 2018, applied to the edge weights). Beside
the correlation it records the mean, median and largest absolute change
of the edge weights. The cases are the sequences of a transition network
and the rows of the data otherwise. It takes a single network, or a
`netobject_group` or an `mcml`, whose networks are assessed one by one.

It complements
[`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md):
that asks whether centrality *rankings* are stable, this asks whether
the *edge-weight structure itself* is. For an `mcml` built from an edge
list, where each row of the data is one transition, this is
case-dropping of **edges**.

The CS-coefficient is the largest proportion in `drop_prop` that passes,
so it can only be one of the values tried. A smaller proportion need not
have passed as well, so look at the correlation at every proportion, not
only at the CS-coefficient:
[`summary()`](https://rdrr.io/r/base/summary.html) gives its mean,
median and quantiles (`q025`, `q975`) per proportion, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it. A
correlation near 1 can hide a large change in a few edges, which the
largest absolute change (`max_abs_dev`) shows. Stability under
case-dropping says the estimate does not rest on a few cases; it gives
no interval for a single edge, which
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
does.

## Usage

``` r
casedrop_reliability(
  x,
  iter = 1000L,
  drop_prop = seq(0.1, 0.9, by = 0.1),
  threshold = 0.7,
  certainty = 0.95,
  method = c("spearman", "pearson", "kendall"),
  include_diag = FALSE,
  seed = NULL,
  target = c("auto", "existing", "all")
)

# S3 method for class 'net_casedrop_reliability'
print(x, digits = 3, ...)

# S3 method for class 'net_casedrop_reliability'
summary(object, ...)

# S3 method for class 'net_casedrop_reliability_group'
print(x, ...)

# S3 method for class 'net_casedrop_reliability_group'
summary(object, drop_prop = NULL, ...)

# S3 method for class 'summary.net_casedrop_reliability_group'
print(x, ...)

# S3 method for class 'net_casedrop_reliability'
plot(x, combined = TRUE, ...)

# S3 method for class 'net_casedrop_reliability_group'
plot(
  x,
  metric = c("correlation", "mean_abs_dev", "median_abs_dev", "max_abs_dev"),
  ...
)
```

## Arguments

- x:

  The network: a `netobject` or `cograph_network` that carries its data,
  or a `netobject_group` or `mcml`, whose networks (for an `mcml`, the
  macro network and each cluster's) are assessed one by one. A network
  without data raises an error, and one that rests on a single sequence,
  or one-hot data with no actor or session, has no cases to drop and
  raises an error of class `nestimate_single_sequence`. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: the
  result of `casedrop_reliability()` (for
  [`print()`](https://rdrr.io/r/base/print.html), also the
  [`summary()`](https://rdrr.io/r/base/summary.html) of a group result).

- iter:

  The number of random subsets drawn for each drop proportion, at
  least 2. Default `1000`.

- drop_prop:

  The proportions of cases to drop, each strictly between 0 and 1, in
  any order. Default `seq(0.1, 0.9, by = 0.1)`. A proportion that drops
  no case of a small sample (`n * p < 1`) is not run (without a warning;
  its metrics are `NA`) and counts as not passing. For the
  [`summary()`](https://rdrr.io/r/base/summary.html) method of a group
  result: the one proportion to report, which must be one of those the
  result was built with; a value outside them raises an error listing
  the available ones. Default `NULL` there, the middle proportion of the
  grid (the lower of the two middle ones for an even number).

- threshold:

  The lowest correlation between the original and the re-estimated edge
  weights for an iteration to count as stable, between 0 and 1. Default
  `0.7`.

- certainty:

  The share of iterations that must reach `threshold` for a drop
  proportion to pass, between 0 and 1. Default `0.95`. It is the share
  of all iterations: one in which the correlation is undefined (the edge
  vector is constant in that subset), or in which the network could not
  be re-estimated, counts as not reaching `threshold`.

- method:

  The correlation: `"spearman"` (default) compares ranks, so it is not
  affected by the scale of the weights, which often span several orders
  of magnitude; `"pearson"` compares the weight magnitudes; `"kendall"`
  is Kendall's tau.

- include_diag:

  Whether self-loops are part of the edge vector. Default `FALSE`.

- seed:

  A number that makes the subsets reproducible. The session's own random
  stream is left as it was. For a `netobject_group` or `mcml`, each
  network is assessed from the same seed. Default `NULL`.

- target:

  The cells compared. `"existing"` compares only the cells that are
  non-zero in the original network; `"all"` compares every cell, zeros
  included, so edges that appear once cases are dropped count. `"auto"`
  (default) uses `"existing"` for unscaled and unthresholded transition
  networks, where an absent transition cannot appear in a subset, and
  `"all"` otherwise. `include_diag` applies to either.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) method of a
  single result: the decimal places of the metric table. Default `3`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) methods:
  the result of `casedrop_reliability()`.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method of a single result: `TRUE` (default) draws the four metrics as
  panels of one plot; `FALSE` returns a named list of four single-panel
  plots, one per metric.

- metric:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method of a group result: the metric drawn, `"correlation"` (default),
  `"mean_abs_dev"`, `"median_abs_dev"` or `"max_abs_dev"`.

## Value

A `net_casedrop_reliability` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the CS-coefficient
and the mean and standard deviation of each metric at each drop
proportion; [`summary()`](https://rdrr.io/r/base/summary.html) returns
that table in full, one row per metric and drop proportion;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
metrics across the drop proportions. A `netobject_group` or `mcml` gives
a `net_casedrop_reliability_group`, a named list of one result per
network: its [`print()`](https://rdrr.io/r/base/print.html) gives each
network's number of cases, number of edges and CS-coefficient, its
[`summary()`](https://rdrr.io/r/base/summary.html) one row per network
at one drop proportion, and its
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) one curve per
network. The fields of a single result:

- `cs`:

  The CS-coefficient: the largest drop proportion at which the
  correlation reaches `threshold` in at least `certainty` of the
  iterations. `0` when no proportion passes; `NA` when the original edge
  weights do not vary (constant up to rounding error, relative to their
  size), so the correlation is undefined.

- `summary`:

  A data frame with one row per metric and drop proportion: `metric`
  (`"mean_abs_dev"`, `"median_abs_dev"`, `"correlation"`,
  `"max_abs_dev"`), `drop_prop`, and the `mean`, `sd`, `median`, `mad`,
  `q025` and `q975` of the metric over the iterations.

- `metrics`:

  The per-iteration values: a named list of four matrices, one per
  metric, with one row per iteration and one column per drop proportion.

- `correlations`:

  The same matrix as the `correlation` entry of `metrics`.

- `drop_prop`, `threshold`, `certainty`, `iter`, `method`,
  `include_diag`:

  The settings used.

- `n_cases`:

  The number of cases resampled from: sequences for a transition
  network, rows of the data otherwise.

- `n_edges`:

  The number of cells compared, set by `target` and `include_diag`.

The resolved `target` is kept as the result's `target` attribute.

When the original edge weights do not vary, a warning is given and the
result has `cs = NA`, a `summary` whose statistics are `NA`, and metric
matrices of `NA`.

The [`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly. The [`summary()`](https://rdrr.io/r/base/summary.html)
method of a single result returns its `summary` field, described above.
The [`summary()`](https://rdrr.io/r/base/summary.html) method of a group
result returns a `summary.net_casedrop_reliability_group` data frame
with one row per network (named in the row names): `n_edges` and the
four metrics (`cor`, `mean_abs_dev`, `median_abs_dev`, `max_abs_dev`) at
the chosen drop proportion, each written as `"mean +/- sd"`; its
[`print()`](https://rdrr.io/r/base/print.html) names that proportion.
The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
return a `ggplot`; with `combined = FALSE`, the named list of plots is
returned invisibly. A single result is drawn as one panel per metric,
the mean across iterations as a line with a band of plus and minus one
standard deviation; a group result draws the chosen metric with one line
and band per network. On the correlation, a dashed line marks
`threshold`. A result whose CS-coefficient is undefined is drawn as a
panel that says so, and in a group such networks are left out and named
in the caption.

## Details

For each drop proportion \\p\\ and each iteration, a subset of \\n -
\lfloor np \rfloor\\ of the \\n\\ cases is drawn **without
replacement**, the network is re-estimated with the input's method,
scaling and threshold, and the compared cells of the new weight matrix
are set against the same cells of the original. For transition networks
the per-sequence transition counts are computed once and summed over
each subset; other networks are re-estimated from the subset of rows. A
re-estimation that fails, or that loses a node (a variable constant in
the subset), leaves the iteration without values; such iterations are
counted in a warning of class `nestimate_resample_failed` and count as
below `threshold`.

With `target = "auto"`, unscaled and unthresholded transition networks
compare only the existing edges: an absent transition cannot appear in a
subset, and the fixed zeros would inflate the correlation. Association
networks and scaled or thresholded networks compare all cells, because a
zero original weight can become non-zero once cases are dropped; set
`target = "existing"` to compare only the original edges there.

Case-dropping drops cases at random, so it treats them as independent
units. When the cases are nested in actors (several sessions per
person), the remaining cases of an actor stand in for the dropped ones,
and the CS-coefficient overstates the stability; drop whole actors
instead (for example with `tna::estimate_cs(actor = )`).

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods* 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## See also

[`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md),
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Examples

``` r
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE),
  V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
es  <- casedrop_reliability(net, iter = 50, drop_prop = c(0.1, 0.3, 0.5),
                      seed = 1)
print(es)
#> Edge-weight Case-dropping Stability
#>   Cases (rows of $data) : 30
#>   Edges assessed        : 12 (diagonal excluded)
#>   Iterations / prop     : 50
#>   Correlation method    : spearman
#>   CS-coefficient (r)    : 0.10  (threshold=0.70, certainty=0.95)
#> 
#> Model-level reliability across iterations (mean +/- sd per drop):
#>   drop_prop      p=0.1        p=0.3        p=0.5      
#>   mean|diff|      0.028+- 0.007   0.062+- 0.014   0.095+- 0.018
#>   MAD             0.021+- 0.008   0.052+- 0.014   0.085+- 0.023
#>   cor             0.927+- 0.045   0.820+- 0.078   0.688+- 0.161
#>   max|diff|       0.081+- 0.023   0.160+- 0.051   0.241+- 0.056
```
