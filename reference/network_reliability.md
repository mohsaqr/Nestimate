# Split-Half Reliability for Network Estimates

`network_reliability()` measures how much a network's weights depend on
the particular sample: it splits the sequences (for an association
network, the rows of its data) at random into two halves, estimates the
network on each half with the original network's method and settings,
compares the two weight matrices, and repeats this `iter` times. Each
comparison gives four numbers: the mean, median and largest absolute
difference between the halves' weights and their Pearson correlation,
taken over all cells of the weight matrix. It takes one network or
several, so the reliability of different methods or groups can be set
side by side; a `netobject_group` contributes one network per group.

Each half holds only part of the sequences, so the agreement is that of
networks estimated from a sample of that size. Differences are on the
scale of the weights: counts of a `"frequency"` network and
probabilities of a `"relative"` one are not comparable without `scale`.
To see how far single edges move, use
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md);
to check centrality rankings,
[`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md).

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

  The networks: one or more `netobject`s from
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  A `netobject_group` is split into its networks, named by group; an
  `mcml` is converted with
  [`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
  first, and a bare `cograph_network` is read as a `netobject`. Each
  network is named by its argument name or, failing that, its method;
  repeated names get a suffix (`"_1"`, `"_2"`, ...). A network without
  its sequence data (which
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  keeps) raises an error, a network resting on one sequence an error of
  class `nestimate_single_sequence`, and networks of different methods
  with `scale = "none"` a warning of class `nestimate_mixed_methods`.
  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- iter:

  The number of split-half iterations, at least 2. Default `1000`.

- split:

  The share of the sequences that goes to the first half, the rest
  forming the second, strictly between 0 and 1. Default `0.5`.

- scale:

  The rescaling applied to both halves' weight matrices before they are
  compared, so that networks on different scales can be set side by
  side: `"none"` (default); `"minmax"`, which maps each matrix to 0–1;
  `"standardize"`, which subtracts the mean of its cells and divides by
  their SD; or `"proportion"`, which divides by the sum of its cells. A
  matrix whose cells are all equal, or sum to 0, is left as it is.

- seed:

  A seed for the random splits, so the result repeats exactly. The
  session's own random stream is left as it was. Default `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_reliability` object.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_reliability` object.

- bins:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the number of histogram bins per panel. Default `60`.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: `TRUE` (default) draws the four metrics as panels of one plot;
  `FALSE` returns a named list of four single-panel plots, one per
  metric, without drawing them.

## Value

A `net_reliability` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the mean and SD of
each metric over the iterations, per network;
[`summary()`](https://rdrr.io/r/base/summary.html) returns them as a
data frame with one row per network and metric and the columns `model`,
`metric` (`mean_dev`, `median_dev`, `cor`, `max_dev`), `mean` and `sd`;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws their
distributions. The object's fields:

- iterations:

  A data frame with one row per iteration and network, and the columns
  `model`, `mean_dev`, `median_dev`, `cor` and `max_dev`.

- summary:

  The table [`summary()`](https://rdrr.io/r/base/summary.html) returns.

- models:

  The networks compared, as a named list.

- iter, split, scale:

  The settings of the call.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws the histograms (see *Methods*) and returns the `ggplot` invisibly;
with `combined = FALSE` it returns, without drawing, a list of four
`ggplot`s named by the metric labels (`"Max Abs. Diff."`,
`"Mean Abs. Diff."`, `"Median Abs. Diff."`, `"Pearson"`).

## Details

The first half holds `round(split * n)` of the `n` sequences (at least
one), the second the rest; a `split` that leaves the second half empty
(few sequences, `split` far from 0.5) raises an error of class
`nestimate_bad_input`. For the transition methods (`"relative"`,
`"frequency"`, `"co_occurrence"`, and one-hot `wtna` networks, which are
counted the same way), each sequence's transition counts are computed
once and each half is summed from them, as in
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).
Other methods are re-estimated on each half of the rows. The network's
own scaling and threshold are applied to both halves. A half on which
the estimator fails, or that loses a variable (a column constant in that
half), gives `NA` for that iteration, counted in a warning of class
`nestimate_resample_failed`. The summary leaves out missing values, such
as a correlation that is undefined because one half's weights are all
equal.

## Methods

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): histograms
  of the split-half metrics over the iterations, one panel per metric,
  the mean marked by a dashed line and a label. With several networks,
  their histograms are overlaid, one colour per network, with a legend
  below.

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
