# Edge-weight Case-dropping Stability

Computes a **CS-coefficient for the edge-weight vector** of a network:
the maximum proportion of cases (rows of `x$data`) that can be dropped
while the flattened edge-weight vector of the re-estimated network still
correlates with the original above `threshold` in at least `certainty`
of iterations.

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
  seed = NULL
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

  A `netobject`, `cograph_network`, `netobject_group`, or `mcml`. For
  the group types this function iterates over each element and returns a
  named list. For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `net_casedrop_reliability` or
  `net_casedrop_reliability_group` (or its
  [`summary()`](https://rdrr.io/r/base/summary.html)).

- iter:

  Integer. Iterations per drop proportion. Default `1000`.

- drop_prop:

  Numeric vector of proportions to evaluate. Each entry must lie
  strictly between 0 and 1. Default `seq(0.1, 0.9, by = 0.1)`. In
  `summary.net_casedrop_reliability_group()`: Drop proportion at which
  to report the four metrics (mean +/- sd per network). Must be one of
  the drop proportions the object was built with. Defaults to the
  object's median grid value (the stored grid is used, not an assumed
  `0.7`); pass an explicit value not in the grid to get an error listing
  the available proportions.

- threshold:

  Numeric in `[0, 1]`. Minimum edge-vector correlation for an iteration
  to count as stable. Default `0.7`.

- certainty:

  Numeric in `[0, 1]`. Required fraction of iterations whose correlation
  must exceed `threshold` for a drop proportion to qualify. Default
  `0.95`.

- method:

  Correlation method: `"pearson"` (weight magnitudes), `"spearman"`
  (ranks, robust to scale), or `"kendall"`. Default `"spearman"` because
  edge weights often span several orders of magnitude and rank stability
  is the typical target.

- include_diag:

  Logical. Include diagonal (self-loop) edges in the edge vector.
  Default `FALSE`.

- seed:

  Optional integer for reproducibility.

- digits:

  Digits to display. Default `3`.

- ...:

  For the S3 methods: further arguments passed to or from other methods.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_casedrop_reliability` or
  `net_casedrop_reliability_group`.

- combined:

  When `TRUE` (default), all four metrics are shown in one ggplot via
  `facet_wrap(~ metric)`. When `FALSE`, returns a named list of four
  single-panel ggplots, one per metric.

- metric:

  Which metric to plot. One of `"correlation"` (default),
  `"mean_abs_dev"`, `"median_abs_dev"`, `"max_abs_dev"`.

## Value

An object of class `net_casedrop_reliability` with:

- `cs`:

  Scalar CS-coefficient - the maximum drop proportion for which the
  edge-vector correlation remains \>= `threshold` in at least
  `certainty` of iterations. Zero if no proportion qualifies.

- `summary`:

  Tidy data frame, one row per metric by drop proportion, with columns
  `metric` (`"mean_abs_dev"`, `"median_abs_dev"`, `"correlation"`,
  `"max_abs_dev"`), `drop_prop`, `mean`, `sd`, `median`, `mad`, `q025`,
  `q975`.

- `metrics`:

  Named list of four `iter` x `length(drop_prop)` matrices, one per
  metric, holding the raw per-iteration values.

- `correlations`:

  `iter` x `length(drop_prop)` matrix of per- iteration correlations
  (the `correlation` entry of `metrics`).

- `drop_prop`, `threshold`, `certainty`, `iter`, `method`,
  `include_diag`:

  Inputs.

- `n_cases`:

  Number of cases resampled from (sequences for transition methods, rows
  of `$data` otherwise).

- `n_edges`:

  Length of the edge vector assessed.

A `netobject_group` or `mcml` input instead returns a
`net_casedrop_reliability_group`: a named list of one result per
constituent network.

When the original edge vector has zero variance a warning is issued and
the object is returned with `cs = 0`, an empty `summary`, and all-`NA`
metric matrices.

In `print.net_casedrop_reliability()`: The input `x` invisibly.

In `summary.net_casedrop_reliability()`: A tidy data frame with columns
`metric`, `drop_prop`, `mean`, `sd` summarising edge-weight stability
across case-dropping iterations.

In `summary.net_casedrop_reliability_group()`: A data frame with one row
per network containing `cor`, `mean_abs_dev`, `median_abs_dev`,
`max_abs_dev` formatted as "mean +/- sd".

In `plot.net_casedrop_reliability()`: A `ggplot` object, or a named list
of four ggplots when `combined = FALSE`.

In `plot.net_casedrop_reliability_group()`: A `ggplot` object.

## Details

Complements
[`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md):
that function asks whether centrality *rankings* are stable; this one
asks whether the *edge-weight structure itself* is stable. For
MCML-derived networks where each row of `$data` is one transition, this
is case-dropping of **edges**.

For each `drop_prop` p and each iteration, a size `n_cases * (1 - p)`
subset of `$data` rows is selected **without replacement**, the network
is re-estimated using the same method/scaling/threshold as the input,
and the upper/lower-triangle (directed: all off-diagonal entries) of the
new weight matrix is flattened and correlated with the corresponding
vector of the original matrix. The correlation method defaults to
Spearman for robustness to the wide dynamic range of transition
probabilities.

Unlike bootstrap CIs, case-dropping does not estimate sampling variance
and so does not rely on the i.i.d. assumption. This makes it the
appropriate robustness check for **edgelist-derived** networks (where
rows of `$data` lack actor grouping), since dropping rows at random is a
well-posed operation regardless of within-actor correlation.

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
