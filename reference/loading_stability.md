# Composite-Weight Stability Under Case Resampling

`loading_stability()` bootstraps the item weights of a
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
fit and returns a `pc_loading_stability` object with each item's
observed signed weight, its bootstrap mean and standard deviation, a
percentile interval and how often its sign flipped. Each replicate
resamples the rows of the raw data with replacement, re-estimates the
node-level network with the fit's estimator (`cor`, `pcor` or `glasso`;
any other node-level method is re-estimated as `cor`) and correlation
type, and reapplies the fitted weighting scheme, data-dependent and
custom weights and the factor-analysis options included. A fit with a
descriptive aggregation (`"average"`, `"escoufier"`, `"cancor"`) has its
item-table weights, the `"strength"` weights, resampled.
**Experimental.**

Wide intervals mean the weighting should not be over-interpreted. Equal
weights keep their magnitudes in every replicate, so for them only the
sign flips carry information; a sign that flips often marks an item
whose reverse-keying is not settled by the data. Supplied numeric
weights are fixed in magnitude and sign, so their intervals have zero
width and resampling them tells nothing. The intervals describe the
weights, not the macro network: edge-weight intervals for it come from
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
on its macro network.

## Usage

``` r
loading_stability(x, iter = 200L, ci_level = 0.05, seed = NULL)

# S3 method for class 'pc_loading_stability'
print(x, digits = 3, ...)

# S3 method for class 'pc_loading_stability'
summary(object, ...)

# S3 method for class 'pc_loading_stability'
plot(x, ...)
```

## Arguments

- x:

  An `mcml_pc` object that carries its raw data; one built from a
  network without data raises an error, and so does an older saved fit
  with custom weights, which must be rebuilt with
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  to keep its weighting. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `pc_loading_stability` object.

- iter:

  The number of bootstrap replicates, at least 2. Default `200`; each
  replicate re-estimates the node-level network, so this is slower than
  a plain bootstrap.

- ci_level:

  The tail level of the percentile intervals: `0.05` (default) gives 95%
  intervals.

- seed:

  An integer that makes the resampling reproducible. Default `NULL`. The
  session's own random stream is left as it was.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: the
  decimal places of the numeric columns. Default `3`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `pc_loading_stability` object.

## Value

A `pc_loading_stability` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the item table and
[`summary()`](https://rdrr.io/r/base/summary.html) returns it;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws each
item's signed weight with its interval, one panel per cluster. The
object's fields:

- summary:

  A data frame with one row per item: `node`, `cluster`, `weight` (the
  observed signed weight), `boot_mean`, `boot_sd`, `ci_lower`,
  `ci_upper` and `sign_flips` (the share of replicates in which the
  item's sign differed from the observed one).

- boot_weights:

  The `iter` x items matrix of resampled signed weights.

- iter, ci_level:

  The settings used.

For the methods: [`print()`](https://rdrr.io/r/base/print.html) returns
`x` invisibly, [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
a `ggplot`.

For the [`summary()`](https://rdrr.io/r/base/summary.html) method: the
item table, a data frame with one row per item and the columns listed
under `summary` above.

## Methods

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws each
  item's observed signed composite weight with its bootstrap percentile
  interval, one panel per cluster, with a dashed line at zero.

## Examples

``` r
set.seed(1)
f <- stats::rnorm(100)
g <- stats::rnorm(100)
df <- data.frame(a1 = f + stats::rnorm(100), a2 = f + stats::rnorm(100),
                 a3 = f + stats::rnorm(100), b1 = g + stats::rnorm(100),
                 b2 = g + stats::rnorm(100), b3 = g + stats::rnorm(100))
cl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
fit <- build_mcml_pc(df, cl, aggregation = "loadings",
                     method = "cor")
stability <- loading_stability(fit, iter = 50, seed = 1)
stability
#> Composite-Weight Stability (case bootstrap, experimental)
#>   50 replicates | 95% percentile CIs
#> 
#>  node cluster weight boot_mean boot_sd ci_lower ci_upper sign_flips
#>    a1       A  0.338     0.332   0.021    0.292    0.370          0
#>    a2       A  0.336     0.339   0.018    0.308    0.375          0
#>    a3       A  0.326     0.329   0.021    0.290    0.367          0
#>    b1       B  0.316     0.317   0.022    0.276    0.348          0
#>    b2       B  0.350     0.346   0.022    0.310    0.393          0
#>    b3       B  0.334     0.337   0.015    0.307    0.360          0
```
