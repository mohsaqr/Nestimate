# Compare Network-Level Statistics of Two Networks

`vertex_compare()` tests whether two observed networks differ in their
network-level statistics, the two-network test Snijders & Borgatti
(1999) proposed the vertex bootstrap for. Each network's statistics get
vertex-bootstrap standard errors from
[`vertex_bootstrap()`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md),
and each difference is tested with \$\$z = (\hat{\theta}\_x -
\hat{\theta}\_y) / \sqrt{SE_x^2 + SE_y^2}\$\$ against a standard normal
distribution, with a normal-approximation confidence interval for the
difference. It covers density, mean weight, centralisation, reciprocity
and any statistic of your own.

The test is about whole-network descriptives, and its standard errors
reflect the resampling of nodes, not of the cases behind the weights.
The p-values are not adjusted for the number of statistics tested; pass
them to [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
when several are read together. When the raw data of both networks are
available,
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
and [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md)
compare them edge by edge.

## Usage

``` r
vertex_compare(
  x,
  y,
  iter = 1000L,
  ci_level = 0.05,
  statistics = NULL,
  statistic_fn = NULL,
  directed = NULL,
  seed = NULL,
  labels = c("x", "y")
)

# S3 method for class 'net_vertex_comparison'
print(x, digits = 3, ...)

# S3 method for class 'net_vertex_comparison'
summary(object, ...)

# S3 method for class 'net_vertex_comparison'
plot(x, ...)
```

## Arguments

- x, y:

  The two networks: `netobject`s, `cograph_network`s, square weight
  matrices, or `net_vertex_bootstrap` results computed beforehand. A
  precomputed result is used as it is, so `iter`, `statistics`,
  `statistic_fn` and `directed` do not apply to it. The two must cover
  the same statistics: a mismatch, such as a directed network's
  `reciprocity` against an undirected network, or two precomputed
  results with different `statistics`, raises an error rather than
  comparing a subset. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods, `x`
  is a `net_vertex_comparison` result.

- iter:

  The number of bootstrap replicates, at least 2. Default `1000`.

- ci_level:

  The significance level of the normal-approximation intervals for the
  differences: `0.05` (default) gives 95% intervals. It is also passed
  to the bootstraps run here; with a precomputed result, the differences
  still use this level.

- statistics:

  The built-in statistics to compare, as named for
  [`vertex_bootstrap()`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md):
  `density`, `mean_weight`, `centralization` and, for directed networks,
  `reciprocity`. Default `NULL`: all that apply. A name not available
  for a network raises an error.

- statistic_fn:

  Statistics of your own, as a named list of functions, each taking the
  weight matrix and returning a single number. They are computed beside
  the built-ins on the whole matrix, diagonal included. Default `NULL`.

- directed:

  Whether the networks are directed, which decides how collisions are
  filled and whether `reciprocity` is available; one value applies to
  both networks, not to a precomputed result. Default `NULL` reads each
  network's own directedness, and for a matrix (or a network that
  records none) treats an asymmetric matrix as directed.

- seed:

  A seed passed to [`set.seed()`](https://rdrr.io/r/base/Random.html)
  once, before both bootstraps run, so the comparison can be repeated;
  the session's random-number stream then continues from it. Default
  `NULL` uses the current stream.

- labels:

  Two names for the networks, used in the print header, the plot axis
  and the column names `observed_<label>` (made syntactic with
  [`make.names()`](https://rdrr.io/r/base/make.names.html)). Default
  `c("x", "y")`.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: the
  number of decimals shown. Default `3`.

- ...:

  Not used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_vertex_comparison` result.

## Value

A `net_vertex_comparison` object. Printing it shows the two networks'
sizes and the first network's number of replicates, then the summary
table with significance stars;
[`summary()`](https://rdrr.io/r/base/summary.html) returns that table as
a data frame with one row per statistic and the columns `statistic`, the
two observed values (`observed_<label>`), `diff` (first minus second),
`se_diff`, `z`, `p_value` (two-sided), `ci_lower` and `ci_upper`;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
differences with their intervals. When a statistic has no bootstrap
variation in either network (both standard errors 0), its `z` and
`p_value` are `NA`. The object holds:

- summary:

  The summary table described above.

- x, y:

  The two `net_vertex_bootstrap` results.

- labels, ci_level:

  The settings used.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly, [`summary()`](https://rdrr.io/r/base/summary.html) the
summary data frame and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) a `ggplot`
object.

## Methods

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws each
statistic's difference as a point with its normal-approximation
interval, beside a dashed line at zero; intervals that exclude zero are
drawn in red, the others in blue.

## References

Snijders, T. A. B., & Borgatti, S. P. (1999). Non-parametric standard
errors and tests for network statistics. *Connections*, 22(2), 161-170.

## See also

[`vertex_bootstrap`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md),
[`nct`](https://pak.dynasite.org/Nestimate/reference/nct.md) for the
permutation-based comparison of edge-level structure when raw data are
available,
[`permutation`](https://pak.dynasite.org/Nestimate/reference/permutation.md).

## Examples

``` r
states <- c("plan", "code", "debug", "test")
s1 <- data.frame(
  T1 = rep(states, 5), T2 = rep(rev(states), 5),
  T3 = rep(states[c(2, 3, 4, 1)], 5)
)
s2 <- data.frame(
  T1 = rep(states[c(3, 1, 4, 2)], 5), T2 = rep(states, 5),
  T3 = rep(states[c(4, 3, 1, 2)], 5)
)
net1 <- build_network(s1, method = "relative")
net2 <- build_network(s2, method = "relative")
cmp <- vertex_compare(net1, net2, iter = 100, seed = 1)
summary(cmp)
#>        statistic observed_x observed_y        diff    se_diff          z
#> 1        density  0.5000000  0.5833333 -0.08333333 0.25478584 -0.3270721
#> 2    mean_weight  0.5000000  0.5714286 -0.07142857 0.09622256 -0.7423267
#> 3 centralization  0.3333333  0.0000000  0.33333333 0.24140627  1.3807981
#> 4    reciprocity  1.0000000  0.2500000  0.75000000 0.35488859  2.1133393
#>      p_value    ci_lower  ci_upper
#> 1 0.74361337 -0.58270441 0.4160377
#> 2 0.45788942 -0.26002131 0.1171642
#> 3 0.16734104 -0.13981427 0.8064809
#> 4 0.03457174  0.05443114 1.4455689
```
