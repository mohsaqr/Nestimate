# Analytic certainty of network edges (Bayesian Dirichlet-Multinomial)

`certainty()` gives each edge of a transition-probability network a
posterior mean, standard deviation, credible interval and a decision on
whether the edge is certain, all in closed form, and returns them in the
layout of
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).
It models the outgoing transitions of each state as a
Dirichlet-Multinomial process: with a Dirichlet prior of concentration
`prior`, the posterior for state \\i\\ is \\\mathrm{Dirichlet}(c_i +
\mathrm{prior})\\, where \\c_i\\ are its transition counts, so each edge
probability is marginally Beta. It reads a `netobject` built with
`method = "relative"`, or a `netobject_group` of them.

With `inference = "stability"`, an edge is certain when the posterior
mass outside a band around its weight (`consistency_range`) is below
`ci_level`; with `"threshold"`, when the mass below `edge_threshold` is.
Certainty says how precisely the observed counts pin an edge down, not
whether the edge survives resampling the sequences
([`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md))
or holds across split halves
([`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)).
Certainty and bootstrap stability agree on homogeneous data. Certainty
treats every transition as independent, so it is over-confident when the
sequences come from a mixture of latent classes, the more so the longer
they are; use
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
there.

## Usage

``` r
certainty(
  x,
  prior = 0.5,
  ci_level = 0.05,
  inference = c("stability", "threshold"),
  consistency_range = c(0.75, 1.25),
  edge_threshold = NULL
)

# S3 method for class 'net_certainty'
print(x, ...)
```

## Arguments

- x:

  The network: a `netobject` built with
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  using `method = "relative"` (or its aliases `"tna"`, `"transition"`),
  without `scaling`, or a `netobject_group` of such networks, each
  handled on its own. A bare `cograph_network` is recounted from the
  sequences it carries, and used only when they reproduce its weights;
  otherwise, and for an `mcml`, an error of class `nestimate_bad_input`
  is raised. Another method or a scaled network raises an error. For the
  [`print()`](https://rdrr.io/r/base/print.html) method: a
  `net_certainty` object.

- prior:

  The Dirichlet concentration added to every cell of the count matrix,
  zero cells and self-loops included; a positive number. Default `0.5`,
  the Jeffreys prior.

- ci_level:

  The tail level, in (0, 1), of the equal-tailed credible intervals and
  of the decision: an edge is certain when its posterior mass (the
  `p_value` column) is below it. Default `0.05`, a 95% interval. Named
  to match
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

- inference:

  The decision rule. `"stability"` (default) takes the posterior mass
  outside the band `consistency_range` times the edge's weight;
  `"threshold"` takes the posterior mass below `edge_threshold`.

- consistency_range:

  For `inference = "stability"`: the two multiplicative bounds of the
  band around each weight, positive numbers in either order. Default
  `c(0.75, 1.25)`.

- edge_threshold:

  For `inference = "threshold"`: the weight an edge must exceed. Default
  `NULL`: the 10th percentile of the non-zero weights. Ignored for
  `"stability"`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) method:
  ignored.

## Value

For a `netobject`: a `net_certainty` object (class
`c("net_certainty", "net_bootstrap")`), so the `net_bootstrap` methods
work on it. [`print()`](https://rdrr.io/r/base/print.html) lists up to
five certain edges with the largest absolute posterior means (with their
intervals), then the prior, the count of certain edges and the settings;
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
non-zero edge with columns `from`, `to`, `weight`, `mean`, `sd`,
`p_value`, `sig`, `ci_lower` and `ci_upper`, plus `cr_lower` and
`cr_upper` (the band) for `inference = "stability"`;
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it. The object has the fields of a
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
result: `original` (the input network), `mean`, `sd`, `p_values`,
`significant` (the weights of the certain edges, zero elsewhere),
`ci_lower`, `ci_upper`, `cr_lower`, `cr_upper` (matrices), `summary`,
`model` (the network of certain edges only), `method`, `params`,
`ci_level`, `inference`, `consistency_range` and `edge_threshold`, plus
`prior`, `ci_method = "analytic"` and `iter = NA` (no iterations).

For a `netobject_group`: a named list of those objects, one per network,
of class `c("net_certainty_group", "net_bootstrap_group", "list")`,
printed and summarised by the `net_bootstrap_group` methods
([`summary()`](https://rdrr.io/r/base/summary.html) stacks the
per-network tables with a `group` column).

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly.

## Details

Nothing is resampled: the posterior quantities are Beta means, variances
and quantiles ([`stats::qbeta()`](https://rdrr.io/r/stats/Beta.html),
[`stats::pbeta()`](https://rdrr.io/r/stats/Beta.html)), so a call takes
milliseconds.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

## See also

[`bootstrap_network`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
[`bayes_compare`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md),
[`network_reliability`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
                   V3 = c("C","A","C","B","A"))
net <- build_network(seqs, method = "relative")
cert <- certainty(net)
cert
#> Certainty (Dirichlet)  [Transition Network (relative) | directed]
#>   Prior      : Dirichlet(0.50)  |  Nodes : 3
#>   Edges      : 0 certain / 3 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
summary(cert)
#>   from to weight      mean        sd   p_value   sig  ci_lower  ci_upper
#> 1    A  B      1 0.7777778 0.1772720 0.3653545 FALSE 0.3485528 0.9927924
#> 2    B  C      1 0.8181818 0.1512819 0.2740159 FALSE 0.4405413 0.9943896
#> 3    C  A      1 0.7777778 0.1772720 0.3653545 FALSE 0.3485528 0.9927924
#>   cr_lower cr_upper
#> 1     0.75     1.25
#> 2     0.75     1.25
#> 3     0.75     1.25
```
