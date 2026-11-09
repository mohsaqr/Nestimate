# Bayesian Dirichlet-Multinomial comparison of two transition networks

`bayes_compare()` compares two transition networks edge by edge and
returns, for every transition, the posterior mean difference in
transition probability, a credible interval, and the probability that
the difference lies on its side of zero. The transitions out of each
state are modelled as a multinomial draw with a Dirichlet prior on the
transition probabilities, so the posterior for the transitions out of
state \\i\\ is \\\mathrm{Dirichlet}(c_i + \alpha)\\, where \\c_i\\ are
the observed outgoing counts and \\\alpha\\ is `prior`. Each edge
probability is then marginally Beta-distributed: the posterior mean
difference has a closed form, and the credible interval comes from Monte
Carlo draws. It takes two networks estimated by
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with method `"relative"` or `"frequency"` (for either, the difference is
in transition probabilities estimated from the counts), a group of
networks, or two edge-betweenness networks from
[`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
(see *Details*). An `mcml` object raises an error of class
`nestimate_bad_input`: its layers hold different states and carry no
transition counts.

The comparison complements
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md):
the permutation test asks whether a difference is more extreme than
chance, and this comparison gives the plausible range of the difference
and how precisely the counts estimate it. An edge whose source state has
few outgoing transitions gets a wide interval even when its
row-normalised probabilities look far apart. An edge is flagged as
credibly different when its interval excludes zero and its mean and
nearest bound pass `mean_threshold` and `bound_threshold`; the
thresholds are conventions, and no correction is made for comparing many
edges at once. The reported `p_value` (the Bayesian p-equivalent \\2(1 -
p_d)\\, with \\p_d\\ the probability of direction) summarises posterior
mass and is not a frequentist p-value. The model treats every transition
as an independent draw, so it does not account for transitions being
nested in sequences, which a permutation of sequences with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
respects.

## Usage

``` r
bayes_compare(
  x,
  y = NULL,
  prior = 0.5,
  draws = 10000L,
  ci = 0.95,
  mean_threshold = 0.01,
  bound_threshold = 0.001,
  seed = NULL
)

# S3 method for class 'net_bayes'
print(x, ...)

# S3 method for class 'net_bayes'
summary(object, ...)

# S3 method for class 'net_bayes'
plot(x, significant_only = TRUE, title = NULL, ...)

# S3 method for class 'net_bayes_group'
print(x, ...)

# S3 method for class 'net_bayes_group'
summary(object, ...)
```

## Arguments

- x:

  The first network: a `netobject` estimated with `"relative"` or
  `"frequency"` or one of their aliases, a `netobject_group`, or a
  `net_edge_betweenness` object. A bare `cograph_network` is recounted
  from the sequences it carries, and used only when they reproduce its
  weights; otherwise, and for an `mcml`, an error of class
  `nestimate_bad_input` is raised. Both networks must use the same
  method and hold the same nodes; any other method raises an error,
  because the model needs transition counts. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_bayes`, or for [`print()`](https://rdrr.io/r/base/print.html) a
  `net_bayes_group`.

- y:

  The second network, of the same kind as `x`. Default `NULL`: with a
  `netobject_group` as `x`, every pair of its groups is compared. With
  two `netobject_group`s, the groups sharing a name are compared, and an
  error is raised when none do.

- prior:

  The Dirichlet concentration added to every cell of the count matrices,
  a positive number. Default `0.5`, the Jeffreys prior; `1` gives a
  uniform (Laplace) prior. Every cell gets it, so a transition observed
  in neither network still has a small posterior probability.

- draws:

  The number of Monte Carlo posterior draws behind the credible
  intervals, at least 2. Default `10000`. For edge-betweenness networks
  the betweenness is recomputed on every draw, so time grows with
  `draws`.

- ci:

  The mass of the equal-tailed credible interval, a number between 0
  and 1. Default `0.95`.

- mean_threshold:

  The absolute posterior mean difference an edge must exceed to be
  flagged as credibly different. Default `0.01`, on the
  transition-probability scale.

- bound_threshold:

  The absolute value the interval bound nearest zero must exceed for an
  edge to be flagged. It guards against differences that are detectable
  but negligibly small. Default `0.001`.

- seed:

  An integer that makes the draws reproducible; the session's own random
  stream is left as it was. Default `NULL`. For a group, it is set again
  before each comparison.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): further
  arguments for
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  when cograph is installed (such as `weight_digits`, the digits the
  drawn differences are rounded to, default 2, or `edge_label_digits`,
  the digits of the edge labels, default 3); otherwise ignored. The
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods ignore it.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a `net_bayes`
  or `net_bayes_group`.

- significant_only:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): whether
  to draw only the edges flagged as credibly different. Default `TRUE`.
  When no edge is left to draw, an error suggests
  `significant_only = FALSE`.

- title:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the plot
  title. Default `NULL`, which gives `"Bayesian network difference"`
  (with cograph) or a title naming the posterior mean difference
  (without it).

## Value

A `net_bayes` object (class
`c("net_bayes", "netdifference", "net_permutation")`), or, for a
`netobject_group` input, a `net_bayes_group`: a list of `net_bayes`
objects named by comparison (`"A vs B"`), or by group name when two
groups are compared. [`print()`](https://rdrr.io/r/base/print.html)
shows the method, the prior, draws, interval and thresholds, and the
numbers of nodes, of edges compared and of edges credibly different (for
a group, it lists the comparisons);
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
edge observed in at least one network (for a group, with a `comparison`
column); and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws the difference network. A `net_bayes` carries the same fields as a
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
result, so it is accepted wherever a `net_permutation` is, and also
carries a difference matrix for the `netdifference` verbs
([`as_netdifference()`](https://pak.dynasite.org/Nestimate/reference/as_netdifference.md),
cograph's difference plots). Its fields:

- x, y:

  The two input networks.

- diff:

  Posterior mean difference matrix (`prob_x - prob_y`); the analogue of
  the permutation test's observed difference.

- difference_matrix:

  The same matrix as `diff`, for cograph's `netdifference` helpers.

- diff_sig:

  The difference where `sig` is `TRUE`, else 0.

- p_values:

  The two-sided Bayesian p-equivalent, in the field a `net_permutation`
  consumer reads as p-values (see `p_bayes`).

- effect_size:

  Posterior mean difference over its posterior standard deviation (0
  where that is 0).

- ci_lower, ci_upper:

  Credible-interval bound matrices.

- p_difference:

  Probability of direction: the share of posterior mass on the dominant
  side of zero, in \\\[0.5, 1\]\\.

- p_bayes:

  The same matrix as `p_values`: the two-sided Bayesian p-*equivalent*
  \\2(1-\mathrm{p\\difference})\\. It summarises posterior mass, not a
  frequentist tail probability, so it is not a p-value and should not be
  reported as one.

- prob_x, prob_y:

  Posterior mean transition-probability matrices (posterior mean
  betweenness for edge-betweenness input).

- sig:

  Logical matrix: the interval excludes zero and the mean and nearest
  bound exceed their thresholds.

- summary:

  The per-edge table that
  [`summary()`](https://rdrr.io/r/base/summary.html) returns, with the
  columns of a `net_permutation` summary (`from`, `to`, `weight_x`,
  `weight_y`, `diff`, `effect_size`, `p_value`, `sig`) plus `count_x`,
  `count_y`, `ci_lower`, `ci_upper`, `ci_width` and `p_difference`.

- method, iter, alpha, paired, adjust:

  Permutation-compatible settings (`iter = draws`, `alpha = 1 - ci`,
  `paired = FALSE`, `adjust = "none"`); `method` is `"edge_betweenness"`
  for edge-betweenness input.

- prior, draws, ci, mean_threshold, bound_threshold:

  The Bayesian settings.

- observed_diff, source_method, edge_betweenness:

  For edge-betweenness input only: the plug-in difference, the source
  networks' method and the betweenness settings.

For the [`print()`](https://rdrr.io/r/base/print.html) methods: the
input, invisibly. For
[`summary()`](https://rdrr.io/r/base/summary.html): the per-edge data
frame described under the `summary` field; for a `net_bayes_group`, the
tables of every comparison stacked, with a `comparison` column first.
For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): with
cograph installed, it draws the difference network with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
each edge labelled with its posterior mean difference and credible
interval, and returns the `cograph_network` invisibly; without cograph,
it returns a `ggplot` chord diagram in which colour gives the direction
of the difference and width its size.

## Details

For two `net_edge_betweenness` objects (source method `"relative"`
only), edge betweenness is a nonlinear function of the whole transition
matrix, so Beta marginals do not apply. Instead the full transition
matrix is drawn from each network's row-wise Dirichlet posterior and
edge betweenness is recomputed on every draw, the Bayesian analogue of
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
on edge-betweenness networks. The result summarises the posterior of
\\EB(x) - EB(y)\\: `diff` is the posterior mean difference, `prob_x` and
`prob_y` the posterior mean betweenness matrices, and `observed_diff`
the plug-in difference of the two input networks. Both inputs must use
the same `invert` setting, and `mean_threshold` and `bound_threshold`
are then on the betweenness scale.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

Jeffreys, H. (1946). An invariant form for the prior probability in
estimation problems. *Proceedings of the Royal Society of London A*,
186(1007), 453-461.

## See also

[`permutation`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
for the frequentist complement;
[`certainty`](https://pak.dynasite.org/Nestimate/reference/certainty.md)
for single-network posterior edge intervals;
[`subtract_networks`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md)
and
[`as_netdifference`](https://pak.dynasite.org/Nestimate/reference/as_netdifference.md)
for the difference verbs;
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
n1 <- build_network(s1, method = "relative")
n2 <- build_network(s2, method = "relative")
bayes_compare(n1, n2, draws = 500, seed = 1)
#> Bayesian Dirichlet-Multinomial Comparison: Transition Network (relative probabilities) 
#>   Prior: Dirichlet(0.50)  |  Draws: 500  |  CI: 95%
#>   Thresholds: |mean diff| > 0.010, nearest CI bound > 0.001
#>   Nodes: 3  |  Edges compared: 6  |  Credibly different: 0
```
