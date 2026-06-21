# Bayesian Dirichlet-Multinomial comparison of two transition networks

Compares two transition networks estimated by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
(method `"relative"` or `"frequency"`) using a Bayesian
Dirichlet-Multinomial model. The outgoing transitions from each source
state are modelled as a Multinomial draw with a Dirichlet prior on the
transition probabilities. With a Jeffreys prior the posterior for the
transitions out of state \\i\\ is \\\mathrm{Dirichlet}(c_i + \alpha)\\,
where \\c_i\\ are the observed outgoing counts. Each edge probability is
then marginally Beta-distributed, so the posterior mean difference
between the two networks is available in closed form and a credible
interval is obtained by Monte Carlo.

This is a complement to
[`permutation`](https://saqr.me/Nestimate/reference/permutation.md): the
permutation test answers "is this difference more extreme than chance?";
the Bayesian comparison answers "what is the plausible range of the true
difference, and how precisely is it estimated given the counts?". An
edge with few outgoing transitions from its source state yields a wide
credible interval even when its row-normalised probability looks
decisive.

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
```

## Arguments

- x:

  A `netobject` (from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)),
  a `netobject_group`, or an `mcml` object. Must use a transition method
  (`"relative"` / `"frequency"` and their aliases).

- y:

  A second object of the same kind as `x`, or `NULL`. When `x` is a
  `netobject_group` and `y` is `NULL`, all pairwise comparisons among
  the groups are returned.

- prior:

  Numeric. Dirichlet prior concentration added to every cell (default
  `0.5`, the Jeffreys prior). Use `1` for a uniform (Laplace) prior.

- draws:

  Integer. Number of Monte Carlo posterior draws used for the credible
  intervals (default `10000`).

- ci:

  Numeric in (0, 1). Credible interval mass (default `0.95`).

- mean_threshold:

  Numeric. An edge is flagged significant only if the absolute posterior
  mean difference exceeds this (default `0.01`).

- bound_threshold:

  Numeric. An edge is flagged significant only if the credible-interval
  bound nearest zero exceeds this in absolute value (default `0.001`).
  Guards against differences that are detectable but negligibly small.

- seed:

  Integer or NULL. RNG seed for reproducible credible intervals.

## Value

An object of class `c("net_bayes", "net_permutation")`. It carries the
same fields as a
[`permutation`](https://saqr.me/Nestimate/reference/permutation.md)
result, so it is a drop-in wherever a `net_permutation` is consumed,
plus Bayesian extras:

- x, y:

  The two input `netobject`s.

- diff:

  Posterior mean difference matrix (`prob_x - prob_y`); the analogue of
  the permutation observed difference.

- diff_sig:

  Difference where `sig`, else 0.

- p_values:

  P-value matrix (the two-sided Bayesian p-equivalent).

- effect_size:

  Posterior mean difference over its posterior SD.

- ci_lower, ci_upper:

  Credible-interval bound matrices.

- pd:

  Probability-of-direction matrix in \\\[0.5, 1\]\\.

- p_bayes:

  Alias of `p_values` (two-sided Bayesian p, \\2(1-pd)\\).

- prob_x, prob_y:

  Posterior mean transition-probability matrices.

- sig:

  Logical significance matrix (CI excludes zero, mean and nearest bound
  exceed their thresholds).

- summary:

  Long-format data frame whose columns are a superset of
  `summary.net_permutation`
  (`from, to, weight_x, weight_y, diff, effect_size, p_value, sig`) plus
  `count_x, count_y, ci_lower, ci_upper, ci_width, pd`.

- method, iter, alpha, paired, adjust:

  Permutation-compatible settings (`iter = draws`, `alpha = 1 - ci`,
  `paired = FALSE`, `adjust = "none"`).

- prior, draws, ci, mean_threshold, bound_threshold:

  Bayesian settings.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

Jeffreys, H. (1946). An invariant form for the prior probability in
estimation problems. *Proceedings of the Royal Society of London A*,
186(1007), 453-461.

## See also

[`permutation`](https://saqr.me/Nestimate/reference/permutation.md),
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)

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
