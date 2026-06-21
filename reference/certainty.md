# Analytic certainty of network edges (Bayesian Dirichlet-Multinomial)

Closed-form alternative to
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
for transition networks. Models the outgoing transitions from each state
as a Dirichlet-Multinomial process: with a Jeffreys prior the posterior
for state \\i\\ is \\\mathrm{Dirichlet}(c_i + \mathrm{prior})\\, so each
edge is marginally Beta and its posterior mean, standard deviation,
credible interval and stability decision are available analytically. No
resampling, so it runs in microseconds.

The return value has the same structure as
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
(same slots and summary columns) and carries class
`c("net_certainty", "net_bootstrap")`, so
[`summary()`](https://rdrr.io/r/base/summary.html) and any code that
consumes a `net_bootstrap` object work unchanged.

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
```

## Arguments

- x:

  A `netobject` from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
  using a transition-probability method (`"relative"` / `"tna"`), or a
  `netobject_group`.

- prior:

  Numeric. Dirichlet prior concentration added to every cell (default
  `0.5`, the Jeffreys prior).

- ci_level:

  Numeric in (0,1). Tail level for credible intervals and the stability
  decision (default `0.05`, i.e. a 95\\ match
  [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md).

- inference:

  Character. `"stability"` (default) tests whether the posterior keeps
  the edge within a multiplicative `consistency_range` of its weight;
  `"threshold"` tests whether the edge exceeds `edge_threshold`.

- consistency_range:

  Numeric vector of length 2. Multiplicative bounds for stability
  inference (default `c(0.75, 1.25)`).

- edge_threshold:

  Numeric or NULL. Fixed threshold for `inference = "threshold"`. If
  NULL, defaults to the 10th percentile of non-zero edge weights.

## Value

An object of class `c("net_certainty", "net_bootstrap")` with the same
fields as
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md):
`original`, `mean`, `sd`, `p_values`, `significant`, `ci_lower`,
`ci_upper`, `cr_lower`, `cr_upper`, `summary`, `model`, `method`,
`params`, `ci_level`, `inference`, `consistency_range`,
`edge_threshold`, plus `prior` and `iter = NA` (no iterations).

## Details

Certainty (this function), stability
([`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md))
and reliability
([`reliability`](https://saqr.me/Nestimate/reference/network_reliability.md))
answer different questions about an edge: how precisely it is pinned
down by the observed counts, whether it survives resampling the
sequences, and whether it is consistent across split-halves. Certainty
and stability agree on homogeneous data; certainty is over-confident
when the data are a mixture of latent classes, because it treats
transitions clustered within a sequence as independent.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

## See also

[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
[`bayes_compare`](https://saqr.me/Nestimate/reference/bayes_compare.md),
[`network_reliability`](https://saqr.me/Nestimate/reference/network_reliability.md)

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
