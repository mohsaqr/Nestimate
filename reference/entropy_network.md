# Transition Entropy Network

`entropy_network()` splits a chain's entropy rate into one term per
transition and returns the terms as a network. The entropy rate \\h(P) =
-\sum\_{ij} \pi_i P\_{ij} \log_b P\_{ij}\\
([`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md);
Krejtz et al., 2015, 2025) is a sum over transitions, so every edge \\i
\to j\\ owns the exact term \\\pi_i P\_{ij} \log_b(1/P\_{ij})\\ of the
chain-level uncertainty. No new quantity is estimated: the network
places the summands of the entropy-rate equation on the transition
graph, which shows where in the process the uncertainty lies. Two other
edge weights are available: each transition's surprisal, and its
contribution to the entropy production rate (irreversibility). The
result is a regular `netobject`, so it prints, summarises and plots
([`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html))
like any other Nestimate network.

Under the default weight, a heavy edge is a transition that is both
frequent and unpredictable; under `"surprisal"`, a heavy edge is a rare
one, whatever its frequency. The weights are computed from one estimated
transition matrix, so an edge that rests on a few observations can carry
a large weight by chance:
[`entropy_bayes()`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
gives each contribution a credible interval and keeps only the edges
that credibly hold a share of the entropy rate.

## Usage

``` r
entropy_network(
  x,
  base = 2,
  weight = c("contribution", "surprisal", "production"),
  scaling = c("none", "share", "chance"),
  normalize = TRUE
)
```

## Arguments

- x:

  The chain: a `netobject` or `cograph_network` with non-negative
  transition weights, a `tna` model, a numeric transition matrix, or
  wide sequence data (a data frame or a character matrix, one row per
  sequence) that holds a numbered run of positions (T1, T2, ...), from
  which a relative transition network is built; other columns are set
  aside with a message. For a frame without such a run, build the
  network with
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  (`state_cols =` or `action =`) and pass that. An association network
  (`cor`, `pcor`, `glasso`, `ising`, `mgm`) or a matrix with negative
  entries raises an error of class `nestimate_not_markov`. A
  `netobject_group` is not accepted; pass one network at a time.

- base:

  The logarithm base, which sets the unit: `2` (default) for bits,
  `exp(1)` for nats, `10` for hartleys. Any positive number other than
  1.

- weight:

  What each edge's weight measures:

  `"contribution"`

  :   (default) \\\pi_i P\_{ij} \log_b(1/P\_{ij})\\, the edge's additive
      share of the entropy rate; all edge weights sum exactly to
      \\h(P)\\. Heavy edges are transitions that are both frequent *and*
      unpredictable.

  `"surprisal"`

  :   \\\log_b(1/P\_{ij})\\, the information content of observing the
      transition, ignoring how often it occurs. Heavy edges are rare,
      surprising transitions; a low surprisal marks an expected
      transition, which is the edge-level reading of predictability.

  `"production"`

  :   \\F\_{ij} \log_b(F\_{ij}/F\_{ji})\\ with \\F\_{ij} = \pi_i
      P\_{ij}\\, the edge's contribution to the entropy production rate
      (irreversibility). It is positive on the dominant direction of a
      pair and negative on the reverse, and the two always sum to
      \\(F\_{ij}-F\_{ji}) \log_b(F\_{ij}/F\_{ji}) \ge 0\\. A pair with
      flow in one direction only has infinite production and is left out
      (weight 0); the number of such pairs is stored with the network's
      parameters. The total over the included pairs is the production
      rate. It is `0` for a reversible chain (detailed balance); when
      one-way pairs were left out, a total of 0 does not show
      reversibility, since such pairs mean the chain is not reversible.

- scaling:

  How the weights are expressed. `"none"` (default) keeps them in `base`
  units. `"share"` (contribution only) gives each edge as a *percentage*
  of the entropy rate: the weights sum to 100, each one reads "this
  transition holds x% of the process's unpredictability", and they
  compare across networks whatever the number of states or the base.
  `"chance"` (surprisal only) divides by \\\log_b n\\, the surprisal of
  a transition when all \\n\\ next states are equally likely: `1` is
  chance level, below 1 an expected transition, above 1 one rarer than
  chance. `"share"` with another weight, or `"chance"` with another
  weight, raises an error.

- normalize:

  Whether rows that do not sum to 1 (a frequency network, a count
  matrix) are divided by their sums, with a warning. With `FALSE`, such
  rows raise an error. Default `TRUE`. A state with no outgoing
  transitions (a row summing to zero) is an error either way.

## Value

A `netobject` (class `c("netobject", "cograph_network")`; for a bare
`cograph_network` input, that network with its own class) with method
`"entropy"`, whose weights are the per-edge entropy quantities. Print it
for an overview; [`summary()`](https://rdrr.io/r/base/summary.html) (see
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md))
gives its network metrics,
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
one row per edge, and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it. When `x` is a network (or sequence data, from which a relative
network is built), the result is that network with the entropy weights
swapped in: its initial probabilities, metadata, node groups and node
coordinates are kept, so it plots with the same TNA styling and layout
as its source. The stored parameters are `base`, `weight`, `scaling`,
the entropy rate `entropy_rate` (whatever the weight) and the stationary
distribution `stationary`, plus `production_rate` and `n_oneway_pairs`
when `weight = "production"`.

The object also declares the entropy house style through its
`meta$splot` producer contract (honoured by cograph \>= 2.4.4): no
minimum-weight pruning (bit values are smaller than probabilities), edge
labels with 2 digits (1 for `"share"`), vermilion edges, node rings
showing the stationary distribution, and TNA rather than psychometric
geometry. The production weights are signed, so they keep cograph's
positive and negative edge colours instead of the vermilion. With
`weight = "surprisal"` and `scaling = "chance"`, the edges are drawn and
sized by the transition probabilities, transitions likelier than chance
in blues (darker for likelier) and those at or below chance dashed grey.
Because the contract states the styling outright, cograph needs no
knowledge of the `"entropy"` method: `cograph::splot(ent)` renders it
with no arguments, and any argument given overrides the contract.

## Details

An impossible transition (\\P\_{ij} = 0\\) gets weight `0` under every
definition (the \\0 \log 0 = 0\\ convention), so the entropy network has
no edge the transition network lacks. Self-loops are kept like any other
edge; their production weight is always 0. A deterministic transition
(\\P\_{ij} = 1\\) also gets contribution and surprisal `0`: observing
the inevitable carries no information.

## References

Cover, T.M. & Thomas, J.A. (2006). *Elements of Information Theory*, 2nd
ed., chapter 4. Wiley.

Krejtz, K., Duchowski, A., Szmidt, T., Krejtz, I., Gonzalez Perilli, F.,
Pires, A., Vilaro, A., & Villalobos, N. (2015). Gaze transition entropy.
*ACM Transactions on Applied Perception*, 13(1), 4:1-4:20.
[doi:10.1145/2834121](https://doi.org/10.1145/2834121)

Krejtz, K., Hughes, C.J., Stasiak, I., Duchowski, A., & Krejtz, I.
(2025). Real-time mobile transition matrix entropy based on eye and head
movements. *Proceedings of ETRA '25*.
[doi:10.1145/3715669.3723128](https://doi.org/10.1145/3715669.3723128)

## See also

[`transition_entropy`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
for the chain- and state-level summary,
[`entropy_bayes`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
for credible intervals on the decomposition,
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

## Examples

``` r
net <- build_network(group_regulation_long,
                     method = "relative",
                     actor = "Actor", action = "Action", time = "Time")
ent <- entropy_network(net)
ent
#> Network (method: entropy) [directed]
#>   Weights: [0.002, 0.132]  |  mean: 0.031
#> 
#>   Weight matrix:
#>              adapt cohesion consensus coregulate discuss emotion monitor  plan
#>   adapt      0.000    0.011     0.011      0.002   0.005   0.008   0.003 0.002
#>   cohesion   0.002    0.009     0.033      0.024   0.016   0.024   0.011 0.027
#>   consensus  0.009    0.023     0.074      0.113   0.113   0.069   0.051 0.132
#>   coregulate 0.008    0.014     0.032      0.010   0.042   0.036   0.025 0.040
#>   discuss    0.041    0.032     0.080      0.046   0.070   0.052   0.019 0.011
#>   emotion    0.002    0.058     0.057      0.018   0.037   0.031   0.019 0.036
#>   monitor    0.004    0.011     0.021      0.012   0.026   0.015   0.005 0.023
#>   plan       0.002    0.033     0.127      0.025   0.065   0.100   0.069 0.130
#>   synthesis  0.013    0.004     0.014      0.005   0.007   0.007   0.002 0.007
#>              synthesis
#>   adapt          0.000
#>   cohesion       0.002
#>   consensus      0.013
#>   coregulate     0.009
#>   discuss        0.060
#>   emotion        0.003
#>   monitor        0.005
#>   plan           0.004
#>   synthesis      0.000 
#> 
#>   Initial probabilities:
#>   consensus     0.214  ████████████████████████████████████████
#>   plan          0.204  ██████████████████████████████████████
#>   discuss       0.175  █████████████████████████████████
#>   emotion       0.151  ████████████████████████████
#>   monitor       0.144  ███████████████████████████
#>   cohesion      0.060  ███████████
#>   synthesis     0.019  ████
#>   coregulate    0.019  ████
#>   adapt         0.011  ██
# \donttest{
if (requireNamespace("cograph", quietly = TRUE)) {
  cograph::splot(ent)
}

# }
```
