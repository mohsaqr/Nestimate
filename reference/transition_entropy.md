# Compute the transition entropy of a Markov chain

`transition_entropy()` measures how uncertain a Markov chain's next step
is, at three levels: the branching entropy of each state (how evenly its
outgoing transitions spread), the entropy of the stationary
distribution, and the chain's entropy rate, the average uncertainty per
step. It takes a transition network, a `tna` model, a transition matrix
or wide sequence data; a `netobject_group` gives one result per group.
The entropy rate \$\$h(P) = -\sum_i \pi_i \sum_j P\_{ij} \log_b
P\_{ij},\$\$ with \\\pi\\ the stationary distribution, is the per-step
uncertainty of the chain's trajectories (Cover & Thomas, 2006, ch. 4).
It was introduced to behavioural research as gaze transition entropy by
Krejtz et al. (2015) and tracked in real time as mobile transition
matrix entropy by Krejtz et al. (2025). The normalised values divide by
\\\log_b n\\, the entropy of uniform branching among \\n\\ states; they
are the scale-free variants those papers report.

An entropy rate near 0 means the current state almost fixes the next
one; near the ceiling \\\log_b n\\ (normalised 1), every next state is
about equally likely. The redundancy \\H(\pi) - h(P)\\ is the
uncertainty that knowing the current state removes, zero when successive
states are independent. All of these are point estimates from one
transition matrix:
[`entropy_bayes()`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
gives their credible intervals,
[`entropy_network()`](https://pak.dynasite.org/Nestimate/reference/entropy_network.md)
shows which transitions carry the entropy, and
[`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
how it changes over time. A chain with more than one closed class has
more than one stationary distribution; check its classes with
[`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md).

## Usage

``` r
transition_entropy(x, base = 2, normalize = TRUE)

# S3 method for class 'net_transition_entropy'
print(x, digits = 3, ...)

# S3 method for class 'net_transition_entropy_group'
print(x, ...)

# S3 method for class 'net_transition_entropy_group'
summary(object, table = c("states", "chain"), ...)

# S3 method for class 'net_transition_entropy_group'
plot(x, title = "Transition Entropy", fill = "#0072B2", ...)

# S3 method for class 'net_transition_entropy'
summary(object, table = c("states", "chain"), ...)

# S3 method for class 'net_transition_entropy'
plot(x, title = "Transition Entropy", fill = "#0072B2", ...)
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
  (`state_cols =` or `action =`) and pass that. A `netobject_group`
  gives one result per group. An association network (`cor`, `pcor`,
  `glasso`, `ising`, `mgm`) or a matrix with negative entries raises an
  error of class `nestimate_not_markov`. A `group_tna` (package tna)
  gives one result per group, as a `netobject_group` does. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `transition_entropy()` result;
  [`print()`](https://rdrr.io/r/base/print.html) also takes a group
  result.

- base:

  The logarithm base, which sets the unit: `2` (default) for bits,
  `exp(1)` for nats, `10` for hartleys. Any finite number greater than
  1; another value raises an error of class `nestimate_bad_input`.

- normalize:

  Whether rows that do not sum to 1 (a frequency network, a count
  matrix) are divided by their sums, with a warning. With `FALSE`, such
  rows raise an error. Default `TRUE`. A state with no outgoing
  transitions (a row summing to zero) is an error either way.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) methods:
  decimal places of the values. Default `3`. The normalised values
  always print with three decimals.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: not
  used.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a
  `transition_entropy()` result.

- table:

  For [`summary()`](https://rdrr.io/r/base/summary.html): `"states"`
  (default), the per-state table, or `"chain"`, the chain-level
  quantities.

- title:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the plot
  title. Default `"Transition Entropy"`.

- fill:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the bar
  fill colour. Default `"#0072B2"`, the Okabe-Ito blue.

## Value

A `net_transition_entropy` object, or for a `netobject_group` a
`net_transition_entropy_group`: a named list of them, one per group.
[`print()`](https://rdrr.io/r/base/print.html) shows the entropy rate,
stationary entropy and redundancy, raw and normalised, and for a group
the entropy rate of each group.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame.
With `table = "states"` (default), one row per state (`state`, the
stationary probability, raw and normalised branching entropy, the
state's additive contribution \\\pi_i H(P\_{i\cdot})\\ to the entropy
rate and that contribution as a percentage of it), sorted by
contribution, largest first; for a group, with a `group` column first.
With `table = "chain"`, one row per chain-level quantity (entropy rate,
stationary entropy, redundancy and the ceiling, raw and normalised, with
the unit); for a group, one row per group (`group`, `n_states`, the raw
and normalised entropy rate and stationary entropy, `redundancy` and
`unit`). For a group,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the bar
chart below in one panel per group.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` bar chart of each state's branching entropy, with bar widths
proportional to the stationary probability (so a bar's area is
proportional to the state's contribution), a dashed line at the entropy
rate \\h(P)\\ and a dotted line at the ceiling \\\log_b n\\ (uniform
branching). The [`print()`](https://rdrr.io/r/base/print.html) methods
return their input invisibly. The object's fields:

- row_entropy:

  Named numeric vector, one value per state: its branching entropy
  \$\$H(P\_{i\cdot}) = -\sum_j P\_{ij} \log_b P\_{ij}.\$\$

- row_entropy_norm:

  `row_entropy` divided by the ceiling \\\log_b n\\, in \\\[0, 1\]\\;
  all zeros when \\n = 1\\.

- stationary:

  Named numeric vector: the stationary distribution \\\pi\\.

- stationary_entropy:

  \\H(\pi) = -\sum_i \pi_i \log_b \pi_i\\, the entropy of \\\pi\\ read
  as an independent draw each step; an upper bound on the entropy rate.

- stationary_entropy_norm:

  `stationary_entropy` divided by the ceiling \\\log_b n\\.

- entropy_rate:

  \\h(P) = \sum_i \pi_i H(P\_{i\cdot})\\, the Shannon-McMillan-Breiman
  entropy rate.

- entropy_rate_norm:

  `entropy_rate` divided by the ceiling \\\log_b n\\.

- redundancy:

  \\H(\pi) - h(P)\\, the entropy deficit due to serial dependence; zero
  for an independent chain (every row of \\P\\ equal to \\\pi\\).

- redundancy_norm:

  The relative redundancy \\(H(\pi) - h(P)) / H(\pi)\\, the share of the
  stationary entropy that order-1 memory removes; **not** `redundancy`
  divided by \\\log_b n\\. `0` when \\H(\pi) = 0\\.

- max_entropy:

  The normalising ceiling \\\log_b n\\.

- base:

  The logarithm base used.

- states:

  Character vector of state names.

## Details

The convention \\0 \log 0 = 0\\ applies, so an absorbing or
deterministic row has zero branching entropy. The chain need not be
irreducible: \\\pi\\ is the eigenvector of \\P^\top\\ for the eigenvalue
closest to 1, as elsewhere in the package. For a chain with more than
one closed class this is one stationary distribution among many; when
some of its probabilities are zero, a warning says so, and they are set
to machine epsilon and the distribution renormalised before the
entropies are computed.

\\h(P) \le H(\pi)\\ holds, with equality if and only if successive
states are independent. The deficit \\H(\pi) - h(P)\\ is reported as
`redundancy`, a measure of how much memory the chain has at order 1.

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

Shannon, C.E. (1948). A mathematical theory of communication. *Bell
System Technical Journal*, 27, 379-423.

## See also

[`entropy_network`](https://pak.dynasite.org/Nestimate/reference/entropy_network.md)
for the edge-level decomposition,
[`entropy_trajectory`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
for the sliding-window version,
[`entropy_bayes`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
for credible intervals;
[`markov_stability`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md),
[`passage_time`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
[`markov_order_test`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md),
[`chain_structure`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)

## Examples

``` r
net <- build_network(as.data.frame(trajectories), method = "relative")
te  <- transition_entropy(net)
print(te)
#> Transition Entropy (3 states, bits; ceiling = 1.585)
#> 
#>                           raw            normalised
#>   Entropy rate    h(P)  = 1.247 bits    0.787
#>   Stationary    H(pi)  = 1.501 bits    0.947
#>   Redundancy   H(pi)-h = 0.255 bits    0.170
#> 
#> Normalised: h(P) and H(pi) are raw / log_2(n_states) (0 = deterministic, 1 = uniform);
#>   redundancy is the relative redundancy (H(pi) - h(P)) / H(pi), not raw / log_2(n_states).
summary(te)
#>        state stationary row_entropy row_entropy_norm contribution
#> 1    Average  0.4431398    1.354486        0.8545857    0.6002268
#> 2     Active  0.3719025    1.040453        0.6564530    0.3869472
#> 3 Disengaged  0.1849577    1.403267        0.8853628    0.2595450
#>   contribution_pct
#> 1         48.14451
#> 2         31.03724
#> 3         20.81825
plot(te)

```
