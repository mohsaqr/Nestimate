# Qualitative structure of a discrete-time Markov chain

`chain_structure()` describes the qualitative structure of a
discrete-time Markov chain: the properties that depend only on which
transitions are possible (the support of the transition matrix), not on
any starting distribution. It classifies every state as absorbing,
recurrent or transient, finds the communicating classes and the period
of each recurrent state, says whether the chain is irreducible,
aperiodic, regular and reversible, and computes the hitting
probabilities and, when transient states exist, the probability that
each ends in each closed class (an absorbing state or a closed class of
several states) and the expected number of steps it takes to enter one.
It reads a transition network, a `tna` model, a transition matrix or
sequences, and analyses a grouped network one network at a time.

It is a diagnostic to run before trusting
[`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
or
[`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md).
Both assume a regular chain (irreducible and aperiodic), whose
stationary distribution is unique; `is_regular` says whether that holds.
The structure follows the estimated matrix: a transition never observed
in the data counts as impossible, so a chain estimated from few
sequences can look reducible because a rare move was not seen.

## Usage

``` r
chain_structure(x, normalize = TRUE, tol = 1e-10)

# S3 method for class 'chain_structure'
print(x, ...)

# S3 method for class 'chain_structure'
plot(x, show_values = TRUE, digits = 2L, ...)

# S3 method for class 'chain_structure'
summary(object, ...)

# S3 method for class 'chain_structure_group'
print(x, ...)

# S3 method for class 'chain_structure_group'
summary(object, ...)

# S3 method for class 'summary_chain_structure'
print(x, ...)
```

## Arguments

- x:

  The chain: a `netobject`, bare `cograph_network` or `tna` model with
  transition weights, a numeric transition matrix, or a wide sequence
  data frame or character matrix (estimated with
  `build_network(method = "relative")`). A `netobject_group` is analysed
  network by network. An association network (`cor`, `pcor`, `glasso`,
  `ising`, `mgm`) or a matrix with negative entries raises an error of
  class `nestimate_not_markov`. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `chain_structure` object
  ([`print()`](https://rdrr.io/r/base/print.html) also takes a
  `chain_structure_group` or a `summary_chain_structure`).

- normalize:

  Whether to rescale rows that do not sum to 1 so that they do, with a
  warning, as
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  does. Default `TRUE`. With `FALSE`, such a matrix raises an error. A
  state with no outgoing transitions (a row of zeros) raises an error
  either way.

- tol:

  The tolerance below which a transition probability counts as zero when
  the support graph is built (it drives `classification`,
  `communicating_classes`, `period` and `hitting_probabilities`), and
  the tolerance of the reversibility check (detailed balance). Default
  `1e-10`. It does **not** govern the absorbing-state test: a state is
  absorbing only when `P[i, i]` equals 1 to a fixed internal tolerance
  of `.Machine$double.eps^0.5`, so raising `tol` to ignore tiny
  transition probabilities never reclassifies a near-deterministic state
  as absorbing.

- ...:

  Ignored, except by the [`print()`](https://rdrr.io/r/base/print.html)
  methods of a `chain_structure_group`, which pass it to each network's
  [`print()`](https://rdrr.io/r/base/print.html), and of a
  `summary_chain_structure`, which pass it to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

- show_values:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: whether to print each cell's probability inside it. Default
  `TRUE`. Values are never printed for more than 20 states; set `FALSE`
  earlier when the labels overlap (from about 10 states).

- digits:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: decimal places of the cell values. Default `2`.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) methods: a
  `chain_structure` or `chain_structure_group` object.

## Value

A `chain_structure` object, or for a `netobject_group` a
`chain_structure_group` (class `c("chain_structure_group", "list")`), a
named list holding one `chain_structure` per network.
[`print()`](https://rdrr.io/r/base/print.html) shows the chain-level
properties (for a group, one block per network);
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
(see *Summary columns*): for a `chain_structure`, one row per state, of
class `c("summary_chain_structure", "data.frame")`, carrying the
chain-level flags `is_regular`, `is_irreducible`, `is_aperiodic`,
`is_reversible`, `n_classes` and `absorbing_states` as attributes that
its [`print()`](https://rdrr.io/r/base/print.html) method shows as a
header above the table; for a `chain_structure_group`, a plain data
frame with one row per group and state.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` heatmap of the hitting probabilities (see *Plot colours*). The
[`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly.

The fields of a `chain_structure`:

- `states`:

  Character vector of state names.

- `classification`:

  Named character vector, one of `"absorbing"`, `"recurrent"` or
  `"transient"` per state.

- `communicating_classes`:

  List of state-name vectors, one per strongly connected component of
  the support graph.

- `recurrent_classes`:

  The communicating classes that are closed (no transition leaves the
  class).

- `transient_classes`:

  The communicating classes that are not closed.

- `absorbing_states`:

  Character vector of states with `P[i, i] = 1`, tested to within
  `.Machine$double.eps^0.5`; `tol` does not relax this.

- `period`:

  Named integer vector: the period of each recurrent state, `NA` for
  transient states.

- `is_irreducible`:

  `TRUE` when there is exactly one communicating class.

- `is_aperiodic`:

  `TRUE` when every recurrent state has period 1.

- `is_regular`:

  `is_irreducible && is_aperiodic`.

- `is_reversible`:

  `TRUE` when the chain satisfies detailed balance against its
  stationary distribution; `NA` for a chain that is not irreducible (no
  unique stationary distribution).

- `hitting_probabilities`:

  `n x n` matrix whose `[i, j]` entry is the probability of ever
  reaching `j` from `i`, computed over the same `tol`-thresholded
  support graph as `classification`, so the two agree: a state in a
  closed class has hitting probability 0 to every state outside it.

- `absorption_probabilities`:

  `n_transient x n_closed` matrix whose `[i, j]` entry is the
  probability that the chain started in transient state `i` ends in
  closed class `j`; each row sums to 1. An absorbing state is a column
  named after it; a closed class of several states is one column named
  by its states joined with `+` (for example `B+C`). `NULL` when the
  chain has no transient state.

- `mean_absorption_time`:

  Named numeric vector of the expected number of steps each transient
  state takes to enter a closed class, whichever it is. `NULL` when
  `absorption_probabilities` is.

- `P`:

  The transition matrix analysed (after any row normalisation).

## Details

The fundamental-matrix absorption analysis follows Kemeny & Snell
(1976), with every closed class lumped into one absorbing state, so each
transient state's absorption probabilities sum to 1; the hitting
probabilities are the minimal non-negative solution of the linear system
in Norris (1997), solved over the states that can reach the target. The
diagonal of the hitting matrix uses the return-time convention, the
probability of returning to the state in one or more steps, as
`markovchain::hittingProbabilities()` does.

## Plot colours

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
draws the hitting-probability matrix as a heatmap, with rows and columns
ordered by communicating class (closed classes first, larger classes
before smaller ones) so the block structure is visible. Cell colour runs
from white (0) to dark blue (1) and encodes
`P(ever reach j | start at i)`; the diagonal uses the return-time
convention, `P(return to j in >= 1 steps)`. The state labels on both
axes name each state's class, as in `Active (recurrent)`, and the
subtitle lists the chain-level properties. A chain that is not
irreducible shows zero entries between classes that cannot reach each
other; when every state ends in one absorbing state, that state's column
is all 1s.

## Summary columns

[`summary()`](https://rdrr.io/r/base/summary.html) of a
`chain_structure` returns one row per state, combining every per-state
result: `state`, `classification`, `period`, `persistence` (the diagonal
of the transition matrix), `return_probability` (the diagonal of the
hitting matrix) and `sojourn_steps` (`1 / (1 - persistence)`, `Inf` for
an absorbing state). An irreducible chain adds `stationary_probability`.
A chain with transient states adds `absorption_probability` (one closed
class) or one `absorbed_in_<class>` column per closed class (named as
the columns of `absorption_probabilities`), and `mean_absorption_time`
(steps to enter a closed class); these are `NA` for states that are not
transient. Probabilities are rounded to 4 decimals and step counts to 2.
Columns are ordered for readability: identifiers first, classification
second, per-state results last.

[`summary()`](https://rdrr.io/r/base/summary.html) of a
`chain_structure_group` stacks these tables into one data frame with a
`group` column first, one row per group and state; a column that only
some groups have is `NA` for the others.

## References

Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*.
Springer-Verlag.

Norris, J. R. (1997). *Markov Chains*. Cambridge University Press.

## See also

[`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
[`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md),
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
net <- build_network(as.data.frame(trajectories), method = "relative")
cs  <- chain_structure(net)
print(cs)
#> Chain structure  [3 states, 1 communicating classes]
#>   irreducible: TRUE   aperiodic: TRUE   regular: TRUE   reversible: FALSE
#>   recurrent classes: 1   transient classes: 0
#> 
#> Use summary(x) for the per-state table, plot(x) for the heatmap.
# \donttest{
summary(cs)
#> Chain structure summary  [3 states, 1 classes]
#>   irreducible: TRUE   aperiodic: TRUE   regular: TRUE   reversible: FALSE
#> 
#>       state classification period persistence return_probability sojourn_steps
#>      Active      recurrent      1      0.6976                  1          3.31
#>     Average      recurrent      1      0.6099                  1          2.56
#>  Disengaged      recurrent      1      0.4831                  1          1.93
#>  stationary_probability
#>                  0.3719
#>                  0.4431
#>                  0.1850
# }
```
