# Qualitative structure of a discrete-time Markov chain

Computes properties that depend only on the transition matrix support,
not on any starting distribution: state classification, communicating
classes, periods, irreducibility / aperiodicity / regularity /
reversibility, hitting probabilities, and absorption analysis when
absorbing states exist.

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

  A `netobject`, `cograph_network`, `tna` model, transition matrix, or
  sequence data.frame (passed through
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  with `method = "relative"`). A `netobject_group` is also accepted and
  analysed constituent by constituent. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `chain_structure`, `chain_structure_group` or
  `summary_chain_structure`.

- normalize:

  Logical. If `TRUE` (default), rows of the transition matrix are
  renormalized to sum to 1 before analysis (see
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  for the same convention).

- tol:

  Numerical tolerance for the reversibility check (detailed balance) and
  for treating near-zero entries as zero when building the support graph
  (which drives `classification`, `communicating_classes`, `period`, and
  `hitting_probabilities`). It does **not** govern the absorbing-state
  test: a state is absorbing only when `P[i, i]` equals 1 to an internal
  fixed tolerance of `.Machine$double.eps^0.5`, independent of `tol` (so
  raising `tol` to ignore tiny transition probabilities never
  reclassifies a near-deterministic state as absorbing). Default
  `1e-10`.

- ...:

  For the S3 methods: further arguments passed to or from other methods.
  In `print.chain_structure_group()`: Forwarded to
  `print.chain_structure`. In `print.summary_chain_structure()`:
  Forwarded to `print.data.frame`.

- show_values:

  Logical. If `TRUE` (default), prints the numeric probability inside
  each cell. Set `FALSE` for large state spaces (n \> 10) where labels
  overlap.

- digits:

  Integer. Decimal places for in-cell labels.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `chain_structure` or `chain_structure_group`.

## Value

For a `netobject_group`, a `c("chain_structure_group", "list")`: a named
list holding one `chain_structure` per constituent network, with its own
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods.

Otherwise a `chain_structure` object: a list with elements

- `states`:

  Character vector of state names.

- `classification`:

  Named character vector. One of `"absorbing"`, `"recurrent"`,
  `"transient"` per state.

- `communicating_classes`:

  List of state-name vectors. Each sublist is a strongly connected
  component of the support graph.

- `recurrent_classes`:

  Subset of `communicating_classes` that are closed (no transitions
  leaving the class).

- `transient_classes`:

  Subset that are not closed.

- `absorbing_states`:

  Character vector of states with `P[i, i] = 1` (tested exactly, to
  within `.Machine$double.eps^0.5`; the user-facing `tol` does not relax
  this).

- `period`:

  Named integer vector. Period of each recurrent state; `NA` for
  transient states.

- `is_irreducible`:

  Logical. `TRUE` iff there is exactly one communicating class.

- `is_aperiodic`:

  Logical. `TRUE` iff every recurrent state has period 1.

- `is_regular`:

  Logical. `is_irreducible && is_aperiodic`.

- `is_reversible`:

  Logical or `NA`. `TRUE` iff the chain satisfies detailed balance
  against its stationary distribution. `NA` for non-irreducible chains
  (no unique stationary).

- `hitting_probabilities`:

  `n x n` matrix. `[i, j] = P(ever reach j starting from i)`, computed
  over the same `tol`-thresholded support graph that drives
  `classification` so the two are mutually consistent (a state
  classified `"absorbing"`/closed never shows hitting probability to
  states outside its class).

- `absorption_probabilities`:

  `n_transient x n_absorbing` matrix or `NULL` if no transient -\>
  absorbing pathway exists.
  `[i, j] = P(eventual absorption in j | start in i)`.

- `mean_absorption_time`:

  Named numeric vector or `NULL`. Expected number of steps until
  absorption from each transient state.

- `P`:

  The (possibly normalized) transition matrix used.

In `print.chain_structure()`, `print.chain_structure_group()` and
`print.summary_chain_structure()`: `x` invisibly.

In `plot.chain_structure()`: A `ggplot` object.

In `summary.chain_structure()`: A `data.frame` with one row per state,
of class `c("summary_chain_structure", "data.frame")`, carrying the
chain-level flags (`is_regular`, `is_irreducible`, `is_aperiodic`,
`is_reversible`, `n_classes`, `absorbing_states`) as attributes, which
its [`print()`](https://rdrr.io/r/base/print.html) method shows as a
header. Columns as described above.

In `summary.chain_structure_group()`: A `data.frame` with columns
`group`, `state`, `classification`, `period`, `persistence`,
`return_probability`, `sojourn_steps`, plus `stationary_probability` if
all groups are irreducible and `mean_absorption_time` if any group has
absorbing states.

## Details

Built specifically as a diagnostic to run *before* trusting the output
of
[`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
or
[`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md).
Both implicitly assume a regular chain (irreducible + aperiodic) so that
the stationary distribution is unique and meaningful. Use `is_regular`
to check.

The fundamental-matrix absorption math follows Kemeny & Snell (1976);
the hitting-probability linear system follows Norris (1997).

## Plot colours

Cell colour encodes `P(ever reach j | start at i)`. The diagonal uses
the return-time convention (`P(return to j in >= 1 steps)`), matching
`markovchain::hittingProbabilities`. A non-irreducible chain shows zero
off-block entries – visual evidence of one-way doors between behavioural
phases. An absorbing chain shows a column of 1's for the absorbing
state.

## Summary columns

Columns are ordered for readability: identifiers first, classification
second, dynamic per-state metrics last.

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
