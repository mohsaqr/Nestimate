# Qualitative structure of a discrete-time Markov chain.

Computes properties that depend only on the transition matrix support,
not on any starting distribution: state classification, communicating
classes, periods, irreducibility / aperiodicity / regularity /
reversibility, hitting probabilities, and absorption analysis when
absorbing states exist.

## Usage

``` r
chain_structure(x, normalize = TRUE, tol = 1e-10)
```

## Arguments

- x:

  A `netobject`, `cograph_network`, `tna` model, transition matrix, or
  sequence data.frame (passed through
  [`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
  with `method = "relative"`).

- normalize:

  Logical. If `TRUE` (default), rows of the transition matrix are
  renormalized to sum to 1 before analysis (see
  [`passage_time()`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  for the same convention).

- tol:

  Numerical tolerance for the reversibility check (detailed balance) and
  for treating near-zero entries as zero when building the support
  graph. Default `1e-10`.

## Value

A `chain_structure` object: a list with elements

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

  Character vector of states with `P[i, i] = 1`.

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

  `n x n` matrix. `[i, j] = P(ever reach j starting from i)`.

- `absorption_probabilities`:

  `n_transient x n_absorbing` matrix or `NULL` if no transient -\>
  absorbing pathway exists.
  `[i, j] = P(eventual absorption in j | start in i)`.

- `mean_absorption_time`:

  Named numeric vector or `NULL`. Expected number of steps until
  absorption from each transient state.

- `P`:

  The (possibly normalized) transition matrix used.

## Details

Built specifically as a diagnostic to run *before* trusting the output
of
[`passage_time()`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
or
[`markov_stability()`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md).
Both implicitly assume a regular chain (irreducible + aperiodic) so that
the stationary distribution is unique and meaningful. Use `is_regular`
to check.

The fundamental-matrix absorption math follows Kemeny & Snell (1976);
the hitting-probability linear system follows Norris (1997).

## References

Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*.
Springer-Verlag.

Norris, J. R. (1997). *Markov Chains*. Cambridge University Press.

## See also

[`passage_time()`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md),
[`markov_stability()`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md),
[`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)

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
