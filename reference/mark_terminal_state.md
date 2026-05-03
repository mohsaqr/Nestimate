# Mark terminal-NA cells with an explicit state label.

Replaces every cell after each row's last observed state with the label
given by `state`, leaving non-terminal NAs untouched. The result, passed
to
[`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
yields a Markov chain in which the marked state is **absorbing** by
construction (`P[state, state] = 1`).

## Usage

``` r
mark_terminal_state(data, state = "End", cols = NULL)
```

## Arguments

- data:

  A wide-format matrix or data.frame (rows = actors, cols = time steps)
  of state labels with `NA` for missing observations.

- state:

  Character. Label to insert in terminal-NA cells. Default `"End"`.

- cols:

  Optional state-column names; otherwise all columns.

## Value

A `data.frame` of the same shape as `data` with terminal NAs filled by
`state`.

## Details

This is the small piece of pre-processing required to turn
right-censored sequence data into an absorbing-chain model. The chain on
the resulting matrix has one extra state (`state`) which is structurally
absorbing because every cell after the actor's last observed step has
been set to `state` — the chain stays there forever once entered.

Use
[`chain_structure()`](https://mohsaqr.github.io/Nestimate/reference/chain_structure.md)
on the result to compute mean absorption time, absorption probabilities,
and per-state classification. Note that
[`markov_stability()`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md)
is *not* the right summary for absorbing chains; its stationary
distribution will collapse to the absorbing state.

## See also

[`actor_endpoints()`](https://mohsaqr.github.io/Nestimate/reference/actor_endpoints.md),
[`chain_structure()`](https://mohsaqr.github.io/Nestimate/reference/chain_structure.md),
[`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)

## Examples

``` r
M <- mark_terminal_state(trajectories, state = "Dropout")
net <- build_network(M, method = "relative")
chain_structure(net)
#> Chain structure  [4 states, 2 communicating classes]
#>   irreducible: FALSE   aperiodic: TRUE   regular: FALSE   reversible: NA
#>   recurrent classes: 1   transient classes: 1
#>   absorbing states: Dropout
#> 
#> Use summary(x) for the per-state table, plot(x) for the heatmap.
```
