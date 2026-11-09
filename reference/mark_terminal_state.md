# Mark terminal-NA cells with an explicit state label

`mark_terminal_state()` fills every missing cell after each sequence's
last observed state with an explicit state, `"End"` by default, and
returns the wide table. This turns right-censored sequences (drop-out,
course end) into data for an absorbing-chain model: in a transition
network built from the result with
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
the marked state is entered once a sequence stops and never left.
Internal and leading missing cells are left as they are, and so are rows
observed to the last position or not at all.

[`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
on that network gives the mean time to absorption, the absorption
probabilities and each state's classification.
[`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
is not the right summary for an absorbing chain: its stationary
distribution collapses onto the absorbing state.

## Usage

``` r
mark_terminal_state(data, state = "End", cols = NULL)
```

## Arguments

- data:

  The sequences: a matrix or data frame with one row per actor and one
  column per position, holding state labels, with `NA` for a missing
  observation. Without `cols`, the frame must be one column or a
  numbered run of positions (`T1`, `T2`, ...); other columns beside them
  raise an error of class `nestimate_bad_input`.

- state:

  The label written into the trailing missing cells, a single non-empty
  string. Default `"End"`. If the label already occurs in `data`, a
  warning of class `nestimate_state_collision` is given and `_1`, `_2`,
  ... is appended until it is unique, so the absorbing marker is never
  confused with an observed state.

- cols:

  The columns holding the positions, in order; only they are returned. A
  name not in `data` raises an error of class `nestimate_bad_input`.
  Default `NULL` (every column).

## Value

A data frame of character columns with the shape of `data` (or of its
`cols`), the trailing missing cells filled with the marker. The label
used is kept in the `"terminal_state"` attribute, which differs from
`state` when it had to be made unique. The columns are marked as
positions, so
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
reads them without `state_cols`.

## Details

The marked state is absorbing because each sequence that reaches it
stays there to the last column, so its only observed move is to itself
(`P[state, state] = 1` in a `"relative"` network). That needs at least
one sequence with two or more trailing missing cells: when every
sequence stops at most one cell early, the state has no outgoing
transitions and its row in the network is zero.

Unlike `build_network(end = )`, which adds an end state only in the cell
after each sequence's last observed state, with no outgoing edges and no
self-loop, this function fills all the trailing cells.

## See also

[`actor_endpoints()`](https://pak.dynasite.org/Nestimate/reference/actor_endpoints.md),
[`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md),
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

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
