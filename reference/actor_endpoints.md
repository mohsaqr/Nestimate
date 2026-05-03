# Tidy per-actor endpoint summary of a wide-format sequence dataset.

For each actor (row), reports the first and last observed states, the
time indices at which they appear, the number of observed steps, and a
`dropped_out` flag that is `TRUE` when the actor has a *terminal-NA*
pattern (after the final observed step, every remaining cell is `NA`).

## Usage

``` r
actor_endpoints(data, cols = NULL)
```

## Arguments

- data:

  A wide-format matrix or data.frame where rows are actors and columns
  are time steps. Cells are state labels; `NA` represents missing
  observations. If `data` is a `data.frame`, non-character/-factor
  columns (e.g. an `id` column) are dropped via the `cols` argument.

- cols:

  Optional character vector of state-column names. If `NULL` (default)
  every column is treated as a state column.

## Value

A tidy `data.frame` with one row per actor and columns:

- `actor`:

  Row number (or row name if present).

- `first_state`:

  First non-NA state.

- `last_state`:

  Last non-NA state.

- `first_step`:

  Column index of the first observed state.

- `last_step`:

  Column index of the last observed state.

- `n_observed`:

  Number of non-NA cells.

- `dropped_out`:

  `TRUE` iff every cell after `last_step` is `NA` and
  `last_step < ncol(data)`.

## See also

[`mark_terminal_state()`](https://mohsaqr.github.io/Nestimate/reference/mark_terminal_state.md),
[`chain_structure()`](https://mohsaqr.github.io/Nestimate/reference/chain_structure.md)

## Examples

``` r
actor_endpoints(trajectories) |> head()
#>   actor first_state last_state first_step last_step n_observed dropped_out
#> 1     1      Active     Active          1        12         12        TRUE
#> 2     2     Average    Average          1        15         15       FALSE
#> 3     3     Average     Active          1        15         15       FALSE
#> 4     4      Active    Average          1        15         15       FALSE
#> 5     5      Active    Average          1        15         15       FALSE
#> 6     6     Average    Average          1        15         15       FALSE
```
