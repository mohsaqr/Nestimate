# Tidy per-actor endpoint summary of a wide-format sequence dataset

`actor_endpoints()` returns one row per sequence (actor) of a
wide-format sequence table, with its first and last observed states, the
positions where they occur, the number of observed cells and whether the
sequence ends early. A sequence ends early (`dropped_out`) when its last
observed state comes before the last column, so every cell after it is
missing: the *terminal-NA* pattern of drop-out, course end or censoring.
Internal gaps do not count.

The flag records only where the observations stop; it cannot say why. A
sequence recorded over fewer positions than the others looks the same as
one that left.
[`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md)
turns the terminal-NA cells into an explicit end state, and
[`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
then describes the chain with that state.

## Usage

``` r
actor_endpoints(data, cols = NULL)
```

## Arguments

- data:

  The sequences: a matrix or data frame with one row per actor and one
  column per position, holding state labels, with `NA` for a missing
  observation. Without `cols`, every column is read as a position, so
  the frame must be one column or a numbered run of positions (`T1`,
  `T2`, ...); a frame with other columns beside them, such as an id,
  raises an error of class `nestimate_bad_input`. A table from
  [`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md)
  or
  [`mark_first_state()`](https://pak.dynasite.org/Nestimate/reference/mark_first_state.md)
  is read as it is, because those functions mark their columns as
  positions.

- cols:

  The columns holding the positions, in order; the other columns are
  ignored. A name not in `data` raises an error of class
  `nestimate_bad_input`. Default `NULL` (every column).

## Value

A data frame with one row per actor and the columns:

- `actor`:

  The row name when `data` has row names other than `1`, `2`, ...;
  otherwise the row number.

- `first_state`, `last_state`:

  The first and the last observed state.

- `first_step`, `last_step`:

  The positions (column numbers among the position columns) of the first
  and the last observed state.

- `n_observed`:

  The number of observed cells.

- `dropped_out`:

  `TRUE` when `last_step` is before the last column, so the row ends in
  missing cells.

A row with no observed state has `NA` states and steps, `n_observed` 0
and `dropped_out` `FALSE`.

## See also

[`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md),
[`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)

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
