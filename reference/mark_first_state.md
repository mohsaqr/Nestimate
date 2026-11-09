# Mark leading-NA cells with an explicit state label

`mark_first_state()` fills every missing cell before each sequence's
first observed state with an explicit state, `"Start"` by default, and
returns the wide table. It is the counterpart of
[`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md)
for *left-censored* data, such as students who join a course at
different times: a network built from the result shows into which states
the late entrants arrive. Internal and trailing missing cells are left
as they are, and so are rows observed from the first position or not at
all.

The marked state is a transient entry point, not an absorbing state:
once a sequence leaves it, it never returns. In a transition network
built from the result, its row holds the moves from it, which go to the
first observed states and, for sequences with two or more leading
missing cells, back to itself.

## Usage

``` r
mark_first_state(data, state = "Start", cols = NULL)
```

## Arguments

- data:

  The sequences: a matrix or data frame with one row per actor and one
  column per position, holding state labels, with `NA` for a missing
  observation. Without `cols`, the frame must be one column or a
  numbered run of positions (`T1`, `T2`, ...); other columns beside them
  raise an error of class `nestimate_bad_input`.

- state:

  The label written into the leading missing cells, a single non-empty
  string. Default `"Start"`. If the label already occurs in `data`, a
  warning of class `nestimate_state_collision` is given and `_1`, `_2`,
  ... is appended until it is unique, so the marker is never confused
  with an observed state.

- cols:

  The columns holding the positions, in order; only they are returned. A
  name not in `data` raises an error of class `nestimate_bad_input`.
  Default `NULL` (every column).

## Value

A data frame of character columns with the shape of `data` (or of its
`cols`), the leading missing cells filled with the marker. The label
used is kept in the `"leading_state"` attribute, which differs from
`state` when it had to be made unique. The columns are marked as
positions, so
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
reads them without `state_cols`.

## See also

[`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md),
[`actor_endpoints()`](https://pak.dynasite.org/Nestimate/reference/actor_endpoints.md)

## Examples

``` r
M <- mark_first_state(trajectories, state = "Start")
# In a chain built from M, "Start" is a transient entry point.
```
