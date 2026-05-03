# Mark leading-NA cells with an explicit state label.

Mirror of
[`mark_terminal_state()`](https://mohsaqr.github.io/Nestimate/reference/mark_terminal_state.md)
for *left-censored* sequence data. Replaces every cell *before* each
row's first observed state with the label given by `state`. The
resulting chain has a structurally *recurrent* "Start" state that
everyone enters from — useful for cohort-entry analyses where students
join at different time points and you want a uniform pre-observation
marker.

## Usage

``` r
mark_first_state(data, state = "Start", cols = NULL)
```

## Arguments

- data:

  A wide-format matrix or data.frame (rows = actors, cols = time steps)
  of state labels with `NA` for missing observations.

- state:

  Character. Label to insert in leading-NA cells. Default `"Start"`.

- cols:

  Optional state-column names; otherwise all columns.

## Value

A `data.frame` of the same shape as `data` with leading NAs filled by
`state`.

## Details

Unlike
[`mark_terminal_state()`](https://mohsaqr.github.io/Nestimate/reference/mark_terminal_state.md),
the marked state is **not** absorbing in the resulting transition matrix
— every transition from "Start" goes to one of the original states (the
actor's first observed state), and the "Start" row is row-stochastic
exactly as the data dictates.

## See also

[`mark_terminal_state()`](https://mohsaqr.github.io/Nestimate/reference/mark_terminal_state.md),
[`actor_endpoints()`](https://mohsaqr.github.io/Nestimate/reference/actor_endpoints.md)

## Examples

``` r
M <- mark_first_state(trajectories, state = "Start")
# In a chain built from M, "Start" is a transient entry point.
```
