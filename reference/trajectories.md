# Student Engagement Trajectories

Wide-format sequences of student engagement over 15 consecutive
observations, one row per student and one column per observation. The
states are `"Active"`, `"Average"` and `"Disengaged"`. A missing value
occurs only as a run at the end of a row, which marks drop-out; two rows
are missing throughout.

## Usage

``` r
trajectories
```

## Format

A character matrix with 138 rows and 15 columns named `"1"` to `"15"`,
in time order. Each entry is `"Active"`, `"Average"`, `"Disengaged"` or
`NA`. The trailing `NA` runs are what the right-censored sequence verbs
([`actor_endpoints()`](https://pak.dynasite.org/Nestimate/reference/actor_endpoints.md),
[`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md))
are built to read.

## Examples

``` r
sequence_plot(trajectories, main = "Engagement trajectories")

sequence_plot(trajectories, k = 3)

sequence_plot(trajectories, type = "distribution")

```
