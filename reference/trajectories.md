# Student engagement trajectories

`trajectories` is a dataset of the engagement states (`"Active"`,
`"Average"` and `"Disengaged"`) of 138 students over the 15 successive
courses of a programme. It is in wide format, one row per student and
one column per course. A missing value occurs only as a run at the end
of a row, which marks drop-out; two rows are missing throughout.

## Usage

``` r
trajectories
```

## Format

A character matrix with 138 rows and 15 columns named `"1"` to `"15"`,
the courses in the order they were taken. Each entry is `"Active"`,
`"Average"`, `"Disengaged"` or `NA`. The trailing `NA` runs are what the
right-censored sequence verbs
([`actor_endpoints()`](https://pak.dynasite.org/Nestimate/reference/actor_endpoints.md),
[`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md))
are built to read.

## Examples

``` r
sequence_plot(trajectories, main = "Engagement trajectories")

sequence_plot(trajectories, k = 3)

sequence_plot(trajectories, type = "distribution")

```
