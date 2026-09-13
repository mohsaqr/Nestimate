# Student Engagement Trajectories

Wide-format state sequences of student engagement over 15 weekly
observations. Each row is one student; columns `1`..`15` hold the
engagement state for that week. States: `"Active"`, `"Average"`,
`"Disengaged"`. Missing weeks are `NA`.

## Usage

``` r
trajectories
```

## Format

A character matrix with 138 rows and 15 columns, one row per student.
Columns are named `"1"`..`"15"` (the week). Entries are one of
`"Active"`, `"Average"`, `"Disengaged"`, or `NA`; `NA` runs at the end
of a row mark drop-out, which the right-censored sequence verbs
([`actor_endpoints`](https://saqr.me/Nestimate/reference/actor_endpoints.md),
[`mark_terminal_state`](https://saqr.me/Nestimate/reference/mark_terminal_state.md))
are built to read.

## Examples

``` r
# \donttest{
sequence_plot(trajectories, main = "Engagement trajectories")

sequence_plot(trajectories, k = 3)

sequence_plot(trajectories, type = "distribution")

# }
```
