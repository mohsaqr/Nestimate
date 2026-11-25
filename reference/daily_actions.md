# Simulated learning actions over seven days

`daily_actions` is a simulated log of the learning actions of 30
students over seven days: `plan`, `read`, `write`, `discuss` and
`review`. Each student has between 51 and 98 actions, spread unevenly
over the days: a student acts at least once and at most 41 times on one
day (median 9), and a student's busiest day has between 2 and 31 times
as many actions as the quietest. Each student favours some actions over
others; the most frequent action of a student accounts for 32% to 78% of
their actions. The rows of each student are in time order. The data suit
sequence plots and networks in which every actor has a different number
of actions per period, as in `sequence_plot(periods = "day")`.

## Usage

``` r
daily_actions
```

## Format

A data frame with 2,280 rows and 3 character columns:

- student:

  Student identifier, `"s01"` to `"s30"`.

- day:

  The day of the action, `"day 1"` to `"day 7"`.

- action:

  The learning action: `plan`, `read`, `write`, `discuss` or `review`.

## Source

Simulated with the script `data-raw/daily_actions.R` (seed 2026).

## Examples

``` r
sequence_plot(daily_actions, actor = "student", action = "action",
              periods = "day")

```
