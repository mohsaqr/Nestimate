# Split each actor's sequence into phases

`mark_phases()` divides the events of each actor into `n` consecutive
phases of nearly equal length and records the phase of every event in a
new column. The result can be passed to
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `group = "phase"` to estimate one network per phase, for example to
compare the first and second halves of each session.

The events of an actor are ordered by `time` (or kept in row order) and
numbered \\k = 1, \dots, L\\. Event \\k\\ belongs to phase \\\lceil k n
/ L \rceil\\, so the phases of one actor differ in length by at most one
event, and with an odd number of events and `n = 2` the middle event
belongs to the second phase. An actor with fewer than `n` events has no
events in some phases.

## Usage

``` r
mark_phases(
  data,
  actor = NULL,
  time = NULL,
  n = 2L,
  labels = NULL,
  name = "phase"
)
```

## Arguments

- data:

  An event log: a data frame with one row per event.

- actor:

  The column identifying whose events form one sequence (a session, a
  student). Default `NULL`, which treats the whole log as one sequence.

- time:

  The column ordering the events within each actor: numeric, a `Date` or
  a `POSIXct`. Default `NULL`, which keeps the row order.

- n:

  The number of phases, a whole number of at least 2. Default `2`.

- labels:

  The labels of the phases, in order, `n` distinct non-empty strings.
  Default `NULL`, which uses `"1"`, `"2"`, ..., `"n"`.

- name:

  The name of the new column. Default `"phase"`. A column of that name
  already in `data` raises an error of class `nestimate_bad_input`.

## Value

`data` with one added character column, `name`, holding each event's
phase label; the rows and the other columns are unchanged and keep their
order. An absent `actor` or `time` column, missing values in either, a
`time` column of another type, or an invalid `n` or `labels` raises an
error of class `nestimate_bad_input`.

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for the per-phase networks,
[`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
and
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
for comparing them.

## Examples

``` r
halves <- mark_phases(human_long, actor = "session_id",
                      time = "order_in_session", n = 2,
                      labels = c("early", "late"))
head(halves)
#>   message_id   project   session_id  timestamp session_date      code
#> 1       3439 Project_7 0086cabebd15 1772661600   2026-03-05   Specify
#> 2       3439 Project_7 0086cabebd15 1772661600   2026-03-05   Command
#> 3       3439 Project_7 0086cabebd15 1772661600   2026-03-05   Specify
#> 4       3440 Project_7 0086cabebd15 1772661600   2026-03-05 Interrupt
#> 5       3442 Project_7 0086cabebd15 1772661600   2026-03-05    Verify
#> 6       3444 Project_7 0086cabebd15 1772661600   2026-03-05   Specify
#>         cluster code_order order_in_session phase
#> 1     Directive          1                1 early
#> 2     Directive          2                2 early
#> 3     Directive          3                3 early
#> 4 Metacognitive          1                4 early
#> 5    Evaluative          1                7 early
#> 6     Directive          1               10 early
net <- build_network(halves, method = "relative", actor = "session_id",
                     action = "cluster", time = "order_in_session",
                     group = "phase")
#> Metadata aggregated per session: ties resolved by first occurrence in 'session_date' (2 sessions), 'code' (125 sessions)
#> Metadata aggregated per session: ties resolved by first occurrence in 'session_date' (2 sessions), 'code' (114 sessions)
net
#> Group Networks (2 groups, group_col: phase)
#> 
#>   Group  Nodes  Edges  Weights
#>   early  3      9      [0.105, 0.637]
#>   late   3      9      [0.131, 0.575]
```
