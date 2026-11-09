# Prepare Event Log Data for Network Estimation

`prepare()` turns an event log (one row per event: the action, and
optionally who did it and when) into sequences, one row per session and
one column per position (`T1`, `T2`, ...). It parses the timestamps,
orders the events, splits each actor's events into sessions by explicit
session columns and by time gaps, and keeps every other column as one
row of metadata per session.
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
runs the same preparation itself when given `action`, with the same
arguments, so `prepare()` is for inspecting the sequences or passing
them to the sequence verbs.

## Usage

``` r
prepare(
  data,
  actor,
  action,
  time = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900,
  custom_format = NULL,
  is_unix_time = FALSE,
  unix_time_unit = c("seconds", "milliseconds", "microseconds"),
  timezone = "UTC"
)

# S3 method for class 'nestimate_data'
print(x, ...)
```

## Arguments

- data:

  The event log: a data frame with one row per event.

- actor:

  The column, or columns, identifying who performed each event (e.g.
  `"student"` or `c("student", "group")`); several columns are combined,
  each observed combination being one actor. Missing or `NULL` treats
  all events as one actor.

- action:

  The column holding each event's action (state or code).

- time:

  The timestamp column. Default `NULL`, in which case the row order (or
  `order`) gives the order of the events and no session is split by
  time. Numeric values are read as Unix timestamps, and POSIXct and
  POSIXlt values are used as they are, in their own time zone. Text may
  be ISO 8601 (with `Z` or an offset), `YYYY-MM-DD HH:MM[:SS]` with `-`,
  `/` or `.`, day-first or month-first dates, dates with month names, or
  compact `YYYYMMDDHHMMSS`, more than 40 layouts in all; text that
  matches none but reads as a number is taken as a Unix timestamp. A
  value must match a layout completely, with a four-digit year, so a
  two-digit year raises an error of class `nestimate_bad_input` instead
  of being misread, as does any value no layout reads. Day-first or
  month-first is decided once for the whole column: any first field
  above 12 means day-first, any second field above 12 means month-first.
  When the column cannot decide, values are read day-first with a
  warning of class `nestimate_ambiguous_dates`; `custom_format` settles
  it. An empty or missing timestamp is not an error: that event is
  sorted after its session's timed events and starts a session of its
  own.

- order:

  The column that orders the events: with `time`, it breaks ties between
  events with the same timestamp; without `time`, it gives the order of
  each session's events. Default `NULL` (the row order).

- session:

  The column, or columns, identifying explicit sessions (e.g. `"course"`
  or `c("course", "semester")`). Each actor's events in one session form
  a sequence, further split by time gaps when `time` is set. Default
  `NULL`.

- time_threshold:

  The longest gap, in seconds, between consecutive events of one
  session; a longer gap starts a new session, numbered `s1`, `s2`, ...
  within its actor and session. Used only with `time`. `FALSE` turns the
  split off, so each actor (or actor-session) forms one sequence however
  long the gaps. Default `900` (15 minutes).

- custom_format:

  A [`strptime()`](https://rdrr.io/r/base/strptime.html) format tried
  before the built-in layouts; values it does not match fall through to
  them. Default `NULL`.

- is_unix_time:

  Whether `time` holds Unix timestamps. `TRUE` reads a text column as
  numbers and never as a date layout (a 12- or 14-digit value would
  otherwise match the compact `YYYYMMDDHHMM[SS]` layouts); a value that
  is not a number then raises an error of class `nestimate_bad_input`.
  With `FALSE` (default), text is matched to the date layouts first and
  read as a Unix timestamp only when it matches none. A numeric `time`
  column is always read as Unix timestamps.

- unix_time_unit:

  The unit of Unix timestamps: `"seconds"` (default), `"milliseconds"`
  or `"microseconds"`. A numeric `time` column that gives dates after
  the year 2500 raises a warning of class `nestimate_time_unit`, the
  usual sign of the wrong unit.

- timezone:

  The Olson time zone (see
  [`OlsonNames()`](https://rdrr.io/r/base/timezones.html)) in which
  timestamps without zone information are read and Unix timestamps are
  expressed. Timestamps that end in `Z`, `UTC`, `GMT` or a numeric
  offset such as `+02:00` are converted from that offset, so parsing
  does not depend on the machine's time zone. Default `"UTC"`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `nestimate_data` object.

- ...:

  Not used.

## Value

A `nestimate_data` object.
[`print()`](https://rdrr.io/r/base/print.html) gives the numbers of
sessions and actions, the longest sequence, the number of actors when
`actor` was given and whether times were parsed, and returns its input
invisibly.
[`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md),
[`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
and
[`sequence_metrics()`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
read it directly. It holds:

- sequence_data:

  The sequences as a data frame: one row per session, one column per
  position (`T1`, `T2`, ...), padded with `NA`.

- long_data:

  The events in their sorted order, with the parsed time and the working
  columns that identify each event's session and position
  (`.session_id`, `.session_label`, `.sequence`, ...).

- meta_data:

  One row per session, in the row order of `sequence_data`:
  `.session_id`, `.session_label`, the actor column (with several, their
  combined label `.actor`, joined by `"-"`, or by the first of `"|"`,
  `" | "`, `" / "`, `" :: "` that keeps distinct actors apart when
  values contain `"-"`), the session column(s) under their own names,
  and every other column aggregated per session.

- time_data:

  The parsed times in the layout of `sequence_data` (columns `time_T1`,
  `time_T2`, ...), in `timezone` (or the zone of a date-time `time`
  column), or `NULL` without `time`.

- statistics:

  A list with `total_sessions`, `total_actions` and
  `max_sequence_length`, plus `unique_actors` (distinct combinations of
  the `actor` columns) only when `actor` was given (without it, every
  event belongs to one actor).

## Details

Within each actor and session, events are sorted by time and then by
`order`. Sessions are identified by the observed combinations of the
`actor` and `session` values, so identifiers containing separator
characters stay distinct and high-cardinality identifiers cannot
overflow. A missing value in any of these columns raises an error of
class `nestimate_bad_input`: drop or relabel those events first.

Every other column is aggregated into the session metadata: numeric
columns by their mean, other columns by their most frequent value. A tie
goes to the tied value that occurs first in the session, and a message
names the columns where that happened.

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
set.seed(1)
df <- data.frame(
  student = rep(1:3, each = 5),
  code = sample(c("read", "write", "test"), 15, replace = TRUE),
  timestamp = seq.POSIXt(as.POSIXct("2024-01-01"), by = "min", length.out = 15)
)
prepared <- prepare(df, actor = "student", action = "code",
                    time = "timestamp")
prepared
#> Prepared Data for Network Estimation
#>   Sessions: 3  |  Actions: 15  |  Max length: 5
#>   Actors: 3
#>   Time data: available

# build_network() runs the same preparation itself
net <- build_network(df, actor = "student", action = "code",
                     time = "timestamp", method = "relative")
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.200, 0.500]  |  mean: 0.333
#> 
#>   Weight matrix:
#>          read  test write
#>   read  0.400 0.200 0.400
#>   test  0.500 0.250 0.250
#>   write 0.333 0.333 0.333 
#> 
#>   Initial probabilities:
#>   test          0.667  ████████████████████████████████████████
#>   read          0.333  ████████████████████
#>   write         0.000  
```
