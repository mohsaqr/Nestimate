# Prepare Event Log Data for Network Estimation

Converts event log data (actor, action, time) into wide sequence format
suitable for
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).
Automatically parses timestamps, detects sessions from time gaps, and
handles tie-breaking.

## Usage

``` r
prepare_data(
  data,
  actor,
  action,
  time = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900,
  custom_format = NULL,
  is_unix_time = FALSE,
  unix_time_unit = c("seconds", "milliseconds", "microseconds")
)
```

## Arguments

- data:

  Data frame with event log columns.

- actor:

  Character or character vector. Column name(s) identifying who
  performed the action (e.g. `"student"` or `c("student", "group")`). If
  missing, all data is treated as one actor.

- action:

  Character. Column name containing the action/state/code.

- time:

  Character or NULL. Column name containing timestamps. Supports
  ISO8601, Unix timestamps (numeric), and 40+ date/time formats. If
  NULL, row order defines the sequence. Default: NULL.

- order:

  Character or NULL. Column name for tie-breaking when timestamps are
  identical. If NULL, original row order is used. Default: NULL.

- session:

  Character, character vector, or NULL. Column name(s) for explicit
  session grouping (e.g. `"course"` or `c("course", "semester")`). When
  combined with `time`, sessions are further split by time gaps.
  Default: NULL.

- time_threshold:

  Numeric. Maximum gap in seconds between consecutive events before a
  new session starts. Only used when `time` is provided. Default: 900
  (15 minutes).

- custom_format:

  Character or NULL. Custom `strptime` format string for parsing
  timestamps. Default: NULL (auto-detect).

- is_unix_time:

  Logical. If TRUE, treat numeric time values as Unix timestamps.
  Default: FALSE (auto-detected for numeric columns).

- unix_time_unit:

  Character. Unit for Unix timestamps: `"seconds"`, `"milliseconds"`, or
  `"microseconds"`. Default: `"seconds"`.

## Value

A list with class `"nestimate_data"` containing:

- sequence_data:

  Data frame in wide format (one row per session, columns T1, T2, ...).

- long_data:

  The processed long-format data with session IDs.

- meta_data:

  Session-level metadata (session ID, actor).

- time_data:

  Parsed time values in wide format (if time provided).

- statistics:

  List with total_sessions, total_actions, max_sequence_length,
  unique_actors, etc.

## See also

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
[`prepare_onehot`](https://mohsaqr.github.io/Nestimate/reference/prepare_onehot.md)

## Examples

``` r
# \donttest{
df <- data.frame(
  student = rep(1:3, each = 5),
  code = sample(c("read", "write", "test"), 15, replace = TRUE),
  timestamp = seq.POSIXt(as.POSIXct("2024-01-01"), by = "min", length.out = 15)
)
prepared <- prepare_data(df, actor = "student", action = "code",
                         time = "timestamp")
net <- build_network(prepared$sequence_data, method = "relative")
# }
```
