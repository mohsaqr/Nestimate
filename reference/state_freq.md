# Print, Plot, and Convert a state_freq Object

[`plot_state_frequencies()`](https://saqr.me/Nestimate/reference/plot_state_frequencies.md)
returns a `state_freq` object holding both the rendered chart and the
tidy frequency table. [`print()`](https://rdrr.io/r/base/print.html)
shows the table in the console *and* draws the chart on the active
graphics device,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the chart
alone, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the tidy table for downstream piping.

## Usage

``` r
# S3 method for class 'state_freq'
print(x, digits = 1, max_states = 20L, ...)

# S3 method for class 'state_freq'
plot(x, ...)

# S3 method for class 'state_freq'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `state_freq` object.

- digits:

  Number of decimal places for proportion / share columns. Default 1.

- max_states:

  Cap on rows shown per group in the per-state table (default 20); the
  surplus is folded into a single `"(+k more)"` row. The full, uncapped
  table is returned by `as.data.frame(x)`.

- ...:

  Unused.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly
(after printing the table and drawing the chart);
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns
`invisible(NULL)` after drawing;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the tidy `data.frame`, one row per (group, state) cell with columns
`group`, `state`, `count`, `proportion`.
