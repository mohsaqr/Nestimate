# Plot method for `net_entropy_trajectory`

Raw per-window entropy as faint lines with a loess-smoothed trend per
group, Okabe-Ito coloured. Declining trend = routinization; level shifts
= phase changes.

## Usage

``` r
# S3 method for class 'net_entropy_trajectory'
plot(
  x,
  normalized = FALSE,
  span = 0.4,
  title = "Transition entropy over time",
  ...
)
```

## Arguments

- x:

  A `net_entropy_trajectory` object.

- normalized:

  Logical. Plot `entropy_norm` instead of raw bits (default `FALSE`).

- span:

  Numeric. Loess span (default `0.4`).

- title:

  Character. Plot title.

- ...:

  Ignored.

## Value

A ggplot object.
