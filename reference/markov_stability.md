# Markov Stability Analysis

`markov_stability()` describes how stable each state of a Markov chain
is and returns a `net_markov_stability` object holding one row per
state: its persistence (the probability of staying, \\P\_{ii}\\),
stationary probability (\\\pi_i\\, the long-run share of steps spent in
it), mean recurrence time, sojourn time (the expected number of
consecutive steps in it, \\1/(1 - P\_{ii})\\), and the mean of its first
passage times to the other states and of theirs to it, taken from
[`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md).
It reads the same inputs as
[`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
and analyses a grouped network one network at a time.

Persistence and sojourn time describe staying in a state; the passage
times describe how long it takes to get from it to the other states, or
from them to it. All of them are computed from the estimated transition
matrix. On an irreducible chain, the mean recurrence time is
\\1/\pi_i\\. On a reducible chain, recurrence times are computed within
each closed class and are `Inf` for transient states; passage times to
states not reached with certainty are `Inf`. The stationary
probabilities are `NA` when there is more than one closed class, because
the stationary distribution is not unique.
[`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
identifies the classes; how far the underlying transition probabilities
move under resampling is shown by
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
markov_stability(x, normalize = TRUE)

# S3 method for class 'net_markov_stability'
print(x, ...)

# S3 method for class 'net_markov_stability_group'
print(x, ...)

# S3 method for class 'net_markov_stability'
summary(object, ...)

# S3 method for class 'net_markov_stability'
plot(
  x,
  metrics = c("persistence", "stationary_prob", "return_time", "sojourn_time",
    "avg_time_to_others", "avg_time_from_others"),
  combined = TRUE,
  ...
)

# S3 method for class 'net_markov_stability_group'
plot(
  x,
  metrics = c("persistence", "stationary_prob", "return_time", "sojourn_time",
    "avg_time_to_others", "avg_time_from_others"),
  ...
)
```

## Arguments

- x:

  The chain, as in
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md):
  a `netobject`, bare `cograph_network` or `tna` model with non-negative
  transition weights, a non-negative numeric transition matrix, or a
  wide sequence data frame or character matrix (a `"relative"` network
  is estimated from it). A `netobject_group` gives one result per group.
  An association network or negative weights raise an error of class
  `nestimate_not_markov`. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_markov_stability` or `net_markov_stability_group`.

- normalize:

  Whether to rescale rows that do not sum to 1 so that they do, with a
  warning, as for a `"frequency"` network. With `FALSE`, such rows raise
  an error. A state with no outgoing transitions (a zero row) raises an
  error either way. Default `TRUE`.

- ...:

  Not used ([`print()`](https://rdrr.io/r/base/print.html) of a group
  passes it on to each group's
  [`print()`](https://rdrr.io/r/base/print.html), which ignores it).

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_markov_stability`.

- metrics:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: the metrics to draw, one panel each, any of `"persistence"`,
  `"stationary_prob"`, `"return_time"`, `"sojourn_time"`,
  `"avg_time_to_others"` and `"avg_time_from_others"`. Default all six.

- combined:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  single network: `TRUE` (default) draws the metrics as panels of one
  plot; `FALSE` returns a named list of one-panel plots, one per metric,
  to print or save separately. Ignored for a group, whose
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) always draws
  one plot.

## Value

A `net_markov_stability` object; for a `netobject_group`, a
`net_markov_stability_group`, a named list holding one such object per
group. [`print()`](https://rdrr.io/r/base/print.html) shows the
per-state table and returns its input invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) prints the most
accessible state (the lowest mean passage time from the other states)
and the most persistent one (the longest sojourn time), and returns the
per-state table as a data frame.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` with one panel of bars per metric and the states ordered by
stationary probability; with `combined = FALSE`, invisibly, a list of
one-panel plots named by metric. For a group,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one panel
per metric with one row per state and one point per group (told apart by
colour and shape) joined by a line, the states ordered by their mean
stationary probability; more than 81 groups (9 colours times 9 shapes)
raise an error of class `nestimate_bad_input`. A `net_markov_stability`
holds:

- stability:

  The per-state table, one row per state, with the columns `state`,
  `persistence` (\\P\_{ii}\\), `stationary_prob` (\\\pi_i\\),
  `return_time` (the mean recurrence time), `sojourn_time`
  (\\1/(1-P\_{ii})\\), `avg_time_to_others` (the mean first passage time
  from the state to the others) and `avg_time_from_others` (the mean
  first passage time from the others to the state). `persistence` and
  `stationary_prob` are rounded to 4 decimals, the times to 2.

- mpt:

  The `net_mpt` from
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  the times come from.

## Details

**Sojourn time** is the expected number of consecutive steps spent in a
state before leaving it, \\1/(1-P\_{ii})\\, the mean of a geometric
distribution. A state with `persistence = 1` never leaves and has
`sojourn_time = Inf`.

**avg_time_to_others** is the mean passage time from the state to each
other state: a high value means the chain takes long to get anywhere
else from it.

**avg_time_from_others** is the mean passage time from each other state
to this one: a low value means the state is reached quickly from
elsewhere.

## References

Kemeny, J.G. and Snell, J.L. (1976). *Finite Markov Chains*.
Springer-Verlag.

## See also

[`passage_time`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)

## Examples

``` r
net <- build_network(as.data.frame(trajectories), method = "relative")
ms  <- markov_stability(net)
print(ms)
#> Markov Stability Analysis
#> 
#>       state persistence stationary_prob return_time sojourn_time
#>      Active      0.6976          0.3719        2.69         3.31
#>     Average      0.6099          0.4431        2.26         2.56
#>  Disengaged      0.4831          0.1850        5.41         1.93
#>  avg_time_to_others avg_time_from_others
#>                6.99                 5.84
#>                6.74                 3.20
#>                4.47                 9.16
# \donttest{
plot(ms)

# }
```
