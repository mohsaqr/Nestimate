# Track transition entropy over time in sliding windows

`entropy_trajectory()` tracks transition entropy over time: it forms the
transitions within each actor's sequence of an event log, pools them in
time order, and slides a window of a fixed number of transitions across
the stream, giving one entropy estimate per window. A `group` column
gives one trajectory per group. This turns transition entropy from a
snapshot into a process measure, as Krejtz et al. (2025) do when they
track gaze transition entropy through the phases of a task.

A falling trajectory means the transitions become more predictable
(routinisation), a rising one less predictable (exploration), and a
level shift marks a change of phase. Each window's estimate rests on
`window` transitions, and with the default `step` neighbouring windows
share four fifths of them, so their estimates are not independent: read
the smoothed trend of
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) rather than
single windows. For the whole-process value with credible intervals, use
[`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
and
[`entropy_bayes()`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md).

## Usage

``` r
entropy_trajectory(
  data,
  action = NULL,
  actor = NULL,
  time = NULL,
  group = NULL,
  window = 500L,
  step = NULL,
  base = 2
)

# S3 method for class 'net_entropy_trajectory'
print(x, digits = 3, ...)

# S3 method for class 'net_entropy_trajectory'
summary(object, ...)

# S3 method for class 'net_entropy_trajectory'
plot(
  x,
  normalized = FALSE,
  span = 0.4,
  title = "Transition entropy over time",
  iter = 200L,
  ci_level = 0.95,
  seed = NULL,
  ...
)
```

## Arguments

- data:

  An event log (a data frame with one row per event, read with
  `action`), or sequences: a network built from sequences (a
  `netobject`, a grouped `netobject_group`, a `tna` model), prepared
  data, or a wide table or character matrix with one row per sequence
  and one column per position. Sequences are read in position order:
  each row is one actor, the position is the time, and the groups of a
  `netobject_group` become `group`. A network whose data are not
  sequences of states (an association or one-hot network) raises an
  error of class `nestimate_bad_input`.

- action:

  The column of an event log holding each event's state. Default `NULL`,
  which reads `data` as sequences.

- actor:

  The column identifying sequences: transitions are formed only between
  consecutive events of the same actor (within the same group). Default
  `NULL`, which reads the data as one sequence.

- time:

  The timestamp column, which orders the events within each actor and
  places the windows on a time axis. Default `NULL`, which keeps the row
  order and uses the row number of the event that ends each transition
  as the axis. The column must be numeric, a `Date` or a `POSIXct`; a
  character timestamp fails, so parse it first.

- group:

  A column splitting the data into parallel trajectories (for example a
  condition or an achievement level), each with its own windows. Default
  `NULL` (one trajectory, labelled `"all"`). A missing value in it
  raises an error of class `nestimate_bad_input`.

- window:

  The number of transitions per window, at least 10. Default `500`.

- step:

  The number of transitions between the starts of consecutive windows.
  Default `NULL`, which uses `window %/% 5` (at least 1).

- base:

  The logarithm base, which sets the unit: `2` (default) for bits. Any
  finite number greater than 1; another value raises an error of class
  `nestimate_bad_input`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  `entropy_trajectory()` result.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html): decimal places of
  the entropies. Default `3`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: not
  used.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): an
  `entropy_trajectory()` result.

- normalized:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): whether
  to draw `entropy_norm` instead of the raw entropy. Default `FALSE`.

- span:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  loess span of the trend line. Default `0.4`. A group gets a trend line
  when `span` times its number of windows is at least 4 (10 windows at
  the default); a group with fewer windows is drawn as its windows
  joined by a line, with a `nestimate_sparse_trajectory` warning naming
  it. A group with enough windows whose loess fit still warns is drawn
  the same way, with a `nestimate_trend_not_fitted` warning giving
  loess's message.

- title:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the plot
  title. Default `"Transition entropy over time"`.

- iter:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  number of bootstrap resamples behind the confidence band of each
  trend. Default `200`; `0` draws the trends without a band.

- ci_level:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  coverage of the band. Default `0.95`.

- seed:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): a seed
  for the resampling. It applies to this call only: the session's random
  state is restored afterwards. Default `NULL`.

## Value

A `net_entropy_trajectory` object.
[`print()`](https://rdrr.io/r/base/print.html) shows each group's number
of windows and its mean and range of entropy, and returns its input
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
data frame with one row per group: `group`, `windows`, the `mean`, `sd`,
`min` and `max` of the window entropies, the entropy of the `first` and
`last` windows, and their difference `change` (negative when the
transitions became more predictable).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` with a loess-smoothed trend per group and its bootstrap
confidence band, each group in an Okabe-Ito colour that no state uses
and named at the end of its trend; more groups than those colours are
drawn one panel each. The object's fields:

- trajectory:

  Data frame, one row per window: `group`, `window` (index within the
  group), `time` (the window's midpoint; a row number when there is no
  `time` column), `time_start` and `time_end` (the time, or row number,
  of the window's first and last transition), `n_transitions`,
  `n_states` (distinct states in the window), `entropy` (in `base` units
  per transition) and `entropy_norm` (divided by \\\log_b\\ of the
  window's number of states, floored at 2 states so the ceiling is never
  zero; in \\\[0, 1\]\\).

- window, step, base:

  The settings used.

- states:

  The states of the whole log.

- stream:

  The transitions in time order within each group (state codes, time,
  group and actor), from which
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) resamples.

- state_palette:

  The colours of the states (the input network's own, or Okabe-Ito dealt
  over the sorted states).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  groups in the Okabe-Ito colours these leave, so a group never takes a
  state's colour.

## Details

Each window's entropy is the empirical conditional entropy \\-\sum\_{ij}
(n\_{ij}/N) \log_b(n\_{ij}/n\_{i\cdot})\\, which weights rows by their
observed occupancy rather than by the stationary distribution from the
eigenvector. Within a short window the chain is often not ergodic
(absorbing fragments, unvisited states), where that eigenvector is
undefined or misleading; the empirical estimator is the standard
windowed choice and converges to the stationary entropy rate over long
stationary stretches.

The transitions of all actors are pooled in time order within each group
(ties broken by actor, then by input order). A window shorter than
`window` at the end of a stream is dropped; a group whose stream is
shorter than `window` gets one window covering all of it, with a
warning. Fewer than two transitions in all raise an error.

The band of [`plot()`](https://rdrr.io/r/graphics/plot.default.html) is
a bootstrap over actors (Davison & Hinkley, 1997): each resample draws a
group's actors with replacement, rebuilds its windows from their
transitions in time order and fits the same loess trend. The band is the
basic bootstrap interval of the trend at each time: the central
`ci_level` of the resampled trends, reflected about the observed trend.
It takes the spread of the resampled trends around the observed one as
the spread of the observed trend around the true one, so a systematic
shift of the resamples is carried to the other side of the trend rather
than kept. The plug-in entropy of a window is biased downward (Miller,
1955), by about \\(c - r) / (2N \ln b)\\ for \\c\\ observed transition
types from \\r\\ states in \\N\\ transitions, and resampled actors
(repeated transitions) shift it further down, so percentile limits tend
to sit below the trend. In a simulation (40 datasets of 100 actors with
20 steps of a stationary 4-state chain, `window = 300`, `step = 100`)
the 95% band covered the true entropy rate at mid-study in 38 of the 40,
missing once on each side. Whole actors are resampled because the
windows overlap (each transition sits in about `window / step` of them)
and an actor's transitions depend on each other, so a band from the
loess standard error, which treats the windows as independent, would be
too narrow. A resample whose trend cannot be fitted (too few windows, or
a loess warning) is left out of the band, and the plot warns how many
were, with class `nestimate_dropped_resamples`; a time reached by fewer
than half of the `iter` resamples has no band. A group with one actor
has no band (warning class `nestimate_single_actor_band`).

## References

Krejtz, K., Hughes, C.J., Stasiak, I., Duchowski, A., & Krejtz, I.
(2025). Real-Time Mobile Transition Matrix Entropy Based on Eye and Head
Movements. *ETRA '25*. doi:10.1145/3715669.3723128

Cover, T.M. & Thomas, J.A. (2006). *Elements of Information Theory*, 2nd
ed. Wiley.

Davison, A.C. & Hinkley, D.V. (1997). *Bootstrap Methods and Their
Application*. Cambridge University Press.

Miller, G.A. (1955). Note on the bias of information estimates. In H.
Quastler (Ed.), *Information Theory in Psychology* (pp. 95–100). Free
Press.

## See also

[`transition_entropy`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
for the whole-process snapshot,
[`entropy_bayes`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
for credible intervals on it.

## Examples

``` r
tr <- entropy_trajectory(group_regulation_long,
                         action = "Action", actor = "Actor",
                         time = "Time", group = "Achiever")
tr
#> Entropy Trajectory (9 states; window 500 transitions, step 100)
#> 
#>   High:        123 windows   mean 2.286   range [2.175, 2.394]
#>   Low:         124 windows   mean 2.294   range [2.143, 2.437]
#> 
#> Use plot() for the trajectory, summary() for per-group statistics.
summary(tr)
#>   group windows     mean         sd      min      max    first     last
#> 1  High     123 2.285872 0.04313978 2.175275 2.394429 2.309764 2.301800
#> 2   Low     124 2.294082 0.05803398 2.142537 2.436502 2.233659 2.347734
#>         change
#> 1 -0.007963838
#> 2  0.114074796
plot(tr)


# sequences: the network of the engagement trajectories, by course
net <- build_network(trajectories, method = "relative")
entropy_trajectory(net, window = 300, step = 50)
#> Entropy Trajectory (3 states; window 300 transitions, step 50)
#> 
#>   all:          29 windows   mean 1.212   range [1.034, 1.352]
#> 
#> Use plot() for the trajectory, summary() for per-group statistics.
```
