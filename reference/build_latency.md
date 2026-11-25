# Estimate a latency network analysis model

`build_network(method = "latency")`, or for short `build_latency()`,
estimates a latency network analysis model from a message log. The edge
from sender A to sender B is the average weight of A's replies to B,
each reply weighted by its speed: 1 for an immediate reply, falling
linearly to `floor` at `cap` (Templeton et al., 2022).

Its edges are the replies of
[`build_reverse()`](https://pak.dynasite.org/Nestimate/reference/build_reverse.md):
within each sequence (an `actor`'s events, split by `session` or
`time_threshold`), each event points from its sender, the responder, to
the sender of the event before it. A reply that takes \\t\\ seconds
weighs \\1 - (1 - f) \min(t / c, 1)\\, with \\f\\ the `floor` and \\c\\
the `cap` in seconds, and the weights of a pair are averaged, or summed
with `aggregate = "sum"`. A sender who posts twice in a row is a
self-loop on the diagonal, weighted by the same rule. Shifting every
time by a constant leaves the weights unchanged, and so does a change of
time unit when the cap is a percentile.
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
adds each edge's number of replies (`n_replies`) and their median
latency in seconds (`median_latency`). `start`, `end` and
`params = list(weighted = TRUE)` raise an error of class
`nestimate_bad_input`.
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
resamples latency networks by their sequences and judges stability on
the speed score; the other resampling verbs, such as
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
raise an error of class `nestimate_latency_resampling` on one.

In
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
the four settings are named `latency_floor`, `latency_cap`,
`latency_aggregate` and `latency_scope`.

## Usage

``` r
build_latency(
  data,
  action,
  time,
  actor = NULL,
  session = NULL,
  order = NULL,
  time_threshold = 900,
  timezone = "UTC",
  floor = 0.5,
  cap = 0.95,
  aggregate = c("mean", "sum"),
  scope = c("global", "actor", "session", "responder"),
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  ...
)
```

## Arguments

- data:

  A message log: one row per message, `action` the sender.

- action:

  The column holding each message's sender.

- time:

  The timestamp column. Events are ordered by it, and a gap longer than
  `time_threshold` starts a new sequence.

- actor:

  The column naming who each row belongs to (a person, a team). In an
  event log, each actor's events form a sequence, split further by
  `session` or `time_threshold`; in one-hot data, each actor's rows.
  Default `NULL`.

- session:

  The column identifying sessions in an event log; each session is its
  own sequence, also without `actor`. Default `NULL`.

- order:

  The column that orders events sharing a timestamp; without it, those
  events keep their row order. Default `NULL`.

- time_threshold:

  The longest gap, in seconds, between two events of the same sequence
  in an event log; a longer gap starts a new sequence. `FALSE` turns the
  split off, so each actor (or actor-session) is one sequence. Default
  `900`.

- timezone:

  The Olson time zone for timestamps without an offset (passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md));
  timestamps with an offset, such as `...Z` or `+02:00`, are converted
  from it. Default `"UTC"`.

- floor:

  The weight of a reply that takes `cap` or longer, a single number from
  0 to 1. `1` ignores time (with `aggregate = "sum"`, the counts of
  [`build_reverse()`](https://pak.dynasite.org/Nestimate/reference/build_reverse.md));
  `0` gives replies at or past the cap no weight. Default `0.5`.

- cap:

  The time at which a reply's weight reaches `floor`. A number in (0,
  1\] is a percentile of the latencies of the replies to another sender,
  pooled over all of `data` (with `group`, over every group, so the
  groups share one cap); a string such as `"120s"` is a number of
  seconds. The resolved cap, in seconds, is stored on the network.
  Default `0.95`.

- aggregate:

  How the replies of a pair of senders make one weight: `"mean"`
  (default) averages their weights, the pair's typical reply speed;
  `"sum"` adds them, so an edge grows with both the speed and the number
  of replies.

- scope:

  With a percentile `cap`: the units that get their own cap. `"global"`
  (default) takes one cap from all of `data`; `"actor"` one per value of
  `actor`; `"session"` one per sequence; `"responder"` one per
  responding sender. A unit with fewer than 5 replies to another sender
  takes the global cap, and a message of class
  `nestimate_latency_fallback` says how many units did. A cap in seconds
  applies to every reply, so a scope other than `"global"` is then
  ignored with a message of class `nestimate_ignored_argument`.

- group:

  A column to split the data by: one network per value, returned as a
  `netobject_group`. Default `NULL`.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` rescales the weights to run from 0 to 1, over every cell
  for the transition methods (tna's rule, so with no zero cell the
  smallest weight becomes 0) and over the non-zero weights, into (0,
  1\], for the other methods, `"max"` divides by the largest absolute
  weight, `"rank"` replaces the non-zero weights by their ranks, and
  `"normalize"` makes each row's absolute weights sum to 1; an
  undirected network whose matrix it makes asymmetric becomes directed,
  so both of a pair's values are kept as edges. Zero weights stay zero
  under every option. Combine them as in `c("rank", "minmax")`. Default
  `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero. Default
  `0` (none).

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.
  Default `NULL`.

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md):
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)'s
  time options `unix_time_unit`, `is_unix_time` and `custom_format`.
  `method` cannot be given (error of class `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it. Beside the weights it keeps `initial`, `n_replies` and
`median_latency` (matrices laid out as the weights),
`latency_cap_seconds` (the resolved global cap), `latency_unit_caps`
(the per-unit caps of a non-global `scope`) and `latency_time`, the
event times in seconds.

## References

Templeton, E. M., Chang, L. J., Reynolds, E. A., Cone LeBeaumont, M. D.,
& Wheatley, T. (2022). Fast response times signal social connection in
conversation. *Proceedings of the National Academy of Sciences*, 119(4),
e2116915119.
[doi:10.1073/pnas.2116915119](https://doi.org/10.1073/pnas.2116915119)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_reverse()`](https://pak.dynasite.org/Nestimate/reference/build_reverse.md)
for the unweighted replies.

## Examples

``` r
build_network(group_regulation_long, method = "latency", actor = "Actor",
              action = "Action", time = "Time")
#> Reply Network (response-latency weighted) [directed]
#>   Weights: [0.694, 0.860]  |  mean: 0.741
#> 
#>   Strongest edges (10 of 78):
#>    from       to        weight
#>    monitor    adapt     0.860 
#>    plan       adapt     0.857 
#>    adapt      cohesion  0.821 
#>    coregulate adapt     0.793 
#>    monitor    synthesis 0.792 
#>    discuss    adapt     0.774 
#>    synthesis  plan      0.770 
#>    adapt      plan      0.769 
#>    coregulate synthesis 0.760 
#>    emotion    synthesis 0.755  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
#> 
#>   Latency: weight 1 at 0 s, falling to 0.5 at the cap of 476 s (0.95 percentile)  |  mean per pair  |  25533 replies
build_latency(group_regulation_long, action = "Action", time = "Time",
              actor = "Actor")
#> Reply Network (response-latency weighted) [directed]
#>   Weights: [0.694, 0.860]  |  mean: 0.741
#> 
#>   Strongest edges (10 of 78):
#>    from       to        weight
#>    monitor    adapt     0.860 
#>    plan       adapt     0.857 
#>    adapt      cohesion  0.821 
#>    coregulate adapt     0.793 
#>    monitor    synthesis 0.792 
#>    discuss    adapt     0.774 
#>    synthesis  plan      0.770 
#>    adapt      plan      0.769 
#>    coregulate synthesis 0.760 
#>    emotion    synthesis 0.755  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
#> 
#>   Latency: weight 1 at 0 s, falling to 0.5 at the cap of 476 s (0.95 percentile)  |  mean per pair  |  25533 replies
```
