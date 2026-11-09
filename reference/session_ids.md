# The session behind each sequence

`session_ids()` names every sequence of a network, or of a clustering
fitted on one, by the actor and session values it was built from, and
returns one row per sequence. A network built from long data with
`actor` and `session` holds one sequence per actor-session; its rows are
ordered by the grouping, not by the input, and a fitted clustering
reports its assignments in that same order. `session_ids()` returns the
key that joins them back to the input data, so no label has to be
parsed. For a clustering it adds each sequence's cluster, and for a
mixed Markov model also the posterior probability of that cluster.

## Usage

``` r
session_ids(x, ...)

# Default S3 method
session_ids(x, ...)

# S3 method for class 'netobject'
session_ids(x, ...)

# S3 method for class 'net_mmm'
session_ids(x, ...)

# S3 method for class 'net_clustering'
session_ids(x, ...)
```

## Arguments

- x:

  A `netobject` built by
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  from long data, or a `net_mmm`
  ([`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md))
  or `net_clustering`
  ([`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md))
  fitted on such a network. Any other object raises an error of class
  `nestimate_no_session_ids` (see *Errors*).

- ...:

  Not used.

## Value

A data frame with one row per sequence, in the row order of the
network's sequence data (and of the fit's assignments):

- sequence:

  The integer row number of the sequence.

- actor and session columns:

  The `actor` and `session` columns given to
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  under their own names and with their own values.

- session_label:

  The readable label of the sequence. With `time`, the sessions split at
  time gaps carry an `" s<n>"` suffix, so this column tells them apart.

- cluster:

  For a `net_mmm` or `net_clustering`: the assigned cluster (integer).

- posterior:

  For a `net_mmm`: the posterior probability of the assigned cluster.

## Errors

Raises an error of class `nestimate_no_session_ids` when `x` carries no
per-sequence metadata (a network built from wide data, a fit on wide
data or on a `tna` model, or a fit made before Nestimate 0.9.6, which
must be refitted) and when `x` is of another class, such as a
`netobject_group`. Raises an error of class
`nestimate_session_ids_misaligned` when the metadata and the sequences
differ in number.

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`build_mmm`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md),
[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)

## Examples

``` r
events <- data.frame(
  student = rep(c("s1", "s2", "s3"), each = 8),
  step    = rep(c("a", "b"), each = 4, times = 3),
  action  = sample(c("read", "write", "test"), 24, replace = TRUE)
)
net <- build_network(events, actor = "student", session = "step",
                     action = "action", method = "relative")
session_ids(net)
#>   sequence student step session_label
#> 1        1      s1    a        s1 | a
#> 2        2      s2    a        s2 | a
#> 3        3      s3    a        s3 | a
#> 4        4      s1    b        s1 | b
#> 5        5      s2    b        s2 | b
#> 6        6      s3    b        s3 | b

# \donttest{
fit <- build_mmm(net, k = 2, n_starts = 2, seed = 1)
session_ids(fit)
#>   sequence student step session_label cluster posterior
#> 1        1      s1    a        s1 | a       2 0.9999906
#> 2        2      s2    a        s2 | a       1 0.9999245
#> 3        3      s3    a        s3 | a       1 0.9999830
#> 4        4      s1    b        s1 | b       2 0.8881216
#> 5        5      s2    b        s2 | b       1 1.0000000
#> 6        6      s3    b        s3 | b       1 0.9999990
# }
```
