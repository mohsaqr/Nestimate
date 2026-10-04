# The session behind each sequence

Names every sequence of a network or of a clustering by the columns it
was built from. A network built from long data with `actor` and
`session` has one sequence per actor-session; its rows are ordered by
the grouping, not by the input, and a fitted clustering reports its
assignments in that same order. `session_ids()` returns the key that
joins them back to the input data, so no label has to be parsed.

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

  A `netobject` from
  [`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  on long data, or a `net_mmm`
  ([`build_mmm`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md))
  or `net_clustering`
  ([`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md))
  fitted on such a network.

- ...:

  Unused.

## Value

A data frame with one row per sequence, in the row order of the
network's `$data` (and of the fit's assignments):

- sequence:

  Integer row number of the sequence.

- actor and session columns:

  The `actor` and `session` columns given to
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  under their own names and with their own values.

- session_label:

  The readable label of the sequence. With `time`, sessions split at
  time gaps carry a `" s<n>"` suffix, so this column separates them.

- cluster:

  For a `net_mmm` or `net_clustering`: the assigned cluster (integer).

- posterior:

  For a `net_mmm`: the posterior probability of the assigned cluster.

## Errors

Raises `nestimate_no_session_ids` when `x` carries no per-sequence
metadata: a network built from wide data, a fit on wide data or on a
`tna` model, or a fit made before Nestimate 0.9.6 (refit it). Raises
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
