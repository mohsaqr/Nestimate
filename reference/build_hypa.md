# Detect Path Anomalies via HYPA

`build_hypa()` tests which paths in a set of sequences occur more or
less often than expected, with HYPA (LaRock et al., 2020). For an order
`k` it builds the order-`k` De Bruijn graph, whose nodes are runs of `k`
states and whose edges are paths of `k + 1` states, and compares each
path's count with a hypergeometric null model that keeps every node's
observed in- and out-strength. Paths that are rarer than the null model
allows are labelled under-represented, paths that are more common
over-represented, and the rest normal. Unobserved paths that the
observed transitions make possible are tested too, so a path that never
occurs can be under-represented. Several orders can be tested in one
call.

A path is over-represented when its last state follows its `k`-state
context more often than the null model predicts from the strengths of
the De Bruijn nodes alone (at order 2, from the first-order
transitions). Each path is tested on its own at `alpha`, so with many
possible paths some are flagged by chance: `p_adjust` corrects for that,
and with small counts the upper-tail test rejects more often than
`alpha` (see `p_adjust`). Whether the sequences need a higher order
overall is a separate question, tested by
[`markov_order_test()`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md);
[`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
selects an order and
[`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
builds a network of variable order.

## Usage

``` r
build_hypa(
  data,
  order = 2L,
  alpha = 0.05,
  min_count = 5L,
  p_adjust = "none",
  k = NULL
)

# S3 method for class 'net_hypa'
print(x, ...)

# S3 method for class 'net_hypa'
summary(
  object,
  n = 10L,
  type = c("all", "over", "under"),
  order_by = c("sig", "freq", "frequency", "ratio", "path"),
  ...
)
```

## Arguments

- data:

  The sequences: a wide data frame (one row per trajectory, position
  columns only, named as a numbered run), a character matrix with the
  same layout, a list of character vectors, a `tna` model, or a
  `netobject` built from sequences. For a `tna` model or a `netobject`,
  numeric state codes are converted to their labels. A frame with other
  columns beside the positions, a frame or matrix of only 0 and 1, a
  `netobject_group`, and a one-hot or association network raise an error
  of class `nestimate_bad_input`.

- order:

  The order or orders of the De Bruijn graph, whole numbers of at least
  1 without duplicates. Default `2`. Order `k` tests paths of `k + 1`
  states (`k` transitions), each conditioned on the `k`-state context
  before its last state. With several orders, one layer is built per
  order, in ascending order whatever the order given; an order with no
  paths (longer than the sequences allow) is dropped without a message,
  and an error is raised when none is left. The network fields
  (`weights`, `edges`, `adjacency`, `xi`, `nodes`, `meta`) describe the
  lowest order built, while `scores` holds every order.

- alpha:

  The significance threshold for the labels, a number in (0, 0.5).
  Default `0.05`. A path is under-represented when its adjusted
  lower-tail p-value (`p_adjusted_under`) is below `alpha`, otherwise
  over-represented when its adjusted upper-tail p-value
  (`p_adjusted_over`) is at most `alpha`, and normal otherwise.

- min_count:

  The smallest count a path needs to be labelled, a whole number of at
  least 1. Default `5`. An over-represented path must be observed at
  least `min_count` times, an under-represented path expected at least
  `min_count` times under the null model, so a path that never occurs
  can still be under-represented. Paths below it are normal.

- p_adjust:

  The correction of the p-values for multiple testing: any of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html), such
  as `"BH"` for Benjamini-Hochberg false discovery rate control. Default
  `"none"`, which compares the raw p-values with `alpha`, as LaRock et
  al. (2020) do. Under- and over-representation p-values are adjusted
  separately, within each order, over every possible path. With small
  counts the paper's upper-tail score `P(X > f)` rejects more often than
  `alpha` under the null, so neither the raw nor the corrected labels
  carry an exact type-I error or false discovery rate guarantee;
  `min_count` limits the problem.

- k:

  The former name of `order`, deprecated: when given, it replaces
  `order`, with a warning of class `deprecatedWarning`. Default `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `net_hypa`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_hypa`.

- n:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: the
  most paths to show and return per direction. Default `10`.

- type:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method:
  which anomalies to show, `"all"` (default), `"over"` or `"under"`.

- order_by:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: how
  the paths are ranked within each direction before the top `n` are
  kept. `"sig"` (default) ranks by the tail probability of that
  direction, `"freq"` (or `"frequency"`) by observed count, `"ratio"` by
  observed / expected (highest first for over-represented paths, lowest
  first for under-represented ones), and `"path"` alphabetically. With
  `type = "all"` the returned rows of both directions are ranked
  together the same way, by how far the ratio is from 1 on the log scale
  (`abs(log(ratio))`) for `"ratio"`.

## Value

A `net_hypa` object (class `c("net_hypa", "cograph_network")`).
[`print()`](https://rdrr.io/r/base/print.html) reports the orders, the
number of paths tested and how many are over- and under-represented (per
order when there are several);
[`summary()`](https://rdrr.io/r/base/summary.html) prints and returns
the strongest anomalies;
[`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
gives the anomalous paths as pathway strings for
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
The object's fields:

- scores:

  A data frame with one row per possible path of each order: every
  observed path and every unobserved one the observed transitions permit
  (`observed = 0`). Within each order the over-represented paths come
  first (by decreasing ratio), then the under-represented (by increasing
  ratio), then the rest. Columns: `path`, the full state sequence
  (`"A -> B -> C"`); `from`, the context (the conditioning states);
  `to`, the next state; `observed`; `expected`, the count under the null
  model, whose propensity matrix keeps the observed in- and out-strength
  of every De Bruijn node (at order 2 this equals the first-order Markov
  expectation \\n(a b \cdot) n(\cdot b c) / n(\cdot b \cdot)\\);
  `ratio`, observed / expected; `p_under`, the hypergeometric lower-tail
  probability `P(X <= observed)`, and `p_value`, the same value under
  its older name; `p_over`, the strict upper tail
  `P(X > observed) = 1 - p_under`, so that `p_over <= alpha` is the
  paper's HYPA score of at least `1 - alpha`; `p_adjusted_under` and
  `p_adjusted_over`, the same after `p_adjust`; `anomaly` (`"over"`,
  `"under"` or `"normal"`); and `order`.

- ho_edges:

  The same table as `scores`, under the name the other higher-order
  builders use.

- over, under:

  The rows of `scores` labelled over- or under-represented.

- adjacency:

  The weighted adjacency matrix (path counts) of the lowest-order De
  Bruijn graph.

- weights:

  The cograph weight matrix of the lowest-order graph.

- xi:

  The fitted propensity matrix of the lowest-order graph.

- edges:

  The cograph edge data frame of the lowest-order graph.

- by_order:

  A list with the results of each order built, named by order.

- order:

  The orders built, ascending.

- k:

  The same as `order`, kept for backward compatibility.

- alpha:

  The significance threshold used.

- p_adjust:

  The correction used.

- n_anomalous, n_over, n_under:

  The number of anomalous, over- and under-represented paths, over all
  orders.

- n_edges:

  The number of possible paths tested, over all orders.

- nodes:

  A data frame (`id`, `label`, `name`) of the nodes of the lowest-order
  De Bruijn graph; `label` and `name` are the node's `k` states joined
  by the internal separator `"\x01"`, as are the dimnames of
  `adjacency`, `weights` and `xi`.

- directed:

  Always `TRUE`.

- meta:

  The cograph meta list of the lowest-order graph.

- node_groups:

  Always `NULL`.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) prints the
settings, the anomaly counts and the top `n` paths of each direction,
and returns invisibly a data frame of those paths, one row per path:
`order` (the De Bruijn order the path was tested at), `path`,
`observed`, `expected`, `ratio`, `p_tail` (the raw tail probability of
its direction, `p_over` for over-represented paths and `p_under` for
under-represented ones) and `direction` (`"over"` or `"under"`). When
nothing is anomalous it is a zero-row data frame with the same columns
except `order`.

## Details

Internal missing states split a trajectory, so no path spans a gap, and
pieces shorter than two states are dropped. State labels must not
contain the control characters `"\x01"` or `"\x02"`, which are used
internally as separators; such labels raise an error of class
`nestimate_hypa_bad_label`.

The possible paths of order `k` link each node with an observed outgoing
path to every node with an observed incoming path that continues it (its
last `k - 1` states are the other's first) through an observed
first-order transition; at order 1 they are the observed transitions.
The propensity matrix of the null model is fitted on them by continuous
iterative proportional fitting; when the fit does not reach its
tolerance, a warning of class `nestimate_hypa_no_converge` is raised.
The original Python implementation rounds the propensities inside every
sweep, so its expected counts can differ from these by a few
thousandths.

As in the original implementation, when an order has only one possible
path, that path's upper-tail p-value is 0, so it is labelled
over-represented once it is observed `min_count` times, even when its
count equals the expected one.

## References

LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
& Schweitzer, F. (2020). HYPA: Efficient Detection of Path Anomalies in
Time Series Data on Networks. *SDM 2020*, 460-468.

## Examples

``` r
seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
hyp <- build_hypa(seqs, order = 2)

# \donttest{
trajs <- list(c("A","B","C"), c("A","B","C"), c("A","B","C"),
              c("A","B","D"), c("C","B","D"), c("C","B","A"))
h <- build_hypa(trajs, order = 2)
print(h)
#> HYPA: Path Anomaly Detection
#>   Order(s):     2
#>   Edges:        6
#>   Anomalous:    0 (alpha=0.05, p_adjust=none)
#>     Over-repr:  0
#>     Under-repr: 0
# }
```
