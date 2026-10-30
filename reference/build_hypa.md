# Detect Path Anomalies via HYPA

Constructs a k-th order De Bruijn graph from sequential trajectory data
and uses a hypergeometric null model to detect paths with anomalous
frequencies. Paths occurring more or less often than expected under the
null model are flagged as over- or under-represented.

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

  A data.frame (rows = trajectories), list of character vectors, `tna`
  object, or `netobject` with sequence data. For `tna`/`netobject`,
  numeric state IDs are automatically converted to label names.

- order:

  Integer scalar or integer vector. Order(s) of the De Bruijn graph
  (default `2L`). An order of `k` tests paths of `k + 1` states (`k`
  transitions), each conditioned on the `k`-state context that precedes
  its last state. Orders must be unique. When a vector is supplied, one
  De Bruijn layer is built per order and the per-order results are
  stored in `$by_order` (named by order). The orders are sorted
  ascending internally, so the `cograph_network` slots (`$weights`,
  `$edges`, `$adjacency`, `$xi`, `$nodes`, `$meta`) always describe the
  network of the *lowest order that produced a layer* (a requested order
  with no edges is dropped from `$by_order`, `$order` and `$k`),
  regardless of the order in which the vector is given; `$scores`
  aggregates every built order.

- alpha:

  Numeric in `(0, 0.5)`. Significance threshold for anomaly
  classification (default 0.05). A path is labelled `"under"` when its
  adjusted lower-tail p-value `p_adjusted_under` is below `alpha`,
  `"over"` when its adjusted upper-tail p-value `p_adjusted_over` is
  below `alpha`, and `"normal"` otherwise (the under-representation test
  is applied first).

- min_count:

  Integer. Minimum count for a path to be classified as anomalous
  (default 5). An over-represented path must be observed at least
  `min_count` times; an under-represented path must be expected at least
  `min_count` times under the null, so a path that never occurs in the
  data can still be flagged as under-represented. Other paths are
  classified as `"normal"`.

- p_adjust:

  Character. Optional multiple testing correction of p-values. Default
  `"none"`: paths are classified against the threshold `alpha` directly,
  as in LaRock et al. (2020), which applies no correction. Accepts any
  method from
  [`p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html), e.g.
  `"BH"` for Benjamini-Hochberg FDR control. Under- and
  over-representation p-values are adjusted separately, within each
  order, over every possible path. With small counts the paper's
  upper-tail score `P(X > f)` rejects more often than `alpha` under the
  null, so neither the uncorrected nor the corrected classification
  carries an exact type-I or FDR guarantee; `min_count` limits the
  problem.

- k:

  Deprecated. Former name of `order`; if supplied it overrides `order`
  and emits a deprecation message. Use `order` instead.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: an
  object of class `net_hypa`.

- ...:

  In `print.net_hypa()` and `summary.net_hypa()`: Additional arguments
  (ignored).

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_hypa`.

- n:

  Integer. Maximum number of paths to display per category (default:
  10).

- type:

  Character. Which anomalies to show: `"all"` (default), `"over"`, or
  `"under"`.

- order_by:

  Character. Ranking used within each anomaly direction: `"sig"`
  (default) ranks by the active tail probability, `"freq"` (or its alias
  `"frequency"`) by observed count, `"ratio"` by observed/expected
  ratio, and `"path"` alphabetically.

## Value

An object of class `c("net_hypa", "cograph_network")` with components:

- scores:

  Data frame with path, from, to, observed, expected, ratio, p_value,
  p_under, p_over, p_adjusted_under, p_adjusted_over, anomaly, order
  columns (one block of rows per requested order). One row per possible
  path of the De Bruijn graph: every observed path plus every unobserved
  path that the lower-order transitions permit (`observed = 0`).
  `expected` comes from a propensity matrix fitted so that the expected
  counts preserve the observed in- and out-strength of every De Bruijn
  node (LaRock et al. 2020); at order 2 this equals the first-order
  Markov expectation \\n(a b \cdot) n(\cdot b c) / n(\cdot b \cdot)\\.
  The `path` column shows the full state sequence (e.g., "A -\> B -\>
  C"); `from` is the context (conditioning states); `to` is the next
  state; `ratio` is observed / expected; `p_value` is retained as an
  alias for `p_under`, the raw lower-tail hypergeometric CDF value;
  `p_over` is the strict upper-tail probability
  `P(X > observed) = 1 - p_under`, so that a path is over-represented
  when its HYPA score is at least `1 - alpha`, i.e. `p_over <= alpha`
  (LaRock et al., 2020); under-represented when `p_under < alpha`;
  `p_adjusted_under` and `p_adjusted_over` are the corrected p-values
  for under- and over-representation tests respectively.

- ho_edges:

  Alias for `scores` (all orders, arrow notation).

- over:

  Subset of `scores` classified as over-represented.

- under:

  Subset of `scores` classified as under-represented.

- adjacency:

  Weighted adjacency matrix of the lowest-order De Bruijn graph.

- weights:

  cograph weight matrix of the lowest-order graph.

- xi:

  Fitted propensity matrix of the lowest-order graph.

- edges:

  cograph edge data.frame of the lowest-order graph.

- by_order:

  Named list of per-order result lists.

- order:

  Integer vector of orders actually built (sorted ascending).

- k:

  Back-compatibility alias for `order`.

- alpha:

  Significance threshold used.

- p_adjust:

  Multiple testing correction method used.

- n_anomalous:

  Number of anomalous paths detected (all orders).

- n_over:

  Number of over-represented paths (all orders).

- n_under:

  Number of under-represented paths (all orders).

- n_edges:

  Total number of tested (possible) paths, all orders.

- nodes:

  data.frame (`id`, `label`, `name`) of the lowest-order De Bruijn graph
  nodes (arrow notation).

- directed:

  Logical. Always `TRUE`.

- meta:

  cograph meta list of the lowest-order graph.

- node_groups:

  Always `NULL`.

In `print.net_hypa()`: The input object, invisibly.

In `summary.net_hypa()`: A data frame of the reported anomalies, at most
`n` rows per direction, with columns `order` (the De Bruijn order the
path was found at), `path`, `observed`, `expected`, `ratio`, `p_tail`
(the raw tail probability in the reported direction: `p_over` for
over-represented paths, `p_under` for under-represented ones) and
`direction` (`"over"`/`"under"`). Returned invisibly; the summary text
and the top-`n` tables are printed as a side effect. When no anomalies
were detected, a zero-row data frame with the same columns except
`order` is returned.

## Details

Internal missing states split a trajectory, so no path spans a gap;
pieces shorter than two states are dropped. State labels must not
contain the control characters `"\x01"` or `"\x02"`, which are used
internally as separators (error class `nestimate_hypa_bad_label`). The
propensity matrix is fitted by continuous iterative proportional
fitting; when the fit does not reach its tolerance a warning of class
`nestimate_hypa_no_converge` is raised. The original Python
implementation rounds the propensities inside every sweep, so its
expected counts can differ from these by a few thousandths.

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
