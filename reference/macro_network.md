# Cluster-Level Network, With One Cluster Expanded

`macro_network()` returns the cluster-level (macro) network of an
`mcml`, optionally with one or more clusters expanded back into their
member states. Every other cluster stays collapsed to a single node, so
the result is a network at mixed resolution: the cluster of interest in
detail, its context in summary, with weights counted from the sequences.
For an `mcml_pc` from
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md),
it returns the cluster-level network that fit already holds.

## Usage

``` r
macro_network(
  x,
  expand = NULL,
  method = "relative",
  scaling = NULL,
  threshold = 0,
  params = list(),
  ...
)
```

## Arguments

- x:

  An `mcml` built from sequence data, or an `mcml_pc` from
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md).
  An `mcml` built from a matrix or an edge list holds no sequences to
  count, so it raises an error of class `nestimate_no_expand_source`,
  even with `expand = NULL`; a name in `expand` that is not a cluster is
  reported first.

- expand:

  The clusters to expand into their member states. `NULL` (default)
  collapses every cluster, which reproduces the macro layer; `"all"` or
  `TRUE` expands every cluster, so each state is its own node. A name
  that is not a cluster raises an error. For an `mcml_pc`, any value
  other than `NULL` raises an error of class `nestimate_no_expand`:
  there are no sequences to count again.

- method:

  The estimator passed to
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  Default `"relative"` (row-normalised transitions). For an `mcml_pc`,
  giving it raises an error: the estimator was set by
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md).

- scaling, threshold, params:

  Passed to
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  for the re-counted network: `scaling` (default `NULL`, none),
  `threshold` (default `0`) and `params` (default
  [`list()`](https://rdrr.io/r/base/list.html)), as there. For an
  `mcml_pc`, giving any raises an error, as for `method`.

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  such as `start` or `end`. For an `mcml_pc`, giving any raises an
  error, as for `method`.

## Value

For an `mcml`: a `netobject` (also a `cograph_network`) whose nodes are
the collapsed clusters plus the member states of any expanded cluster,
with weights counted by
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
It is printed, resampled (for example with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md))
and drawn like any other `netobject`; its node grouping is read by
cograph's group-aware plots, such as
[`cograph::plot_mtna()`](https://sonsoles.me/cograph/reference/plot_mtna.html).
Beside the usual `netobject` fields it holds `node_groups`, a data frame
with the columns `node` and `group` that maps each node to its cluster
(an expanded cluster's states map to that cluster, a collapsed cluster
maps to itself), the same labels as a factor in the `groups` column of
`nodes`, and `expanded`, the names of the clusters expanded (`NULL` when
none were).

For an `mcml_pc`: its cluster-level network, the `netobject` estimated
by
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md),
unchanged.

## Details

The weights are counted again from the sequences, after each state is
recoded to its cluster (or kept, for an expanded cluster), because a
cluster-level matrix cannot be split back into its states.

## See also

[`build_mcml`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md),
[`as_tna`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)

## Examples

``` r
seqs <- data.frame(
  t1 = c("A", "C", "A", "B"), t2 = c("B", "D", "C", "A"),
  t3 = c("C", "A", "D", "C"), stringsAsFactors = FALSE
)
mc <- build_mcml(seqs, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))
macro_network(mc)                    # every cluster collapsed
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.333, 0.667]  |  mean: 0.500
#> 
#>   Weight matrix:
#>         G1    G2
#>   G1 0.400 0.600
#>   G2 0.333 0.667 
#> 
#>   Initial probabilities:
#>   G1            0.750  ████████████████████████████████████████
#>   G2            0.250  █████████████
macro_network(mc, expand = "G2")     # G2 shown as C and D
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.400, 1.000]  |  mean: 0.750
#> 
#>   Weight matrix:
#>        C D  G1
#>   C  0.0 1 0.0
#>   D  0.0 0 1.0
#>   G1 0.6 0 0.4 
#> 
#>   Initial probabilities:
#>   G1            0.750  ████████████████████████████████████████
#>   C             0.250  █████████████
#>   D             0.000  
```
