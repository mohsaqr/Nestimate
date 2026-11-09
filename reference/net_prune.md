# Prune a Network's Edges

`net_prune()` removes weak or non-significant edges from a network and
returns the pruned network with a record of what it removed, so the
pruning can be undone and re-applied. It is Nestimate's counterpart of
[`tna::prune()`](https://sonsoles.me/tna/reference/prune.html); the
`net_` prefix avoids a name clash with it. It takes a `netobject` or a
`netobject_group`, whose networks are pruned one by one with the same
settings, and offers four methods (see *Methods*).

Pruning is non-destructive: the network carries a `"pruning"` attribute
with the original weights, the pruned weights, the settings used and a
table of the removed edges.
[`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
restores the original weights and
[`net_reprune()`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md)
re-applies the pruning, both without recomputing anything;
[`net_pruning_details()`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
lists what was removed.

`"threshold"` and `"lowest"` remove small weights; they do not show that
those edges are absent, and the cut-off is the analyst's choice.
`"disparity"` and `"bootstrap"` keep an edge on a statistical criterion.
Before reading a pruned network, see which edges went with
[`net_pruning_details()`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
and whether the conclusions hold at another cut-off;
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
gives the uncertainty of each edge.

## Usage

``` r
net_prune(
  x,
  method = "threshold",
  threshold = 0.1,
  lowest = 0.05,
  level = 0.5,
  boot = NULL,
  ...
)

# S3 method for class 'netobject'
net_prune(
  x,
  method = "threshold",
  threshold = 0.1,
  lowest = 0.05,
  level = 0.5,
  boot = NULL,
  ...
)

# S3 method for class 'netobject_group'
net_prune(
  x,
  method = "threshold",
  threshold = 0.1,
  lowest = 0.05,
  level = 0.5,
  boot = NULL,
  ...
)

# Default S3 method
net_prune(
  x,
  method = "threshold",
  threshold = 0.1,
  lowest = 0.05,
  level = 0.5,
  boot = NULL,
  ...
)
```

## Arguments

- x:

  The network: a `netobject` or a `netobject_group`. Another class
  raises an error. A network whose pruning is active raises an error of
  class `nestimate_bad_input`; after
  [`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
  it can be pruned again with other settings, which replaces its pruning
  record.

- method:

  The pruning method: `"threshold"` (default), `"lowest"`, `"disparity"`
  or `"bootstrap"` (see *Methods*).

- threshold:

  For `method = "threshold"`: the largest weight that is removed.
  Default `0.1`. Ignored by the other methods, but must still be a
  single number.

- lowest:

  For `method = "lowest"`: the quantile, between 0 and 1, of the
  positive weights that sets the cut-off. Default `0.05`. Ignored by the
  other methods, but checked whatever the method: a value outside
  `[0, 1]` raises an error.

- level:

  For `method = "disparity"`: the significance level, between 0 and 1;
  an edge is kept when its p-value is below it. Default `0.5`. Ignored
  by the other methods, but checked whatever the method: a value outside
  `[0, 1]` raises an error.

- boot:

  For `method = "bootstrap"`: a `net_bootstrap` from
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  of this network. With the default `NULL`, the bootstrap is run here on
  `x`. For a `netobject_group`, pass the `net_bootstrap_group` of
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  on the group, matched to the networks by name, or `NULL` to bootstrap
  each network on its own. Another object, a bootstrap of another
  network (its weights differ), a single `net_bootstrap` for a group, or
  a group bootstrap missing one of the networks raises an error of class
  `nestimate_bad_input`.

- ...:

  For `method = "bootstrap"` with `boot = NULL`: arguments passed to
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
  such as `iter`. Ignored otherwise.

## Value

The input, of the same class, with the pruned weights, its edge table
and edge count updated to them, and the `"pruning"` attribute described
above.
[`net_pruning_details()`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
reads the record,
[`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
and
[`net_reprune()`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md)
switch between the original and the pruned weights. For a
`netobject_group`, each network is pruned and carries its own record.

## Methods

- `"threshold"`:

  Removes positive edges with weight at most `threshold`.

- `"lowest"`:

  Removes positive edges with weight at most the `lowest` quantile of
  the positive weights (self-loops included in the quantile; each pair
  counted once on an undirected network).

- `"disparity"`:

  Keeps the backbone found by the disparity filter of Serrano, Boguna &
  Vespignani (2009): an edge stays when, at its source's outgoing or its
  target's incoming edges, its share of the strength is significant at
  `level` against a uniform split of that strength; the filter is
  defined for non-negative weights.

- `"bootstrap"`:

  Removes the edges that
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  finds non-significant: pass a computed result as `boot`, or bootstrap
  arguments through `...`.

For `"threshold"` and `"lowest"`, the candidate edges are visited one at
a time, column by column and down each column, and each is removed only
when the network stays weakly connected without it, so a network that is
not weakly connected to begin with loses no edge to these two methods.
`"disparity"` and `"bootstrap"` apply no connectivity check. Negative
edges are never candidates for `"threshold"` and `"lowest"`, whatever
their size. On a signed network, `"disparity"` treats negative edges as
not significant and removes them at the default `level`; at a `level`
above 0.5, a negative edge leaving a state with no positive outgoing
edge (or entering one with no positive incoming edge) can stay.
`"bootstrap"` removes the non-significant edges of either sign.

On an undirected network an edge is a pair of cells: it is removed in
both directions or kept in both, so the pruned matrix stays symmetric,
and the record lists and counts it once.

Self-loops (the diagonal) are observed data: they count towards the
`"lowest"` quantile but are never removed by any method. This is a
deliberate difference from
[`tna::prune()`](https://sonsoles.me/tna/reference/prune.html), which
prunes self-loops like any other edge.

## References

Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the
multiscale backbone of complex weighted networks. *Proceedings of the
National Academy of Sciences*, 106(16), 6483–6488.
[doi:10.1073/pnas.0808904106](https://doi.org/10.1073/pnas.0808904106)

## See also

[`net_deprune`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md),
[`net_reprune`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md),
[`net_pruning_details`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)

## Examples

``` r
seqs <- data.frame(
  V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
  V3 = c("C","A","C","B","A"))
net    <- build_network(seqs, method = "relative")
pruned <- net_prune(net, method = "threshold", threshold = 0.2)
net_pruning_details(pruned)
#> Pruning details
#>   Method:  user-specified threshold (0.2)
#>   Removed: 0 edges
#>   Retained: 3 edges
#> 
#> (no edges removed)
```
