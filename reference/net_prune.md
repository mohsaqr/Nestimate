# Prune a Network's Edges

Removes weak or non-significant edges from a network, keeping a record
so the operation can be reversed. This is Nestimate's counterpart of
[`tna::prune()`](http://sonsoles.me/tna/reference/prune.md); the `net_`
prefix avoids a name clash with
[`tna::prune()`](http://sonsoles.me/tna/reference/prune.md).

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

  A `netobject` or `netobject_group`.

- method:

  One of `"threshold"`, `"lowest"`, `"disparity"`, `"bootstrap"`.
  Default `"threshold"`.

- threshold:

  Numeric cut-off for `method = "threshold"`. Default `0.1`.

- lowest:

  Quantile (0-1) for `method = "lowest"`. Default `0.05`.

- level:

  Significance level (0-1) for `method = "disparity"`. Default `0.5`.

- boot:

  Optional precomputed `net_bootstrap` for `method = "bootstrap"`.

- ...:

  Passed to
  [`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
  when `method = "bootstrap"`.

## Value

The input network (or group) with pruned `$weights` and a `"pruning"`
attribute. Class is unchanged.

## Details

Pruning is non-destructive: the pruned network carries a `"pruning"`
attribute holding the original weights, the pruned weights, the
parameters used, and a tidy table of removed edges. Use
[`net_deprune`](https://saqr.me/Nestimate/reference/net_deprune.md) to
restore the original weights and
[`net_reprune`](https://saqr.me/Nestimate/reference/net_reprune.md) to
re-apply the pruning, both without recomputation.
[`net_pruning_details`](https://saqr.me/Nestimate/reference/net_pruning_details.md)
reports what was removed.

Methods:

- `"threshold"`:

  Remove edges with weight \\\le\\ `threshold`.

- `"lowest"`:

  Remove the lowest `lowest` quantile of non-zero edges.

- `"disparity"`:

  Serrano disparity-filter backbone at significance `level`.

- `"bootstrap"`:

  Remove edges deemed non-significant by
  [`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
  (pass a precomputed result via `boot`, or extra bootstrap arguments
  via `...`).

For `"threshold"`, `"lowest"`, and `"disparity"` an edge is dropped only
when its removal leaves the network weakly connected.

Diagonal self-loops (self-transitions) are observed data: they are
counted equally when computing the cut-off but are never removed by any
method. (This is a deliberate divergence from
[`tna::prune()`](http://sonsoles.me/tna/reference/prune.md), which
prunes self-loops like any other edge.)

## See also

[`net_deprune`](https://saqr.me/Nestimate/reference/net_deprune.md),
[`net_reprune`](https://saqr.me/Nestimate/reference/net_reprune.md),
[`net_pruning_details`](https://saqr.me/Nestimate/reference/net_pruning_details.md)

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
