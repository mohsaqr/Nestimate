# Estimate a Network (Deprecated)

`estimate_network()` is deprecated: use
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
instead. It passes its arguments to
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
unchanged and returns what that returns, after a warning of class
`deprecatedWarning`.

## Usage

``` r
estimate_network(
  data,
  method = "relative",
  params = list(),
  scaling = NULL,
  threshold = 0,
  level = NULL,
  ...
)
```

## Arguments

- data:

  The data to estimate from: a data frame of sequences (one row per
  sequence, one column per position), an event log (one row per event,
  with `action` set), one-hot code columns, numeric variables (one row
  per observation) for the association methods, or a square symmetric
  correlation or covariance matrix. A fitted `net_clustering` or
  `net_mmm` object also works: its per-cluster networks are built, with
  the `scaling`, `threshold`, `start`, `end` and `params` given here,
  and returned as a `netobject_group`. The relative networks of a
  `net_mmm` are the fitted chains: `scaling` and `threshold` are applied
  to them, and any other argument raises an error of class
  `nestimate_bad_input`.

- method:

  The estimator, as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  Default `"relative"`, kept for backward compatibility;
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  has no default.

- params:

  Settings passed to the estimator, as a named list, such as
  `list(gamma = 0.5)` for `glasso` or `list(n_gram = 3)` for `ngram`.
  The list is stored on the network, so
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  and the other resampling verbs re-estimate with exactly the same
  settings. The transition methods take tna's sequence options
  `weighted` and `concat`, and the low-level `begin_state` / `end_state`
  that `start` and `end` set. `window_size` (with `window_span`, or
  `windowed = TRUE` for a size of 2) counts `relative`, `frequency` and
  `co_occurrence` between consecutive non-overlapping blocks of that
  many positions, as
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
  does; other transition methods raise an error of class
  `nestimate_bad_input`, and `mode = "overlapping"` a warning of class
  `nestimate_ignored_mode`. Column roles given here (`action`, `id`,
  `id_col`, `actor`, `time`, `session`, `order`, `cols`, `codes`,
  `group`) must name existing columns, and must agree with the same role
  given directly.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` maps the non-zero weights into (0, 1\], `"max"` divides by
  the largest absolute weight, `"rank"` replaces the non-zero weights by
  their ranks, and `"normalize"` makes each row's absolute weights sum
  to 1. Zero weights stay zero under every option. Combine them as in
  `c("rank", "minmax")`. Default `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero. Default
  `0` (none). For `method = "mgm"`, `"LW"` or `"none"` instead selects
  mgm's coefficient-thresholding rule (a `threshold` in `params` takes
  precedence); any other non-numeric value raises an error.

- level:

  Multilevel decomposition for the undirected association methods
  (`cor`, `pcor`, `glasso`, `mgm`): `"between"` estimates the network of
  person means, `"within"` the network of deviations from them, and
  `"both"` returns the two as a `netobject_ml`. It needs an id column,
  given as `actor` or as `params$id`; a directed method raises an error.
  Default `NULL` (no decomposition). For `method = "mgm"`, a numeric
  vector instead gives one category count per modelled variable (`1` for
  continuous ones), after group, actor and metadata columns are set
  aside; it is checked and kept for compatibility, and the fit reads the
  categories from the data.

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  such as `actor`, `action`, `time` or `group`.

## Value

A `netobject`, exactly as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it (a `netobject_group` or `netobject_ml` when `...` or `level`
asks for one).

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
data <- data.frame(T1 = c("x","y","z","x"), T2 = c("y","x","z","y"))
net <- estimate_network(data, method = "relative")
#> Warning: 'estimate_network' is deprecated.
#> Use 'build_network' instead.
#> See help("Deprecated")
```
