# Build a Correlation Network

`build_cor()` estimates a correlation network: an undirected graph of
variables in which each edge is the correlation of two variables. It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "cor"`. The correlation is Pearson's by default;
`params = list(cor_method = "spearman")` (or `"kendall"`) passes another
method to [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

A correlation edge is a marginal association: it can come from a shared
link with a third variable. To separate direct from indirect links, use
[`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md)
or
[`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md).
Check the stability of the edges with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
build_cor(data, actor = NULL, level = NULL, ...)
```

## Arguments

- data:

  The variables: a data frame with one row per observation and one
  column per variable, or a square symmetric correlation or covariance
  matrix, which needs the sample size in `params` (as in
  `params = list(n = 250)`). Non-numeric columns are set aside; a frame
  of two or more character or factor columns and fewer than two numeric
  ones is read as sequences and turned into per-row state counts first.
  A column named `rid` is set aside; columns with non-syntactic names,
  all-missing and constant columns, and rows with a missing value are
  dropped, each with a message; fewer than 2 variables or 3 complete
  rows raise an error.

- actor:

  The id column naming whose observation each row is (a person). With
  `level`, the decomposition groups by it; without, it is set aside, not
  modelled. Default `NULL`.

- level:

  Multilevel decomposition, with `actor`: `"between"` estimates the
  network of person means, `"within"` the network of deviations from
  them, and `"both"` returns the two as a `netobject_ml`. Default `NULL`
  (no decomposition).

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  such as `params` (`cor_method`, `n`), `threshold`, `scaling` or
  `group`. `method` cannot be given: it is always `"cor"` (error of
  class `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`, a `netobject_ml` with
`level = "both"`), exactly as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it: print it for an overview, read the edges with
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md),
and draw it with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
data(srl_strategies)
net <- build_cor(srl_strategies)
```
