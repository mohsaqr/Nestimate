# Build a Partial Correlation Network

`build_pcor()` estimates a partial correlation network: an undirected
graph of variables in which each edge is the correlation of two
variables after the others are held constant, read from the inverse of
the correlation matrix. It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "pcor"`, with no regularisation.

Unregularised partial correlations are noisy when there are many
variables for the number of observations, and they cannot be estimated
when the correlation matrix is singular (more variables than
observations, or collinear ones): that raises an error, and a
near-singular matrix a warning.
[`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)
estimates a sparse network in those cases. Check the stability of the
edges with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
build_pcor(data, actor = NULL, level = NULL, ...)
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
  `group`. `method` cannot be given: it is always `"pcor"` (error of
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
net <- build_pcor(srl_strategies)
```
