# Build a Graphical Lasso Network (EBICglasso)

`build_glasso()` estimates a regularised partial correlation network
with the graphical lasso: an L1 penalty shrinks small partial
correlations to exactly zero, and the penalty is chosen along a path of
`nlambda` values by the extended Bayesian information criterion (EBIC;
Foygel & Drton, 2010). It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "glasso"`. The settings go in `params`: `gamma` (the EBIC
hyperparameter, default `0.5`; larger values give sparser networks),
`nlambda` (default `100`), `lambda.min.ratio` (default `0.01`),
`penalize.diagonal` (default `FALSE`), `refit` (default `FALSE`; `TRUE`
re-estimates the selected edges without the penalty) and `cor_method`.

An absent edge means the penalty set it to zero at the selected lambda,
not that the two variables are shown to be independent. Check the
stability of the edges and centralities with
[`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md).

## Usage

``` r
build_glasso(data, actor = NULL, level = NULL, ...)
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
  such as `params` (`gamma`, `nlambda`, `refit`, `n` and the others
  above), `threshold`, `scaling` or `group`. `method` cannot be given:
  it is always `"glasso"` (error of class `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`, a `netobject_ml` with
`level = "both"`), exactly as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it: print it for an overview, read the edges with
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md),
and draw it with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604–612.

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
data(srl_strategies)
net <- build_glasso(srl_strategies)
```
