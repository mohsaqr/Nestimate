# Build an Ising Network

`build_ising()` estimates an Ising network for binary (0/1) variables:
each variable is regressed on all the others with L1-penalised logistic
regression, the penalty is chosen by EBIC, and the two estimates of each
edge are combined into one undirected weight (van Borkulo et al., 2014).
It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "ising"`, and needs the glmnet package. The settings go
in `params`: `gamma` (the EBIC hyperparameter, default `0.25`),
`nlambda` (default `100`) and `rule`, `"AND"` (default: an edge only
when both regressions select it) or `"OR"` (when either does).

An edge is a conditional association between two variables given all the
others; an absent edge means the penalty removed it. Check the stability
of the edges with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
build_ising(data, actor = NULL, ...)
```

## Arguments

- data:

  The binary variables: a data frame with one row per observation and
  one column per variable, holding only 0 and 1. Non-numeric columns and
  a column named `rid` are set aside; columns with non-syntactic names,
  all-missing and constant columns, and rows with a missing value are
  dropped, each with a message; fewer than 2 variables or 3 complete
  rows, or any value other than 0 and 1, raise an error.

- actor:

  An id column naming whose observation each row is; it is set aside,
  not modelled. Default `NULL`.

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  such as `params` (`gamma`, `nlambda`, `rule`), `threshold`, `scaling`
  or `group`. `method` cannot be given: it is always `"ising"` (error of
  class `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`), exactly as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it: print it for an overview, read the edges with
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md),
and draw it with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## References

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
constructing networks from binary data. *Scientific Reports*, 4, 5918.
[doi:10.1038/srep05918](https://doi.org/10.1038/srep05918)

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(1)
  bin_data <- data.frame(matrix(rbinom(200, 1, 0.5), ncol = 5))
  net <- build_ising(bin_data)
}
```
