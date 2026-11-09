# Compute Node Predictability

`predictability()` computes how well each node of an association network
is predicted by the others: the proportion of its variance explained
(R\\^2\\), with the root mean squared error (RMSE) of that prediction,
following Haslbeck & Waldorp (2018). It reads networks estimated with
`method = "glasso"`, `"pcor"` or `"cor"`, single, multilevel
(`netobject_ml`) or grouped (`netobject_group`).
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
stores the same R\\^2\\ values on these networks when
`predictability = TRUE` (the default), and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws them as a ring around each node.

R\\^2\\ says how much of a node's variance its connections in the
network account for; a low value leaves most of it to variables outside
the network or to noise. It is computed from the matrices estimated on
the sample, so it describes that sample. For `"glasso"` it reflects the
regularised estimate, and for `"cor"` it depends on which edges survive
`threshold`.

## Usage

``` r
predictability(object, ...)

# S3 method for class 'netobject'
predictability(object, data = NULL, ...)

# S3 method for class 'netobject_ml'
predictability(object, ...)

# S3 method for class 'netobject_group'
predictability(object, ...)
```

## Arguments

- object:

  The network: a `netobject` estimated with `method = "glasso"`,
  `"pcor"` or `"cor"`, or a `netobject_ml` or `netobject_group` of them.
  A network from any other method raises an error.

- ...:

  Not used.

- data:

  The variables the network was estimated from, used only for the `RMSE`
  column; R\\^2\\ comes from the matrices stored on the network. Its
  columns are matched by node label or, failing that, by node name.
  Default `NULL`, which uses the data stored on the network; without
  either, or when the columns do not match, `RMSE` is `NA`. The
  `netobject_ml` and `netobject_group` methods take no `data` (one
  passed to them is ignored) and use each network's own.

## Value

A data frame with one row per node and the columns `node` (the node
label), `R2` (between 0 and 1) and `RMSE` (`NA` when no data is
available). A `netobject_ml` gives a list of two such data frames,
`between` and `within`, and a `netobject_group` a named list of them,
one per group.

## Details

For `"glasso"` and `"pcor"`, R\\^2\\ comes from the estimated precision
matrix \\\Omega\\ (the inverse of the correlation matrix, regularised
for `"glasso"`): \$\$R^2_j = 1 - 1 / \Omega\_{jj}\$\$ For `"cor"`, it is
the multiple R\\^2\\ of each node regressed on its neighbours (the nodes
it shares a non-zero edge with), computed from the stored correlation
matrix: 0 for a node without neighbours, the squared correlation for a
node with one. Values are clipped to \[0, 1\].

RMSE predicts each standardised variable from the others (for `"cor"`,
from its neighbours) with the regression weights implied by the same
matrix, on the complete rows of the data, and reports the error on the
variable's original scale. For `"cor"`, a node without neighbours has an
RMSE equal to its standard deviation (prediction by the mean), and a
singular correlation block of its neighbours gives an R\\^2\\ of 0 and
an RMSE of `NA`.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
predict observations? On the importance of predictability in network
models. *Behavior Research Methods*, 50(2), 853–861.
[doi:10.3758/s13428-017-0910-x](https://doi.org/10.3758/s13428-017-0910-x)

## Examples

``` r
set.seed(42)
mat <- matrix(rnorm(60), ncol = 4)
colnames(mat) <- LETTERS[1:4]
net <- build_network(as.data.frame(mat), method = "glasso")
predictability(net)
#>   node R2      RMSE
#> 1    A  0 0.9904791
#> 2    B  0 1.3109578
#> 3    C  0 0.9630332
#> 4    D  0 1.0523124
```
