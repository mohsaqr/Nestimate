# Network Comparison Test

Tests whether two networks estimated from independent samples differ at
three levels: **global strength** (M-statistic), **network structure**
(S-statistic, max absolute edge difference), and **individual edges**
(E-statistic per edge). Inference is via permutation of group labels.

## Usage

``` r
nct(
  data1,
  data2,
  iter = 1000L,
  gamma = 0.5,
  paired = FALSE,
  abs = TRUE,
  weighted = TRUE,
  p_adjust = "none"
)
```

## Arguments

- data1:

  A numeric matrix or data.frame of observations from group 1.

- data2:

  A numeric matrix or data.frame of observations from group 2. Same
  number of columns as `data1`.

- iter:

  Integer. Number of permutation iterations. Default 1000.

- gamma:

  EBIC tuning parameter for glasso. Default 0.5.

- paired:

  Logical. If `TRUE`, perform a paired permutation (within-subject
  swap). Default `FALSE`.

- abs:

  Logical. If `TRUE`, compute global strength on absolute edge weights.
  Default `TRUE`.

- weighted:

  Logical. If `TRUE`, use weighted networks for the tests. If `FALSE`,
  binarize before computing statistics. Default `TRUE`.

- p_adjust:

  P-value adjustment method for the per-edge tests (any method in
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html)).
  Default `"none"`.

## Value

A list of class `net_nct` with elements:

- nw1, nw2:

  Estimated weighted adjacency matrices.

- M:

  List with `observed`, `perm`, `p_value` for the global strength test.
  P-values are permutation p-values,
  `(sum(perm >= observed) + 1) / (iter + 1)`.

- S:

  Same structure for the maximum absolute edge difference.

- E:

  Same structure for the per-edge tests (`observed` and `p_value` are
  one value per upper-triangle edge, `perm` an `iter` by edges matrix),
  plus `edge_names`, a two-column data frame of the node pairs (`NULL`
  when `data1` has no column names).

- n_iter:

  Number of permutations.

- paired:

  Whether a paired test was used.

- params:

  List of the settings used: `gamma`, `abs`, `weighted`, `p_adjust`.

## Details

Follows `NetworkComparisonTest::NCT()` with defaults `abs = TRUE`,
`weighted = TRUE`, `paired = FALSE`. The network estimator is
EBIC-selected glasso applied to a Pearson correlation matrix, with
[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
symmetrization (matching NCT's `NCT_estimator_GGM` default). The glasso
solver is not the Fortran one NCT wraps, so results agree to
independent-solver precision (of the order of `1e-4` on the test
statistics) rather than bit-for-bit, even under the same seed.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
x1 <- matrix(rnorm(200 * 5), 200, 5)
x2 <- matrix(rnorm(200 * 5), 200, 5)
colnames(x1) <- colnames(x2) <- paste0("V", 1:5)
res <- nct(x1, x2, iter = 100)
res
summary(res)
} # }
```
