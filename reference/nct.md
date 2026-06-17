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

- S:

  Same structure for the maximum absolute edge difference.

- E:

  Same structure for per-edge tests.

- n_iter:

  Number of permutations.

- paired:

  Whether a paired test was used.

## Details

Implementation matches
[`NetworkComparisonTest::NCT()`](https://rdrr.io/pkg/NetworkComparisonTest/man/NCT.html)
with defaults `abs = TRUE`, `weighted = TRUE`, `paired = FALSE` at
machine precision when the same seed is used. The network estimator is
EBIC-selected glasso applied to a Pearson correlation matrix, with
[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
symmetrization (matching NCT's `NCT_estimator_GGM` default).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
x1 <- matrix(rnorm(200 * 5), 200, 5)
x2 <- matrix(rnorm(200 * 5), 200, 5)
colnames(x1) <- colnames(x2) <- paste0("V", 1:5)
res <- nct(x1, x2, iter = 100)
res$M$p_value
res$S$p_value
} # }
```
