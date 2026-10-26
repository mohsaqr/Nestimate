# Network Comparison Test

Tests whether two networks estimated from independent samples differ at
three levels (van Borkulo et al. 2022): **network structure**
(M-statistic, the largest absolute edge difference), **global strength**
(S-statistic, the difference in summed absolute edge weights; with
`abs = FALSE`, summed signed weights, i.e. global expected influence),
and **individual edges** (E-statistic per edge). Inference is via
permutation of group labels.

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
  p_adjust = "none",
  na_method = c("pairwise", "listwise")
)

# S3 method for class 'net_nct'
print(x, n = Inf, ...)

# S3 method for class 'net_nct'
summary(object, ...)

# S3 method for class 'net_nct'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- data1:

  A numeric matrix or data.frame of observations from group 1.

- data2:

  A numeric matrix or data.frame of observations from group 2, with the
  same variables as `data1`.

- iter:

  Integer. Number of permutation iterations. Default 1000.

- gamma:

  EBIC tuning parameter for glasso. Default 0.5.

- paired:

  Logical. If `TRUE`, perform a paired permutation (within-subject swap;
  row `i` of both groups is the same unit). Default `FALSE`.

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

- na_method:

  `"pairwise"` (default) or `"listwise"` handling of missing values; see
  Details.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods: an object of class `net_nct`.

- n:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: number
  of edges to list, largest difference first. Default `Inf` (all).

- ...:

  In `print.net_nct()`, `summary.net_nct()` and
  `as.data.frame.net_nct()`: Ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_nct`.

- row.names, optional:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
  ignored (base generic compatibility).

## Value

A list of class `net_nct` with elements:

- nw1, nw2:

  Estimated weighted adjacency matrices.

- structure:

  The network-structure test (statistic M, the largest absolute edge
  difference): list with `observed`, `perm`, `p_value`. P-values are
  permutation p-values, `(sum(perm >= observed) + 1) / (iter + 1)`.

- strength:

  The global-strength test (statistic S): same structure, plus
  `by_group`, each network's global strength (NCT's `glstrinv.sep`).

- E:

  Same structure for the per-edge tests (`observed` and `p_value` are
  one value per upper-triangle edge, `perm` an `iter` by edges matrix),
  plus `edge_names`, a data frame with columns `from` and `to` naming
  each node pair.

- n_iter:

  Number of permutations.

- groups:

  Display names of the two groups: the argument names when they are
  plain variables, otherwise "group 1" / "group 2".

- n:

  Rows used in each group.

- paired:

  Whether a paired test was used.

- params:

  List of the settings used: `gamma`, `abs`, `weighted`, `p_adjust`,
  `na_method`.

The letters follow van Borkulo et al. (2022) and
`NetworkComparisonTest`. Before Nestimate 0.9.27 the fields `$M` and
`$S` held the global-strength and network-structure tests the other way
round; accessing them now raises a `nestimate_renamed` error. Errors of
class `nestimate_bad_input` for malformed arguments or groups whose
columns differ (in name or order); of class
`nestimate_degenerate_sample` when (pairwise) a variable pair shares
fewer than two observed rows in a group or a permuted split. A warning
of class `nestimate_constant_variable` when a variable is constant in an
observed group.

In `print.net_nct()`: The input object, invisibly.

In `summary.net_nct()`: a data frame with one row per global test,
columns `test` (`"network structure"`; `"global strength"`, or
`"global expected influence"` when `abs = FALSE`), `statistic` (`"M"`,
`"S"`), `observed`, `p_value`, `significant` (`p_value < 0.05`).

In `as.data.frame.net_nct()`: a data frame with one row per node pair,
columns `from`, `to`, `weight_1`, `weight_2` (the edge in each group's
network), `difference` (statistic E), `p_value`, `significant`.

## Details

Follows `NetworkComparisonTest::NCT()` with defaults `abs = TRUE`,
`weighted = TRUE`, `paired = FALSE`. The network estimator is
EBIC-selected glasso
([`psychnets::ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.html))
applied to a Pearson correlation matrix, with
[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
symmetrization (matching NCT's `NCT_estimator_GGM` default). Under the
same seed both functions permute the same rows, so permutation p-values
agree; the glasso solver is not the Fortran one NCT wraps, so the test
statistics agree to independent-solver precision (of the order of
`1e-4`) rather than bit-for-bit.

The two groups must hold the same variables in the same order
(non-numeric columns such as an id or a group label are dropped first).

Missing data follow `na_method`. `"pairwise"` (default) keeps every row,
permutes rows as usual, and estimates each (sub)sample from its
pairwise-complete correlation matrix, with the median pairwise
co-observation count as the sample size. `"listwise"` drops incomplete
rows first (in a paired design, a pair is dropped when either member is
incomplete). Complete data are unaffected by the choice.

A variable without variation in a group, or in a permuted split (a rare
binary item can be constant in a split), enters that network as an
isolated node. The statistic is then defined on every split, so the
permutation test stays exact. A constant variable in an observed group
raises a `nestimate_constant_variable` warning.

## Methods

- Printing the object shows every result. `summary.net_nct()` returns
  the M (network structure) and S (global strength) tests;
  `as.data.frame.net_nct()` returns one row per edge.

## References

van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022). Comparing
network structures on three aspects: A permutation test. *Psychological
Methods*. [doi:10.1037/met0000476](https://doi.org/10.1037/met0000476)

## Examples

``` r
set.seed(1)
x1 <- matrix(rnorm(100 * 4), 100, 4)
x2 <- matrix(rnorm(100 * 4), 100, 4)
colnames(x1) <- colnames(x2) <- paste0("V", 1:4)
# iter = 20 keeps the example fast; a real analysis uses 1000 or more.
res <- nct(x1, x2, iter = 20)
res
#> Network Comparison Test: x1 vs x2
#> 20 permutations, independent groups, n = 100 and 100
#> 
#>   Test                  Statistic       p
#>   M  network structure      0.000   1.000 
#>   S  global strength        0.000   1.000 
#>        x1 0.000 vs x2 0.000
#>   E  edges differing       0 of 6
#> 
#>   Edges, largest difference first
#>   Edge         x1      x2  |diff|       p
#>   V1 - V2   0.000   0.000   0.000   1.000 
#>   V1 - V3   0.000   0.000   0.000   1.000 
#>   V1 - V4   0.000   0.000   0.000   1.000 
#>   V2 - V3   0.000   0.000   0.000   1.000 
#>   V2 - V4   0.000   0.000   0.000   1.000 
#>   V3 - V4   0.000   0.000   0.000   1.000 
#> 
#>   * p < 0.05
summary(res)
#>                test statistic observed p_value significant
#> 1 network structure         M        0       1       FALSE
#> 2   global strength         S        0       1       FALSE
```
