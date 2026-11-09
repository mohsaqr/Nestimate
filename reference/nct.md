# Network Comparison Test

`nct()` tests whether two networks estimated from two samples of the
same variables differ, at the three levels of van Borkulo et al. (2022):
**network structure** (statistic M, the largest absolute edge
difference), **global strength** (statistic S, the difference in the
summed absolute edge weights, or with `abs = FALSE` in the summed signed
weights, the global expected influence) and **individual edges**
(statistic E, the absolute difference of each edge). Each network is an
EBIC graphical lasso network of the Pearson correlations, and each
p-value comes from re-estimating both networks on permuted group labels
(within-pair swaps with `paired = TRUE`).

A p-value counts the permutations whose statistic is at least the
observed one, with the observed split added:
`(sum(perm >= observed) + 1) / (iter + 1)`, so the smallest attainable
value is `1 / (iter + 1)`. A non-significant test does not show that the
networks are equal; with small samples the test has little power (van
Borkulo et al. 2022). The edge tests are many tests at once; set
`p_adjust` before reading single edges. Both networks are regularised
estimates, so check how stable each one is with
[`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md).

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

  The observations of group 1: a numeric matrix or data frame, one row
  per observation. Non-numeric columns are dropped; at least two numeric
  variables and three rows are needed.

- data2:

  The observations of group 2, with the same variables as `data1` in the
  same order; columns that differ in name or order raise an error of
  class `nestimate_bad_input`. An unnamed matrix is matched by position
  and takes the other group's names, or `V1`, `V2`, ... when neither
  group has names.

- iter:

  The number of permutations, a whole number of at least 1. Default
  `1000`.

- gamma:

  The EBIC tuning parameter of the glasso, a non-negative number; larger
  values give sparser networks. Default `0.5`.

- paired:

  Whether row `i` of both groups is the same unit (such as a person
  measured twice). Each permutation then swaps the two rows of a unit at
  random instead of shuffling all rows, and the groups need the same
  number of rows, otherwise an error of class `nestimate_bad_input` is
  raised. Default `FALSE`.

- abs:

  Whether global strength sums the absolute edge weights. With `FALSE`,
  it sums the signed weights (global expected influence). Default
  `TRUE`.

- weighted:

  Whether the tests compare the weighted networks. With `FALSE`, both
  networks are binarised (an edge is present or absent) before the
  statistics are computed. Default `TRUE`.

- p_adjust:

  The multiplicity correction for the per-edge p-values, any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"none"`.

- na_method:

  How missing values are handled. `"pairwise"` (default) keeps every
  row, permutes rows as usual, and estimates each sample or permuted
  split from its pairwise-complete correlation matrix, with the median
  pairwise co-observation count as the sample size; a variable pair
  sharing fewer than two observed rows there raises an error of class
  `nestimate_degenerate_sample`. `"listwise"` drops incomplete rows
  first (with `paired = TRUE`, a pair is dropped when either member is
  incomplete). Complete data are unaffected by the choice.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods: a `net_nct` object.

- n:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: the
  number of edges to list, largest difference first. Default `Inf`
  (all).

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods: ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_nct` object.

- row.names, optional:

  For the [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  method: ignored, kept for compatibility with the base generic.

## Value

A `net_nct` object. [`print()`](https://rdrr.io/r/base/print.html)
reports the three tests, each group's global strength and the edges,
largest difference first, and returns the object invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the M and S
tests as a data frame with one row per test. The object's fields:

- nw1, nw2:

  The estimated weighted adjacency matrices of the two groups (binarised
  with `weighted = FALSE`).

- structure:

  The network-structure test (statistic M, the largest absolute edge
  difference): `observed`, `perm` (the statistic in each permutation)
  and `p_value`.

- strength:

  The global-strength test (statistic S): the same, plus `by_group`,
  each network's global strength (NCT's `glstrinv.sep`).

- E:

  The per-edge tests: `observed` and `p_value` hold one value per node
  pair (the upper triangle), `perm` is an `iter` by node-pair matrix,
  and `edge_names` is a data frame with the columns `from` and `to`
  naming each node pair.

- n_iter:

  The number of permutations.

- groups:

  The display names of the two groups: the argument names when they are
  plain variables, otherwise `"group 1"` and `"group 2"`.

- n:

  The rows used in each group.

- paired:

  Whether a paired test was used.

- params:

  The settings used: `gamma`, `abs`, `weighted`, `p_adjust` and
  `na_method`.

The letters follow van Borkulo et al. (2022) and
`NetworkComparisonTest`. Before Nestimate 0.9.27 the fields `M` and `S`
held the global-strength and network-structure tests the other way
round; reading them now raises an error of class `nestimate_renamed`.

For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
data frame with one row per global test and the columns `test`
(`"network structure"`, and `"global strength"`, or
`"global expected influence"` when `abs = FALSE`), `statistic` (`"M"`,
`"S"`), `observed`, `p_value` and `significant` (`p_value < 0.05`, a
fixed level). The
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
gives one row per node pair, with the columns `from`, `to`, `weight_1`
and `weight_2` (the edge in each group's network), `difference`
(statistic E), `p_value` and `significant` (`p_value < 0.05`).

## Details

`nct()` follows `NetworkComparisonTest::NCT()` with its defaults
`abs = TRUE`, `weighted = TRUE` and `paired = FALSE`. The network
estimator is EBIC-selected glasso
([`psychnets::ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.html))
applied to a Pearson correlation matrix made positive definite with
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) (NCT's
`NCT_estimator_GGM` default). Under the same seed both functions permute
the same rows, so the permutation p-values agree; the glasso solver is
not the Fortran one NCT wraps, so the test statistics agree to
independent-solver precision (of the order of `1e-4`) rather than bit
for bit. The Matrix package (in Suggests) must be installed.

The two groups must hold the same variables in the same order;
non-numeric columns, such as an id or a group label, are dropped first.
Infinite values are refused; use finite observations or `NA` for missing
data.

A variable without variation in a group, or in a permuted split (a rare
binary item can be constant in a split), enters that network as an
isolated node. The statistic is then defined on every split, so the
permutation test stays exact. A constant variable in an observed group
raises a warning of class `nestimate_constant_variable`.

Malformed `iter`, `gamma`, `paired`, `abs`, `weighted` or data raise an
error of class `nestimate_bad_input`.

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
