# Bootstrap a Network Estimate

Non-parametric bootstrap for any network estimated by
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).
Works with all built-in methods (transition and association) as well as
custom registered estimators.

For transition methods (`"relative"`, `"frequency"`, `"co_occurrence"`),
uses a fast pre-computation strategy: per-sequence count matrices are
computed once, and each bootstrap iteration only resamples sequences via
`colSums` (C-level) plus lightweight post-processing. Data must be in
wide format for transition bootstrap; use
[`convert_sequence_format`](https://mohsaqr.github.io/Nestimate/reference/convert_sequence_format.md)
to convert long-format data first.

For association methods (`"cor"`, `"pcor"`, `"glasso"`, and custom
estimators), the full estimator is called on resampled rows each
iteration.

## Usage

``` r
bootstrap_network(
  x,
  iter = 1000L,
  ci_level = 0.05,
  inference = "stability",
  consistency_range = c(0.75, 1.25),
  edge_threshold = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A `netobject` from
  [`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).
  The data, method, params, scaling, threshold, and level are all
  extracted from this object.

- iter:

  Integer. Number of bootstrap iterations (default: 1000).

- ci_level:

  Numeric. Significance level for CIs and p-values (default: 0.05).

- inference:

  Character. `"stability"` (default) tests whether bootstrap replicates
  fall within a multiplicative consistency range around the original
  weight. `"threshold"` tests whether replicates exceed a fixed edge
  threshold.

- consistency_range:

  Numeric vector of length 2. Multiplicative bounds for stability
  inference (default: `c(0.75, 1.25)`).

- edge_threshold:

  Numeric or NULL. Fixed threshold for `inference = "threshold"`. If
  NULL, defaults to the 10th percentile of absolute original edge
  weights.

- seed:

  Integer or NULL. RNG seed for reproducibility.

## Value

An object of class `"net_bootstrap"` containing:

- original:

  The original `netobject`.

- mean:

  Bootstrap mean weight matrix.

- sd:

  Bootstrap SD matrix.

- p_values:

  P-value matrix.

- significant:

  Original weights where p \< ci_level, else 0.

- ci_lower:

  Lower CI bound matrix.

- ci_upper:

  Upper CI bound matrix.

- cr_lower:

  Consistency range lower bound (stability only).

- cr_upper:

  Consistency range upper bound (stability only).

- summary:

  Long-format data frame of edge-level statistics.

- model:

  Pruned `netobject` (non-significant edges zeroed).

- method, params, iter, ci_level, inference:

  Bootstrap config.

- consistency_range, edge_threshold:

  Inference parameters.

## See also

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
[`print.net_bootstrap`](https://mohsaqr.github.io/Nestimate/reference/print.net_bootstrap.md),
[`summary.net_bootstrap`](https://mohsaqr.github.io/Nestimate/reference/summary.net_bootstrap.md)

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
boot <- bootstrap_network(net, iter = 100)
print(boot)
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 100  |  Nodes : 4
#>   Edges      : 0 significant / 12 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
summary(boot)
#>    from to     weight       mean         sd   p_value   sig   ci_lower
#> 1     A  A 0.27272727 0.24366540 0.10695384 0.6039604 FALSE 0.05025253
#> 2     A  B 0.22727273 0.22934223 0.09023840 0.5544554 FALSE 0.06875000
#> 3     A  C 0.13636364 0.15377496 0.08801395 0.6039604 FALSE 0.00000000
#> 4     A  D 0.36363636 0.37321741 0.08903158 0.2772277 FALSE 0.20937500
#> 5     B  A 0.13043478 0.13950739 0.06991792 0.5643564 FALSE 0.00000000
#> 6     B  B 0.26086957 0.25706101 0.07198383 0.3168317 FALSE 0.08693182
#> 7     B  C 0.21739130 0.21055293 0.08965278 0.5940594 FALSE 0.04441700
#> 8     B  D 0.39130435 0.39287866 0.10391133 0.3861386 FALSE 0.20395833
#> 9     C  A 0.19047619 0.17665609 0.07578627 0.6039604 FALSE 0.04761905
#> 10    C  B 0.19047619 0.19566986 0.09939174 0.6435644 FALSE 0.03919231
#> 11    C  C 0.28571429 0.28402163 0.08666471 0.4059406 FALSE 0.09375000
#> 12    C  D 0.33333333 0.34365242 0.11159900 0.4257426 FALSE 0.13944805
#> 13    D  A 0.33333333 0.33808647 0.09367446 0.3564356 FALSE 0.16316667
#> 14    D  B 0.37500000 0.36456841 0.08298671 0.2871287 FALSE 0.20826087
#> 15    D  C 0.20833333 0.21139596 0.08015006 0.5049505 FALSE 0.07612045
#> 16    D  D 0.08333333 0.08594916 0.05569708 0.6930693 FALSE 0.00000000
#>     ci_upper   cr_lower  cr_upper
#> 1  0.4372115 0.20454545 0.3409091
#> 2  0.4104947 0.17045455 0.2840909
#> 3  0.3500000 0.10227273 0.1704545
#> 4  0.5478409 0.27272727 0.4545455
#> 5  0.2766667 0.09782609 0.1630435
#> 6  0.3893058 0.19565217 0.3260870
#> 7  0.3515441 0.16304348 0.2717391
#> 8  0.5956818 0.29347826 0.4891304
#> 9  0.3012443 0.14285714 0.2380952
#> 10 0.3800481 0.14285714 0.2380952
#> 11 0.4342949 0.21428571 0.3571429
#> 12 0.5812500 0.25000000 0.4166667
#> 13 0.5000000 0.25000000 0.4166667
#> 14 0.5201923 0.28125000 0.4687500
#> 15 0.3552500 0.15625000 0.2604167
#> 16 0.1864005 0.06250000 0.1041667
# }
```
