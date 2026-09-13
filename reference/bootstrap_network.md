# Bootstrap a Network Estimate

Non-parametric bootstrap for any network estimated by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
Works with all built-in methods (transition and association) as well as
custom registered estimators.

For transition methods (`"relative"`, `"frequency"`, `"co_occurrence"`),
uses a fast pre-computation strategy: per-sequence count matrices are
computed once, and each bootstrap iteration only resamples sequences via
`colSums` (C-level) plus lightweight post-processing. Data must be in
wide format for transition bootstrap; use
[`convert_sequence_format`](https://saqr.me/Nestimate/reference/convert_sequence_format.md)
to convert long-format data first.

For association methods (`"cor"`, `"pcor"`, `"glasso"`, and custom
estimators), the full estimator is called on resampled rows each
iteration.

If a transition network contains only one sequence, the function warns
that such a network is not recommended for bootstrap or other
confirmatory testing.

## Usage

``` r
bootstrap_network(
  x,
  iter = 1000L,
  ci_level = 0.05,
  inference = "stability",
  consistency_range = c(0.75, 1.25),
  edge_threshold = NULL,
  seed = NULL,
  boundary = c("inclusive", "strict"),
  ci_method = c("percentile", "basic")
)
```

## Arguments

- x:

  A `netobject` from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
  The data, method, params, scaling, threshold, and level are all
  extracted from this object. A `cograph_network` is coerced first; a
  `netobject_group` or `mcml` bootstraps every constituent network, and
  a `wtna_mixed` bootstraps both of its components (see **Value**).

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

- boundary:

  Character. Comparison rule when computing the consistency-range
  p-value. `"inclusive"` (default, tna-compatible) counts iterations
  that meet the bound (\\\le\\ / \\\ge\\); `"strict"` counts only
  iterations strictly outside (\\\<\\ / \\\>\\).

- ci_method:

  Character. Method for the edge-weight confidence intervals.
  `"percentile"` (default) uses the empirical bootstrap quantiles
  (Efron). `"basic"` reflects those quantiles around the observed
  weight, \\(2\hat{\theta} - q\_{1-\alpha/2}, 2\hat{\theta} -
  q\_{\alpha/2})\\ (Davison & Hinkley 1997, eq. 5.6), which corrects
  first-order bootstrap bias but can produce bounds outside the natural
  weight range near boundaries (e.g., below 0 for transition
  probabilities close to 0).

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

  Long-format data frame, one row per non-zero original edge (undirected
  networks keep one row per unordered pair), with columns `from`, `to`,
  `weight`, `mean`, `sd`, `p_value`, `sig`, `ci_lower`, `ci_upper`, plus
  `cr_lower` and `cr_upper` when `inference = "stability"`.

- model:

  Pruned `netobject` (non-significant edges zeroed).

- method, params, iter, ci_level, inference, ci_method:

  Bootstrap config.

- consistency_range, edge_threshold:

  Inference parameters.

A `netobject_group` or `mcml` input returns a `"net_bootstrap_group"`
(named list of `net_bootstrap` results); a `wtna_mixed` input returns a
`"wtna_boot_mixed"` with `$transition` and `$cooccurrence` results.

## See also

[`certainty`](https://saqr.me/Nestimate/reference/certainty.md) for the
closed-form Bayesian counterpart (same result layout, no resampling);
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`print.net_bootstrap`](https://saqr.me/Nestimate/reference/print.net_bootstrap.md),
[`summary.net_bootstrap`](https://saqr.me/Nestimate/reference/summary.net_bootstrap.md)

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
  method = "relative")
boot <- bootstrap_network(net, iter = 10)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
boot <- bootstrap_network(net, iter = 100)
print(boot)
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 100  |  Nodes : 4
#>   Edges      : 0 significant / 16 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
summary(boot)
#>    from to    weight      mean         sd   p_value   sig   ci_lower  ci_upper
#> 1     A  A 0.1304348 0.1274665 0.08074773 0.6831683 FALSE 0.00000000 0.3022826
#> 2     A  B 0.3043478 0.2997861 0.08278592 0.3663366 FALSE 0.14284091 0.4582168
#> 3     A  C 0.1304348 0.1418228 0.07482973 0.6732673 FALSE 0.00000000 0.3000000
#> 4     A  D 0.4347826 0.4309246 0.11390253 0.3267327 FALSE 0.24720280 0.6570000
#> 5     B  A 0.1724138 0.1737516 0.06474250 0.5346535 FALSE 0.07821429 0.3180000
#> 6     B  B 0.4137931 0.4205944 0.09336440 0.2970297 FALSE 0.23541667 0.5884848
#> 7     B  C 0.2068966 0.2050824 0.07346017 0.4752475 FALSE 0.06981818 0.3420833
#> 8     B  D 0.2068966 0.2005716 0.08216875 0.5049505 FALSE 0.05967023 0.3751667
#> 9     C  A 0.4000000 0.3857941 0.10757329 0.3861386 FALSE 0.17766798 0.6044118
#> 10    C  B 0.2000000 0.1884783 0.09391798 0.6534653 FALSE 0.04761905 0.3828755
#> 11    C  C 0.1500000 0.1571277 0.07338105 0.6039604 FALSE 0.00000000 0.2994885
#> 12    C  D 0.2500000 0.2686000 0.09868840 0.5544554 FALSE 0.10360963 0.4642308
#> 13    D  A 0.1666667 0.1651025 0.08394913 0.6336634 FALSE 0.02261905 0.3436275
#> 14    D  B 0.2222222 0.2362980 0.11192342 0.6930693 FALSE 0.00000000 0.4332589
#> 15    D  C 0.4444444 0.4507200 0.14203236 0.3861386 FALSE 0.18312325 0.7498039
#> 16    D  D 0.1666667 0.1478795 0.11004755 0.8118812 FALSE 0.00000000 0.3750000
#>      cr_lower  cr_upper
#> 1  0.09782609 0.1630435
#> 2  0.22826087 0.3804348
#> 3  0.09782609 0.1630435
#> 4  0.32608696 0.5434783
#> 5  0.12931034 0.2155172
#> 6  0.31034483 0.5172414
#> 7  0.15517241 0.2586207
#> 8  0.15517241 0.2586207
#> 9  0.30000000 0.5000000
#> 10 0.15000000 0.2500000
#> 11 0.11250000 0.1875000
#> 12 0.18750000 0.3125000
#> 13 0.12500000 0.2083333
#> 14 0.16666667 0.2777778
#> 15 0.33333333 0.5555556
#> 16 0.12500000 0.2083333
# }
```
