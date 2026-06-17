# Compare Network-Level Statistics of Two Networks

Snijders & Borgatti (1999) two-network test: each network's statistics
get vertex-bootstrap standard errors, and each difference is tested with
\$\$z = (\hat{\theta}\_x - \hat{\theta}\_y) / \sqrt{SE_x^2 + SE_y^2}\$\$
against a standard normal reference. This is the comparison the vertex
bootstrap was originally proposed for: deciding whether two observed
networks differ in density, centralization, reciprocity, or any other
whole-network descriptive.

## Usage

``` r
vertex_compare(
  x,
  y,
  iter = 1000L,
  ci_level = 0.05,
  statistics = NULL,
  statistic_fn = NULL,
  directed = NULL,
  seed = NULL,
  labels = c("x", "y")
)
```

## Arguments

- x, y:

  The two networks: `netobject`s, `cograph_network`s, square weight
  matrices, or precomputed `net_vertex_bootstrap` objects (then `iter`,
  `statistics`, `statistic_fn`, `directed`, and `seed` are ignored for
  that argument). The two sides must cover exactly the same statistics;
  a mismatch (e.g., a directed network's `reciprocity` against an
  undirected one, or precomputed objects built with different
  `statistics` selections) is an error, never a silent subset.

- iter:

  Integer. Number of bootstrap replicates (default 1000).

- ci_level:

  Numeric. Significance level for the confidence intervals (default 0.05
  for 95% CIs).

- statistics:

  Character vector selecting built-in statistics (see Details). Default:
  all applicable to the network's directedness.

- statistic_fn:

  Optional named list of functions, each taking the weight matrix and
  returning a single numeric value. Computed alongside the built-ins.

- directed:

  Logical or NULL. Directedness of the network. NULL (default) reads
  `x$directed` when available, otherwise falls back to matrix symmetry.

- seed:

  Integer or NULL. RNG seed for reproducibility.

- labels:

  Character vector of length 2 naming the networks in the output
  (default `c("x", "y")`).

## Value

An object of class `"net_vertex_comparison"` containing:

- summary:

  Tidy data frame, one row per statistic: `statistic`, the two observed
  values, `diff`, `se_diff`, `z`, `p_value`, and a normal-approximation
  confidence interval for the difference.

- x, y:

  The two `net_vertex_bootstrap` results.

- labels, ci_level:

  Configuration.

When both bootstrap SEs are zero (a statistic with no resampling
variation in either network) `z` and `p_value` are `NA`.

## References

Snijders, T. A. B., & Borgatti, S. P. (1999). Non-parametric standard
errors and tests for network statistics. *Connections*, 22(2), 161-170.

## See also

[`vertex_bootstrap`](https://saqr.me/Nestimate/reference/vertex_bootstrap.md),
[`nct`](https://saqr.me/Nestimate/reference/nct.md) for the
permutation-based comparison of edge-level structure when raw data are
available, `permutation_test`.

## Examples

``` r
states <- c("plan", "code", "debug", "test")
s1 <- data.frame(
  T1 = rep(states, 5), T2 = rep(rev(states), 5),
  T3 = rep(states[c(2, 3, 4, 1)], 5)
)
s2 <- data.frame(
  T1 = rep(states[c(3, 1, 4, 2)], 5), T2 = rep(states, 5),
  T3 = rep(states[c(4, 3, 1, 2)], 5)
)
net1 <- build_network(s1, method = "relative")
net2 <- build_network(s2, method = "relative")
cmp <- vertex_compare(net1, net2, iter = 100, seed = 1)
cmp$summary
#>        statistic observed_x observed_y        diff    se_diff          z
#> 1        density  0.5000000  0.5833333 -0.08333333 0.25478584 -0.3270721
#> 2    mean_weight  0.5000000  0.5714286 -0.07142857 0.09622256 -0.7423267
#> 3 centralization  0.3333333  0.0000000  0.33333333 0.24140627  1.3807981
#> 4    reciprocity  1.0000000  0.2500000  0.75000000 0.35488859  2.1133393
#>      p_value    ci_lower  ci_upper
#> 1 0.74361337 -0.58270441 0.4160377
#> 2 0.45788942 -0.26002131 0.1171642
#> 3 0.16734104 -0.13981427 0.8064809
#> 4 0.03457174  0.05443114 1.4455689
```
