# Composite-Weight Stability Under Case Resampling

**Experimental.** Bootstraps the item weights of a
[`build_mcml_pc`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
fit: rows of the raw data are resampled, the node-level network is
re-estimated each time, and the connectivity-based composite weights are
recomputed. Wide intervals mean the weighting (and therefore the
`"loadings"` macro network) should not be over-interpreted.

## Usage

``` r
loading_stability(x, iter = 200L, ci_level = 0.05, seed = NULL)
```

## Arguments

- x:

  An `mcml_pc` object that carries raw data.

- iter:

  Integer. Bootstrap replicates (default 200; node-level re-estimation
  makes this heavier than a plain bootstrap).

- ci_level:

  Numeric. Significance level for percentile CIs (default 0.05).

- seed:

  Integer or NULL. RNG seed.

## Value

An object of class `"pc_loading_stability"`: a list with `summary` (tidy
data frame: `node`, `cluster`, `weight`, `boot_mean`, `boot_sd`,
`ci_lower`, `ci_upper`, `sign_flips` — the proportion of replicates in
which the item's sign differed from the observed one), `boot_weights`
(iter x n_nodes matrix), `iter`, and `ci_level`. Has print and plot
methods.

## Examples

``` r
# \donttest{
set.seed(1)
df <- as.data.frame(matrix(rnorm(600), 100, 6))
names(df) <- c("a1", "a2", "a3", "b1", "b2", "b3")
cl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
fit <- build_mcml_pc(df, cl, aggregation = "loadings",
                     method = "cor")
#> Warning: Item(s) more strongly connected to another cluster than their own (possible misassignment): a1, a2, a3, b1, b2. See $loadings (misfit, cross_cluster).
#> Warning: Reverse-keyed item(s) flipped in composites: a2, b3. See $loadings (sign).
ls <- loading_stability(fit, iter = 50, seed = 1)
ls$summary
#>   node cluster     weight  boot_mean   boot_sd   ci_lower  ci_upper sign_flips
#> 1   a1       A  0.1405873  0.2119349 0.2931353 -0.4594061 0.4760502       0.18
#> 2   a2       A -0.3666270  0.1041219 0.3064605 -0.4215860 0.4384606       0.68
#> 3   a3       A  0.4927857  0.1710647 0.3171502 -0.4661878 0.4662406       0.24
#> 4   b1       B  0.2856282  0.2469936 0.2441808 -0.4553946 0.4688604       0.10
#> 5   b2       B  0.2708179  0.2257656 0.2509415 -0.4394226 0.4624575       0.12
#> 6   b3       B -0.4435540 -0.1735737 0.3268693 -0.4691528 0.4829590       0.26
# }
```
