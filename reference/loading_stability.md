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
                     estimator = "cor")
#> Error: 'estimator' in ... is the lavaan estimator (e.g. "WLSMV") and applies only to fa_method = "cfa". The network estimator is selected by 'method'.
ls <- loading_stability(fit, iter = 50, seed = 1)
#> Error: object 'fit' not found
ls$summary
#> Error in ls$summary: object of type 'closure' is not subsettable
# }
```
