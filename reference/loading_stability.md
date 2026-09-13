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
`ci_lower`, `ci_upper`, `sign_flips` - the proportion of replicates in
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
stability <- loading_stability(fit, iter = 50, seed = 1)
stability
#> Composite-Weight Stability (case bootstrap, experimental)
#>   50 replicates | 95% percentile CIs
#> 
#>  node cluster weight boot_mean boot_sd ci_lower ci_upper sign_flips
#>    a1       A  0.141     0.212   0.293   -0.459    0.476       0.18
#>    a2       A -0.367     0.104   0.306   -0.422    0.438       0.68
#>    a3       A  0.493     0.171   0.317   -0.466    0.466       0.24
#>    b1       B  0.286     0.247   0.244   -0.455    0.469       0.10
#>    b2       B  0.271     0.226   0.251   -0.439    0.462       0.12
#>    b3       B -0.444    -0.174   0.327   -0.469    0.483       0.26
# }
```
