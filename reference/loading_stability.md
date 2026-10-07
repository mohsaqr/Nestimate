# Composite-Weight Stability Under Case Resampling

**Experimental.** Bootstraps the item weights of a
[`build_mcml_pc`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
fit: rows of the raw data are resampled, the node-level network is
re-estimated each time, and the connectivity-based composite weights are
recomputed. Wide intervals mean the weighting (and therefore the
`"loadings"` macro network) should not be over-interpreted.

## Usage

``` r
loading_stability(x, iter = 200L, ci_level = 0.05, seed = NULL)

# S3 method for class 'pc_loading_stability'
print(x, digits = 3, ...)

# S3 method for class 'pc_loading_stability'
plot(x, ...)
```

## Arguments

- x:

  An `mcml_pc` object that carries raw data. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `pc_loading_stability`.

- iter:

  Integer. Bootstrap replicates (default 200; node-level re-estimation
  makes this heavier than a plain bootstrap).

- ci_level:

  Numeric. Significance level for percentile CIs (default 0.05).

- seed:

  Integer or NULL. RNG seed.

- digits:

  Number of digits to display (default 3).

- ...:

  For the S3 methods: further arguments passed to or from other methods.

## Value

An object of class `"pc_loading_stability"`: a list with `summary` (tidy
data frame: `node`, `cluster`, `weight`, `boot_mean`, `boot_sd`,
`ci_lower`, `ci_upper`, `sign_flips` - the proportion of replicates in
which the item's sign differed from the observed one), `boot_weights`
(iter x n_nodes matrix), `iter`, and `ci_level`. Has print and plot
methods.

In `print.pc_loading_stability()`: `x`, invisibly.

In `plot.pc_loading_stability()`: A ggplot object.

## Examples

``` r
set.seed(1)
f <- stats::rnorm(100)
g <- stats::rnorm(100)
df <- data.frame(a1 = f + stats::rnorm(100), a2 = f + stats::rnorm(100),
                 a3 = f + stats::rnorm(100), b1 = g + stats::rnorm(100),
                 b2 = g + stats::rnorm(100), b3 = g + stats::rnorm(100))
cl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
fit <- build_mcml_pc(df, cl, aggregation = "loadings",
                     method = "cor")
stability <- loading_stability(fit, iter = 50, seed = 1)
stability
#> Composite-Weight Stability (case bootstrap, experimental)
#>   50 replicates | 95% percentile CIs
#> 
#>  node cluster weight boot_mean boot_sd ci_lower ci_upper sign_flips
#>    a1       A  0.338     0.332   0.021    0.292    0.370          0
#>    a2       A  0.336     0.339   0.018    0.308    0.375          0
#>    a3       A  0.326     0.329   0.021    0.290    0.367          0
#>    b1       B  0.316     0.317   0.022    0.276    0.348          0
#>    b2       B  0.350     0.346   0.022    0.310    0.393          0
#>    b3       B  0.334     0.337   0.015    0.307    0.360          0
```
