# Tidy coefficients from a fitted mlvar model

`coefs()` returns the temporal coefficients of a
[`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
fit as a data frame, one row per `(outcome, predictor)` pair: the
fixed-effect lagged regression coefficients that make up the temporal
network, with their standard errors and tests.

Only the within-person (temporal) coefficients are tabulated. The
between-subjects effects that enter the between-subjects network through
the `D (I - Gamma)` transformation are not given as a table of their
own. The p-values are not adjusted for the number of coefficients, and
`significant` is `p < 0.05` for each coefficient on its own.

## Usage

``` r
coefs(x, ...)

# S3 method for class 'net_mlvar'
coefs(x, ...)

# Default S3 method
coefs(x, ...)
```

## Arguments

- x:

  A fitted model; only a `net_mlvar` from
  [`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  has a method, and any other object raises an error.

- ...:

  Not used.

## Value

A data frame with one row per `(outcome, predictor)` pair and the
columns `outcome` (the variable predicted), `predictor` (the variable at
the earlier occasion), `beta` (the fixed-effect coefficient), `se` (its
standard error), `t` (`beta / se`), `p` (the two-sided p-value from the
normal distribution, `2 * pnorm(-abs(t))`, which stays above zero for
large `t`), `ci_lower` and `ci_upper` (the 95% Wald interval, `beta`
plus or minus `qnorm(0.975)` standard errors) and `significant`
(`p < 0.05`).

## Examples

``` r
# A three-variable ESM panel: 20 people x 20 beeps. `tired` is driven by
# `happy` one beep earlier, so the temporal network should recover it.
if (requireNamespace("lme4", quietly = TRUE)) {
  set.seed(1)
  n_beep <- 20
  ar1 <- function(n, phi) as.numeric(stats::filter(stats::rnorm(n), phi,
                                                   method = "recursive"))
  panel <- do.call(rbind, lapply(seq_len(20), function(i) {
    happy <- ar1(n_beep, 0.4)
    data.frame(
      id    = i,
      beep  = seq_len(n_beep),
      happy = happy + stats::rnorm(1),
      calm  = ar1(n_beep, 0.3) + stats::rnorm(1),
      tired = 0.5 * c(0, happy[-n_beep]) + stats::rnorm(n_beep) +
              stats::rnorm(1)
    )
  }))
  fit <- build_mlvar(panel, vars = c("happy", "calm", "tired"),
                     id = "id", beep = "beep")
  fit
  coefs(fit)
  summary(fit)
}
#> === mlVAR Summary ===
#> Subjects: 20  | Observations: 380  | Variables: 3  | Lag: 1 
#> Standardized: FALSE 
#> 
#> --- Temporal Network (rows: predictor at t - lag; columns: outcome) ---
#>         happy    calm   tired
#> happy  0.3034 -0.0699  0.4950
#> calm  -0.0521  0.2252 -0.0236
#> tired  0.0025 -0.0511 -0.0898
#> 
#> Significant temporal edges (p < 0.05):
#>   outcome predictor   beta     se     t        p
#> 1   happy     happy 0.3034 0.0509 5.961 2.50e-09
#> 2    calm      calm 0.2252 0.0507 4.446 8.76e-06
#> 3   tired     happy 0.4950 0.0513 9.656 4.62e-22
#> 
#> --- Contemporaneous Network ---
#>         happy    calm   tired
#> happy  0.0000 -0.1005  0.0356
#> calm  -0.1005  0.0000 -0.0230
#> tired  0.0356 -0.0230  0.0000
#> 
#> --- Between-Subjects Network ---
#>         happy    calm  tired
#> happy  0.0000 -0.2812 0.0967
#> calm  -0.2812  0.0000 0.0838
#> tired  0.0967  0.0838 0.0000
#> 
#>   outcome predictor         beta         se           t            p
#> 1   happy     happy  0.303439780 0.05090043  5.96143847 2.500271e-09
#> 2   happy      calm -0.052058528 0.05405686 -0.96303275 3.355311e-01
#> 3   happy     tired  0.002524451 0.04659591  0.05417753 9.567937e-01
#> 4    calm     happy -0.069875366 0.04770161 -1.46484298 1.429638e-01
#> 5    calm      calm  0.225222041 0.05065968  4.44578533 8.757136e-06
#> 6    calm     tired -0.051100146 0.04366760 -1.17020726 2.419176e-01
#> 7   tired     happy  0.494960751 0.05125695  9.65646054 4.615350e-22
#> 8   tired      calm -0.023644879 0.05443549 -0.43436510 6.640233e-01
#> 9   tired     tired -0.089803421 0.04692228 -1.91387584 5.563602e-02
#>      ci_lower    ci_upper significant
#> 1  0.20367677 0.403202790        TRUE
#> 2 -0.15800803 0.053890974       FALSE
#> 3 -0.08880185 0.093850752       FALSE
#> 4 -0.16336880 0.023618068       FALSE
#> 5  0.12593090 0.324513180        TRUE
#> 6 -0.13668707 0.034486783       FALSE
#> 7  0.39449897 0.595422534        TRUE
#> 8 -0.13033649 0.083046729       FALSE
#> 9 -0.18176940 0.002162561       FALSE
```
