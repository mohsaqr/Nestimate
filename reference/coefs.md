# Tidy coefficients from a fitted mlvar model

Generic accessor for the tidy coefficient table stored on a
[`build_mlvar()`](https://saqr.me/Nestimate/reference/build_mlvar.md)
result. Returns a `data.frame` with one row per `(outcome, predictor)`
pair and columns `outcome`, `predictor`, `beta`, `se`, `t`, `p`,
`ci_lower`, `ci_upper`, `significant`.

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

  A fitted model object - currently only `net_mlvar` is supported.

- ...:

  Unused.

## Value

A tidy `data.frame` of coefficient estimates.

## Details

Only the within-person (temporal) coefficients are tabulated - these are
the lagged fixed effects that populate `fit$temporal`. The
between-subjects effects that go into `fit$between` are handled via the
`D (I - Gamma)` transformation and are not exposed as a separate tidy
table.

## Examples

``` r
# \donttest{
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
#>           network n_nodes n_edges density mean_abs_weight n_positive n_negative
#> 1        temporal       3       6       1      0.11569402          2          4
#> 2 contemporaneous       3       3       1      0.05304465          1          2
#> 3         between       3       3       1      0.15389239          2          1
# }
```
