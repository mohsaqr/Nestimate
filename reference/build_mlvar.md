# Build a Multilevel Vector Autoregression (mlVAR) network

Estimates three networks from ESM/EMA panel data, matching
`mlVAR::mlVAR()` with `estimator = "lmer"`, `temporal = "fixed"`,
`contemporaneous = "fixed"` at machine precision: (1) a directed
temporal network of fixed-effect lagged regression coefficients, (2) an
undirected contemporaneous network of partial correlations among
residuals, and (3) an undirected between-subjects network of partial
correlations derived from the person-mean fixed effects.

## Usage

``` r
build_mlvar(
  data,
  vars,
  id,
  day = NULL,
  beep = NULL,
  lag = 1L,
  standardize = FALSE
)
```

## Arguments

- data:

  A `data.frame` containing the panel data.

- vars:

  Character vector of variable column names to model.

- id:

  Character string naming the person-ID column.

- day:

  Character string naming the day/session column, or `NULL`. When
  provided, lag pairs are only formed within the same day.

- beep:

  Character string naming the measurement-occasion column, or `NULL`.
  When `NULL`, row position within each (id, day) is used.

- lag:

  Integer. The lag order (default 1).

- standardize:

  Logical. If `TRUE`, each variable is grand-mean centered and divided
  by its pooled SD *before* augmentation. Default `FALSE`, matching
  `mlVAR::mlVAR(scale = FALSE)` - the only setting for which numerical
  equivalence has been validated.

## Value

A dual-class `c("net_mlvar", "netobject_group")` object - a named list
of three full netobjects, one per network, plus model-level metadata
stored as attributes. Each element is a standard
`c("netobject", "cograph_network")` weight-matrix wrapper (no raw
`$data`), so [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coefs()`](https://saqr.me/Nestimate/reference/coefs.md), and
`cograph::splot(fit$temporal)` work directly. See **Dispatch
limitation** for the verbs that do *not* work on this object
([`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
`centrality()`, reliability/stability). Structure:

- `fit$temporal`:

  Directed netobject for the `d x d` matrix of fixed-effect lagged
  coefficients. `$weights[i, j]` is the effect of variable j at t-lag on
  variable i at t. `method = "mlvar_temporal"`, `directed = TRUE`.

- `fit$contemporaneous`:

  Undirected netobject for the `d x d` partial-correlation network of
  within-person lmer residuals. `method = "mlvar_contemporaneous"`,
  `directed = FALSE`.

- `fit$between`:

  Undirected netobject for the `d x d` partial-correlation network of
  person means, derived from `D (I - Gamma)`.
  `method = "mlvar_between"`, `directed = FALSE`.

- `attr(fit, "coefs")` /
  [`coefs()`](https://saqr.me/Nestimate/reference/coefs.md):

  Tidy `data.frame` with one row per `(outcome, predictor)` pair and
  columns `outcome`, `predictor`, `beta`, `se`, `t`, `p`, `ci_lower`,
  `ci_upper`, `significant`. Filter, sort, or plot with base R or the
  tidyverse. Retrieve with `coefs(fit)`.

- `attr(fit, "n_obs")`:

  Number of rows in the augmented panel after na.omit.

- `attr(fit, "n_subjects")`:

  Number of unique subjects remaining.

- `attr(fit, "lag")`:

  Lag order used.

- `attr(fit, "standardize")`:

  Logical; whether pre-augmentation standardization was applied.

## Details

Estimation is delegated to
[`idiographic::fit_mlvar()`](https://pak.dynasite.org/idiographic/reference/fit_mlvar.html)
(the clean-room home of the temporal idiographic estimators), called
with `estimator = "lmer"`, `temporal = "fixed"`,
`contemporaneous = "fixed"`. The pipeline follows mlVAR's lmer path
exactly:

1.  Drop rows with NA in id/day/beep and optionally grand-mean
    standardize each variable.

2.  Expand the per-(id, day) beep grid and right-join original values,
    producing the augmented panel (`augData`).

3.  Add within-person lagged predictors (`L1_*`) and person-mean
    predictors (`PM_*`).

4.  For each outcome variable fit
    `lmer(y ~ within + between-except-own-PM + (1 | id))` with
    `REML = FALSE`. Collect the fixed-effect temporal matrix `B`,
    between-effect matrix `Gamma`, random-intercept SDs (`mu_SD`), and
    lmer residual SDs.

5.  Contemporaneous network:
    `cor2pcor(D %*% cov2cor(cor(resid)) %*% D)`.

6.  Between-subjects network:
    `cor2pcor(pseudoinverse(forcePositive(D (I - Gamma))))`.

Validated to machine precision (max_diff \< 1e-10) against
`mlVAR::mlVAR()` on 25 real ESM datasets from `openesm` and 20 simulated
configurations, and to exact equality (max_diff == 0) against the
pre-delegation Nestimate implementation on all layers, coefficients, and
observation counts across lag/standardize/day/beep configurations.

When the data carry no between-person variance (a random-intercept SD of
zero), the between-subjects network is not estimable. Since 0.9.0 this
raises a warning from idiographic and still returns the zero matrix by
convention; before 0.9.0 the zero matrix was returned silently.

## Dispatch limitation

There is no [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method for `net_mlvar` - plot a single constituent
(`cograph::splot(fit$temporal)`) instead. The three constituents are
matrix-wrapped and carry no `$data`, so the data-resampling and
data-reading verbs do **not** work on the fitted object or its parts:
[`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
[`certainty()`](https://saqr.me/Nestimate/reference/certainty.md),
[`network_reliability()`](https://saqr.me/Nestimate/reference/network_reliability.md),
[`centrality_stability()`](https://saqr.me/Nestimate/reference/centrality_stability.md)
and `centrality()` all need the source panel. Extract a constituent and
rebuild it through
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
if you need those. Use
[`coefs()`](https://saqr.me/Nestimate/reference/coefs.md) for the tidy
model output.

## See also

[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)

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
