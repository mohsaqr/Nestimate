# Build a Multilevel Vector Autoregression (mlVAR) network

`build_mlvar()` fits a multilevel vector autoregression to intensive
longitudinal panel data (ESM or EMA: many people, each measured on many
occasions) and returns its three networks. The *temporal* network is
directed and holds the fixed-effect lagged regression coefficients: the
weight in row `i`, column `j` is the coefficient of variable `i` at an
earlier occasion on variable `j` now, so each edge runs from the earlier
variable to the one it predicts. The *contemporaneous* network is
undirected: the partial correlations of the within-person residuals at
the same occasion. The *between-subjects* network is undirected too: the
partial correlations of the person means. Each variable is regressed, in
a mixed model, on every variable `lag` occasions earlier and on the
other variables' person means; the temporal and contemporaneous effects
are fixed and only the intercepts are random.

The temporal and contemporaneous edges describe within-person relations
averaged over people; the between-subjects edges describe how people who
differ in one variable on average differ in another (Epskamp et al.,
2018). Only random intercepts are estimated, so the networks do not show
how people differ in their dynamics. A temporal edge is a lagged
predictive association, and a contemporaneous edge has no temporal
order.
[`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md) gives
each temporal coefficient with its standard error, p-value and
confidence interval; check them before interpreting any single temporal
edge.

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

# S3 method for class 'net_mlvar'
print(x, ...)

# S3 method for class 'net_mlvar'
summary(object, ...)

# S3 method for class 'net_mlvar'
plot(x, layer = NULL, ...)
```

## Arguments

- data:

  The panel data: a data frame with one row per measurement occasion. A
  column named in `vars`, `id`, `day` or `beep` that is not in it raises
  an error.

- vars:

  The columns to model, at least two.

- id:

  The column identifying the person each row belongs to.

- day:

  The column identifying the day (or session) of each row. When given,
  lagged pairs are formed only within the same day. Default `NULL`.

- beep:

  The column numbering the measurement occasions. When `NULL` (default),
  the row position within each person (and day) is used.

- lag:

  How many occasions back the predictors of the temporal network are, a
  whole number of at least 1: `lag = 2` regresses each variable on every
  variable two occasions earlier (that lag only, not lags 1 and 2).
  Default `1`.

- standardize:

  Whether to grand-mean centre each variable and divide it by its pooled
  SD before the panel is augmented, as `mlVAR::mlVAR(scale = TRUE)`
  does. Default `FALSE`, matching `mlVAR::mlVAR(scale = FALSE)` (mlVAR's
  own default is `TRUE`); `FALSE` is the only setting for which
  agreement with mlVAR has been validated.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `net_mlvar`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: not used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_mlvar`.

- layer:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: `"temporal"`, `"contemporaneous"` or `"between"` to draw that
  network alone. Default `NULL` draws all three side by side.

## Value

A `net_mlvar` object (class `c("net_mlvar", "netobject_group")`): a
named list of three networks, with the model-level results stored as
attributes. [`print()`](https://rdrr.io/r/base/print.html) gives the
number of subjects, observations and variables, the lag, and how many
temporal edges are significant;
[`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md)
returns the temporal coefficients, one row per outcome and predictor;
[`summary()`](https://rdrr.io/r/base/summary.html) prints the three
weight matrices and the significant temporal edges and returns the same
table as
[`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md),
visibly, so at the console the matrices are followed by the table.
idiographic registers its own
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods for
`net_mlvar`; Nestimate loads it first and registers its methods after,
so these two are the ones that run. Each network is a `netobject` (class
`c("netobject", "cograph_network")`) built from its weight matrix,
without data; the *Dispatch limitation* section lists the verbs that do
not work on it. The parts:

- `temporal`:

  Directed network of the `d x d` matrix of fixed-effect lagged
  coefficients: the weight in row `i`, column `j` is the effect of
  variable `i` at `t - lag` on variable `j` at `t`, so edges run from
  predictor to outcome. Method `"mlvar_temporal"`.

- `contemporaneous`:

  Undirected `d x d` network of the partial correlations of the
  within-person lmer residuals. Method `"mlvar_contemporaneous"`.

- `between`:

  Undirected `d x d` network of the partial correlations of the person
  means, derived from `D (I - Gamma)`. Method `"mlvar_between"`.

- coefficient table (attribute `coefs`):

  The data frame
  [`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md)
  returns: one row per `(outcome, predictor)` pair, with columns
  `outcome`, `predictor`, `beta`, `se`, `t`, `p`, `ci_lower`, `ci_upper`
  and `significant`.

- attribute `n_obs`:

  The number of rows of the augmented panel after rows with missing
  values are removed.

- attribute `n_subjects`:

  The number of subjects remaining.

- attribute `lag`:

  The lag used.

- attribute `standardize`:

  Whether the variables were standardised before augmentation.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly.

## Details

Estimation is done by
[`idiographic::fit_mlvar()`](https://pak.dynasite.org/idiographic/reference/fit_mlvar.html),
called with `estimator = "lmer"`, `temporal = "fixed"` and
`contemporaneous = "fixed"`. It follows mlVAR's lmer path:

1.  Drop rows with NA in id/day/beep and, with `standardize = TRUE`,
    grand-mean standardise each variable.

2.  Expand the per-(id, day) beep grid and right-join the original
    values, giving the augmented panel (`augData`).

3.  Add within-person lagged predictors (`L1_*`) and person-mean
    predictors (`PM_*`).

4.  For each outcome variable fit
    `lmer(y ~ within + between-except-own-PM + (1 | id))` with
    `REML = FALSE`. Collect the fixed-effect temporal matrix `B`, the
    between-effect matrix `Gamma`, the random-intercept SDs (`mu_SD`)
    and the lmer residual SDs.

5.  Contemporaneous network:
    `cor2pcor(D %*% cov2cor(cor(resid)) %*% D)`.

6.  Between-subjects network:
    `cor2pcor(pseudoinverse(forcePositive(D (I - Gamma))))`.

The result agrees with `mlVAR::mlVAR()` to machine precision (largest
difference below 1e-10) on 25 real ESM datasets from `openesm` and 20
simulated configurations, and exactly with the implementation Nestimate
carried before 0.9.0 on all layers, coefficients and observation counts,
across settings of `lag`, `standardize`, `day` and `beep`.

The temporal weight matrix has predictors in rows and outcomes in
columns, the transpose of mlVAR's `Beta` and the same matrix as
`mlVAR::getNet(, "temporal")`. A Nestimate network reads its weight
matrix from row to column, so the edges of the temporal network (its
edge list and
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md))
run from the predictor at `t - lag` to the outcome at `t`. Nestimate
1.0.9 and earlier did not transpose it, so every temporal edge pointed
from the outcome to its predictor.

When the data carry no between-person variance (a random-intercept SD of
zero), the between-subjects network cannot be estimated. Since 0.9.0
idiographic then raises a warning and the zero matrix is returned by
convention; before 0.9.0 the zero matrix was returned silently.

## Dispatch limitation

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the three
networks side by side, or one with `layer`, with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
when cograph is installed and in base graphics otherwise; each
constituent can also be drawn on its own
([`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
on the temporal network). The constituents are wrappers around their
weight matrices and carry no data, so the verbs that resample or reread
the data,
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
[`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
and
[`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md),
raise an error on the fit and on each constituent, and
[`certainty()`](https://pak.dynasite.org/Nestimate/reference/certainty.md)
refuses them as not transition networks.
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) read the weights
and work on each constituent.
[`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md) gives
the model output.

## Methods

- `summary.net_mlvar()` prints the three weight matrices and the
  significant temporal edges, then returns the coefficient table of
  [`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md).

## References

Epskamp, S., Waldorp, L. J., Mõttus, R., & Borsboom, D. (2018). The
Gaussian graphical model in cross-sectional and time-series data.
*Multivariate Behavioral Research*, 53(4), 453–480.
[doi:10.1080/00273171.2018.1454823](https://doi.org/10.1080/00273171.2018.1454823)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

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
