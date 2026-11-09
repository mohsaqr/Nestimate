# GIMME: Group Iterative Multiple Model Estimation

`build_gimme()` estimates a directed network for each person from
intensive longitudinal data with GIMME (group iterative multiple model
estimation; Gates & Molenaar, 2012), in the unified structural equation
model (uSEM): each variable at time \\t\\ is regressed on the variables
at \\t - 1\\ (temporal, lagged paths) and on the other variables at
\\t\\ (contemporaneous paths). The search runs in two stages:

1.  **Group-level paths**: directed paths added for everyone because
    they would improve the model for a large share of the people (more
    than the share `groupcutoff`, by default `0.75`, of the people whose
    models converged).

2.  **Individual-level paths**: further paths specific to each person,
    added after the group paths until that person's model fits well (see
    `n_excellent`).

It takes one long data frame with a person-ID column, not a directory of
CSV files.

The temporal and contemporaneous count matrices say how many people have
each path, which shows how common a path is in the sample, not how
strong it is; the coefficients are estimated per person. Group paths are
shared by construction, while individual paths are what the search adds
for one person's series: before reading them, check each person's number
of time points and fit statistics, which the result carries.

## Usage

``` r
build_gimme(
  data,
  vars,
  id,
  time = NULL,
  ar = TRUE,
  standardize = FALSE,
  groupcutoff = 0.75,
  subcutoff = 0.5,
  paths = NULL,
  exogenous = NULL,
  hybrid = FALSE,
  rmsea_cutoff = 0.05,
  srmr_cutoff = 0.05,
  nnfi_cutoff = 0.95,
  cfi_cutoff = 0.95,
  n_excellent = 2L,
  seed = NULL
)
```

## Arguments

- data:

  A data frame in long format: one row per person and time point, with
  the person-ID column, the variables and, optionally, a time or beep
  column.

- vars:

  The names of the variables to model, at least two, all columns of
  `data`; a missing name raises an error.

- id:

  The name of the person-ID column. The data must hold at least two
  people; a person with fewer than three time points is left out of the
  analysis without a warning.

- time:

  The name of the time or order column. When given, the rows are sorted
  by `id` and then `time` before the lagged variables are formed; with
  `NULL` (default), the row order is used.

- ar:

  Whether each variable's autoregressive path (itself at lag 1) is fixed
  into every model. Default `TRUE`.

- standardize:

  Whether the variables are standardised within each person before
  estimation. Default `FALSE`. The coefficients come from lavaan's
  standardised solution either way, so this changes only the
  scale-dependent fit statistics (the log-likelihood, AIC and BIC); the
  chi-square, RMSEA, SRMR, NNFI and CFI do not change.

- groupcutoff:

  A path is added at the group level when it is significant for more
  than this share of the people whose models converged, a number in (0,
  1\]. No group path is added when half or fewer of the models converge.
  Default `0.75`.

- subcutoff:

  Not used: reserved for subgrouping, which neither Nestimate nor
  [`idiographic::fit_gimme()`](https://pak.dynasite.org/idiographic/reference/fit_gimme.html)
  implements, and accepted for compatibility with gimme; it is not
  passed on. Default `0.50`.

- paths:

  Paths, in lavaan syntax, to fix into every model, such as
  `"V2~V1lag"`. Default `NULL`.

- exogenous:

  The variables treated as exogenous (predictors only), all among
  `vars`; at least one variable must remain endogenous, and either rule
  broken raises an error. Default `NULL`.

- hybrid:

  Whether the search also considers residual covariances between the
  variables at the same time point. Default `FALSE`.

- rmsea_cutoff:

  The RMSEA at or below which a person's fit counts as excellent.
  Default `0.05`.

- srmr_cutoff:

  The SRMR at or below which a person's fit counts as excellent. Default
  `0.05`.

- nnfi_cutoff:

  The NNFI (TLI) at or above which a person's fit counts as excellent.
  Default `0.95`.

- cfi_cutoff:

  The CFI at or above which a person's fit counts as excellent. Default
  `0.95`.

- n_excellent:

  How many of the four fit indices must be excellent for a person's
  individual search to stop, a whole number from 1 to 4. Default `2`.

- seed:

  An integer that makes the estimation reproducible; the session's own
  random stream is left as it was. Default `NULL`.

## Value

The `net_gimme` object (class
`c("net_gimme", "cograph_network", "list")`) that
[`idiographic::fit_gimme()`](https://pak.dynasite.org/idiographic/reference/fit_gimme.html)
returns; its [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods are
idiographic's, not Nestimate's.
[`print()`](https://rdrr.io/r/base/print.html) shows the group-level
paths, the number of individual paths and the share of people with each
path; [`summary()`](https://rdrr.io/r/base/summary.html) returns a data
frame with one row per network layer (temporal, contemporaneous) and its
size, density and weights;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
network.
[`idiographic::edges()`](https://pak.dynasite.org/idiographic/reference/edges.html)
gives one row per path with its level (group or individual), and
[`idiographic::coefs()`](https://pak.dynasite.org/idiographic/reference/coefs.html)
one row per person and path with its coefficient. The object is a
*superset* of the pre-0.9.0 in-package field contract: every field below
is present, alongside idiographic's additions (`contemp_cov`,
`contemp_cov_avg`, `contemp_is_cov`). Its fields:

- `temporal`:

  p x p matrix of temporal (lagged) path counts: entry `[i, j]` is the
  number of people with the path j(t-1) -\> i(t), group and individual
  paths alike.

- `contemporaneous`:

  p x p matrix of contemporaneous path counts: entry `[i, j]` is the
  number of people with the path j(t) -\> i(t).

- `temporal_avg`, `contemporaneous_avg`:

  p x p group-average coefficient matrices.

- `coefs`:

  List of per-person p x 2p coefficient matrices (rows the endogenous
  variables, columns the lagged and then the contemporaneous variables).

- `psi`:

  List of per-person residual covariance matrices.

- `fit`:

  Data frame of per-person fit indices (`file`, the person, then chisq,
  df, pvalue, rmsea, srmr, nnfi, cfi, bic, aic, logl and status).

- `path_counts`:

  p x 2p matrix: how many people have each path.

- `paths`:

  List of per-person character vectors of lavaan path syntax.

- `group_paths`:

  Character vector of the group-level paths found.

- `individual_paths`:

  List of per-person character vectors of the individual-level paths
  (beyond the group paths).

- `syntax`:

  List of per-person full lavaan syntax strings.

- `labels`:

  Character vector of variable names.

- `n_subjects`:

  Number of people.

- `n_obs`:

  Integer vector of the lagged time points per person (its rows minus
  one).

- `config`:

  List of the settings used, including the corrections and cutoffs that
  are not arguments here.

The object also carries idiographic's network fields (`nodes`, `edges`,
`directed`, `data`, `meta`, `node_groups`, and a `weights` entry that is
NULL) for rendering with cograph.

## Details

Estimation is delegated to
[`idiographic::fit_gimme()`](https://pak.dynasite.org/idiographic/reference/fit_gimme.html),
the clean-room home of the temporal idiographic estimators, whose search
reproduces the upstream gimme package (\>= 10.0) exactly (verified at
tolerance 0 on path counts and per-person coefficient matrices in
idiographic's own parity suite). The models are fitted with lavaan,
which also gives the modification indices the search uses;
`build_gimme()` raises an error when lavaan is not installed. The
settings `build_gimme()` does not expose (such as the
multiple-comparison corrections and `alpha`) keep
[`idiographic::fit_gimme()`](https://pak.dynasite.org/idiographic/reference/fit_gimme.html)'s
defaults.

## Results changed in 0.9.0

Before 0.9.0 the search ran in an in-package implementation that was
*not* upstream-gimme-exact. Delegating to
[`idiographic::fit_gimme()`](https://pak.dynasite.org/idiographic/reference/fit_gimme.html)
changed which paths the search selects on the same data –
individual-level paths in particular – so numeric results are not
comparable with Nestimate \<= 0.8.5. The returned object also gained
fields (see **Value**); nothing was removed.

## References

Gates, K. M., & Molenaar, P. C. M. (2012). Group search algorithm
recovers effective connectivity maps for individuals in homogeneous and
heterogeneous samples. *NeuroImage*, 63(1), 310–319.
[doi:10.1016/j.neuroimage.2012.06.026](https://doi.org/10.1016/j.neuroimage.2012.06.026)

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
# \donttest{
# Create simple panel data (3 subjects, 4 variables, 30 time points).
set.seed(42)
n_sub <- 3; n_t <- 30; vars <- paste0("V", 1:4)
rows <- lapply(seq_len(n_sub), function(i) {
  d <- as.data.frame(matrix(rnorm(n_t * 4), ncol = 4))
  names(d) <- vars; d$id <- i; d
})
panel <- do.call(rbind, rows)
res <- build_gimme(panel, vars = vars, id = "id")
print(res)
#> GIMME Network Analysis
#> ------------------------------ 
#> Subjects:   3 
#> Variables:  4  ( V1, V2, V3, V4 )
#> AR paths:   yes 
#> Hybrid:     no 
#> 
#> Group-level paths found: 0 
#> 
#> Individual-level paths:  mean 0.3, range 0-1
#> 
#> Proportion of subjects with each path:
#> 
#>   Temporal [directed]
#>     weights [1.000, 1.000]  |  +4 / -0 edges
#>        V1 V2 V3 V4
#>     V1  1  0  0  0
#>     V2  0  1  0  0
#>     V3  0  0  1  0
#>     V4  0  0  0  1
#> 
#>   Contemporaneous [directed]
#>     weights [0.333, 0.333]  |  +1 / -0 edges
#>        V1 V2 V3   V4
#>     V1  0  0  0 0.33
#>     V2  0  0  0 0.00
#>     V3  0  0  0 0.00
#>     V4  0  0  0 0.00
#> 
#>   plot(x)  (faithful gimme-style mixed network) | plot(x, layer = "temporal") 
#>   edges(x) | nodes(x) | summary(x) | coefs(x) | matrices(x)
# }
```
