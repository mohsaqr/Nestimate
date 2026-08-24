# GIMME: Group Iterative Multiple Model Estimation

Estimates person-specific directed networks from intensive longitudinal
data using the unified Structural Equation Modeling (uSEM) framework.
Implements a data-driven search that identifies:

1.  **Group-level paths**: Directed edges present for a majority
    (default 75\\

2.  **Individual-level paths**: Additional edges specific to each
    person, found after group paths are established.

Estimation is delegated to
[`idiographic::fit_gimme()`](https://pak.dynasite.org/idiographic/reference/fit_gimme.html),
the clean-room home of the temporal idiographic estimators, whose search
reproduces the upstream gimme package (\>= 10.0) exactly (verified at
tolerance 0 on path counts and per-person coefficient matrices in
idiographic's own parity suite). Uses `lavaan` for SEM estimation and
modification indices. Accepts a single data frame with an ID column (not
CSV directories).

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

  A `data.frame` in long format with columns for person ID, time-varying
  variables, and optionally a time/beep column.

- vars:

  Character vector of variable names to model.

- id:

  Character string naming the person-ID column.

- time:

  Character string naming the time/order column, or `NULL`. When
  provided, data is sorted by `id` then `time` before lagging.

- ar:

  Logical. If `TRUE` (default), autoregressive paths (each variable
  predicting itself at lag 1) are included as fixed paths.

- standardize:

  Logical. If `TRUE` (default `FALSE`), variables are standardized per
  person before estimation.

- groupcutoff:

  Numeric between 0 and 1. Proportion of individuals for whom a path
  must be significant to be added at group level. Default `0.75`.

- subcutoff:

  Numeric. Not used (reserved for future subgrouping).

- paths:

  Character vector of lavaan-syntax paths to force into the model (e.g.,
  `"V2~V1lag"`). Default `NULL`.

- exogenous:

  Character vector of variable names to treat as exogenous. Default
  `NULL`.

- hybrid:

  Logical. If `TRUE`, also searches residual covariances. Default
  `FALSE`.

- rmsea_cutoff:

  Numeric. RMSEA threshold for excellent fit (default 0.05).

- srmr_cutoff:

  Numeric. SRMR threshold for excellent fit (default 0.05).

- nnfi_cutoff:

  Numeric. NNFI/TLI threshold for excellent fit (default 0.95).

- cfi_cutoff:

  Numeric. CFI threshold for excellent fit (default 0.95).

- n_excellent:

  Integer. Number of fit indices that must be excellent to stop
  individual search. Default `2`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

An S3 object of class `"net_gimme"` containing:

- `temporal`:

  p x p matrix of group-level temporal (lagged) path counts – entry
  `[i,j]` = number of individuals with path j(t-1)-\>i(t).

- `contemporaneous`:

  p x p matrix of group-level contemporaneous path counts – entry
  `[i,j]` = number of individuals with path j(t)-\>i(t).

- `temporal_avg`, `contemporaneous_avg`:

  p x p group-average coefficient matrices.

- `coefs`:

  List of per-person p x 2p coefficient matrices (rows = endogenous,
  cols = `[lagged, contemporaneous]`).

- `psi`:

  List of per-person residual covariance matrices.

- `fit`:

  Data frame of per-person fit indices (chisq, df, pvalue, rmsea, srmr,
  nnfi, cfi, bic, aic, logl, status).

- `path_counts`:

  p x 2p matrix: how many individuals have each path.

- `paths`:

  List of per-person character vectors of lavaan path syntax.

- `group_paths`:

  Character vector of group-level paths found.

- `individual_paths`:

  List of per-person character vectors of individual-level paths (beyond
  group).

- `syntax`:

  List of per-person full lavaan syntax strings.

- `labels`:

  Character vector of variable names.

- `n_subjects`:

  Integer. Number of individuals.

- `n_obs`:

  Integer vector. Time points per individual.

- `config`:

  List of configuration parameters.

The object additionally carries idiographic's netobject fields
(`weights`, `nodes`, `edges`, ...) so it renders directly with cograph,
and dispatches to idiographic's `print`, `summary`, and `plot` methods.

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)

## Examples

``` r
# \donttest{
# Create simple panel data (3 subjects, 4 variables, 50 time points).
set.seed(42)
n_sub <- 3; n_t <- 50; vars <- paste0("V", 1:4)
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
#> Individual-level paths:  mean 1.3, range 0-3
#> 
#> Proportion of subjects with each path:
#> 
#>   Temporal [directed]
#>     weights [0.333, 1.000]  |  +6 / -0 edges
#>        V1 V2   V3   V4
#>     V1  1  0 0.00 0.00
#>     V2  0  1 0.33 0.33
#>     V3  0  0 1.00 0.00
#>     V4  0  0 0.00 1.00
#> 
#>   Contemporaneous [directed]
#>     weights [0.333, 0.333]  |  +2 / -0 edges
#>        V1   V2 V3   V4
#>     V1  0 0.00  0 0.00
#>     V2  0 0.00  0 0.00
#>     V3  0 0.00  0 0.33
#>     V4  0 0.33  0 0.00
#> 
#>   plot(x)  (faithful gimme-style mixed network) | plot(x, layer = "temporal") 
#>   edges(x) | nodes(x) | summary(x) | coefs(x) | matrices(x)
# }
```
