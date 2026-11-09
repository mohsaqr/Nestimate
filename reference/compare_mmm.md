# Compare MMM fits across different k

`compare_mmm()` fits a mixed Markov model with
[`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
for each number of subgroups in `k` and returns one row per fit with its
log-likelihood, information criteria (AIC, BIC, ICL), overall average
posterior probability (AvePP), entropy and convergence. Every fit uses
the same data and the same
[`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
settings.

Lower AIC, BIC and ICL indicate more support relative to the other fits.
ICL adds a penalty for uncertain classification to BIC, so it favours
subgroups that separate cleanly (Biernacki, Celeux & Govaert, 2000);
AvePP and entropy describe that separation within each fit. The criteria
do not settle `k` on their own: check that every fit converged, that the
subgroups' transition matrices of the chosen fit differ in interpretable
ways ([`summary()`](https://rdrr.io/r/base/summary.html) of that fit
prints them), and, when the criteria keep falling at the largest `k`,
extend the range.

## Usage

``` r
compare_mmm(
  data,
  k = 2:5,
  return_fits = FALSE,
  n_starts = 50L,
  max_iter = 200L,
  tol = 1e-06,
  smooth = 0.01,
  seed = NULL,
  covariates = NULL,
  covariate_effect = c("em", "posthoc"),
  estimator = c("auto", "firth", "multinom", "chisq")
)

# S3 method for class 'mmm_compare'
print(x, ...)

# S3 method for class 'mmm_compare'
summary(object, ...)

# S3 method for class 'mmm_compare'
plot(x, ...)
```

## Arguments

- data:

  The sequences, in any form
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  accepts: a wide data frame, a `netobject` or `cograph_network` built
  from sequences, or a `tna` model.

- k:

  The numbers of subgroups to fit, a non-empty vector of whole numbers
  of at least 2; one model is fitted per value. Default `2:5`.

- return_fits:

  Whether to keep the fitted models. With `TRUE`, the `net_mmm` fits are
  attached to the result as its `"fits"` attribute, a list named by `k`,
  so the chosen model needs no second run of the EM. Default `FALSE`,
  which returns only the comparison table.

- n_starts:

  How many starting configurations to try, a positive whole number. More
  starts make it likelier that the best local maximum is found. Default
  `50`.

- max_iter:

  The most EM iterations a start may run, a positive whole number.
  Default `200`.

- tol:

  The convergence tolerance: EM stops when the log-likelihood changes by
  less than this between iterations (with `covariate_effect = "em"`,
  also when its relative change does). A finite positive number; default
  `1e-6`.

- smooth:

  A pseudo-count added to every transition count before the matrices are
  estimated, so a transition unseen in a subgroup gets a small
  probability instead of zero. A finite non-negative number; default
  `0.01`.

- seed:

  A seed for reproducible starts, applied afresh to each fit (every `k`
  starts from the same random state). Default `NULL`.

- covariates:

  Variables describing each sequence (such as a course or a grade),
  given as a column name, a character vector of names, a one-sided
  formula or a data frame, as in
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md).
  For a `netobject`, names are looked up in its metadata first, so a
  typical call is `build_mmm(net, k = 3, covariates = "session_label")`.
  How they enter the model is set by `covariate_effect`. Default `NULL`.

- covariate_effect:

  How `covariates` enter the model. `"em"` (default) makes the mixing
  proportions depend on them during EM, so they shape the subgroups
  found and add `(k - 1) * p` parameters for `p` covariate columns; rows
  with a missing covariate are dropped first, with a warning.
  `"posthoc"` fits a plain mixture on every sequence and relates the
  covariates to the subgroups only afterwards, with a multinomial logit,
  so neither their values nor their missingness change the subgroups.
  Ignored without `covariates`.

- estimator:

  How the after-fit multinomial logit is estimated; it never affects the
  EM. `"auto"` (default) uses Firth's penalised likelihood
  ([`brglm2::brmultinom()`](https://rdrr.io/pkg/brglm2/man/brmultinom.html),
  finite under separation) when a cell of the subgroup-by-covariate
  table has fewer than 5 sequences, and
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
  otherwise. `"firth"` and `"multinom"` force one of them (`"multinom"`
  warns about separation risk), and `"chisq"` runs descriptive tests
  instead of a logit.
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  has the details. Ignored without `covariates`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  `mmm_compare` table.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: not
  used; any extra argument raises an error.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  `mmm_compare` table.

## Value

An `mmm_compare` object, a data frame with one row per value of `k`.
[`print()`](https://rdrr.io/r/base/print.html) shows the table with a
`best` column marking the lowest BIC (`"<-- BIC"`) and, when it falls on
a different row, the lowest ICL (`"<-- ICL"`), and returns that marked
table invisibly. [`summary()`](https://rdrr.io/r/base/summary.html)
returns the same table as a plain data frame with a `best` column
holding `"BIC"`, `"ICL"` or `""`; when one row has both minima it is
marked `"BIC"` only.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws BIC, AIC
and ICL against `k`, and returns the `ggplot` invisibly after printing
it. The table's columns are `k`, `log_likelihood`, `AIC`, `BIC`, `ICL`,
`AvePP`, `Entropy` and `converged`; `return_fits = TRUE` adds the
`"fits"` attribute.

## References

Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture
model for clustering with the integrated completed likelihood. *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 22(7),
719–725. [doi:10.1109/34.865189](https://doi.org/10.1109/34.865189)

## See also

[`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md),
[`cluster_mmm()`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md)

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
comp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
comp
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -62.32972      158.6594 182.4798 184.6518 0.9647836 0.2121121  TRUE    
#>  3 -62.32987      176.6597 213.0909 217.2174 0.9343265 0.2487206 FALSE    
#>  best   
#>  <-- BIC
#>         
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
comp <- compare_mmm(seqs, k = 2:3, seed = 42)
print(comp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy    converged
#>  2 -120.8229      275.6458 299.4661 299.9897 0.9920817 0.03537724 TRUE     
#>  3 -115.1023      282.2045 318.6357 319.2157 0.9908021 0.03304310 TRUE     
#>  best   
#>  <-- BIC
#>         

# Retain the fits so the chosen model needs no re-run; summary() marks
# the minimum-BIC and minimum-ICL rows in its `best` column.
comp_with_fits <- compare_mmm(seqs, k = 2:3, seed = 42, return_fits = TRUE)
summary(comp_with_fits)
#>   k log_likelihood      AIC      BIC      ICL     AvePP    Entropy converged
#> 1 2      -120.8229 275.6458 299.4661 299.9897 0.9920817 0.03537724      TRUE
#> 2 3      -115.1023 282.2045 318.6357 319.2157 0.9908021 0.03304310      TRUE
#>   best
#> 1  BIC
#> 2     
# }
```
