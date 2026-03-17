# Design: Post-hoc Covariate Analysis in cluster_data()

## Summary

Add a `covariates` parameter to `cluster_data()` that runs multinomial logistic regression (`nnet::multinom`) on cluster assignments after clustering completes. Produces descriptive profiles per cluster and formal odds-ratio tests of covariate-membership associations. Covariates are post-hoc -- they do not influence cluster assignments.

## API

```r
cluster_data(data, k, ..., covariates = NULL)
```

### `covariates` parameter

Accepts five forms (all equivalent):

| Form | Example |
|------|---------|
| Formula | `~ Age + Gender` |
| Character vector | `c("Age", "Gender")` |
| String | `"Age + Gender"` |
| Data frame | `demographics` (all columns used) |
| `NULL` (default) | No covariate analysis |

### Covariate source resolution

When `covariates` is a formula, character vector, or string, column values are looked up in order:

1. `netobject$metadata` -- if input is a netobject
2. Non-sequence columns in the input data.frame itself
3. Error if columns not found

For `tna` and `cograph_network` inputs, covariate column lookup is not supported -- user must pass covariates as a data.frame. An informative error is thrown if formula/character/string covariates are used with these input types.

When `covariates` is a data.frame, it's used directly (must have `nrow(cov_df) == nrow(sequence_data)` after extraction).

### Input coercion

- Character columns are coerced to factors before fitting.
- Ordered factors are converted to unordered (to avoid polynomial contrasts).
- Intercept-only formulas (`~ 1`) are rejected with an error.
- Interaction terms (`Age * Gender`) are supported -- they pass through to `multinom()` and appear as coefficient names in output. A note in documentation warns that OR interpretation for interaction terms requires care.

## What runs

1. Clustering proceeds exactly as before (no change to distance computation or assignment).
2. If `covariates` is non-NULL:
   a. Resolve covariate data into a data.frame `cov_df`.
   b. Validate: `nrow(cov_df) == n`, columns exist, at least one non-constant column.
   c. Handle NAs: if any rows have NA in covariates, subset both `cov_df` and `assignments` to complete cases. Warn with the number of rows dropped.
   d. Guard: if `min(table(assignments_subset)) < ncol(model_matrix)`, warn that some clusters have fewer observations than parameters.
   e. Compute descriptive statistics per cluster per covariate.
   f. Fit `nnet::multinom(cluster ~ ..., data, trace = FALSE)`.
   g. Extract coefficients, SEs, compute z-values, p-values, CIs, odds ratios.
   h. Compute model fit: AIC, BIC, McFadden R-squared.

## Dependencies

`nnet` -- recommended base R package, always installed. Add to `Suggests` in DESCRIPTION and guard with `requireNamespace("nnet", quietly = TRUE)` + informative error, consistent with how `stringdist` and `cograph` are handled.

## Output: `$covariates` field in `net_clustering`

When covariates are provided, the returned `net_clustering` gains a `$covariates` list:

```r
cl$covariates$profiles     # list with $numeric (data.frame) and $categorical (data.frame)
cl$covariates$coefficients # data.frame: cluster, variable, estimate, std_error, odds_ratio,
                           #   ci_lower, ci_upper, z, p, sig
cl$covariates$fit          # list: aic, bic, deviance, mcfadden_r2, reference_cluster
cl$covariates$model        # raw nnet::multinom object
```

### profiles$numeric

Per-cluster descriptives for numeric covariates:

| cluster | n | pct | variable | mean | sd | median |
|---------|---|-----|----------|------|----|--------|
| 1 | 45 | 30.0 | Age | 22.4 | 3.1 | 22.0 |
| 1 | 45 | 30.0 | GPA | 3.21 | 0.45 | 3.30 |
| 2 | 58 | 38.7 | Age | 25.1 | 4.2 | 24.5 |

The `n` and `pct` columns repeat per variable within a cluster (intentional denormalization for easy printing/export).

### profiles$categorical

Per-cluster counts and percentages for categorical/factor covariates:

| cluster | n | variable | level | count | pct |
|---------|---|----------|-------|-------|-----|
| 1 | 45 | Gender | M | 18 | 40.0 |
| 1 | 45 | Gender | F | 27 | 60.0 |

### coefficients

One row per cluster x variable combination (reference cluster excluded):

| cluster | variable | estimate | std_error | odds_ratio | ci_lower | ci_upper | z | p | sig |
|---------|----------|----------|-----------|------------|----------|----------|---|---|-----|
| 2 | Age | 0.077 | 0.029 | 1.08 | 1.02 | 1.15 | 2.66 | 0.008 | ** |
| 2 | GenderM | -0.80 | 0.39 | 0.45 | 0.21 | 0.97 | -2.05 | 0.041 | * |

- `estimate`: log-odds coefficient
- `odds_ratio`: exp(estimate)
- `ci_lower`, `ci_upper`: 95% CI on odds ratio scale
- `sig`: `***` < 0.001, `**` < 0.01, `*` < 0.05, `.` < 0.1, `` otherwise

Reference cluster is Cluster 1 by default (first level of factor). Documented in output and print methods.

### fit

```r
list(aic, bic, deviance, mcfadden_r2, reference_cluster)
```

McFadden R-squared = 1 - (deviance / null_deviance). Note: McFadden R-squared values are typically much lower than OLS R-squared; values above 0.2 indicate excellent fit.

## S3 method updates

### print.net_clustering

When `$covariates` exists, append one line:

```
  Covariates:     Age, Gender (post-hoc, 2 predictors)
```

### summary.net_clustering

Return value changes when covariates are present: returns a list with `$cluster_stats` (existing data.frame) and `$covariates` (the new output), class `"summary.net_clustering"`. Without covariates, returns the existing data.frame (backwards compatible).

When `$covariates` exists, append after the per-cluster table:

```
Post-hoc Covariate Analysis (does not influence cluster membership)

Cluster Profiles:
  Numeric covariates:
    Cluster   N (%)         Age                    GPA
                            Mean (SD)   Median     Mean (SD)   Median
    1         45 (30%)      22.4 (3.1)  22.0       3.21 (0.45) 3.30
    2         58 (39%)      25.1 (4.2)  24.5       2.87 (0.62) 2.90
    3         47 (31%)      21.8 (2.8)  21.0       3.54 (0.38) 3.60

  Categorical covariates:
    Cluster   GenderM      Nationality_FI
              N (%)        N (%)
    1         18 (40%)     30 (67%)
    2         16 (28%)     42 (72%)
    3         29 (62%)     25 (53%)

Predictors of Membership (reference: Cluster 1):
  Cluster   Variable    OR   95% CI          p      Sig
  2         Age        1.08  [1.02, 1.15]   0.008   **
  2         GenderM    0.45  [0.21, 0.97]   0.041   *
  3         Age        0.97  [0.91, 1.04]   0.382
  3         GenderM    2.31  [1.05, 5.08]   0.037   *

Model: AIC = 142.3 | BIC = 158.7 | McFadden R-squared = 0.12

Note: Covariates are post-hoc and do not influence cluster assignments.
```

### plot.net_clustering

New type `"predictors"` -- odds ratio dot plot with CI whiskers, faceted by cluster, colored by significance. Horizontal line at OR = 1. If `type = "predictors"` is requested but no covariates were run, throw an informative error.

## Verification strategy

### k=2 cross-validation against glm(binomial)

For 2 clusters, multinomial logistic regression is equivalent to binary logistic regression:

```r
# Our implementation
cl <- cluster_data(df, k = 2, covariates = c("Age", "Gender"))
our <- cl$covariates$coefficients

# glm reference
fit_glm <- glm(factor(assignments) ~ Age + Gender, family = binomial, data = cov_df)
ref <- summary(fit_glm)$coefficients
```

Tolerances: coefficients and SEs at 1e-4, p-values at 1e-3 (different optimization: BFGS vs IRLS).

### k>2 verification against manual computation

Extract raw `summary(multinom_fit)$coefficients` and `$standard.errors`, compute z = coef/se, p = 2*(1-pnorm(|z|)), OR = exp(coef), CI = exp(coef +/- 1.96*se) by hand, compare to our output table.

### Descriptive statistics verification

Compare `$profiles$numeric` against `aggregate()` + manual `sd()`/`median()` per cluster. Compare `$profiles$categorical` against `table()` per cluster.

### Input form equivalence

All five input forms (formula, character vector, string, data.frame, NULL) must produce identical `$covariates` output when given the same underlying data.

### S3 method tests

- `print()` output contains "Covariates:" line when covariates present
- `summary()` output contains "Post-hoc Covariate Analysis" and "Predictors of Membership"
- `plot(cl, type = "predictors")` returns ggplot when covariates present
- `plot(cl, type = "predictors")` errors informatively when no covariates

### build_network round-trip

`build_network(cl)` on a clustering with covariates preserves `$covariates` in the clustering metadata via `attr(grp, "clustering")`.

## Files to modify

| File | Changes |
|------|---------|
| `R/cluster_data.R` | Add `covariates` param, `.resolve_covariates()`, `.run_covariate_analysis()`, `.compute_cluster_profiles()`, update `print`/`summary`/`plot` S3 methods |
| `tests/testthat/test-cluster_data.R` | Add test blocks: covariate input forms, profiles, coefficients, k=2 glm cross-validation, k>2 manual verification, S3 methods, edge cases |
| `DESCRIPTION` | Add `nnet` to Suggests |
| `NAMESPACE` | No changes (no new exports) |

## Edge cases

- Single covariate: formula with one predictor
- All-categorical covariates: no numeric profile section printed
- All-numeric covariates: no categorical profile section printed
- Covariate with NAs: subset both `cov_df` and `assignments` to complete cases, warn with count dropped
- Constant covariate: error before fitting
- Intercept-only formula (`~ 1`): error
- k=2: verify binary equivalence
- Small cluster: warn if any cluster has fewer observations than number of parameters
- Perfect separation: capture and relay multinom convergence warning
- `tna`/`cograph_network` input with non-dataframe covariates: informative error
- Covariate data.frame row count mismatch: clear error stating expected N
