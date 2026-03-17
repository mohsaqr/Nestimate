# Design: Covariate-Integrated Mixture of Markov Models (build_mmm)

## Summary

Extend `build_mmm()` to support covariates that influence cluster membership probabilities within the EM algorithm. Instead of fixed mixing proportions shared across all sequences, each sequence gets covariate-dependent mixing proportions via multinomial logistic regression estimated jointly with transition parameters. This is a proper generative model, not post-hoc analysis (Dayton & Macready 1988; Vermunt 2010).

## Mathematical Model

### Current (no covariates)

Mixing proportions are fixed across all sequences:

```
P(component k) = pi_k    (same for all i)
pi_k = (1/N) * sum_i posterior[i, k]
```

### With covariates

Mixing proportions depend on covariates:

```
P(component k | x_i) = exp(beta_k' x_i) / sum_j exp(beta_j' x_i)
```

Where `x_i` is the covariate vector for sequence i, `beta_1 = 0` (reference class), and `nnet::multinom()` adds the intercept automatically.

## API

```r
build_mmm(data, k, covariates = NULL, ...)
```

### `covariates` parameter

Same 4 forms as `cluster_data()`:

| Form | Example |
|------|---------|
| Formula | `~ Age + Gender` |
| Character vector | `c("Age", "Gender")` |
| String | `"Age + Gender"` |
| Data frame | `demographics` (all columns used) |
| `NULL` (default) | No covariates, same behavior as current |

Covariate resolution uses `.resolve_covariates()` from `cluster_data.R` (shared infrastructure). Called with the *original* `data` argument (before sequence extraction), so covariate columns are found in `netobject$metadata`, raw data.frame extra columns, or the provided data.frame.

## EM Algorithm Changes

### E-step changes

Currently `log_pi` is a length-k vector broadcast to all sequences:

```r
log_lik <- log_lik + rep(log_pi, each = N)
```

With covariates, `log_pi` becomes an N x k matrix:

```r
# log_pi_mat[i, k] = log P(component k | x_i)
log_lik <- log_lik + log_pi_mat
```

### M-step changes

Only the mixing proportion update changes. Currently (line 69 of mmm.R):

```r
pi_mix <- .colMeans(post, N, n_comp)
log_pi <- log(pi_mix + 1e-300)
```

With covariates, the M-step maximizes:

```
Q(beta) = sum_i sum_k gamma_{ik} * log pi_k(x_i; beta)
```

where `gamma_{ik}` is the full N x k posterior matrix from the E-step.

**Implementation:** Pass the posterior matrix directly as the response to `nnet::multinom()`. Each row of the posterior is treated as multinomial proportions -- this is the standard soft-EM approach. No `weights` argument is used; the posteriors already encode the soft assignments.

```r
# post is N x k matrix, cov_df is N x p data.frame
fit <- nnet::multinom(post ~ Age + Gender, data = cov_df, trace = FALSE,
                       maxit = 200)
log_pi_mat <- log(fitted(fit) + 1e-300)  # N x k
```

`nnet::multinom()` adds the intercept automatically. No manual intercept column needed in `cov_df`.

### What stays unchanged

- Transition matrix M-step (weighted counts via `crossprod(post, counts)`)
- Initial state distribution M-step
- Log-sum-exp numerical stability in E-step
- Convergence criterion (absolute log-likelihood change)
- Screen-and-refine multi-start strategy
- K-means initialization (for transition parameters)

### Beta initialization

First iteration uses existing posterior initialization (K-means or random). Beta coefficients are first estimated in the M-step of iteration 1. No explicit beta initialization needed -- the first M-step multinomial regression uses the initial posteriors as soft labels.

### Inner multinom convergence

The `nnet::multinom()` call inside the M-step uses `maxit = 200`. If it warns about non-convergence, the warning is captured and suppressed during EM iterations. After final convergence, if the last M-step's `multinom` did not converge, a warning is issued to the user.

### Component collapse protection

If any component's total posterior mass drops below 1e-6 * N (effectively empty), the component's mixing is floored to a small positive value and a warning is issued. This prevents degenerate `multinom` fits on near-zero posteriors.

## Implementation Structure

### Modified functions

**`.mmm_em()`** -- core EM loop. Append `cov_mat = NULL` to existing signature (no other changes to parameter names or order):

```r
.mmm_em <- function(counts, init_state, n_comp, max_iter, tol, smooth, K,
                     init_posterior, from_ind, cov_mat = NULL)
```

When `cov_mat` is NULL, behavior is identical to current. When provided, the M-step for mixing calls `.mmm_softmax_mstep()`.

**`build_mmm()`** -- adds `covariates` parameter. Resolves covariates before EM, passes covariate data.frame to `.mmm_em()`. Pre-EM: if covariate rows have NAs, subset `counts`, `init_state`, and `cov_mat` to complete cases and warn about dropped rows.

### New internal functions

**`.mmm_softmax_mstep()`** -- fits `nnet::multinom()` on posteriors and covariates, returns N x k log-probability matrix.

```r
.mmm_softmax_mstep <- function(post, cov_df) {
  # post: N x k posterior matrix (used as matrix response)
  # cov_df: N x p data.frame (no intercept -- nnet adds it)
  # Returns: list(log_pi_mat = N x k, beta = coefficient matrix, fit = multinom object)
  # Handles k=2 case: summary(fit)$coefficients returns vector, wrap to matrix
}
```

**`.mmm_predict_log_pi()`** -- given a fitted multinom object and covariate data.frame, compute per-sequence log mixing proportions. Used to avoid refitting when only prediction is needed.

```r
.mmm_predict_log_pi <- function(fit, cov_df) {
  # Returns: N x k matrix of log(pi_k(x_i))
  probs <- predict(fit, newdata = cov_df, type = "probs")
  # Handle k=2: predict returns vector, wrap to N x 2 matrix
  log(probs + 1e-300)
}
```

### Return object changes

`net_mmm` gains a `$covariates` field when covariates are used:

```r
$covariates$coefficients  # data.frame: cluster, variable, estimate, std_error,
                          #   odds_ratio, ci_lower, ci_upper, z, p, sig
$covariates$profiles      # list: $numeric, $categorical (same as cluster_data)
$covariates$fit           # list: aic, bic, deviance, mcfadden_r2, reference_cluster
$covariates$beta          # raw (k-1) x (p+1) coefficient matrix
$covariates$model         # final nnet::multinom object from last M-step
```

The `$mixing` field stores marginal mixing proportions: `colMeans(fitted(final_multinom))` -- the average of covariate-predicted mixing proportions across all sequences.

### Parameter count

Current:
```
n_params = k * K * (K-1) + k * (K-1) + (k-1)
         = transitions   + initials   + mixing
```

With covariates (p predictors, intercept added by nnet):
```
n_params = k * K * (K-1) + k * (K-1) + (k-1) * (p+1)
         = transitions   + initials   + covariate coefficients
```

This affects BIC, AIC, and ICL computation.

## Pre-EM Covariate Handling

Before the EM loop starts:

1. Resolve covariates via `.resolve_covariates(covariates, data, n)` using the *original* `data` argument.
2. Coerce characters to factors, unorder ordered factors (same as `cluster_data`).
3. Check for NAs: if any, subset `counts`, `init_state`, and `cov_df` to complete cases. Warn with count dropped.
4. Check for constant columns: error if any covariate has < 2 unique values.
5. Pass `cov_df` to `.mmm_em()`.

## Verification Strategy

### 1. Intercept-only equivalence

With no covariates provided vs. the integrated model, the mixing proportion M-step should produce the same result. Test by running both and comparing:

```r
set.seed(42)
mmm_plain <- build_mmm(data, k = 2, n_starts = 3, seed = 42)

# Run integrated model with NO covariates at all -- should be identical
# (this tests that the cov_mat=NULL path is truly unchanged)
```

Additionally, verify that the per-sequence mixing proportions from the integrated model with an uninformative covariate (random noise) are approximately constant across sequences (low variance), matching the behavior of the covariate-free model.

### 2. Coefficient direction recovery on synthetic data

Generate data with known ground truth where higher Age increases probability of component 2:

```r
# True: P(comp 2 | Age) = logistic(1.5 * Age)
# After fitting, verify Age coefficient for component 2 is POSITIVE
# Check both label orderings to handle permutation
```

This is more robust than exact value matching. For N=200, verify sign; for N=1000, verify magnitude within 0.3 of true value.

### 3. Integrated vs post-hoc classification

On synthetic data with known true labels:

```r
# Integrated: build_mmm(data, k = 2, covariates = cov_df)
# Post-hoc: build_mmm(data, k = 2) then post-hoc multinomial regression
# Adjusted Rand index: integrated >= post-hoc when covariates are informative
```

### 4. Log-likelihood monotonicity

Track log-likelihood at each EM iteration. Verify `ll[t+1] >= ll[t] - tol` for all t. This tests that the M-step (including the inner `multinom` optimization) produces valid updates.

### 5. Nested model comparison

```r
# Final log-likelihood: covariate model >= covariate-free model
# (more parameters, strictly nested)
expect_true(mmm_cov$log_likelihood >= mmm_plain$log_likelihood - 1e-4)
```

BIC may favor the simpler model if covariates are noise -- this is correct behavior.

### 6. BIC model selection with noise covariates

Fit with informative covariates vs pure noise covariates. Verify BIC favors the informative model or that the noise model's BIC is worse than the covariate-free model.

## Files to Modify

| File | Changes |
|------|---------|
| `R/mmm.R` | Add `covariates` param to `build_mmm()`, append `cov_mat` to `.mmm_em()`, add `.mmm_softmax_mstep()`, `.mmm_predict_log_pi()`, update parameter count, update return object, update print/summary/plot S3 methods |
| `R/cluster_data.R` | No changes -- `.resolve_covariates()` and `.compute_cluster_profiles()` already available as internal functions |
| `tests/testthat/test-mmm.R` | Add test blocks: intercept-only equivalence, coefficient direction recovery, integrated vs post-hoc, LL monotonicity, nested model comparison, BIC noise test |
| `DESCRIPTION` | No changes (`nnet` already in Suggests) |

## Edge Cases

- Single covariate: works with 1-predictor data.frame
- All categorical covariates: automatic factor conversion, dummy coding
- Covariate NAs: pre-EM subsetting of counts + init_state + cov_df to complete cases, warn
- k=2: `nnet::multinom()` returns vector coefficients, not matrix -- wrap to 1-row matrix in `.mmm_softmax_mstep()`
- Perfect separation by covariate: multinom may warn -- capture during EM, relay after convergence
- Component collapse: if any component posterior mass < 1e-6 * N, floor mixing and warn
- Highly correlated covariates: may cause unstable multinom -- warn if design matrix is near-singular
- Covariates + parallel starts: covariate data.frame passed to each parallel worker (automatic via mclapply fork)
- `compare_mmm()`: covariates pass through via `...`. BIC/AIC comparisons across k values are only valid when the same covariates are used for all k
- Few sequences relative to parameters: warn if N < (k-1) * (p+1) + k * K * (K-1)

## S3 Method Updates

### print.net_mmm

When `$covariates` exists, append:

```
  Covariates:     Age, Gender (integrated, 2 predictors)
```

### summary.net_mmm

Same profile + coefficient table format as `cluster_data()` summary, but with note:

```
Covariate Analysis (integrated into EM -- influences cluster membership)
```

Instead of the post-hoc disclaimer.

### plot.net_mmm

New `type = "covariates"` for the odds ratio dot plot (same as `cluster_data` predictor plot).
