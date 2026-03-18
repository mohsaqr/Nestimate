# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

Nestimate is an R package for network estimation from sequence and panel data. Split from Saqrlab v0.3.0 — contains all computation code, no simulation or benchmarking. Those live in Saqrlab, which depends on Nestimate.

## Build & Test Commands

```r
devtools::test()                                           # Run all tests
testthat::test_file("tests/testthat/test-build_network.R") # Run single test file
devtools::document()                                       # Rebuild roxygen2 → man/ + NAMESPACE
devtools::check()                                          # Full R CMD check
devtools::install()                                        # Install locally
```

Tests use testthat edition 3. Test helpers in `tests/testthat/helper-simulate.R` provide standalone `simulate_mtna()` and `simulate_data("mlvar")` — no Saqrlab dependency needed for tests.

## Architecture

See `docs/ARCHITECTURE.md` for the full module map.

### Estimator Registry Pattern

`build_network()` is the **single entry point** for all network estimation. It dispatches to estimators via `.estimator_registry` environment (`R/estimator_registry.R`). Eight built-in estimators registered in `.onLoad()`:

| Name | Function | Type |
|------|----------|------|
| `"relative"` | `.estimator_relative` | Row-normalized transitions (directed) |
| `"frequency"` | `.estimator_frequency` | Raw transition counts (directed) |
| `"co_occurrence"` | `.estimator_co_occurrence` | Co-occurrence (undirected) |
| `"cor"` | `.estimator_cor` | Pearson correlation (undirected) |
| `"pcor"` | `.estimator_pcor` | Partial correlations (undirected) |
| `"glasso"` | `.estimator_glasso` | EBICglasso (undirected) |
| `"ising"` | `.estimator_ising` | L1-penalized logistic (undirected) |
| `"attention"` | `.estimator_attention` | Decay-weighted attention (directed) |

**Method aliases** resolved in `.resolve_method_alias()`:
- `"tna"`, `"transition"` → `"relative"`
- `"ftna"`, `"counts"` → `"frequency"`
- `"cna"` → `"co_occurrence"`
- `"corr"`, `"correlation"` → `"cor"`
- `"partial"` → `"pcor"`
- `"ebicglasso"`, `"regularized"` → `"glasso"`
- `"atna"` → `"attention"`

Custom estimators added via `register_estimator(name, fn, description, directed)`. `estimate_network()` is a thin wrapper around `build_network()`.

Returns dual-class `c("netobject", "cograph_network")` S3 objects. The `netobject` class provides Nestimate-specific methods; the `cograph_network` class enables direct use with cograph functions (`splot()`, `communities()`, etc.). Key fields: `$weights` (matrix), `$nodes` (data.frame with id/label/name/x/y), `$edges` (integer from/to + weight), `$meta`, `$node_groups`, plus Nestimate extras (`$method`, `$params`, `$scaling`, `$threshold`, `$level`, `$n_nodes`, `$n_edges`, `$data`, `$metadata`). Supports multilevel decomposition (`level = "between"/"within"/"both"` → `netobject_ml`).

### S3 Class System

All classes have `print`, `summary`, `plot` methods. Naming convention:

| Class | Source | Returns from |
|-------|--------|-------------|
| `c("netobject", "cograph_network")` / `netobject_ml` | `build_network.R` | `build_network()` |
| `netobject_group` | `build_network.R`, `cluster_data.R`, `mcml.R` | `build_network()` (grouped), `cluster_data()`, `build_mcml()` |
| `net_bootstrap` | `bootstrap_network.R` | `bootstrap_network()` |
| `boot_glasso` | `boot_glasso.R` | `boot_glasso()` |
| `net_permutation` | `permutation_test.R` | `permutation_test()` |
| `net_gimme` | `gimme.R` | `build_gimme()` |
| `mcml_network` | `mcml.R` | `build_mcml()` |
| `mlvar_result` | `mlvar.R` | `mlvar()` |
| `net_hon` | `hon.R` | `build_hon()` |
| `net_honem` | `honem.R` | `build_honem()` |
| `net_hypa` | `hypa.R` | `build_hypa()` |
| `net_mogen` | `mogen.R` | `build_mogen()` |
| `net_mmm` | `mmm.R` | `build_mmm()` |
| `mmm_compare` | `mmm.R` | `compare_mmm()` |
| `net_clustering` | `cluster_data.R` | `cluster_data()` |
| `gvar_result` | `graphical_var.R` | `graphical_var()` |
| `ml_graphical_var_result` | `graphical_var.R` | `ml_graphical_var()` |
| `net_reliability` | `reliability.R` | `reliability_network()` |
| `net_stability` | `centrality_stability.R` | `centrality_stability()` |
| `temporal_network` | `temporal_network.R` | `temporal_network()` |
| `tna_velocity` | `velocity_tna.R` | `velocity_tna()` |

### Window-based TNA

`wtna()` — computes networks from one-hot binary matrices using temporal windowing. Supports transition (directed), co-occurrence (undirected), or both. Uses `crossprod()` for efficient computation. Tumbling/sliding window modes. Per-actor grouping.

`import_onehot()` — converts binary indicator data to wide sequence format with optional windowed aggregation.

### Bootstrap & Inference

- `bootstrap_network()` — universal bootstrap with fast precomputed path for transitions (per-sequence count matrices, single `tabulate()` per iteration, 2.8x speedup)
- `permutation_test()` — edge-level network comparison (paired/unpaired, 8 p-value corrections). Shares precomputed infrastructure with bootstrap
- `boot_glasso()` — specialized EBICglasso bootstrap (edge/centrality CIs, CS-coefficient, difference tests)

### Higher-Order Networks

HON, HONEM, HYPA, MOGen share infrastructure:
- `.HON_SEP` (`\x01`) separator for tuple encoding — avoids collision with state names
- `.hon_encode()` / `.hon_decode()` for tuple ↔ string conversion
- `.hon_parse_input()` for input normalization
- `.mogen_count_kgrams()` shared by HYPA and MOGen
- Unified edge columns: `path`, `from`, `to`, `count`, `probability`
- Arrow notation for node names: `"A -> B -> C"`

### Mixed Markov Models

`build_mmm()` — EM-based mixture of Markov chains (soft assignments, per-component transition matrices as netobjects, BIC/AIC/ICL, optional covariate regression). Returns `net_mmm`. `compare_mmm()` compares BIC/AIC/ICL across k values → `mmm_compare`.

### Clustering & Grouping

`cluster_data()` — clusters sequences (ward/k-means/latent-class), then calls `build_network()` per cluster → `netobject_group`. Underlying clustering results returned as `net_clustering`.

### Graphical VAR

`graphical_var()` / `ml_graphical_var()` — idiographic and multilevel VAR network estimation. Return `gvar_result` / `ml_graphical_var_result`.

### Reliability & Stability

`reliability_network()` — split-half reliability of network edges → `net_reliability`.
`centrality_stability()` — CS-coefficient via case-dropping subsets → `net_stability`.

### Network Utilities

`pathways()` — generic dispatching on `net_hon`, `net_hypa`, `net_mogen`; returns character vector of pathway strings in arrow notation. No dedicated return class.
Network comparison functions in `network_comparison.R` (correlation, RMSE, MAE, cosine between two matrices).

### Temporal Networks

`temporal_network()` — 21 snapshot metrics, temporal BFS (multi-pass for non-DAG), proximity timeline plot.
`velocity_tna()` — edge velocity/acceleration (regression, GLLA, finite difference).

## Key Design Decisions

- MCML uses `build_network()` pipeline, NOT `cograph::cluster_summary()`
- MCML `$between` and `$within` are plain matrices, not tna objects
- MCML between-cluster: does NOT collapse consecutive same-cluster states
- HON unified output: columns `path`, `from`, `to`, `count`, `probability`
- Bootstrap precomputes per-sequence counts — resample via `colSums()` not re-counting
- Proximity timeline: strength-based positioning, interleaved slots, spline interpolation
- Plotting uses `cograph::splot()` for static networks

## Internal Conventions

- **Private functions**: `.name()` prefix (e.g., `.estimator_relative()`, `.count_transitions()`)
- **Input validation**: `stopifnot()` at function entry + `match.arg()` for enum params
- **Format detection**: `format = "auto"` with fallback to "wide" or "long"
- **ID column**: NULL means single sequence; character vector for grouping
- **Fast transitions**: `tabulate()` on integer-encoded pairs (vectorized, no loops)
- **Error handling**: `tryCatch()` with `stop(..., call. = FALSE)` for informative messages

## Known Open Issues (Pre-existing)

- HON/HONEM/HYPA/MOGen classes still use their own `$matrix`/`$nodes` fields — they are NOT `cograph_network` compatible
- `print.mcml` S3 method conflict between Nestimate and cograph
- `test-mmm.R` fails during `R CMD check` due to parallel fork restriction in the check sandbox — not a regression

## Testing Gotchas

- `expect_s3_class()` and `expect_null()` do NOT accept an `info` parameter — use `expect_true(inherits(...), info=)` instead
- tna v1.2.1 lowercases `hclust` method names, so `"ward.D2"` becomes `"ward.d2"` and errors — skip ward.D2 in tna cross-validation tests
- `build_network()` classifies columns as state vs metadata by checking if values are in node names — numeric columns always go to `$metadata`; character columns with letters overlapping node names can be misclassified
- `nnet::multinom` with k=2: `summary()$coefficients` returns a named vector, not a matrix — must detect and wrap into a 1-row matrix
- Formula environment: build `as.formula()` inside the fitting function so variables resolve in the correct environment

## Dependencies

**Imports (4):** ggplot2, glasso, data.table, cluster
**Suggests (source, 7):** tna, glmnet, lavaan, lme4, stringdist, nnet, cograph
**Suggests (test-only, 8):** testthat, igraph, IsingFit, bootnet, gimme, mlVAR, qgraph, reticulate

Nestimate is a computation engine — no plot methods delegate to cograph. Users call `cograph::splot(net)` directly. The `centrality_fn` parameter in `centrality_stability()` and `boot_glasso()` accepts an external centrality function (igraph is not required).
