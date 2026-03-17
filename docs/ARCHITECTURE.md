# Nestimate Architecture

## Overview

Nestimate is an R package for network estimation from sequence and panel data. It provides estimation, inference, and temporal analysis — no simulation, no benchmarking, no grid search. Those live in Saqrlab.

## Core Pipeline

```
Input data (sequences, matrices, long-format panel)
  → format conversion (frequencies.R, data_conversion.R)
  → estimation (build_network.R → estimators via registry)
  → inference (bootstrap_network.R, permutation_test.R, boot_glasso.R)
  → analysis/plotting
```

## Modules

### Estimation (`build_network.R`, `estimators.R`, `estimator_registry.R`)

`build_network()` is the single entry point. It dispatches to registered estimators via the registry pattern in `.estimator_registry` environment. Seven built-in estimators: relative, frequency, co_occurrence, cor, pcor, glasso, ising. New estimators added via `register_estimator()` without modifying core code. Estimators are registered in `.onLoad()`.

Returns `netobject` S3 class. Supports multilevel decomposition (`level = "between"/"within"/"both"` → `netobject_ml`).

### Bootstrap & Inference

- `bootstrap_network.R` — universal bootstrap for all estimators. Fast precomputed path for transitions (per-sequence count matrices, single `tabulate()` per iteration). Returns `saqr_bootstrap`.
- `permutation_test.R` — edge-level network comparison (paired/unpaired, 8 p-value corrections). Shares precomputed infrastructure with bootstrap.
- `boot_glasso.R` — specialized bootstrap for EBICglasso (edge/centrality CIs, CS-coefficient, difference tests). Returns `boot_glasso`.

### Advanced Methods

- `gimme.R` — GIMME: lavaan-based uSEM with group + individual path search via modification indices. Returns `saqr_gimme`.
- `mcml.R` — Multi-cluster multi-layer networks. Uses `build_network()` for within/between. Returns `mcml_network`.
- `mlvar.R` — Multilevel VAR for ESM/panel data (temporal + contemporaneous + between-subjects networks).

### Higher-Order Networks

- `hon.R` — BuildHON/BuildHON+ variable-order network construction. Returns `saqr_hon`.
- `honem.R` — HONEM embedding (exponentially-decayed neighborhood matrices, truncated SVD). Returns `saqr_honem`.
- `hypa.R` — HYPA path anomaly detection (hypergeometric null model, IPF). Returns `saqr_hypa`.
- `mogen.R` — Multi-Order Generative Model (De Bruijn graphs, order selection via AIC/BIC/LRT). Returns `saqr_mogen`.

All HON functions use unified edge columns: `path`, `from`, `to`, `count`, `probability`. Node names use arrow notation (`"A -> B -> C"`).

### Temporal Networks

- `temporal_network.R` — dynamic network analysis with temporal BFS (multi-pass convergence for non-DAG edges), 21 snapshot metrics, proximity timeline plot with variable-width lines. Returns `temporal_network`.
- `velocity_tna.R` — edge velocity/acceleration (regression, GLLA, finite difference). Returns `tna_velocity`.

### Data Utilities

- `frequencies.R` — transition counting and format conversion (wide/long, frequency/onehot/edgelist)
- `data_conversion.R` — long↔wide, sequence normalization
- `extraction.R` — edge/weight/centrality extraction from fitted models
- `network_comparison.R` — compare two networks (correlation, RMSE, MAE, cosine)

## Dependencies

**Imports:** tna, igraph, ggplot2, glasso, data.table, dplyr, tidyr
**Suggests:** testthat, lavaan, cograph, glmnet, network, sna

## Key Design Decisions

- `build_network()` is the ONLY entry point for estimation — no alternative paths
- Bootstrap/permutation share precomputed infrastructure for speed
- MCML uses `build_network()` pipeline, NOT `cograph::cluster_summary()`
- MCML `$between` and `$within` are plain matrices, not tna objects
- MCML between-cluster: does NOT collapse consecutive same-cluster states
- All S3 classes have print/summary/plot methods
