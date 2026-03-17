# Nestimate Features

## Network Estimation
- **build_network()** — unified entry point for 7 estimation methods
- Transition networks: relative, frequency, co_occurrence
- Association networks: cor, pcor, glasso, ising
- Extensible estimator registry (`register_estimator()`)
- Multilevel decomposition (between/within/both)

## Bootstrap & Inference
- **bootstrap_network()** — universal bootstrap with fast precomputed path
- **permutation_test()** — edge-level comparison (paired/unpaired, 8 corrections)
- **boot_glasso()** — EBICglasso bootstrap (CIs, CS-coefficient, difference tests)

## Higher-Order Networks
- **build_hon()** — BuildHON/BuildHON+ variable-order networks
- **build_honem()** — HONEM embedding (truncated SVD)
- **build_hypa()** — HYPA anomaly detection (hypergeometric null model)
- **mogen_transitions()** — Multi-Order Generative Model (AIC/BIC/LRT)

## Advanced Methods
- **build_gimme()** — GIMME (lavaan SEM, group + individual path search)
- **build_mcml()** — multi-cluster multi-layer networks
- **build_mlvar()** — multilevel VAR for ESM/panel data

## Temporal Networks
- **temporal_network()** — temporal BFS, 21 snapshot metrics, 12 plot types
- **velocity_tna()** — edge velocity/acceleration (regression, GLLA, finite difference)
- Proximity timeline with variable-width lines and highlight mode

## Data Utilities
- Transition frequency matrices from sequences
- Long↔wide format conversion
- Edge/weight/centrality extraction
- Network comparison (correlation, RMSE, MAE, cosine)
