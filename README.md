# Nestimate

Network estimation for dynamic, probabilistic, and multi-layer networks in R.

## Installation

```r
# From GitHub
devtools::install_github("mohsaqr/Nestimate")
```

## Overview

Nestimate provides a unified framework for building and validating networks from sequence and panel data. It supports two complementary paradigms:

**Dynamic Network Analysis (DNA)** — estimate networks from sequential behavioral data (event logs, coded interactions, learning traces). Methods: transition (relative frequency), frequency counts, co-occurrence.

**Probabilistic Network Analysis (PNA)** — estimate networks from cross-sectional or repeated-measures data using statistical associations. Methods: correlation, partial correlation, graphical lasso (EBIC), Ising model.

Both paradigms share the same `build_network()` interface, validation engine, and output format.

## Core Functions

### Network Estimation

| Function | Description |
|---|---|
| `build_network()` | Universal network builder (DNA and PNA methods) |
| `prepare_data()` | Event log to sequence format conversion |
| `boot_glasso()` | Bootstrapped graphical lasso with edge inclusion |

### Validation

| Function | Description |
|---|---|
| `bootstrap_network()` | Bootstrap confidence intervals and significance |
| `permutation_test()` | Permutation-based edge significance |
| `centrality_stability()` | Centrality stability analysis (CS coefficient) |
| `reliability()` | Split-half and test-retest reliability |

### Higher-Order Networks

| Function | Description |
|---|---|
| `build_hon()` | Higher-Order Network (variable-length memory) |
| `build_hypa()` | Hyper-Path Anomaly detection |
| `build_mogen()` | Multi-Order Generative Model |
| `build_honem()` | Higher-Order Network Embedding |

### Simplicial Analysis

| Function | Description |
|---|---|
| `build_simplicial()` | Simplicial complex construction (clique/pathway/VR) |
| `betti_numbers()` | Topological invariants |
| `persistent_homology()` | Filtration-based feature persistence |
| `q_analysis()` | Q-connectivity and structure vectors |

### Multi-Layer and Advanced

| Function | Description |
|---|---|
| `build_mcml()` | Multi-Cluster Multi-Layer networks |
| `mlvar()` | Multilevel vector autoregression |
| `graphical_var()` | Graphical VAR (temporal + contemporaneous) |
| `wtna()` | Windowed temporal network analysis |

### Data Utilities

| Function | Description |
|---|---|
| `prepare_data()` | Parse event logs with auto-detected timestamps and sessions |
| `cluster_sequences()` | Cluster sequences by similarity |
| `frequencies()` | State and transition frequency tables |
| `path_counts()` | Higher-order path enumeration |

## Quick Example

```r
library(Nestimate)

# DNA: Build a transition network from sequence data
sequences <- data.frame(
  student = rep(c("s1", "s2", "s3"), each = 10),
  action = sample(c("Read", "Write", "Discuss"), 30, replace = TRUE)
)
net <- build_network(sequences, method = "relative",
                     actor = "student", action = "action")

# PNA: Build a partial correlation network
net_pna <- build_network(iris[, 1:4], method = "pcor")

# Validate with bootstrap
boot <- bootstrap_network(net, iter = 1000)

# Plot with cograph
cograph::splot(net)
cograph::splot(boot)
```

## Plotting

Nestimate objects integrate directly with [cograph](https://github.com/sonsoleslp/cograph). `splot()` dispatches on `netobject`, `net_bootstrap`, `net_permutation`, `boot_glasso`, `netobject_group`, and `netobject_ml` — no conversion needed.

## License

MIT
