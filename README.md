# Nestimate <img src="man/figures/logo.png" align="right" height="139" alt="" />

> Estimate, validate, and compare networks from sequential and cross-sectional data in R.

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/mohsaqr/Nestimate)
[![CRAN status](https://www.r-pkg.org/badges/version/Nestimate)](https://CRAN.R-project.org/package=Nestimate)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Nestimate is a unified framework for building and validating networks from sequence and panel data. It implements two complementary paradigms through a single `build_network()` interface:

- **Transition Network Analysis (TNA)** models the dynamics of temporal processes as weighted directed networks using stochastic Markov models — capturing how events follow one another.
- **Psychological Network Analysis (PNA)** estimates the conditional dependency structure among variables using regularized partial correlations and graphical models — capturing how variables relate to each other.

Both paradigms share the same validation engine (bootstrap, permutation, centrality stability), clustering methods, and output format.

### What Sets Nestimate Apart

- **Dynamic networks from binary event data.** Most network packages require either sequential event logs or continuous variables. Nestimate also builds dynamic networks directly from binary indicator data (multiple states active/inactive simultaneously) through co-occurrence networks and windowed TNA — producing directed, undirected, or *mixed* networks that capture both contemporaneous co-occurrence and temporal transitions in a single model.

- **Self-contained implementations, minimal dependencies.** Nestimate implements its own EBICglasso estimation, Floyd-Warshall shortest paths, betweenness/closeness centrality, and coordinate descent regularization — all from scratch. The entire package requires only 4 imports (ggplot2, glasso, data.table, cluster). No dependency on igraph, bootnet, qgraph, or graphicalVAR for core functionality.

- **Numerically validated against reference packages.** Every estimator and validation method has been cross-validated against established R packages: transition networks and bootstrap produce byte-identical results to `tna`; permutation tests match to Monte Carlo precision; EBICglasso yields correlations &ge;0.98 with `graphicalVAR`. These equivalence tests use identical synthetic datasets to verify output value by value.

- **Multi-Cluster Multi-Layer (MCML) decomposition.** A unique analysis that clusters sequences, then decomposes dynamics into within-cluster transition patterns and between-cluster flow — revealing both what happens inside each behavioral subgroup and how subgroups connect.

- **Full validation pipeline built in.** Bootstrap, permutation, reliability, centrality stability, and difference tests — all working identically across TNA and PNA paradigms, not bolted on as separate packages.

## Installation

```r
# From GitHub
devtools::install_github("mohsaqr/Nestimate")
```

## Quick Start

```r
library(Nestimate)

# --- Transition network from event-log data ---
data(human_cat)
net <- build_network(human_cat, method = "tna",
                     action = "category", actor = "session_id",
                     time = "timestamp")

# --- Psychological network from cross-sectional data ---
data(srl_strategies)
net_pna <- build_network(srl_strategies, method = "glasso",
                         params = list(gamma = 0.5))

# --- Validate with bootstrap ---
boot <- bootstrap_network(net, iter = 1000)

# --- Plot with cograph ---
# install.packages("cograph", repos = "https://mohsaqr.r-universe.dev")
cograph::splot(net)
cograph::splot(boot)
```

## Network Estimation Methods

All methods are accessed through `build_network()` with a `method` argument:

### Transition Networks (Sequential Data)

| Method | Aliases | Description |
|--------|---------|-------------|
| `"relative"` | `"tna"`, `"transition"` | Row-normalized transition probabilities (directed) |
| `"frequency"` | `"ftna"`, `"counts"` | Raw transition counts (directed) |
| `"co_occurrence"` | `"cna"` | Co-occurrence counts from binary data (undirected) |
| `"attention"` | `"atna"` | Decay-weighted transitions emphasizing recent events (directed) |

### Psychological Networks (Cross-Sectional/Panel Data)

| Method | Aliases | Description |
|--------|---------|-------------|
| `"cor"` | `"corr"`, `"correlation"` | Pearson correlations (undirected) |
| `"pcor"` | `"partial"` | Partial correlations controlling for all other variables (undirected) |
| `"glasso"` | `"ebicglasso"`, `"regularized"` | L1-regularized precision matrix with EBIC selection (undirected, sparse) |
| `"ising"` | — | L1-regularized logistic regression for binary variables (undirected, sparse) |

All PNA estimators are implemented from scratch within Nestimate — including EBICglasso with coordinate descent regularization, partial correlations via precision matrix inversion, and EBIC model selection. This eliminates hard dependencies on external network packages while producing numerically equivalent results (validated against `graphicalVAR` and `bootnet` across multiple synthetic and real datasets).

Custom estimators can be added via `register_estimator()`.

## Statistical Validation

Nestimate provides a full statistical validation toolkit at the edge, node, and network level:

```r
# Split-half reliability
reliability(net)

# Bootstrap confidence intervals and significance
boot <- bootstrap_network(net, iter = 1000)

# Centrality stability (CS-coefficient)
centrality_stability(net)

# Permutation-based group comparison
perm <- permutation_test(net_group1, net_group2)

# Specialized glasso bootstrap (edge CIs, centrality stability, difference tests)
boot_gl <- boot_glasso(net_pna, iter = 1000,
                       centrality = c("strength", "expected_influence"))
```

| Function | Purpose |
|----------|---------|
| `reliability()` | Split-half reliability of edge weights |
| `bootstrap_network()` | Bootstrap CIs, p-values, and significance for each edge |
| `centrality_stability()` | CS-coefficient via case-dropping subsets |
| `permutation_test()` | Edge-level comparison between two networks (paired/unpaired) |
| `boot_glasso()` | Edge inclusion, centrality stability, and difference tests for glasso networks |

## Centrality

Centrality measures quantify how important each node (state or variable) is within the network:

```r
centrality(net)
```

Computes InStrength, OutStrength, and Betweenness for directed networks; Strength for undirected. Node predictability (how well each variable is explained by its neighbors) is available for psychological networks:

```r
predictability(net_pna)
```

## Clustering

Nestimate supports two approaches to discover subgroups with distinct network structures:

### Dissimilarity-Based Clustering

`cluster_data()` computes pairwise sequence distances and partitions into `k` groups. Supports 9 distance metrics (Hamming, Levenshtein, LCS, cosine, Jaccard, and more), 8 clustering methods (PAM, Ward, complete/average/single linkage), and optional temporal weighting:

```r
clust <- cluster_data(net, k = 3, dissimilarity = "hamming", method = "ward.D2")
plot(clust, type = "silhouette")
plot(clust, type = "mds")

# Build per-cluster networks
cluster_nets <- build_network(clust, method = "tna")
```

### Mixed Markov Models

`build_mmm()` fits a mixture of Markov chains via EM, clustering sequences by their transition dynamics rather than sequence similarity. Supports soft assignments, BIC/AIC/ICL model selection, and covariate regression:

```r
mmm <- build_mmm(net, k = 3, covariates = c("project"))
compare_mmm(net, k = 2:6)  # Compare model fit across k values
```

### Multi-Cluster Multi-Layer (MCML)

MCML is a two-stage analysis that first clusters sequences into behavioral subgroups, then decomposes the dynamics into two layers: **within-cluster** networks (the transition patterns internal to each subgroup) and a **between-cluster** network (how the system moves from one subgroup's behavior to another's). This reveals structure that neither clustering alone nor a single aggregate network can capture — for example, that learners in a "strategic" cluster cycle between planning and monitoring internally, while the system-level flow shows transitions from "disengaged" behavior toward "strategic" behavior over time.

```r
mcml <- build_mcml(net, k = 3)
mcml$within    # List of per-cluster transition matrices
mcml$between   # Between-cluster flow matrix
summary(mcml)
```

## Dynamic Networks from Binary Data

Many real-world datasets are binary: at each time point, multiple states are either active (1) or inactive (0) — learning activities, coded behaviors, sensor signals. Nestimate builds dynamic networks directly from this data type:

- **Co-occurrence networks** (`method = "co_occurrence"`) count how often pairs of states are both active simultaneously. Windowed co-occurrence aggregates across temporal neighborhoods rather than exact time points.
- **Window-based TNA** (`wtna()`) applies temporal windowing to binary matrices, computing directed transitions between consecutive windows, undirected co-occurrence within windows, or both.

The `method = "both"` mode produces a **mixed network** — a single model that combines directed edges (state A leads to state B across windows) with undirected edges (states A and B co-occur within the same window). This captures both the temporal sequencing and the contemporaneous structure that neither a purely directed nor a purely undirected network can represent alone:

```r
data(learning_activities)

# Co-occurrence network
net_co <- build_network(learning_activities, method = "cna", actor = "student")

# Mixed network: transitions + co-occurrence in one model
net_mixed <- wtna(learning_activities, actor = "student",
                  method = "both", type = "relative")
```

For conditional dependency estimation on binary variables (controlling for all others), see the **Ising model** (`method = "ising"`) under Psychological Networks above.

## Higher-Order Networks

Capture dependencies beyond first-order (pairwise) transitions:

| Function | Method |
|----------|--------|
| `build_hon()` | Higher-Order Network — variable-length memory dependencies |
| `build_honem()` | Higher-Order Network Embedding |
| `build_hypa()` | Hyper-Path Anomaly detection — identifies statistically anomalous pathways |
| `build_mogen()` | Multi-Order Generative model — optimal Markov order per node |

```r
hon <- build_hon(sequences, k = 2)
pathways(hon)  # Extract pathway strings
```

## Simplicial Complex Analysis

Topological analysis of network structure:

```r
sc <- build_simplicial(net, method = "clique")
betti_numbers(sc)
euler_characteristic(sc)

ph <- persistent_homology(net)   # Filtration-based persistence
qa <- q_analysis(net)            # Atkin's Q-analysis
```

## Data Preparation

Nestimate accepts data in long format (event logs), wide format (sequences as columns), or one-hot encoded binary matrices. `build_network()` auto-detects the format, or you can prepare data explicitly:

```r
# Long → wide conversion
prepared <- prepare_data(event_log, action = "code",
                         actor = "student", time = "timestamp")

# Format conversion utilities
wide_to_long(wide_data)
long_to_wide(long_data, action = "action", actor = "id", time = "time")
action_to_onehot(long_data, action = "action", actor = "id", time = "time")
```

## Bundled Datasets

| Dataset | Description | N |
|---------|-------------|---|
| `human_cat` | Human interactions in AI pair programming (9 categories) | 10,796 events, 429 sessions |
| `human_detailed` | Same interactions at fine-grained code level | 10,796 events |
| `learning_activities` | Binary learning activity indicators | 6,000 obs (200 students x 30 timepoints) |
| `srl_strategies` | Self-regulated learning strategy frequencies | 250 students, 9 strategies |
| `group_regulation_long` | Group regulation sequences with covariates | Long format with Actor, Action, Time |
| `human_ai_edges` | Pre-computed edge list | — |

See `?vibcoding-data` for the full family of human-AI coding datasets at three granularity levels.

## Visualization

Nestimate objects integrate directly with [cograph](https://github.com/sonsoleslp/cograph) — no conversion needed:

```r
# install.packages("cograph", repos = "https://mohsaqr.r-universe.dev")
library(cograph)

splot(net)                # Single network
splot(boot)               # Bootstrap-filtered network
splot(perm)               # Permutation difference heatmap
splot(cluster_nets)       # Per-cluster networks
```

## Documentation

- [Transition Networks](https://mohsaqr.github.io/Nestimate/articles/transition-networks.html) — TNA, FTNA, ATNA, co-occurrence, WTNA, validation, clustering, MMM
- [Co-occurrence & Ising Networks](https://mohsaqr.github.io/Nestimate/articles/co-occurrence-networks.html) — binary data analysis
- [Psychological Networks](https://mohsaqr.github.io/Nestimate/articles/psychological-networks.html) — correlation, partial correlation, glasso, bootstrap
- [Clustering](https://mohsaqr.github.io/Nestimate/articles/clustering.html) — dissimilarity-based clustering, MMM, per-cluster networks
- [Full Reference](https://mohsaqr.github.io/Nestimate/reference/)

## Citation

If you use Nestimate in your research, please cite:

> Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition Network Analysis: A Novel Framework for Modeling, Visualizing, and Identifying the Temporal Patterns of Learners and Learning. *Proceedings of the 15th Learning Analytics and Knowledge Conference*. doi: [10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

> Saqr, M., Beck, E., & López-Pernas, S. (2024). Psychological Networks. In M. Saqr & S. López-Pernas (Eds.), *Learning Analytics Methods and Tutorials* (pp. 513–546). Springer. doi: [10.1007/978-3-031-54464-4_19](https://doi.org/10.1007/978-3-031-54464-4_19)

## License

MIT
