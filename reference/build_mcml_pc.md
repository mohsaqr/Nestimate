# Multi-Cluster Multi-Level Aggregation for Psychometric Networks

**Experimental.** Aggregates a node-level psychometric network
(correlation, partial correlation, or EBICglasso) into a cluster-level
macro network plus per-cluster within networks — the MCML view that
[`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md)
provides for transition networks, adapted to the statistics of
undirected association networks. The API and the exact aggregation
formulas may change between releases.

Unlike transition counts, partial correlations do not aggregate by
arithmetic: the submatrix of a pcor matrix is *not* the pcor network of
the subsystem (the conditioning set changes), and averaging pcor entries
across blocks is descriptive only. The five `aggregation` methods
therefore have explicitly different statuses:

- `"average"` (descriptive):

  Macro edge A–B = mean of the signed node-level weights between members
  of A and members of B; macro diagonal = mean within-block off-diagonal
  weight. Needs only the weight matrix (works on data-less netobjects).
  Caveat: signed averaging can cancel opposite-sign edges; interpret as
  net average association, not connectivity strength.

- `"composite"` (re-estimated):

  Per-observation cluster scores = (optionally standardized) mean of
  member variables; the chosen `method` is then re-fit on the k
  composite columns, so the macro network is a genuine correlation /
  pcor / EBICglasso network among clusters. Requires raw data.

- `"loadings"` (re-estimated, connectivity-weighted):

  As `"composite"`, but member variables are weighted by their
  within-cluster connection strength in the node-level network (the
  absolute summed weight to the other members of their own cluster,
  normalized to sum to 1 per cluster). Nodes that anchor their cluster
  contribute more to its composite. This is Nestimate's own weighting —
  related in spirit to network loadings (Christensen & Golino 2021) but
  not a reimplementation of any EGA-family estimator. Requires raw data.

- `"rv"` (descriptive, multivariate):

  Macro edge A–B = Escoufier's RV coefficient between the member blocks
  — a matrix correlation in `[0, 1]` computed from the block covariance
  structure. No composites are formed and no estimator is re-fit, so
  nothing is lost to averaging; signs are not represented. Requires raw
  data.

- `"canonical"` (descriptive, multivariate):

  Macro edge A–B = the first canonical correlation between the member
  blocks: the strongest linear relationship any weighting of A's items
  can have with any weighting of B's items — an upper bound on what
  composite methods can recover. Requires raw data.

## Usage

``` r
build_mcml_pc(
  x,
  clusters,
  aggregation = c("average", "composite", "loadings", "rv", "canonical"),
  method = c("pcor", "glasso", "cor"),
  within = c("reestimate", "subnetwork"),
  weighting = c("equal", "strength", "eigen", "closeness", "betweenness",
    "expected_influence", "specificity", "pca", "factor", "item_total"),
  scale = TRUE,
  cor_method = c("pearson", "spearman", "polychoric"),
  signed = TRUE,
  id_col = NULL,
  fa_method = c("ml", "paf", "minres", "cfa"),
  ...
)
```

## Arguments

- x:

  A `netobject` estimated with an undirected association method (`cor`,
  `pcor`, `glasso`; aliases accepted) — its `$data` and `$method` are
  reused — or a numeric data.frame of raw observations (then `method`
  decides the estimator).

- clusters:

  Cluster membership: a named list of node-label vectors (names =
  cluster labels), or a vector of cluster labels named by node. Every
  node must be assigned to exactly one cluster.

- aggregation:

  Character. `"average"`, `"composite"`, `"loadings"`, `"rv"`, or
  `"canonical"` (see Description). Default `"average"`.

- method:

  Character. Network estimator for the re-estimation paths and for
  data.frame input: `"pcor"` (default), `"glasso"`, or `"cor"` — the
  same vocabulary as
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
  Ignored (with the netobject's own method used instead) when `x` is a
  netobject and `method` is not given.

- within:

  Character. `"reestimate"` (default) or `"subnetwork"` (see Details).

- weighting:

  Character. How items are weighted inside their cluster composite
  (`aggregation = "composite"` only):

  `"equal"`

  :   (default) 1/m per item — the scale as scored.

  `"strength"`

  :   Mean absolute connection to the other own-cluster members in the
      node-level network (what `aggregation = "loadings"` selects).

  `"eigen"`

  :   Leading-eigenvector weights of the within-cluster block of the
      node-level network — like strength but giving extra weight to
      items connected to other well-connected items.

  `"pca"`

  :   First principal component of the member items' correlation matrix
      — a data-statistical weighting, blind to the estimated network.

  `"factor"`

  :   Standardized loadings of a one-factor model per cluster — the
      classical latent-variable weighting. The extraction method is
      chosen by `fa_method` and the correlation input respects
      `cor_method` (so polychoric factor analysis of ordinal items is
      one call). Clusters with fewer than 3 items (or non-converging
      fits) fall back to `"pca"` with a warning.

  `"closeness"`, `"betweenness"`

  :   Within-cluster closeness / betweenness centrality of the item in
      the node-level network block (absolute weights). Betweenness can
      be all-zero in densely connected clusters; equal weights are used
      then, with a warning.

  `"expected_influence"`

  :   Mean *signed* connection to the other own-cluster members
      (Robinaugh et al. 2016) — like strength, but opposite-sign
      connections subtract, and negative expected influence marks
      reverse-keyed items.

  `"specificity"`

  :   The misfit margin: own-cluster strength minus the strongest
      cross-cluster strength, floored at 0. Items that belong as much to
      another cluster contribute nothing to their composite — the
      weighting twin of the `misfit` diagnostic.

  `"item_total"`

  :   Corrected item-total correlation: each item against the mean of
      the other (standardized) members — the classical
      scale-construction weighting.

  **Custom weightings:** `weighting` also accepts a *named numeric
  vector* (one entry per item; absolute values are normalized within
  each cluster, signs flip items when `signed = TRUE`) or a *function*
  `function(W_block, data_block, nodes)` returning one numeric weight
  per item, evaluated per cluster.

  Network-based and data-based weightings answer different questions;
  comparing them is informative — divergence means the network's view of
  the cluster differs from its latent-variable view. Schemes with
  inherently non-negative weights (equal, closeness, betweenness,
  specificity) keep the eigenvector-based item signs; sign-carrying
  schemes (eigen, pca, factor, expected_influence, item_total, custom)
  use their own.

- scale:

  Logical. Standardize member variables before forming composites
  (default `TRUE`). Ignored for `"average"`, `"rv"`, and `"canonical"`.

- cor_method:

  Character. Correlation type for estimation from raw data: `"pearson"`
  (default), `"spearman"`, or `"polychoric"` (needs lavaan; see
  Details).

- signed:

  Logical. Flip reverse-keyed items in composites (default `TRUE`; see
  Details).

- id_col:

  Character vector or NULL. Identifier column(s) to drop from data.frame
  input before analysis (e.g., the `rid`/actor columns produced by
  [`convert_sequence_format`](https://saqr.me/Nestimate/reference/convert_sequence_format.md)`(format = "frequency")`,
  whose output otherwise feeds this function directly as per-actor
  behavior profiles). Same convention as the association estimators.

- fa_method:

  Character. Extraction method for `weighting = "factor"`:

  `"ml"`

  :   (default) Maximum likelihood
      ([`stats::factanal`](https://rdrr.io/r/stats/factanal.html) on the
      `cor_method`-consistent correlation matrix).

  `"paf"`

  :   Iterated principal axis factoring (SMC start, communalities
      iterated on the reduced correlation matrix).

  `"minres"`

  :   Minimum residual / unweighted least squares (uniquenesses
      optimized to minimize squared off-diagonal residuals).

  `"cfa"`

  :   One-factor confirmatory model in lavaan (standardized loadings);
      with `cor_method = "polychoric"` the items are declared ordered,
      giving the categorical (DWLS) factor model.

  On well-behaved unidimensional clusters the four agree closely;
  divergence indicates Heywood-prone or non-unidimensional clusters.

- ...:

  Further arguments forwarded directly to the `fa_method` backend,
  exactly as that backend spells them (requires `weighting = "factor"`):
  lavaan arguments for `"cfa"` (`estimator = "WLSMV"`,
  `missing = "fiml"`, `se = "robust"`, ...),
  [`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) arguments
  for `"ml"`, and `max_iter` / `tol` for `"paf"`. Arguments managed
  internally (`model`, `data`, `covmat`, `n.obs`, `factors`) are ignored
  with a warning — the model is always the one-factor model per cluster,
  because the composite needs exactly one weight per item.

## Value

An object of class `"mcml_pc"` containing:

- macro:

  Cluster-level netobject (k x k, undirected).

- clusters:

  Named list of within-cluster netobjects (`NULL` for singleton
  clusters).

- cluster_members:

  Named list of member node labels.

- loadings:

  Tidy item-diagnostic data frame (one row per node): `node`, `cluster`,
  `loading` (signed own-cluster), `weight`, `sign`, `max_cross`,
  `cross_cluster`, `misfit`. `NULL` only when no node-level network is
  available.

- node_network:

  The node-level netobject the aggregation was based on (kept for
  diagnostics and `loading_stability`).

- data:

  The raw data (data.frame) when available, else `NULL`.

- meta:

  List: `aggregation`, `method`, `within`, `scale`, `cor_method`,
  `signed`, `n_nodes`, `n_clusters`, `cluster_sizes`, `n_misfit`,
  `n_flipped`, `directed = FALSE`, `source = "pc"`,
  `experimental = TRUE`.

## Details

**Item diagnostics.** Whenever raw data or a node-level network is
available, every item's connection strength to *every* cluster is
computed. The `$loadings` table reports, per item: its signed
own-cluster loading, its composite weight, its strongest cross-cluster
loading, and a `misfit` flag set when the cross-cluster loading exceeds
the own-cluster loading — evidence the item is assigned to the wrong
cluster. Misfit items trigger a warning; every aggregation silently
inherits a bad membership, so fix the assignment rather than ignoring
the flag.

**Reverse-keyed items.** With `signed = TRUE` (default), items whose
summed within-cluster association is negative are flipped (their
standardized values enter composites with weight sign -1), so a
reverse-keyed item reinforces its cluster composite instead of
cancelling it. Flips are reported in the `sign` column and via a
warning.

**Missing data.** Composites are per-row weighted means over the
*observed* members (weights renormalized per row); rows with no observed
member yield `NA` and are dropped by the estimator with a message.
Node-level estimation applies the estimators' own complete-case
handling.

**Ordinal items.** `cor_method = "polychoric"` (requires the lavaan
package) estimates the node-level and within-cluster networks from
polychoric correlations — appropriate for Likert items. Composites are
continuous sums, so the macro re-estimation uses Pearson correlations of
the composites regardless.

**Within-cluster networks** follow `within`: `"reestimate"` (default)
re-fits the estimator on the member columns alone — the honest
conditional structure of the subsystem; `"subnetwork"` slices the
node-level weight matrix and is descriptive (for pcor/glasso it retains
conditioning on out-of-cluster nodes). Modes without raw data force
`"subnetwork"`. Singleton clusters get no within network (`NULL`).

**Uncertainty.** The composite/loadings macro network is a full
netobject carrying its composite data, so
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)`(fit$macro)`
(edge-weight CIs) and
[`vertex_bootstrap`](https://saqr.me/Nestimate/reference/vertex_bootstrap.md)`(fit$macro)`
(network-level CIs) work directly;
[`vertex_compare`](https://saqr.me/Nestimate/reference/vertex_compare.md)`(fit1$macro, fit2$macro)`
compares two groups.
[`loading_stability`](https://saqr.me/Nestimate/reference/loading_stability.md)
quantifies how stable the composite weights themselves are under case
resampling.

All constituent networks are undirected (`meta$directed = FALSE`), so
renderers that auto-detect directedness (e.g.
[`cograph::plot_mcml()`](https://sonsoles.me/cograph/reference/plot_mcml.html))
draw the result without arrowheads.

## References

Escoufier, Y. (1973). Le traitement des variables vectorielles.
*Biometrics*, 29(4), 751-760.

Hotelling, H. (1936). Relations between two sets of variates.
*Biometrika*, 28(3/4), 321-377.

Robinaugh, D. J., Millner, A. J., & McNally, R. J. (2016). Identifying
highly influential nodes in the complicated grief network. *Journal of
Abnormal Psychology*, 125(6), 747-757.

Christensen, A. P., & Golino, H. (2021). On the equivalency of factor
and network loadings. *Behavior Research Methods*, 53, 1563-1580.

Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
correlation networks. *Psychological Methods*, 23(4), 617-634.

## See also

[`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md) for
transition networks,
[`loading_stability`](https://saqr.me/Nestimate/reference/loading_stability.md)
for composite-weight stability,
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
and
[`vertex_bootstrap`](https://saqr.me/Nestimate/reference/vertex_bootstrap.md)
for uncertainty on the macro network.

## Examples

``` r
set.seed(1)
n <- 200
sigma <- matrix(0.15, 6, 6)
sigma[1:3, 1:3] <- 0.5
sigma[4:6, 4:6] <- 0.5
diag(sigma) <- 1
z <- matrix(rnorm(n * 6), n, 6) %*% chol(sigma)
df <- as.data.frame(z)
names(df) <- c("a1", "a2", "a3", "b1", "b2", "b3")
clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))

fit <- build_mcml_pc(df, clusters, aggregation = "composite",
                     method = "cor")
fit$macro$weights
#>          A        B
#> A 0.000000 0.235001
#> B 0.235001 0.000000
fit$loadings
#>   node cluster   loading    weight sign max_cross cross_cluster misfit
#> 1   a1       A 0.4738902 0.3333333    1 0.1671718             B  FALSE
#> 2   a2       A 0.4537919 0.3333333    1 0.1602877             B  FALSE
#> 3   a3       A 0.4736153 0.3333333    1 0.1400550             B  FALSE
#> 4   b1       B 0.5368416 0.3333333    1 0.1618540             A  FALSE
#> 5   b2       B 0.5156674 0.3333333    1 0.1811768             A  FALSE
#> 6   b3       B 0.5168036 0.3333333    1 0.1244837             A  FALSE
```
