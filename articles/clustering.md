# Sequence Clustering

## Introduction

Clusters represent subgroups within the data that share similar
patterns. Such patterns may reflect similar temporal dynamics (when we
are analyzing sequence data, for example) or relationships between
variables (as is the case in psychological networks). Units within the
same cluster are more similar to each other, while units in different
clusters differ more substantially. In this vignette, we demonstrate how
to perform clustering on sequence data using `Nestimate`.

## Data

To illustrate clustering, we will use the `human_long` dataset, which
contains 10,796 coded human interactions from 429 human-AI pair
programming sessions across 34 projects. Each row represents a single
interaction event with a timestamp, session identifier, and action
label.

``` r

library(Nestimate)
data("human_long")

# Subsample for vignette speed (CRAN build-time limit). Keep only
# sessions of at least 5 turns -- some metrics are undefined on very
# short sequences and downstream sequence_plot row-ordering uses an
# LCS-distance hclust that needs >= 2 observed turns per row.
set.seed(1)
session_lengths <- table(human_long$session_id)
long_enough     <- names(session_lengths)[session_lengths >= 5]
keep            <- sample(long_enough, 80)
human_sub       <- human_long[human_long$session_id %in% keep, ]

head(human_sub)
#>     message_id   project   session_id  timestamp session_date      code
#> 395       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Specify
#> 396       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Command
#> 397       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Request
#> 398       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Specify
#> 399       2903 Project_7 0605767ae57f 1772229600   2026-02-28 Interrupt
#> 400       2905 Project_7 0605767ae57f 1772229600   2026-02-28   Command
#>           cluster code_order order_in_session
#> 395     Directive          1                1
#> 396     Directive          2                2
#> 397     Directive          3                3
#> 398     Directive          4                4
#> 399 Metacognitive          1                5
#> 400     Directive          1                8
```

We can build a transition network using this dataset using
`build_network`. We need to determine the `actor` (`session_id`), the
`action` (`cluster`), and the `time` (`timestamp`). We will use the
overall network object as the starting point to find subgroups since it
structures the raw data into the appropriate units of analysis to
perform clustering.

``` r

net <- build_network(human_sub,
                     method = "tna",
                     action = "cluster",
                     actor  = "session_id",
                     time   = "timestamp")
#> Metadata aggregated per session: ties resolved by first occurrence in 'code' (23 sessions)
```

## Dissimilarity-based Clustering

Dissimilarity-based clustering groups units of analysis (in our case,
sessions, since that is what we provided as `actor`) by directly
comparing their observed sequences. In our case, each session is
represented by its sequence of actions, and similarity between sessions
is defined using a distance metric that quantifies how different two
sequences are.

To implement this method using `Nestimate`, we can use the
[`build_clusters()`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md)
function, which takes either raw sequence data or a network object such
as the `net` object that we estimated (which also contains the original
sequences in `$data`):

``` r

clust <- build_clusters(net, k = 3)

clust
#> Sequence Clustering [pam]
#>   Sequences: 96  |  Clusters: 3
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.289
#> 
#>   Cluster  N           Mean within-dist  Medoid
#>   1        43 (44.8%)  7.378             15
#>   2        39 (40.6%)  18.088            52
#>   3        14 (14.6%)  42.484            42
```

The default clustering mechanism uses **Hamming distance** (number of
positions where sequences differ) with **PAM** (Partitioning Around
Medoids).

The result contains the cluster `assignments` (which cluster each
session belongs to), the cluster `sizes`, and a `silhouette` score that
reflects the quality of the clustering (higher values indicate better
separation between clusters), among other useful information.

``` r

# Cluster assignments (first 20 sessions)
head(clust$assignments, 20)
#>  [1] 1 1 2 2 2 1 1 1 2 1 2 1 1 2 1 2 1 2 1 2

# Cluster sizes
clust$sizes
#>  1  2  3 
#> 43 39 14

# Silhouette score (clustering quality: higher is better)
clust$silhouette
#> [1] 0.2892913
```

### Visualizing Clusters

The silhouette plot shows how well each sequence fits its assigned
cluster. Values near 1 indicate good fit; values near 0 suggest the
sequence is between clusters; negative values indicate possible
misclassification.

``` r

plot(clust, type = "silhouette")
```

![Silhouette plot showing cluster
quality](clustering_files/figure-html/cluster-plot-1.png) The MDS
(multidimensional scaling) plot projects the distance matrix to 2D,
showing cluster separation.

``` r

plot(clust, type = "mds")
```

![MDS plot showing cluster
separation](clustering_files/figure-html/cluster-mds-1.png)

### Distance Metrics

A distance metric defines how (dis)similarity between sequences is
measured. In other words, it quantifies how different two sequences are
from each other. `Nestimate` currently supports 9 distance metrics for
comparing sequences:

| Metric | Description | Best for |
|----|----|----|
| `hamming` | Positions where sequences differ | Equal-length sequences |
| `lv` | Levenshtein (edit distance) | Variable-length, insertions/deletions |
| `osa` | Optimal string alignment | Edit distance + transpositions |
| `dl` | Damerau-Levenshtein | Full edit + adjacent transpositions |
| `lcs` | Longest common subsequence | Preserving order, ignoring gaps |
| `qgram` | Q-gram frequency difference | Pattern-based similarity |
| `cosine` | Cosine of q-gram vectors | Normalized pattern similarity |
| `jaccard` | Jaccard index of q-grams | Set-based pattern overlap |
| `jw` | Jaro-Winkler | Short strings, typo detection |

Different metrics may produce different clustering results. You need to
choose this based on your research question:

- **Hamming**: compares sequences position by position (best for aligned
  sequences of equal length).
- **Edit distances** (lv, osa, dl): allow insertions and deletions (best
  when sequences may be shifted or vary in length).
- **LCS**: captures shared subsequences (best when overall patterns
  matter more than exact alignment).

We can specify which distance metric we want to use through the
`dissimilarity` argument:

``` r

# Levenshtein distance (allows insertions/deletions)
clust_lv <- build_clusters(net, k = 3, dissimilarity = "lv")
clust_lv$silhouette
#> [1] 0.2284073

# Longest common subsequence
clust_lcs <- build_clusters(net, k = 3, dissimilarity = "lcs")
clust_lcs$silhouette
#> [1] 0.2526863
```

Some distance metrics may take additional arguments. For example, the
Hamming distance accepts **temporal weighting** to emphasize earlier or
later positions:

``` r

# Emphasize earlier positions (higher lambda = faster decay)
clust_weighted <- build_clusters(net, 
                               k = 3,
                               dissimilarity = "hamming",
                               weighted = TRUE,
                               lambda = 0.5)
clust_weighted$silhouette
#> [1] 0.2653452
```

### Clustering Methods

By default, Nestimate uses PAM (Partitioning Around Medoids) to form
clusters, which assigns each sequence to the cluster represented by the
most central sequence (the medoid). Besides PAM, `Nestimate` supports
hierarchical clustering methods, which build clusters step by step by
progressively merging similar units into a tree-like structure (a
dendrogram):

- `ward.D2` (“Ward’s Method, Squared Distances”): Minimizes the increase
  in within-cluster variance using squared distances. Typically produces
  compact, well-separated clusters.
- `ward.D` (“Ward’s Method”): An alternative implementation of Ward’s
  approach using a different distance formulation. Similar behavior, but
  results may vary slightly.
- `complete` (“Complete Linkage”): Defines the distance between clusters
  as the maximum distance between their members. Produces tight, compact
  clusters.
- `average` (“Average Linkage”): Uses the average distance between all
  pairs of points across clusters. Provides a balance between
  compactness and flexibility.
- `single` (“Single Linkage”): Uses the minimum distance between points
  in two clusters. Can capture chain-like structures but may lead to -
  loosely connected clusters.
- `mcquitty` (“McQuitty’s Method” / “WPGMA”): A weighted version of
  average linkage that gives equal weight to clusters regardless of
  size.
- `centroid` (“Centroid Linkage”): Defines cluster distance based on the
  distance between cluster centroids (means). Can produce intuitive
  groupings but may introduce inconsistencies in the hierarchy.

To use any of these methods instead of PAM, we need to provide the
`method` argument to `build_clusters`.

``` r

# Ward's method (minimizes within-cluster variance)
clust_ward <- build_clusters(net, k = 3, method = "ward.D2")
clust_ward$silhouette
#> [1] 0.3132877

# Complete linkage
clust_complete <- build_clusters(net, k = 3, method = "complete")
clust_complete$silhouette
#> [1] 0.55901
```

### Choosing k, dissimilarity, and method

[`cluster_choice()`](https://mohsaqr.github.io/Nestimate/reference/cluster_choice.md)
sweeps any combination of `k`, `dissimilarity`, and `method` in a single
call and returns one row per configuration with silhouette, mean
within-cluster distance, and cluster-size balance. Pass a vector to any
axis to sweep it; `dissimilarity = "all"` and `method = "all"` expand to
every supported option.

``` r

ch <- cluster_choice(net, k = 2:4,
                      method = c("pam", "ward.D2", "complete", "average"))
ch
#> Cluster Choice (sweep: k x method)
#> 
#>  k method   silhouette within_dist sizes    ratio  best    
#>  2 pam      0.459      19.321      [27, 69]  2.556         
#>  3 pam      0.289      16.848      [14, 43]  3.071         
#>  4 pam      0.273      15.941      [9, 43]   4.778         
#>  2 ward.D2  0.600      19.937      [11, 85]  7.727         
#>  3 ward.D2  0.313      16.770      [11, 46]  4.182         
#>  4 ward.D2  0.323      16.028      [3, 46]  15.333         
#>  2 complete 0.666      22.034      [3, 93]  31.000 <-- best
#>  3 complete 0.559      19.194      [3, 85]  28.333         
#>  4 complete 0.555      18.652      [1, 85]  85.000         
#>  2 average  0.666      22.034      [3, 93]  31.000         
#>  3 average  0.642      21.493      [1, 93]  93.000         
#>  4 average  0.555      18.652      [1, 85]  85.000
```

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) accepts an
explicit `type =` so you can pick the chart that suits the sweep shape
(`"lines"`, `"bars"`, `"heatmap"`, `"tradeoff"`, `"facet"`, or the
default `"auto"`).

``` r

plot(ch, type = "lines")
```

![Silhouette across k for each clustering
method](clustering_files/figure-html/choose-k-plot-1.png)

To compare several dissimilarities at a fixed `k`, pass a vector and use
the `bars` plot type. `dissimilarity = "all"` expands to every supported
metric (`hamming`, `osa`, `lv`, `dl`, `lcs`, `qgram`, `cosine`,
`jaccard`, `jw`); some metrics require denser sequences than others, so
for this short-sequence sample we name a curated subset directly.
`abbrev = TRUE` shortens the metric names on the tick labels.

``` r

ch_d <- cluster_choice(net, k = 2,
                        dissimilarity = c("hamming", "lv", "lcs"),
                        method = "ward.D2")
ch_d
#> Cluster Choice (sweep: dissimilarity)
#> 
#>  dissimilarity silhouette within_dist sizes    ratio  best    
#>  hamming       0.600      19.937      [11, 85]  7.727         
#>  lv            0.621      15.793      [11, 85]  7.727 <-- best
#>  lcs           0.620      19.112      [6, 90]  15.000
```

``` r

plot(ch_d, type = "bars", abbrev = TRUE)
```

![Silhouette per dissimilarity at k =
2](clustering_files/figure-html/choose-d-plot-1.png)

A high silhouette can come from one big cluster and a tiny outlier
cluster. The `tradeoff` type plots silhouette against the cluster-size
ratio so you can see both axes at once.

``` r

plot(ch_d, type = "tradeoff", abbrev = TRUE)
```

![Quality vs cluster-size
balance](clustering_files/figure-html/choose-tradeoff-1.png)

Higher silhouette indicates better-defined clusters. The `<-- best`
marker in the printed table is the silhouette-max row;
`summary(ch)$best` returns it programmatically. We select `ward.D2` with
2 clusters here.

``` r

clust <- build_clusters(net, k = 2, method = "ward.D2", seed = 42)
summary(clust)
#> Sequence Clustering Summary
#>   Method:        ward.D2 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.6003 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1   85         16.92829
#>        2   11         43.18182
#>   cluster size mean_within_dist
#> 1       1   85         16.92829
#> 2       2   11         43.18182
```

### Validating the choice with `cluster_diagnostics()`

Once a clustering is fit,
[`cluster_diagnostics()`](https://mohsaqr.github.io/Nestimate/reference/cluster_diagnostics.md)
returns a uniform diagnostic surface that works for both distance-based
fits (`net_clustering`) and model-based fits (`net_mmm`,
`net_mmm_clustering`). The print method shows per-cluster size, mean
within-distance, and silhouette for distance-based fits; AvePP, mixing
share, and per-cluster classification error for MMM fits.

``` r

diag <- cluster_diagnostics(clust)
diag
#> Cluster Diagnostics (distance) [ward.D2 / hamming]
#>   Sequences: 96  |  Clusters: 2
#>   Quality: silhouette = 0.600
#> 
#>   Cluster  N           Mean within-dist  Silhouette
#>   1        85 (88.5%)  16.928            0.661
#>   2        11 (11.5%)  43.182            0.133
```

``` r

plot(diag, type = "silhouette")
```

![Per-observation silhouette by
cluster](clustering_files/figure-html/cluster-diagnostics-plot-1.png)

`as.data.frame(diag)` returns the per-cluster table for downstream use.

## Mixture Markov Models

Instead of clustering sequences based on how similar they are to one
another, we can cluster them together based on their transition
dynamics. Mixture Markov models (MMM) fit separate Markov models, and
sequences are assigned to the cluster whose transition structure best
matches their observed behavior.

To implement MMM, we can use
[`build_mmm()`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md),
which returns a `net_mmm` object with full model details:

``` r

mmm_fit <- build_mmm(net, k = 2)
summary(mmm_fit)
#> Mixed Markov Model
#>   Sequences: 96  |  Clusters: 2  |  States: 3
#>   ICs: LL = -1912.188  |  BIC = 3901.969  |  AIC = 3858.375  |  ICL = 3935.706
#>   Quality: AvePP = 0.852  |  Entropy = 0.475  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        22 (22.9%)  29.1%  0.813
#>   2        74 (77.1%)  70.9%  0.863
#> 
#> --- Cluster 1 (29.1%, n=22) ---
#>               Directive Evaluative Metacognitive
#> Directive         0.672      0.114         0.214
#> Evaluative        0.484      0.244         0.272
#> Metacognitive     0.314      0.195         0.491
#> 
#> --- Cluster 2 (70.9%, n=74) ---
#>               Directive Evaluative Metacognitive
#> Directive         0.614      0.242         0.145
#> Evaluative        0.421      0.473         0.106
#> Metacognitive     0.511      0.333         0.156
#>   component     prior n_assigned mean_posterior     avepp
#> 1         1 0.2913827         22      0.8126370 0.8126370
#> 2         2 0.7086173         74      0.8634383 0.8634383
```

The `net_mmm` object contains posterior probabilities, model fit
statistics (BIC, AIC, ICL), and per-cluster transition matrices in
`$models`.
[`cluster_diagnostics()`](https://mohsaqr.github.io/Nestimate/reference/cluster_diagnostics.md)
works on it too, with MMM-specific columns (mixing share, average
posterior probability, per-cluster classification error):

``` r

diag_mmm <- cluster_diagnostics(mmm_fit)
diag_mmm
#> Cluster Diagnostics (mmm) [k = 2]
#>   Sequences: 96  |  Clusters: 2  |  States: 3
#>   Quality: AvePP = 0.852  |  Entropy = 0.475  |  Class.Err = 0.0%
#>   ICs: LL = -1912.188  |  BIC = 3901.969  |  AIC = 3858.375  |  ICL = 3935.706
#> 
#>   Cluster  N           Mix%   AvePP  Class.Err%
#>   1        22 (22.9%)  29.1%  0.813   0.0%
#>   2        74 (77.1%)  70.9%  0.863   0.0%
```

The `posterior` plot type shows the distribution of each sequence’s max
posterior probability, faceted by cluster. Tight distributions near 1
indicate clean assignment; spread-out distributions indicate fuzzy
clusters.

``` r

plot(diag_mmm, type = "posterior")
```

![Posterior certainty per MMM
cluster](clustering_files/figure-html/mmm-diagnostics-plot-1.png)

## Building Networks per Cluster

Nestimate provides two consistent functions for clustering + network
building that both return `netobject_group`:

| Function | Clustering method | Returns |
|----|----|----|
| [`cluster_network()`](https://mohsaqr.github.io/Nestimate/reference/cluster_network.md) | Distance-based (Hamming, LCS, etc.) | `netobject_group` |
| [`cluster_mmm()`](https://mohsaqr.github.io/Nestimate/reference/cluster_mmm.md) | Model-based (MMM) | `netobject_group` |

### Using `cluster_network()` (distance-based)

``` r

# One-shot: cluster + build networks
grp_dist <- cluster_network(net, k = 2, cluster_by = "ward.D2")
grp_dist
#> Group Networks (2 clusters via ward.D2 / hamming)
#> 
#>   Group      Nodes  Edges  Weights         N
#>   Cluster 1  3      9      [0.143, 0.631]  85 (88.5%)
#>   Cluster 2  3      9      [0.108, 0.627]  11 (11.5%)
```

### Using `cluster_mmm()` (model-based)

``` r

# One-shot: MMM cluster + networks
grp_mmm <- cluster_mmm(net, k = 2)
grp_mmm
#> Group Networks (2 clusters from MMM)
#> 
#>   Group      Nodes  Edges  Weights         N
#>   Cluster 1  3      9      [0.114, 0.672]  22 (22.9%)
#>   Cluster 2  3      9      [0.106, 0.614]  74 (77.1%)
```

Both return a `netobject_group` — a list of `netobject`s with clustering
info accessible via `attr(, "clustering")`:

``` r

# Access cluster assignments
attr(grp_dist, "clustering")$assignments[1:10]
#>  [1] 1 1 1 1 1 1 1 1 1 1
attr(grp_mmm, "clustering")$assignments[1:10]
#>  [1] 2 2 1 2 1 2 2 2 2 2

# Access individual cluster networks
grp_dist[[1]]$weights[1:3, 1:3]
#>               Directive Evaluative Metacognitive
#> Directive     0.6310811  0.1918919     0.1770270
#> Evaluative    0.4129213  0.4438202     0.1432584
#> Metacognitive 0.4142259  0.2510460     0.3347280
```

### Visualizing Clustered Sequences

Both `netobject_group` results can be passed directly to
[`sequence_plot()`](https://mohsaqr.github.io/Nestimate/reference/sequence_plot.md)
for visualization:

``` r

sequence_plot(grp_dist, type = "index", main = "Distance-based clusters")
```

![Sequence plot by
cluster](clustering_files/figure-html/cluster-seqplot-1.png)

``` r

sequence_plot(grp_mmm, type = "index", main = "MMM clusters")
```

![Sequence plot by MMM
cluster](clustering_files/figure-html/cluster-seqplot-mmm-1.png)

### Converting `build_clusters()` to networks

If you used
[`build_clusters()`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md)
directly, you can convert to `netobject_group` via
[`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md):

``` r

# build_clusters() returns clustering only (no networks)
clust <- build_clusters(net, k = 2, method = "ward.D2")

# Convert to netobject_group
cluster_net <- build_network(clust)
```

### Comparing Clusters

We can compare which transition probabilities differ significantly among
clusters using permutation testing:

``` r

comparison <- permutation(grp_dist, iter = 100)
```

## Workflow Summary

The end-to-end path from raw sequence data to validated per-cluster
networks:

``` r

# 1. Build the network from long-format data
net <- build_network(human_long, method = "tna",
                     actor = "session_id",
                     action = "cluster",
                     time   = "timestamp")

# 2. Sweep the parameter space
ch <- cluster_choice(net, k = 2:5,
                     dissimilarity = c("hamming", "lcs", "cosine"),
                     method = c("pam", "ward.D2"))
plot(ch, type = "facet", abbrev = TRUE)

# 3. Pick a configuration and fit
clust <- build_clusters(net, k = 2,
                         dissimilarity = "hamming",
                         method = "ward.D2")

# 4. Validate
diag <- cluster_diagnostics(clust)
diag
plot(diag, type = "silhouette")

# 5. Build per-cluster networks
grp <- build_network(clust)

# 6. Optional: model-based second opinion
mmm <- build_mmm(net, k = 2)
plot(cluster_diagnostics(mmm), type = "posterior")
```
