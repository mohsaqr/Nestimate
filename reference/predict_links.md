# Predict Missing or Future Links in a Network

Computes link prediction scores for all node pairs using one or more
structural similarity methods. Accepts `netobject`, `mcml`,
`cograph_network`, or a raw weight matrix.

All methods are fully vectorized using matrix operations - no loops.
Supports both weighted and binary adjacency, directed and undirected
networks.

## Usage

``` r
predict_links(
  x,
  methods = c("common_neighbors", "resource_allocation", "adamic_adar", "jaccard",
    "preferential_attachment", "katz"),
  weighted = TRUE,
  top_n = NULL,
  exclude_existing = TRUE,
  include_self = FALSE,
  katz_damping = NULL
)

# S3 method for class 'net_link_prediction'
print(x, ...)

# S3 method for class 'net_link_prediction'
summary(object, ...)
```

## Arguments

- x:

  A `netobject`, `mcml`, `cograph_network`, or numeric square matrix.
  For the [`print()`](https://rdrr.io/r/base/print.html) method: an
  object of class `net_link_prediction`.

- methods:

  Character vector. One or more of: `"common_neighbors"`,
  `"resource_allocation"`, `"adamic_adar"`, `"jaccard"`,
  `"preferential_attachment"`, `"katz"`. Default: all six methods.

- weighted:

  Logical. If `TRUE`, use the weight matrix directly instead of
  binarizing. Default: `TRUE`.

- top_n:

  Integer or NULL. Return only the top N predictions per method.
  Default: `NULL` (all pairs).

- exclude_existing:

  Logical. If `TRUE`, exclude node pairs that already have an edge.
  Default: `TRUE`.

- include_self:

  Logical. If `TRUE`, include self-loop predictions. Default: `FALSE`.

- katz_damping:

  Numeric or NULL. Attenuation factor for Katz index. If NULL,
  auto-computed as `0.9 / spectral_radius(A)`. Default: `NULL`.

- ...:

  In `print.net_link_prediction()` and `summary.net_link_prediction()`:
  Additional arguments (ignored).

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_link_prediction`.

## Value

An object of class `"net_link_prediction"` containing:

- predictions:

  Data frame, one row per (node pair, method), with columns `from`,
  `to`, `method`, `score`, `existing` (was the pair already an edge?)
  and `rank`. Sorted by score (descending) within each method.

- consensus:

  Data frame, one row per node pair, with columns `from`, `to`,
  `avg_rank`, `n_methods` and `consensus_rank`, ordered by `avg_rank`.
  `NULL` when only one method was requested.

- scores:

  Named list of score matrices (one per method).

- adjacency:

  Integer 0/1 adjacency matrix of the input network.

- methods:

  Character vector of methods used.

- nodes:

  Character vector of node names.

- directed:

  Logical.

- weighted:

  Logical.

- n_nodes:

  Integer.

- n_existing:

  Integer. Number of existing edges.

In `print.net_link_prediction()`: The input object, invisibly.

In `summary.net_link_prediction()`: A data frame, one row per method,
with columns `method`, `n_predictions`, `score_mean`, `score_sd`,
`score_max` and `score_min`. A method with no predictions (every
possible link already exists) has `n_predictions = 0` and `NA` scores.

## Details

### Methods

- common_neighbors:

  Number of shared neighbors. For directed graphs, sums shared
  out-neighbors and shared in-neighbors. Vectorized as
  `A %*% t(A) + t(A) %*% A`.

- resource_allocation:

  Zhou et al. (2009). Like common neighbors but weights each shared
  neighbor z by `1/degree(z)`. Penalizes hubs, rewards rare shared
  connections.

- adamic_adar:

  Adamic & Adar (2003). Like resource allocation but weights by
  `1/log(degree(z))`. Less aggressive penalty than RA.

- jaccard:

  Ratio of shared neighbors to total neighbors. For directed graphs,
  computed on combined (out+in) neighbor sets.

- preferential_attachment:

  Product of source out-degree and target in-degree. Captures the
  "rich-get-richer" effect.

- katz:

  Katz (1953). Weighted sum of all paths between nodes, exponentially
  damped by path length. Computed via matrix inversion:
  `(I - beta * A)^{-1} - I`. Captures global structure.

## References

Liben-Nowell, D. & Kleinberg, J. (2007). The link-prediction problem for
social networks. *JASIST*, 58(7), 1019–1031.

Zhou, T., Lu, L. & Zhang, Y.-C. (2009). Network topology and link
prediction. *European Physical Journal B*, 71, 623–630.

Adamic, L. A. & Adar, E. (2003). Friends and neighbors on the Web.
*Social Networks*, 25(3), 211–230.

Katz, L. (1953). A new status index derived from sociometric analysis.
*Psychometrika*, 18(1), 39–43.

Jaccard, P. (1901). Etude comparative de la distribution florale dans
une portion des Alpes et des Jura. *Bulletin de la Societe Vaudoise des
Sciences Naturelles*, 37, 547–579.

Barabasi, A.-L. & Albert, R. (1999). Emergence of scaling in random
networks. *Science*, 286(5439), 509–512.

## See also

[`evaluate_links`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md)
for prediction evaluation,
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for network estimation.

## Examples

``` r
seqs <- data.frame(
  V1 = c("A", "B", "C", "D", "A", "C", "E", "B"),
  V2 = c("B", "C", "D", "E", "C", "E", "A", "D"),
  V3 = c("C", "D", "E", "A", "D", "A", "B", "E")
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net)
print(pred)
#> Link Prediction  [directed | weighted | 5 nodes | 8 existing edges]
#>   Methods: common_neighbors, resource_allocation, adamic_adar, jaccard, preferential_attachment, katz
#> 
#>   Top predicted links (consensus across 6 methods):
#>     1. C -> B  (avg rank: 4.5, agreed: 6/6)
#>     2. B -> A  (avg rank: 4.8, agreed: 6/6)
#>     3. C -> A  (avg rank: 5.0, agreed: 6/6)
#>     4. D -> C  (avg rank: 5.5, agreed: 6/6)
#>     5. A -> D  (avg rank: 5.8, agreed: 6/6)
#>     6. D -> A  (avg rank: 6.2, agreed: 6/6)
#>     7. E -> D  (avg rank: 6.5, agreed: 6/6)
#>     8. E -> C  (avg rank: 7.5, agreed: 6/6)
#>     9. B -> E  (avg rank: 7.7, agreed: 6/6)
#>     10. A -> E  (avg rank: 8.0, agreed: 6/6)
#>     ... and 2 more predictions
summary(pred)
#>                    method n_predictions score_mean   score_sd score_max
#> 1        common_neighbors            12 0.11284722 0.18530385 0.4722222
#> 2     resource_allocation            12 0.03477045 0.06055055 0.1574074
#> 3             adamic_adar            12 0.09626854 0.16577482 0.4298352
#> 4                 jaccard            12 0.14555751 0.11315182 0.3251534
#> 5 preferential_attachment            12 2.25000000 1.13818037 4.0000000
#> 6                    katz            12 1.71583160 0.85510259 3.3856655
#>    score_min
#> 1 0.00000000
#> 2 0.00000000
#> 3 0.00000000
#> 4 0.01408451
#> 5 1.00000000
#> 6 0.64215374
```
