# Predict Missing or Future Links in a Network

Computes link prediction scores for all node pairs using one or more
structural similarity methods. Accepts `netobject`, `mcml`,
`cograph_network`, or a raw weight matrix.

All methods are fully vectorized using matrix operations — no loops.
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
```

## Arguments

- x:

  A `netobject`, `mcml`, `cograph_network`, or numeric square matrix.

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

## Value

An object of class `"net_link_prediction"` containing:

- predictions:

  Data frame with columns: from, to, method, score, rank. Sorted by
  score (descending) within each method.

- scores:

  Named list of score matrices (one per method).

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

## See also

[`evaluate_links`](https://mohsaqr.github.io/Nestimate/reference/evaluate_links.md)
for prediction evaluation,
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
for network estimation.

## Examples

``` r
seqs <- data.frame(
  V1 = sample(LETTERS[1:5], 50, TRUE),
  V2 = sample(LETTERS[1:5], 50, TRUE),
  V3 = sample(LETTERS[1:5], 50, TRUE)
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net)
print(pred)
#> Link Prediction  [directed | weighted | 5 nodes | 19 existing edges]
#>   Methods: common_neighbors, resource_allocation, adamic_adar, jaccard, preferential_attachment, katz
#> 
#>   Top predicted links (consensus across 6 methods):
#>     1. E -> C  (avg rank: 1.0, agreed: 6/6)
summary(pred)
#>                    method n_predictions  score_mean score_sd   score_max
#> 1        common_neighbors             1  0.39998149       NA  0.39998149
#> 2     resource_allocation             1  0.04079016       NA  0.04079016
#> 3             adamic_adar             1  0.17519419       NA  0.17519419
#> 4                 jaccard             1  0.09383030       NA  0.09383030
#> 5 preferential_attachment             1 16.00000000       NA 16.00000000
#> 6                    katz             1  1.10171908       NA  1.10171908
#>     score_min
#> 1  0.39998149
#> 2  0.04079016
#> 3  0.17519419
#> 4  0.09383030
#> 5 16.00000000
#> 6  1.10171908
```
