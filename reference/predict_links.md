# Predict Missing or Future Links in a Network

`predict_links()` scores every pair of nodes for how strongly the
network's structure suggests a link between them, with one or more
similarity indices: common neighbours, resource allocation, Adamic-Adar,
Jaccard, preferential attachment and Katz (Liben-Nowell & Kleinberg,
2007). It returns one row per pair and index, ranked within each index,
and a consensus ranking across the indices. It takes a `netobject`, a
`cograph_network`, a `tna` model or a square weight matrix, weighted or
binary, directed or undirected.

A high score says that a pair resembles linked pairs in the network's
structure (shared neighbours, busy endpoints, many connecting paths),
not that the link exists in the population or will form. The indices are
heuristics and often disagree; the consensus rank shows the pairs on
which they agree. How well an index recovers links is measured by
scoring held-out edges with
[`evaluate_links()`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md);
[`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
writes the top predictions as pathway strings.

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

  The network: a `netobject`, a `cograph_network`, a `tna` model, or a
  numeric square weight matrix. A network keeps its own direction; a
  matrix or `tna` model is directed unless it is symmetric. A
  `netobject_group` or an `mcml` raises an error of class
  `nestimate_bad_input`: predict each network separately. For the
  [`print()`](https://rdrr.io/r/base/print.html) method: a
  `predict_links()` result.

- methods:

  The indices to compute, one or more of `"common_neighbors"`,
  `"resource_allocation"`, `"adamic_adar"`, `"jaccard"`,
  `"preferential_attachment"` and `"katz"` (defined under *Details*). By
  default all six.

- weighted:

  Whether common neighbours, resource allocation, Adamic-Adar and
  Jaccard use the edge weights (`TRUE`, default) or the binary adjacency
  (`FALSE`). Preferential attachment, Katz and the degrees that resource
  allocation and Adamic-Adar divide by always use the binary adjacency.

- top_n:

  The number of predictions kept per index, the highest scores. Default
  `NULL` (every pair). The ranks are computed before the cut, so they
  stay the same.

- exclude_existing:

  Whether pairs that already have an edge are left out of the
  predictions. Default `TRUE`.

- include_self:

  Whether self-pairs (a node with itself) are listed for a directed
  network; an undirected network never lists them. All indices score
  self-pairs 0. It also decides whether self-loops count in the number
  of existing edges. Default `FALSE`.

- katz_damping:

  The attenuation factor \\\beta\\ of the Katz index. Default `NULL`,
  which uses \\0.9 / \rho(A)\\ for the spectral radius \\\rho(A)\\ of
  the adjacency (`0.1` when it is 0). A value at or above \\1 /
  \rho(A)\\, where the series diverges, is replaced by \\0.9 / \rho(A)\\
  with a warning.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: not used.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a
  `predict_links()` result.

## Value

A `net_link_prediction` object.
[`print()`](https://rdrr.io/r/base/print.html) lists the top ten
predicted links (by consensus rank with several indices, by score with
one) and returns its input invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per index: `method`, `n_predictions`, and the `score_mean`,
`score_sd`, `score_max` and `score_min` of its predictions; an index
with no predictions (every possible link already exists) has
`n_predictions = 0` and `NA` scores.
[`evaluate_links()`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md)
scores the predictions against known edges, and
[`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
writes them as pathway strings. The object's fields:

- predictions:

  Data frame, one row per node pair and index: `from`, `to`, `method`,
  `score`, `existing` (whether the pair already has an edge) and `rank`,
  sorted by score (highest first) within each index. An undirected
  network lists each pair once.

- consensus:

  Data frame, one row per node pair: `from`, `to`, `avg_rank` (the
  pair's mean rank over the indices that list it), `n_methods` (how many
  do) and `consensus_rank`, ordered by `avg_rank`. `NULL` when only one
  index was computed.

- scores:

  Named list of score matrices, one per index.

- adjacency:

  The 0/1 adjacency matrix of the input network.

- methods:

  The indices computed.

- nodes:

  The node names.

- directed, weighted:

  Whether the network was read as directed, and whether weights were
  used.

- n_nodes:

  The number of nodes.

- n_existing:

  The number of existing edges: non-zero cells of the adjacency matrix,
  each undirected edge counted once, self-loops counted only with
  `include_self = TRUE`.

## Details

### Methods

- common_neighbors:

  The number of shared neighbours. For a directed network, shared
  out-neighbours plus shared in-neighbours, computed as
  `A %*% t(A) + t(A) %*% A`.

- resource_allocation:

  Zhou et al. (2009). Like common neighbours, but each shared neighbour
  z is weighted by `1/degree(z)` (in- plus out-degree for a directed
  network). It penalises hubs and rewards rare shared connections.

- adamic_adar:

  Adamic & Adar (2003). Like resource allocation, but weighted by
  `1/log(degree(z))`, a milder penalty; a neighbour of degree 1 adds
  nothing.

- jaccard:

  Shared neighbours divided by all neighbours of the two nodes. With
  `weighted = TRUE` it is the weighted Jaccard (Ruzicka) index, \\\sum_z
  \min(\|w\_{iz}\|, \|w\_{jz}\|) / \sum_z \max(\|w\_{iz}\|,
  \|w\_{jz}\|)\\, which equals the binary index on 0/1 weights and lies
  between 0 and 1. For a directed network, computed on the combined
  (out + in) neighbourhoods, each neighbour weighted by the larger of
  its two directions.

- preferential_attachment:

  The source's out-degree times the target's in-degree (Barabasi &
  Albert, 1999): the "rich get richer" effect.

- katz:

  Katz (1953). A sum over all paths between the two nodes, damped
  exponentially by path length, computed by matrix inversion as
  `(I - beta * A)^{-1} - I`. It reads the global structure. If the
  inversion fails, the series is truncated after six terms.

Every index is computed with matrix operations on the whole network at
once.

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

Ioffe, S. (2010). Improved consistent sampling, weighted minhash and L1
sketching. *IEEE International Conference on Data Mining*, 246–255.

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
#>     1. C -> B  (avg rank: 4.3, agreed: 6/6)
#>     2. B -> A  (avg rank: 4.8, agreed: 6/6)
#>     3. C -> A  (avg rank: 5.0, agreed: 6/6)
#>     4. D -> C  (avg rank: 5.5, agreed: 6/6)
#>     5. A -> D  (avg rank: 5.8, agreed: 6/6)
#>     6. D -> A  (avg rank: 6.2, agreed: 6/6)
#>     7. E -> D  (avg rank: 6.7, agreed: 6/6)
#>     8. E -> C  (avg rank: 7.5, agreed: 6/6)
#>     9. B -> E  (avg rank: 7.7, agreed: 6/6)
#>     10. A -> E  (avg rank: 7.8, agreed: 6/6)
#>     ... and 2 more predictions
summary(pred)
#>                    method n_predictions score_mean   score_sd score_max
#> 1        common_neighbors            12 0.11284722 0.18530385 0.4722222
#> 2     resource_allocation            12 0.03477045 0.06055055 0.1574074
#> 3             adamic_adar            12 0.09626854 0.16577482 0.4298352
#> 4                 jaccard            12 0.31542021 0.22234475 0.6896552
#> 5 preferential_attachment            12 2.25000000 1.13818037 4.0000000
#> 6                    katz            12 1.71583160 0.85510259 3.3856655
#>    score_min
#> 1 0.00000000
#> 2 0.00000000
#> 3 0.00000000
#> 4 0.06122449
#> 5 1.00000000
#> 6 0.64215374
```
