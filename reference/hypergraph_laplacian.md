# Normalised hypergraph Laplacian

`hypergraph_laplacian()` computes the normalised Laplacian of a
hypergraph and returns it as a symmetric node-by-node matrix, the
operator that
[`hypergraph_cluster()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
and
[`hypergraph_transduction()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
are built on. It comes in two forms. `type = "zhou"` is the Laplacian of
Zhou, Huang and Schölkopf (2006), computed on the binary incidence
pattern, so the values in the incidence cells play no part.
`type = "random_walk"` is the Laplacian of Hayashi et al. (2020) for a
random walk with edge-dependent vertex weights: the weighted incidence
cells (for example the summed weights
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
produces with `weight`) set where a walker that enters a hyperedge
lands, and the walk, which is not reversible, is symmetrised through its
stationary distribution (Chung 2005).

Both Laplacians are symmetric positive semi-definite with eigenvalues in
\[0, 2\]. On a binary incidence with unit hyperedge weights the two
forms coincide. The hypergraph must be connected, so the eigenvalue 0
occurs once; further eigenvalues close to 0 indicate groups of nodes
that the walk seldom leaves, which is what
[`hypergraph_cluster()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
partitions by.

## Usage

``` r
hypergraph_laplacian(hg, type = c("zhou", "random_walk"), edge_weights = NULL)
```

## Arguments

- hg:

  The hypergraph: a `net_hypergraph` from
  [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md),
  with at least two nodes and one hyperedge. A hypergraph that is not
  connected has no unique stationary distribution and raises an error of
  class `nestimate_hypergraph_disconnected`; analyse its components
  separately.

- type:

  The Laplacian. `"zhou"` (default) is the Zhou et al. (2006) Laplacian
  on the binary incidence pattern; `"random_walk"` is the Hayashi et
  al. (2020) random-walk Laplacian on the weighted incidence.

- edge_weights:

  Positive hyperedge weights, one per hyperedge in the order of the
  incidence columns. Default `NULL`, which gives unit weights for
  `"zhou"` and, for `"random_walk"`, the heuristic of Hayashi et al.:
  the population standard deviation of each hyperedge's non-zero
  incidence cells plus one, which is 1 on a binary incidence. A vector
  of the wrong length, or with a value that is not finite and positive,
  raises an error.

## Value

A symmetric numeric matrix with one row and one column per node (the
node names as dimnames). Its attributes record the `type`, `pi` (the
stationary distribution of the walk, named by node) and `edge_weights`
(the hyperedge weights used).

## Details

For `"zhou"`, with \\H\\ the binary incidence, \\W\\ the diagonal of
hyperedge weights, \\D_e\\ the hyperedge sizes and \\D_v\\ the node
degrees (the summed weights of the hyperedges a node belongs to), \$\$L
= I - D_v^{-1/2} H W D_e^{-1} H^T D_v^{-1/2},\$\$ and the stationary
distribution is proportional to \\D_v\\. For `"random_walk"`, a walker
at node \\v\\ picks a hyperedge containing it with probability
proportional to its weight, then a node \\u\\ of that hyperedge with
probability proportional to the incidence cell \\\gamma_e(u)\\. The
stationary distribution \\\pi\\ is the leading left eigenvector of this
transition matrix \\P\\, and \\L = I - (\Pi^{1/2} P \Pi^{-1/2} +
\Pi^{-1/2} P^T \Pi^{1/2}) / 2\\. The random-walk Laplacian, its
stationary distribution and its default hyperedge weights follow the
reference implementation in HyperNetX (`laplacians_clustering.py`).

## References

Zhou, D., Huang, J., & Schölkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 19*.

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Chung, F. (2005). Laplacians and the Cheeger inequality for directed
graphs. *Annals of Combinatorics*, 9(1), 1-19.

## Examples

``` r
events <- data.frame(
  person = c("a", "b", "c", "a", "b", "d", "c", "d", "e", "e", "a"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
              "m4", "m4"),
  hours = c(2, 1, 1, 3, 2, 1, 2, 2, 4, 1, 1)
)
hg <- bipartite_groups(events, player = "person", group = "meeting",
                       weight = "hours")
L <- hypergraph_laplacian(hg, type = "random_walk")
range(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
#> [1] -1.803454e-16  1.001279e+00
```
