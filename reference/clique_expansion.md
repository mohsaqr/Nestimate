# Clique expansion of a hypergraph

`clique_expansion()` projects a hypergraph onto an ordinary pairwise
network, its *clique expansion* (also called the "downgrade" of a
hypergraph to a dyadic graph): each hyperedge links every pair of its
members, and the weight between two nodes sums these links over the
hyperedges. It takes a `net_hypergraph` from
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
or
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
and returns an undirected
[netobject](https://pak.dynasite.org/Nestimate/reference/build_network.md).
With binary incidence, the weight `W[i, j]` is the number of hyperedges
containing both `i` and `j`; with weighted incidence, it is the sum over
hyperedges of the product of the two members' incidence values.

The projection keeps which pairs co-occur and how often, and drops the
groups themselves: three hyperedges `{A, B}`, `{B, C}` and `{A, C}` give
the same network as one hyperedge `{A, B, C}`. Read its edges as
pairwise co-occurrence;
[`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md)
and
[`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md)
work on the groups.

## Usage

``` r
clique_expansion(hg, weighted = TRUE)
```

## Arguments

- hg:

  The hypergraph: a `net_hypergraph` from
  [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md).

- weighted:

  Whether to use the hypergraph's incidence values. With `TRUE`
  (default), a weighted hypergraph (from
  [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
  with `weight`) gives a weighted projection; with `FALSE`, the
  incidence is binarised first, so each weight counts the hyperedges the
  pair shares. For an unweighted hypergraph the two agree.

## Value

A `netobject` (class `c("netobject", "cograph_network")`) with method
`"clique_expansion"`, undirected.
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md),
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
read it as any network from
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
It carries no data, so the resampling verbs, such as
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
and sequence clustering refuse it. Its weight matrix is
`W = incidence %*% t(incidence)` with a zero diagonal, and it holds the
standard `netobject` fields (`weights`, `nodes`, `edges`, `n_nodes`,
`n_edges`, `meta`). `edges` has one row per non-zero upper-triangle
cell, with integer `from`/`to` node indices and `weight`. `params`
records `source`, `weighted`, `n_hyperedges` and
`hypergraph_size_distribution`.

## Details

The weight matrix is computed in one BLAS call, `tcrossprod(incidence)`,
in `O(n_nodes^2 * n_hyperedges)` time. With binary incidence (an
unweighted hypergraph, or `weighted = FALSE`), its off-diagonal counts
are those of
[`cooccurrence()`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md)
on the same groups, which also keeps the diagonal.

## Note

Experimental. Validated against `tcrossprod(incidence)` with a zero
diagonal. No external R package exposes clique expansion as a primitive;
the implementation restates the definition in one line.

## References

Tian, H., & Zafarani, R. (2024). Higher-order networks representation
and learning: A survey. *ACM SIGKDD Explorations Newsletter* 26(1),
1-18.

## See also

[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md),
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md),
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

## Examples

``` r
df <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
)
hg  <- bipartite_groups(df, player = "player", group = "session")
net <- clique_expansion(hg)
extract_edges(net, threshold = 1)
#>   from to weight
#> 1    A  B      2
```
