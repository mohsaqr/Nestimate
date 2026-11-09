# Build a Hypergraph from Group Memberships

`bipartite_groups()` builds a hypergraph from long-format membership
data, in which each row records that a `player` took part in a `group`
(a session, team, project, transaction or any other group context), and
returns a `net_hypergraph`. Each distinct group becomes one hyperedge
spanning the players that appeared in it. The incidence matrix of
players by groups is 0/1, or holds summed weights when `weight` is
given.

## Usage

``` r
bipartite_groups(data, player, group, weight = NULL)
```

## Arguments

- data:

  A data frame with one row per membership (a player in a group),
  holding the `player`, `group` and, when used, `weight` columns. Rows
  with a missing value in any of these columns are dropped without a
  message; when no row is left, an error is raised.

- player:

  The name of the column whose values become the nodes (players,
  participants, actors). Values are read as character and the nodes are
  sorted. A name that is not a single column of `data` raises an error,
  as for `group` and `weight`.

- group:

  The name of the column whose values become the hyperedges (groups,
  sessions, teams). Values are read as character, and the hyperedges
  follow their sorted order (as text, so `"10"` comes before `"2"`).

- weight:

  The name of a numeric column summed per player and group to give a
  weighted incidence matrix. With the default `NULL`, the matrix is 0/1,
  and a player listed more than once in a group counts once. With
  weights, a player belongs to a group when the summed weight is
  positive; a group with no such member is dropped, and a player left
  without a group stays a node.

## Value

A `net_hypergraph`, the class
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
returns. Print it for the numbers of nodes and hyperedges and the
hyperedge sizes; [`summary()`](https://rdrr.io/r/base/summary.html)
prints the counts with the mean and largest hyperedge size, and returns
a data frame of node degrees (`node`, and `degree`, the number of
hyperedges a player belongs to).
[`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md),
[`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md),
[`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md)
and
[`hypergraph_cluster()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
read it. The fields:

- hyperedges:

  An unnamed list with one vector per hyperedge, in the sorted order of
  the groups: the indices of its players in `nodes`, sorted and named by
  player.

- incidence:

  The players-by-groups matrix, one column per hyperedge, with player
  and group names as dimnames: integer 0/1, or the summed weights (a
  player whose sum is zero or negative keeps that sum here but is not a
  member of the hyperedge).

- nodes:

  The player names, sorted.

- n_nodes, n_hyperedges:

  The numbers of players and of hyperedges (after dropping empty
  groups).

- size_distribution:

  A named integer vector with the number of hyperedges of each size
  (`size_1`, `size_2`, ...).

- params:

  A list with `source = "bipartite_groups"`, the column names given as
  `player`, `group` and `weight`, and `n_observations`, the number of
  rows used.

## Details

The bipartite representation preserves the full group structure without
projecting to a pairwise network. A group of three players A, B, C
produces a single 3-hyperedge containing all three, not three pairwise
edges AB, AC, BC. This avoids information loss when group interactions
are the primary unit of analysis (Perc et al. 2013).

Unlike
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
(which derives hyperedges from a network's clique structure),
`bipartite_groups()` takes group memberships directly. The two functions
are complementary:

- `bipartite_groups()` - when group membership is observed (sessions,
  transactions, co-authorships).

- [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md) -
  when only pairwise interactions are observed and group structure must
  be inferred from cliques (triangles by default).

## Note

(experimental) Validated against a hand-computed
[`table()`](https://rdrr.io/r/base/table.html) incidence reference only;
no independent R package exposes the long-format-to-binary-incidence
primitive, because the operation is definitionally
[`table()`](https://rdrr.io/r/base/table.html). The code path is a
direct one-to-one restatement of its definition.

## References

Perc, M., Gomez-Gardenes, J., Szolnoki, A., Floria, L. M., & Moreno, Y.
(2013). Evolutionary dynamics of group interactions on structured
populations: a review. *Journal of the Royal Society Interface* 10(80),
20120997.
[doi:10.1098/rsif.2012.0997](https://doi.org/10.1098/rsif.2012.0997)

## See also

[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
for the clique-based constructor.

## Examples

``` r
df <- data.frame(
  player = c("Alice", "Bob", "Carol", "Alice", "Bob",
             "Dave", "Carol", "Dave", "Eve"),
  session = c("S1", "S1", "S1", "S2", "S2",
              "S3", "S3", "S3", "S3")
)
hg <- bipartite_groups(df, player = "player", group = "session")
print(hg)
#> Hypergraph: 5 nodes, 3 hyperedges
#> Size distribution:
#>   size_2   : 1
#>   size_3   : 2
summary(hg)
#> Hypergraph summary
#>   Nodes:         5
#>   Hyperedges:    3
#>   Mean size:     2.67
#>   Max size:      3
#>    node degree
#> 1 Alice      2
#> 2   Bob      2
#> 3 Carol      2
#> 4  Dave      1
#> 5   Eve      1
```
