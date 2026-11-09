# Build a Simplicial Complex

`build_simplicial()` builds a simplicial complex, the higher-order
structure of a network: besides nodes (0-simplices) and edges
(1-simplices) it holds filled triangles (2-simplices), tetrahedra and
larger groups of nodes. Three constructions are available through
`type`:

- **Clique complex** (`"clique"`): every clique of the network whose
  edges reach `threshold` becomes a simplex, the standard bridge from a
  graph to algebraic topology.

- **Pathway complex** (`"pathway"`): the states on each higher-order
  pathway of a `net_hon`, `net_hypa` or `net_mogen` (or of one built
  from a `tna` model or `netobject`) become a simplex.

- **Vietoris-Rips complex** (`"vr"`, alias `"rips"`): the input is a
  non-negative distance matrix, every group of nodes whose pairwise
  distances are all at most `max_scale` becomes a simplex, and each
  simplex records the scale at which it enters (its largest pairwise
  distance), so
  [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  can read the filtration.

A clique complex is fixed by the thresholded graph, since every clique
is filled in; its Betti numbers count the holes that remain, such as
cycles of edges that no triangles fill. Its shape depends on
`threshold`, so read
[`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md)
at more than one threshold, or follow the features across all of them
with
[`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md).
A pathway complex is built from paths observed in the sequences, so its
simplices say which states occur together on a path.
[`simplicial_degree()`](https://pak.dynasite.org/Nestimate/reference/simplicial_degree.md)
and
[`q_analysis()`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
describe how the nodes take part in the complex.

## Usage

``` r
build_simplicial(
  x,
  type = "clique",
  threshold = 0,
  max_dim = 10L,
  max_pathways = NULL,
  anomaly = c("all", "over", "under"),
  max_scale = NULL,
  ...
)

# S3 method for class 'simplicial_complex'
print(x, ...)

# S3 method for class 'simplicial_complex'
plot(x, combined = TRUE, ...)
```

## Arguments

- x:

  The network: a square weight matrix, a `tna` model, a `netobject` or
  `cograph_network`, or a `net_hon`, `net_hypa` or `net_mogen`. With
  `type = "vr"`, the matrix holds non-negative distances. The clique and
  Vietoris-Rips constructions read the weight matrix of any of these
  networks; with `type = "pathway"`, a bare `cograph_network` or a
  matrix raises an error. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `simplicial_complex` object.

- type:

  The construction: `"clique"` (default), `"pathway"` or `"vr"` (alias
  `"rips"`, stored as `"vr"`).

- threshold:

  For `type = "clique"`: the smallest absolute edge weight an edge needs
  to enter (edges at exactly `threshold` are included). Zero weights are
  never edges. Default `0`. Ignored for `"vr"`, which uses `max_scale`,
  and for `"pathway"`.

- max_dim:

  The highest simplex dimension built (a k-simplex has k + 1 nodes), a
  non-negative whole number. Default `10`. When a larger clique is cut
  at this dimension, the clique or Vietoris-Rips complex is marked
  truncated;
  [`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md)
  then warns with class `nestimate_truncated_complex`, as does
  [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  on a Vietoris-Rips complex. A pathway with more than `max_dim + 1`
  distinct states is cut to its faces of that size, without the
  truncated mark.

- max_pathways:

  For `type = "pathway"`, how many pathways to keep: the most frequent
  higher-order paths for a `net_hon`, the most frequent transitions of
  [`mogen_transitions()`](https://pak.dynasite.org/Nestimate/reference/mogen_transitions.md)
  at the optimal order for a `net_mogen` (the highest order fitted when
  the optimal one is 0), and for a `net_hypa` the anomalous paths ranked
  by their observed / expected ratio (see `anomaly`). Default `NULL`
  (all).

- anomaly:

  For `type = "pathway"` on a `net_hypa`: which anomalous paths become
  simplices: `"all"` (default), `"over"` or `"under"`. Over-represented
  and `"all"` paths are ranked by the largest ratio, under-represented
  ones by the smallest. Possible paths that never occur are left out,
  since a pathway simplex is a set of states observed together on a
  path.

- max_scale:

  For `type = "vr"`: the largest distance at which an edge enters;
  longer edges are left out. Default `NULL`, the largest finite
  distance.

- ...:

  For `type = "pathway"` on a `tna` model or `netobject`: `method`
  (`"hon"`, the default, `"hypa"` or `"mogen"`) chooses the higher-order
  model built from it, and the other arguments go to
  [`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md),
  [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  or
  [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md).
  Any other construction or input does not read `...` and warns with
  class `nestimate_ignored_argument` when it is given. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: `TRUE` (default) draws the four panels (f-vector, Betti
  numbers, simplicial degree per node and degree by dimension) as one 2
  x 2 figure, which needs gridExtra; `FALSE` returns the four ggplots in
  a named list (`f_vector`, `betti`, `degree`, `degree_heatmap`) to
  print, save or arrange separately. A complex without edges
  (dimension 0) has no simplicial degree to draw and raises an error of
  class `nestimate_bad_input`.

## Value

A `simplicial_complex` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the counts, the
density, the f-vector (simplices per dimension), the Euler
characteristic and the non-zero Betti numbers;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
f-vector, the Betti numbers and the simplicial degree of each node;
[`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md),
[`euler_characteristic()`](https://pak.dynasite.org/Nestimate/reference/euler_characteristic.md),
[`simplicial_degree()`](https://pak.dynasite.org/Nestimate/reference/simplicial_degree.md),
[`q_analysis()`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
and
[`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
read it. The object's fields:

- simplices:

  List of integer vectors, one per simplex, each holding the sorted
  indices of the nodes it spans. Every face of every simplex is present,
  as are all 0-simplices (isolated nodes included).

- nodes:

  Character vector of node labels, which `simplices` index.

- n_nodes, n_simplices:

  Counts of nodes and simplices.

- dimension:

  Highest simplex dimension present.

- f_vector:

  Named integer vector `dim_0`, `dim_1`, ...: the number of simplices of
  each dimension.

- density:

  `n_simplices` divided by the number of simplices a complete complex on
  `n_nodes` nodes would have up to `dimension`.

- mean_dim:

  Mean simplex dimension.

- type:

  `"clique"`, `"pathway"` or `"vr"`.

- truncated:

  For the clique and Vietoris-Rips complexes: whether `max_dim` cut a
  larger clique.

A Vietoris-Rips complex also holds `filtration` (the scale at which each
simplex enters, in the order of `simplices`) and `max_scale` (the cap
used), which
[`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
reads.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
drawn grid grob invisibly with `combined = TRUE`, and the named list of
four ggplots invisibly with `combined = FALSE`.

## Details

For the clique and Vietoris-Rips complexes the weight matrix is made
symmetric first (the larger absolute weight of the two directions for a
clique complex, the larger distance for a Vietoris-Rips complex) and its
diagonal is not read. Cliques are found with
[`igraph::cliques()`](https://r.igraph.org/reference/cliques.html) when
igraph is installed, and by Bron-Kerbosch otherwise.

## See also

[`betti_numbers`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md),
[`persistent_homology`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md),
[`simplicial_degree`](https://pak.dynasite.org/Nestimate/reference/simplicial_degree.md),
[`q_analysis`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
print(sc)
#> Clique Complex 
#>   3 nodes, 7 simplices, dimension 2
#>   Density: 100.0%  |  Mean dim: 0.71  |  Euler: 1
#>   f-vector: (f0=3 f1=3 f2=1)
#>   Betti: b0=1
#>   Nodes: A, B, C 
betti_numbers(sc)
#> b0 b1 b2 
#>  1  0  0 

# Vietoris-Rips on a distance matrix:
d <- 1 - mat
diag(d) <- 0
sc_vr <- build_simplicial(d, type = "vr", max_scale = 0.6)
```
