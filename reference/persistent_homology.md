# Persistent Homology

`persistent_homology()` follows the holes of a network's simplicial
complex as the complex is grown edge by edge: each component, loop or
void is born when it first appears and dies when it is filled in, and
the result lists every such feature with its birth, death and
persistence (how long it lasts). Two filtrations are available:

- `type = "clique"`:

  A weighted clique filtration on similarities: edges enter from the
  strongest to the weakest, and a simplex enters at the smallest
  absolute weight among its edges, so the thresholds run from high to
  low.

- `type = "vr"`:

  A Vietoris-Rips filtration on a non-negative distance matrix: a
  simplex enters at the largest distance among its nodes, so the
  thresholds run from low to high. `max_scale` caps the distances
  considered.

A `simplicial_complex` from `build_simplicial(type = "vr")`, which
stores the scale at which each simplex enters, is also read directly.

Features with long persistence last over a wide range of thresholds;
features that die soon after they are born are commonly read as noise.
Essential features never die: the component that remains at the end, and
any hole that no simplex in the filtration fills.
[`bottleneck_distance()`](https://pak.dynasite.org/Nestimate/reference/bottleneck_distance.md)
compares the diagrams of two networks, and
[`persistence_landscape()`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
turns a diagram into functions that can be averaged.

## Usage

``` r
persistent_homology(
  x,
  n_steps = 20L,
  max_dim = 3L,
  type = "clique",
  max_scale = NULL
)

# S3 method for class 'persistent_homology'
print(x, ...)

# S3 method for class 'persistent_homology'
plot(x, combined = TRUE, ...)
```

## Arguments

- x:

  The network: a square weight matrix (similarities for
  `type = "clique"`, non-negative distances for `type = "vr"`), a `tna`
  model, a `netobject` or `cograph_network` (a `net_hon` among them), or
  a `simplicial_complex` from `build_simplicial(type = "vr")`. The
  complex's stored entry values are reduced as they are, so its own type
  is used and `type` is ignored; a complex without them (a clique or
  pathway complex) raises an error. A complex built with a `max_dim`
  below its largest clique warns with class
  `nestimate_truncated_complex`, since its top-dimension features cannot
  die, when the complex's dimension is at most `max_dim` and the
  truncated part lies within `max_scale`. All-zero weights raise an
  error in the clique filtration, and negative distances in the
  Vietoris-Rips one. So does a Vietoris-Rips input in which no edge can
  enter: every off-diagonal distance missing or infinite, or
  `max_scale = 0` with all distances positive. An all-zero distance
  matrix gives a single component, and distances that all exceed
  `max_scale` give one essential component per node. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `persistent_homology` object.

- n_steps:

  The number of thresholds on which the Betti curve is reported, a
  positive whole number. Default `20`. The persistence diagram does not
  depend on it.

- max_dim:

  The highest feature dimension reported, a non-negative whole number.
  Default `3`. A `simplicial_complex` input supplies only the simplices
  it holds.

- type:

  The filtration: `"clique"` (default), on similarity weights, or `"vr"`
  (Vietoris-Rips), on distances.

- max_scale:

  For `type = "vr"`, and for a Vietoris-Rips `simplicial_complex`: the
  largest distance at which an edge enters; longer edges and the
  simplices that hold them are left out. Default `NULL`, the largest
  finite distance (for a complex, every simplex it holds).

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: `TRUE` (default) draws the Betti curve and the persistence
  diagram side by side, which needs gridExtra; `FALSE` returns them as a
  named list of ggplots (`betti_curve`, `persistence`).

## Value

A `persistent_homology` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the threshold
range, the number of features per dimension and the three longest-lived;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the Betti
curve and the persistence diagram;
[`bottleneck_distance()`](https://pak.dynasite.org/Nestimate/reference/bottleneck_distance.md)
and
[`persistence_landscape()`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
read it. The object's fields:

- betti_curve:

  Data frame with one row per threshold and dimension: `threshold`,
  `dimension` and `betti`, the number of features alive at that
  threshold.

- persistence:

  Data frame with one row per feature, sorted by decreasing persistence:
  `dimension`, `birth`, `death` and `persistence`. Essential features
  are included, with `death = 0` in the clique filtration and
  `death = Inf` in the Vietoris-Rips one.

- thresholds:

  The thresholds of the Betti curve.

- mode:

  `"clique"` or `"vr"`.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
drawn grid grob invisibly with `combined = TRUE`, and the named list of
two ggplots invisibly with `combined = FALSE`. The diagram plots each
feature's birth against its death; essential features are drawn at the
end of the filtration, at `death = 0` in the clique filtration and, for
the Vietoris-Rips one, at the largest finite value in the diagram or on
the threshold grid.

## Details

The persistence pairs come from the full reduction of the boundary
matrix over \\\mathbb{Z}/2\\ (Edelsbrunner, Letscher & Zomorodian,
2000): each feature is paired with the simplex whose entry creates it
(birth) and the simplex whose entry destroys it (death). Pairs born and
killed at the same value are dropped. The diagram is exact; `n_steps`
only sets the grid on which the Betti curve is reported.

In the clique filtration the weight matrix is made symmetric with the
larger absolute weight of the two directions, its diagonal is not read,
and every node enters at the largest weight. An essential feature has
`death = 0` and `persistence = birth`; the Betti curve is read on
`n_steps` thresholds from the largest weight down to 1% of it. In the
Vietoris-Rips filtration the matrix is made symmetric with the larger
distance, every node enters at 0, and an essential feature has
`death = Inf` and `persistence = Inf`; the curve runs from 0 to the
largest entry value. The filtration holds simplices one dimension above
`max_dim`, so a feature of dimension `max_dim` can die (the Ripser
convention).

## References

Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2000). Topological
persistence and simplification. In *Proceedings of the 41st Annual
Symposium on Foundations of Computer Science*, 454-463. Journal version:
*Discrete & Computational Geometry* (2002) **28**, 511-533.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 10)
print(ph)
#> Persistent Homology
#>   10 filtration steps [0.6000 → 0.0060]
#>   Features: b0: 2 (1 persistent) 
#>   Longest-lived:
#>     b0: 0.6000 → 0.0000 (life: 0.6000)
#>     b0: 0.6000 → 0.5000 (life: 0.1000)
```
