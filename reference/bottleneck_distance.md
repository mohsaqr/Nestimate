# Compute the bottleneck distance between persistence diagrams

`bottleneck_distance()` returns the bottleneck distance between two
persistence diagrams, one value per homology dimension. It takes two
[`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
results, or two data frames of persistence pairs. For the finite pairs,
the distance is \$\$W\_\infty(D_1, D_2) = \inf\_{\gamma} \sup\_{p \in
D_1} \\p - \gamma(p)\\\_\infty,\$\$ where \\\gamma\\ ranges over
bijections \\D_1 \cup \Delta \to D_2 \cup \Delta\\ and \\\Delta =
\\(x,x)\\\\ is the diagonal. Each point is matched to a point of the
other diagram or to its projection onto the diagonal, at cost \\\|d -
b\|/2\\.

The distance is the largest shift needed to match the two diagrams'
features: a small value means every feature of either diagram has a
close counterpart in the other or is short-lived enough to sit near the
diagonal. It is a distance between two observed diagrams and carries no
sampling uncertainty of its own; to compare diagrams on a common grid,
or to average several, use
[`persistence_landscape()`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md).

## Usage

``` r
bottleneck_distance(d1, d2, dimension = NULL, tol = .Machine$double.eps^0.5)
```

## Arguments

- d1, d2:

  The two diagrams:
  [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  results, or data frames with the columns `dimension`, `birth` and
  `death` (other columns are ignored). Any other input, or a data frame
  missing one of those columns, raises an error.

- dimension:

  The homology dimensions to compare, as whole numbers. Default `NULL`:
  every dimension that appears in either diagram. A dimension absent
  from both gives 0.

- tol:

  The tolerance added to each candidate cost when checking whether a
  matching exists, so that floating-point error does not miss an exact
  match. Default `.Machine$double.eps ^ 0.5`.

## Value

A named numeric vector with one distance per dimension, named
`"dim_<k>"` (empty when there is no dimension to compare). `Inf` marks a
structural mismatch: different numbers of essential classes in that
dimension. The distance of a diagram to itself is 0.

## Details

The smallest feasible \\\varepsilon\\ is found by binary search over the
candidate costs (every pairwise \\L\_\infty\\ distance and every
diagonal projection), checking each with a Kuhn bipartite matching.

Essential classes are matched one-to-one within each dimension by their
sorted births. A class is essential when it never dies: death `Inf` in a
Vietoris-Rips filtration, or death 0 with a positive birth in a clique
filtration, as in
[`persistence_landscape()`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md).
A `(0, 0)` pair is a finite point of length zero. When the two diagrams
have different numbers of essential classes in a dimension, the distance
for that dimension is `Inf`.

## References

Edelsbrunner, H. & Harer, J. (2010). *Computational Topology: An
Introduction*. AMS. Section VIII.

## Examples

``` r
mat1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
rownames(mat1) <- colnames(mat1) <- c("A","B","C")
ph1 <- persistent_homology(mat1, n_steps = 5)
bottleneck_distance(ph1, ph1)  # self-distance is 0
#> dim_0 
#>     0 
```
