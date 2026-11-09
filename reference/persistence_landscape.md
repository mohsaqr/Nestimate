# Persistence Landscape

`persistence_landscape()` computes the persistence landscape (Bubenik
2015) of one homology dimension of a persistence diagram, evaluated on a
grid. Each (birth, death) pair contributes a tent function
\$\$\Lambda\_{(b,d)}(t) = \max(0, \min(t - b, d - t)),\$\$ and the
\\k\\-th landscape function \\\lambda^{(k)}(t)\\ is the \\k\\-th largest
of \\\\\Lambda\_{(b_i,d_i)}(t)\\\_i\\ at each \\t\\. It takes a
[`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
result or a data frame of persistence pairs.

At each filtration value \\t\\, a tent is the distance from \\t\\ to the
nearer end of a feature's interval (0 outside it), so the first
landscape function follows the feature that \\t\\ lies deepest inside,
the second the feature with the second largest tent, and so on. This
need not be the most persistent feature alive: a long-lived feature near
its birth or death can sit below a shorter one near its midpoint.
Landscapes are stable: the sup-norm distance between two landscapes is
at most the
[`bottleneck_distance()`](https://pak.dynasite.org/Nestimate/reference/bottleneck_distance.md)
between their diagrams. They also live in a Banach space, where they can
be averaged (Bubenik 2015).

Essential classes are left out, because the tent of an infinitely lived
class is undefined on a finite grid: pairs with `death = Inf` (VR mode)
and pairs with `death = 0` but `birth > 0` (the clique-mode encoding of
an essential class) are dropped before the landscape is built. When no
finite pair remains in the requested dimension, every landscape function
is zero on the grid.

## Usage

``` r
persistence_landscape(ph, k_max = 5L, dimension = 1L, t_grid = NULL)

# S3 method for class 'persistence_landscape'
print(x, ...)

# S3 method for class 'persistence_landscape'
plot(x, ...)
```

## Arguments

- ph:

  The diagram: a
  [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  result, or a data frame with the columns `dimension`, `birth` and
  `death`.

- k_max:

  The number of landscape functions to compute, a positive whole number.
  Functions beyond the number of pairs are zero. Default `5`.

- dimension:

  The homology dimension, a single non-negative whole number. Default
  `1`.

- t_grid:

  The filtration values to evaluate the landscape at. Default `NULL`:
  200 evenly spaced points from the smallest birth to the largest death
  of the remaining pairs (each pair is ordered so that birth comes
  first, since clique filtrations descend), or from 0 to 1 when no pair
  remains. Pass the same grid to landscapes you want to compare.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `persistence_landscape` object.

- ...:

  Not used.

## Value

A `persistence_landscape` object. Printing it shows the dimension,
`k_max`, the grid range and the sup norm (largest value) of each
landscape function;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
functions as lines over the grid. The object holds:

- landscape:

  A data frame with one row per landscape function and grid point, and
  the columns `k`, `t` and `value`.

- dimension:

  The homology dimension.

- k_max:

  The number of landscape functions.

- t_grid:

  The evaluation grid.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) a `ggplot`
object with one line per landscape function.

## References

Bubenik, P. (2015). Statistical topological data analysis using
persistence landscapes. *Journal of Machine Learning Research* **16**,
77-102.

## Examples

``` r
mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
rownames(mat) <- colnames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 5)
pl <- persistence_landscape(ph, k_max = 3, dimension = 0)
```
