# Q-Analysis

`q_analysis()` describes how the maximal simplices of a simplicial
complex (those not contained in a larger one) hang together, after Atkin
(1974). Two maximal simplices are q-near when they share a face of
dimension at least q (q + 1 shared nodes), and q-connected when a chain
of q-near simplices joins them. The result has two parts:

- **Q-vector**: \\Q_q\\, the number of q-connected components among the
  maximal simplices of dimension at least q, for each q from the largest
  simplex dimension down to 0. Simplices of lower dimension do not exist
  at level q and are not counted.

- **Structure vector**: here the dimension of the largest simplex each
  node belongs to. (Atkin uses "structure vector" for the Q-vector
  itself; the name is kept for compatibility.)

A Q-vector of 1 at every level means the maximal simplices stay in one
piece however many shared nodes are required. Requiring more shared
nodes splits the complex more, but fewer simplices take part at higher
q, so the count can rise or fall with q;
[`print()`](https://rdrr.io/r/base/print.html) reports the lowest level
at which there is more than one. Like the complex itself, the vectors
depend on the threshold or pathways
[`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
was given.

## Usage

``` r
q_analysis(sc)

# S3 method for class 'q_analysis'
print(x, ...)

# S3 method for class 'q_analysis'
plot(x, combined = TRUE, ...)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md).

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `q_analysis` object.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: `TRUE` (default) draws the Q-vector and the structure vector
  side by side, which needs gridExtra; `FALSE` returns them as a named
  list of ggplots (`q_vector`, `structure_vector`).

## Value

A `q_analysis` object. [`print()`](https://rdrr.io/r/base/print.html)
shows the components at each level, the level at which the complex first
splits and the structure vector;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the two
vectors. The object's fields:

- q_vector:

  Named integer vector, one entry per level from the highest down
  (`q_2`, `q_1`, `q_0` for a complex of dimension 2): the number of
  q-connected components at that level.

- max_q:

  The largest simplex dimension, the highest level.

- structure_vector:

  Named integer vector, one entry per node: the dimension of the largest
  simplex containing it.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
drawn grid grob invisibly with `combined = TRUE`, and the named list of
two ggplots invisibly with `combined = FALSE`.

## References

Atkin, R. H. (1974). *Mathematical Structure in Human Affairs*.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
q_analysis(sc)
#> Q-Analysis (max q = 2)
#>   Components: q2:1 q1:1 q0:1
#>   Fully connected at all q levels
#>   Structure: A:2 B:2 C:2
```
