# Q-Analysis

Computes Q-connectivity structure (Atkin 1974). Two maximal simplices
are q-connected if they share a face of dimension \\\geq q\\. Reports:

- **Q-vector**: number of connected components at each q-level

- **Structure vector**: highest simplex dimension per node

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

  A `simplicial_complex` object.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `q_analysis`.

- ...:

  For the S3 methods: further arguments passed to or from other methods.

- combined:

  When `TRUE` (default), the two panels are stitched side-by-side via
  [`gridExtra::arrangeGrob`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html).
  When `FALSE`, returns a named list (`q_vector`, `structure_vector`) of
  ggplots.

## Value

A `q_analysis` object with `$q_vector`, `$structure_vector`, and
`$max_q`.

In `print.q_analysis()`: The input object, invisibly.

In `plot.q_analysis()`: A grid grob (invisibly) when `combined = TRUE`;
a named list of two ggplots when `combined = FALSE`.

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
