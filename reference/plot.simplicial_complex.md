# Plot a Simplicial Complex

Produces a multi-panel summary: f-vector, simplicial degree ranking, and
degree-by-dimension heatmap.

## Usage

``` r
# S3 method for class 'simplicial_complex'
plot(x, ...)
```

## Arguments

- x:

  A `simplicial_complex` object.

- ...:

  Ignored.

## Value

A grid grob (invisibly).

## Examples

``` r
# \donttest{
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
if (requireNamespace("gridExtra", quietly = TRUE)) plot(sc)

# }
```
