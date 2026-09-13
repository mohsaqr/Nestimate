# Print Method for an mcml Layer

Compact view of one mcml layer (the macro layer or a single
within-cluster network): a header line with the node and non-zero edge
counts and the weight range, the rounded weight matrix, the initial
probabilities as a bar chart, and the dimensions of any attached data –
rather than the raw list contents.

## Usage

``` r
# S3 method for class 'mcml_layer'
print(x, ...)
```

## Arguments

- x:

  An `mcml_layer`, as held in `$macro` and in each element of
  `$clusters` of an `mcml`.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input `mcml_layer`, invisibly.
