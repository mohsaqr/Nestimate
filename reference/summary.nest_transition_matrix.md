# Summary Method for Transition Matrices

Summary Method for Transition Matrices

## Usage

``` r
# S3 method for class 'nest_transition_matrix'
summary(object, ...)
```

## Arguments

- object:

  A `nest_transition_matrix` returned by
  [`extract_transition_matrix`](https://mohsaqr.github.io/Nestimate/reference/extract_transition_matrix.md).

- ...:

  Additional arguments (ignored).

## Value

A tidy data frame with columns `from`, `to`, `weight`, with one row per
non-zero entry.
