# Summary Method for Transition Count Matrices

Summary Method for Transition Count Matrices

## Usage

``` r
# S3 method for class 'nest_transition_counts'
summary(object, ...)
```

## Arguments

- object:

  A `nest_transition_counts` matrix returned by
  [`frequencies()`](https://mohsaqr.github.io/Nestimate/reference/frequencies.md).

- ...:

  Additional arguments (ignored).

## Value

A tidy data frame with columns `from`, `to`, `count`, with one row per
non-zero transition.
