# Summary Method for Initial Probability Vectors

Summary Method for Initial Probability Vectors

## Usage

``` r
# S3 method for class 'nest_initial_probs'
summary(object, ...)
```

## Arguments

- object:

  A `nest_initial_probs` named numeric vector returned by
  [`extract_initial_probs`](https://mohsaqr.github.io/Nestimate/reference/extract_initial_probs.md).

- ...:

  Additional arguments (ignored).

## Value

A tidy data frame with columns `state` and `prob`, sorted by decreasing
probability.
