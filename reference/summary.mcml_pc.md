# Summarize an MCML-PC Object

Summarize an MCML-PC Object

## Usage

``` r
# S3 method for class 'mcml_pc'
summary(object, ...)
```

## Arguments

- object:

  An `mcml_pc` object.

- ...:

  Additional arguments (ignored).

## Value

Tidy data frame with one row per macro edge (upper triangle), columns
`from`, `to`, `weight`.
