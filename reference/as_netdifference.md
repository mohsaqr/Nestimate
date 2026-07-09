# Coerce an inferential comparison to a network difference

Coerce an inferential comparison to a network difference

## Usage

``` r
as_netdifference(x, ...)

# S3 method for class 'net_bayes'
as_netdifference(x, significant_only = TRUE, ...)

# S3 method for class 'netdifference'
as_netdifference(x, ...)

# Default S3 method
as_netdifference(x, ...)
```

## Arguments

- x:

  An object with network-difference fields.

- ...:

  Additional arguments passed to methods.

- significant_only:

  Logical. For inferential objects, keep only supported differences in
  the plotted weight matrix while retaining the full difference and
  interval matrices. Default `TRUE`.

## Value

A `netdifference` object suitable for
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
