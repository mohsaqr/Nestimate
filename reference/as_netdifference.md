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

## Examples

``` r
s1 <- data.frame(V1 = c("A", "B", "C"), V2 = c("B", "C", "A"))
s2 <- data.frame(V1 = c("A", "C", "B"), V2 = c("C", "B", "A"))
b <- bayes_compare(build_network(s1, method = "relative"),
                   build_network(s2, method = "relative"),
                   draws = 500, seed = 1)
as_netdifference(b, significant_only = FALSE)
#> Network difference (x - y): 3 nodes, 6 differing edges
#> Plot: cograph::splot(d) or cograph::plot_difference(d)
#> 
#>  from to   x   y difference
#>     B  A 0.2 0.6       -0.4
#>     C  A 0.6 0.2        0.4
#>     A  B 0.6 0.2        0.4
#>     C  B 0.2 0.6       -0.4
#>     A  C 0.2 0.6       -0.4
#>     B  C 0.6 0.2        0.4
#>     A  A 0.2 0.2        0.0
#>     B  B 0.2 0.2        0.0
#>     C  C 0.2 0.2        0.0
```
