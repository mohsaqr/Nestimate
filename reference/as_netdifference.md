# Coerce an inferential comparison to a network difference

`as_netdifference()` turns the result of a network comparison into a
plain `netdifference` object, the signed difference network that
[`subtract_networks()`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md)
returns, so that cograph's difference renderer draws it. It reads a
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
result (`net_bayes`), and returns a `netdifference` unchanged. The drawn
weights are the posterior mean differences, by default only those of the
credible edges, while the full difference, its credible intervals and
its probabilities of direction are kept beside them.

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

  The comparison: a `net_bayes` from
  [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md),
  or a `netdifference`, which is returned as it is. Any other object
  raises an error.

- ...:

  Not used.

- significant_only:

  For a `net_bayes`: whether the drawn weights keep only the credible
  differences (those
  [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  marks in `sig`), with every other edge set to 0. The full difference
  and interval matrices are kept either way. Default `TRUE`.

## Value

A `netdifference` object, which is also a `netobject` (method
`"difference"`); [`print()`](https://rdrr.io/r/base/print.html) lists
its edges and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it. From a `net_bayes`, its `weights` are the posterior mean
differences (all of them, or only the credible ones, by
`significant_only`), and it keeps the full posterior mean difference
(`difference_matrix`), the two posterior mean matrices (`x`, `y`), the
interval bounds (`ci_lower`, `ci_upper`), the Bayesian p-equivalents
(`p_values`), the probabilities of direction (`p_difference`), the
credibility matrix (`sig`), `significant_only`, and
`netdifference_type = "bayes_compare"`.

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
