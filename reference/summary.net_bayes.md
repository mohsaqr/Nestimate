# Summary method for net_bayes

Summary method for net_bayes

## Usage

``` r
# S3 method for class 'net_bayes'
summary(object, ...)
```

## Arguments

- object:

  A `net_bayes` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with edge-level posterior differences and intervals.

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
b <- bayes_compare(build_network(s1, method = "relative"),
                   build_network(s2, method = "relative"),
                   draws = 500, seed = 1)
summary(b)
#>   from to weight_x weight_y count_x count_y diff effect_size   ci_lower
#> 1    A  B      0.6      0.2       1       0  0.4    1.139724 -0.3669372
#> 2    A  C      0.2      0.6       0       1 -0.4   -1.228861 -0.9348496
#> 3    B  A      0.2      0.6       0       1 -0.4   -1.177673 -0.9260328
#> 4    B  C      0.6      0.2       1       0  0.4    1.148881 -0.3028972
#> 5    C  A      0.6      0.2       1       0  0.4    1.240274 -0.2947713
#> 6    C  B      0.2      0.6       0       1 -0.4   -1.136604 -0.9355968
#>    ci_upper ci_width    pd p_value   sig
#> 1 0.9320787 1.299016 0.866   0.268 FALSE
#> 2 0.2913513 1.226201 0.886   0.228 FALSE
#> 3 0.3991889 1.325222 0.864   0.272 FALSE
#> 4 0.9419299 1.244827 0.848   0.304 FALSE
#> 5 0.9360172 1.230789 0.902   0.196 FALSE
#> 6 0.4434302 1.379027 0.864   0.272 FALSE
```
