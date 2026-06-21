# Summary method for net_bayes_group

Summary method for net_bayes_group

## Usage

``` r
# S3 method for class 'net_bayes_group'
summary(object, ...)
```

## Arguments

- object:

  A `net_bayes_group` object.

- ...:

  Additional arguments (ignored).

## Value

A combined data frame with a `comparison` column.

## Examples

``` r
s <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
                grp = c("X","X","Y","Y"))
nets <- build_network(s, method = "relative", group = "grp")
summary(bayes_compare(nets, draws = 200, seed = 1))
#>   comparison from to  weight_x  weight_y count_x count_y       diff effect_size
#> 1     X vs Y    A  B 0.6000000 0.6000000       1       1  0.0000000   0.0000000
#> 2     X vs Y    B  C 0.6000000 0.3333333       1       0  0.2666667   0.7130862
#> 3     X vs Y    C  A 0.3333333 0.6000000       0       1 -0.2666667  -0.6867313
#>     ci_lower  ci_upper ci_width    pd p_value   sig
#> 1 -0.7377704 0.7185043 1.456275 0.505    0.99 FALSE
#> 2 -0.5104689 0.9081864 1.418655 0.785    0.43 FALSE
#> 3 -0.9225347 0.4859351 1.408470 0.740    0.52 FALSE
```
