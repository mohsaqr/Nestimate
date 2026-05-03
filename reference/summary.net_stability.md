# Summary Method for net_stability

Returns the mean correlation at each drop proportion for each measure.

## Usage

``` r
# S3 method for class 'net_stability'
summary(object, ...)
```

## Arguments

- object:

  A `net_stability` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with columns `measure`, `drop_prop`, `mean_cor`, `sd_cor`,
`prop_above`.

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#> Warning: All centrality measures have zero variance. No stability can be assessed.
summary(cs)
#>       measure drop_prop mean_cor sd_cor prop_above
#> 1  InStrength       0.3      NaN     NA        NaN
#> 2 OutStrength       0.3      NaN     NA        NaN
#> 3 Betweenness       0.3      NaN     NA        NaN
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
net <- build_network(seqs, method = "relative")
stab <- centrality_stability(net, measures = c("InStrength","OutStrength"),
                              iter = 10)
summary(stab)
#>        measure drop_prop   mean_cor     sd_cor prop_above
#> 1   InStrength       0.1  0.9222246 0.08733393  1.0000000
#> 2   InStrength       0.2  0.5991203 0.40939988  0.5000000
#> 3   InStrength       0.3  0.8122252 0.19362303  0.7000000
#> 4   InStrength       0.4  0.4092091 0.78743599  0.6000000
#> 5   InStrength       0.5  0.1002185 0.76460340  0.3000000
#> 6   InStrength       0.6  0.7769158 0.14463237  0.7000000
#> 7   InStrength       0.7 -0.3303545 0.83617633  0.2000000
#> 8   InStrength       0.8  0.2274100 0.74756117  0.3000000
#> 9   InStrength       0.9  0.4906391 0.58239756  0.4444444
#> 10 OutStrength       0.1  0.9778194 0.01855856  1.0000000
#> 11 OutStrength       0.2  0.8194689 0.24134040  0.8000000
#> 12 OutStrength       0.3  0.8599518 0.15567073  0.8000000
#> 13 OutStrength       0.4  0.8882647 0.15233279  0.9000000
#> 14 OutStrength       0.5  0.6882411 0.59774050  0.8000000
#> 15 OutStrength       0.6  0.1968403 0.74154752  0.4000000
#> 16 OutStrength       0.7  0.4776833 0.62724172  0.5000000
#> 17 OutStrength       0.8  0.3197964 0.44307581  0.2000000
#> 18 OutStrength       0.9  0.4301138 0.41200379  0.2222222
# }
```
