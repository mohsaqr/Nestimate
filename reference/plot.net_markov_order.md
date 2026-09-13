# Plot Method for net_markov_order

Two-panel professional visualization:

- Panel A: log-likelihood, AIC, BIC across tested orders with the
  selected order highlighted (both the permutation-selected order and
  the BIC-minimizing order are marked).

- Panel B: permutation null density per order with the observed \\G^2\\
  as a vertical marker; colored by rejection at `alpha`.

Uses the Okabe-Ito colorblind-safe palette.

## Usage

``` r
# S3 method for class 'net_markov_order'
plot(x, panel = c("both", "ic", "permutation"), combined = TRUE, ...)
```

## Arguments

- x:

  A `net_markov_order` object.

- panel:

  Which panel(s) to render: `"both"`, `"ic"`, or `"permutation"`.
  Default `"both"`.

- combined:

  When `panel = "both"` and `combined = TRUE` (default), the two panels
  are drawn side-by-side. If gridExtra is installed they are arranged
  into a single drawable/saveable gtable (returned); otherwise base
  `grid` viewports draw both panels and a named list of the two ggplots
  is returned invisibly. When `FALSE`, returns that named list (`ic`,
  `permutation`) without drawing. Ignored when `panel != "both"`.

- ...:

  Ignored.

## Value

A ggplot (single panel); for `panel = "both"`, either a `gridExtra`
gtable (when gridExtra is installed) or a named list of two ggplots
(`ic`, `permutation`) drawn side-by-side and returned invisibly.

## Examples

``` r
# \donttest{
# Is one previous state enough to predict the next one?
res <- markov_order_test(as.data.frame(trajectories),
                         max_order = 2, n_perm = 99, seed = 1)
res
#> Markov Order Test  [within-w permutation, n_perm = 99, alpha = 0.050]
#>   131 sequences / 1865 observations / 3 states
#> 
#>   Selected order  BIC: 2   AIC: 2   permutation-LRT: 2
#> 
#>  order   loglik     AIC     BIC df     g2 p_permutation  p_asymptotic
#>      0 -1945.82 3895.63 3906.69 NA     NA            NA            NA
#>      1 -1636.59 3289.19 3333.43  4 618.31          0.01 1.691520e-132
#>      2 -1557.51 3167.02 3310.83 12 154.97          0.01  5.554253e-27
#>  significant
#>           NA
#>         TRUE
#>         TRUE
summary(res)
#>   order    loglik      AIC      BIC df       g2 p_permutation  p_asymptotic
#> 1     0 -1945.815 3895.630 3906.692 NA       NA            NA            NA
#> 2     1 -1636.593 3289.185 3333.434  4 618.3053          0.01 1.691520e-132
#> 3     2 -1557.511 3167.023 3310.829 12 154.9677          0.01  5.554253e-27
#>   significant
#> 1          NA
#> 2        TRUE
#> 3        TRUE
plot(res)

# }
```
