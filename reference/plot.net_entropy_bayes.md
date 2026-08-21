# Plot method for `net_entropy_bayes`

Forest plot of the per-edge entropy contributions: posterior mean and
credible interval, credible edges in Okabe-Ito blue, unstable
(non-credible) edges in grey. The dashed line marks `min_share` of the
posterior-mean entropy rate - the stability criterion.

## Usage

``` r
# S3 method for class 'net_entropy_bayes'
plot(x, top = 25, title = "Bayesian edge entropy contributions", ...)
```

## Arguments

- x:

  A `net_entropy_bayes` object.

- top:

  Integer. Show at most this many edges, by posterior mean contribution
  (default `25`).

- title:

  Character. Plot title.

- ...:

  Ignored.

## Value

A ggplot object.
