# Plot method for `net_transition_entropy`

Bar chart of per-state row entropy with overlaid horizontal lines at the
entropy rate \\h(P)\\ (chain-level summary) and the maximum row entropy
\\\log_b n\\ (uniform branching). Bar widths are proportional to the
stationary probability so the visual area sums to the entropy rate.

## Usage

``` r
# S3 method for class 'net_transition_entropy'
plot(x, title = "Transition Entropy", fill = "#0072B2", ...)
```

## Arguments

- x:

  A `net_transition_entropy` object.

- title:

  Character. Plot title.

- fill:

  Character. Bar fill colour. Default Okabe-Ito blue.

- ...:

  Ignored.

## Value

A ggplot object.
