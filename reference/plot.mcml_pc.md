# Plot an MCML-PC Object

Heatmap of the macro (cluster-level) weights with the package's
diverging palette. For the two-layer network rendering use
[`cograph::plot_mcml()`](https://sonsoles.me/cograph/reference/plot_mcml.html),
which accepts `mcml_pc` objects and draws them undirected.

## Usage

``` r
# S3 method for class 'mcml_pc'
plot(x, digits = 2, ...)
```

## Arguments

- x:

  An `mcml_pc` object.

- digits:

  Number of digits for tile labels (default 2).

- ...:

  Additional arguments (ignored).

## Value

A ggplot object.
