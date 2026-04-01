# Plot Method for mcml

Plots an MCML network. When cograph is available, delegates to
[`cograph::plot_mcml()`](http://sonsoles.me/cograph/reference/plot_mcml.md)
which renders a two-layer visualization (macro summary on top,
within-cluster detail on bottom). Otherwise, converts to a
`netobject_group` and plots each layer as a separate panel.

## Usage

``` r
# S3 method for class 'mcml'
plot(x, ...)
```

## Arguments

- x:

  An `mcml` object.

- ...:

  Additional arguments passed to
  [`cograph::plot_mcml()`](http://sonsoles.me/cograph/reference/plot_mcml.md)
  (e.g., `colors`, `edge_labels`, `mode`).

## Value

The input object, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
seqs <- data.frame(
  T1 = sample(LETTERS[1:6], 30, TRUE),
  T2 = sample(LETTERS[1:6], 30, TRUE),
  T3 = sample(LETTERS[1:6], 30, TRUE)
)
clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
cs <- build_mcml(seqs, clusters)
plot(cs)
} # }
```
