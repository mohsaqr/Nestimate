# Summary Method for net_sequence_comparison

Summary Method for net_sequence_comparison

## Usage

``` r
# S3 method for class 'net_sequence_comparison'
summary(object, ...)
```

## Arguments

- object:

  A `net_sequence_comparison` object.

- ...:

  Additional arguments (ignored).

## Value

The patterns data.frame (tidy: one row per k-gram pattern, per group,
with frequency and proportion columns; includes p-values when a
permutation test was run).

## Examples

``` r
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 60, TRUE),
  V2 = sample(LETTERS[1:4], 60, TRUE),
  V3 = sample(LETTERS[1:4], 60, TRUE),
  V4 = sample(LETTERS[1:4], 60, TRUE)
)
grp <- rep(c("X", "Y"), 30)
net <- build_network(seqs, method = "relative")
res <- sequence_compare(net, group = grp, sub = 2:3, test = "chisq")
```
