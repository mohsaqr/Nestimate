# Summary Method for net_nct

Returns a tidy data frame with one row per edge test. The global M
(strength) and S (structure) statistics are attached as attributes.

## Usage

``` r
# S3 method for class 'net_nct'
summary(object, ...)
```

## Arguments

- object:

  A `net_nct` object.

- ...:

  Ignored.

## Value

A data frame with columns `from`, `to`, `diff_observed`, `p_value`,
`significant`. Attributes `m_stat` and `s_stat` each hold a one-row data
frame with `observed` and `p_value`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
x1 <- matrix(rnorm(200 * 5), 200, 5)
x2 <- matrix(rnorm(200 * 5), 200, 5)
colnames(x1) <- colnames(x2) <- paste0("V", 1:5)
res <- nct(x1, x2, iter = 100)
res$M$p_value
res$S$p_value
} # }
```
