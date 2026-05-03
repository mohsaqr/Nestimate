# Print Method for net_nct

Print Method for net_nct

## Usage

``` r
# S3 method for class 'net_nct'
print(x, ...)
```

## Arguments

- x:

  A `net_nct` object.

- ...:

  Ignored.

## Value

The input object, invisibly.

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
