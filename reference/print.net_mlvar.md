# Print method for net_mlvar

Print method for net_mlvar

## Usage

``` r
# S3 method for class 'net_mlvar'
print(x, ...)
```

## Arguments

- x:

  A `net_mlvar` object returned by
  [`build_mlvar()`](https://mohsaqr.github.io/Nestimate/reference/build_mlvar.md).

- ...:

  Unused; present for S3 consistency.

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
d <- simulate_data("mlvar", seed = 1)
fit <- build_mlvar(d, vars = attr(d, "vars"),
                   id = "id", day = "day", beep = "beep")
print(fit)
summary(fit)
} # }
```
