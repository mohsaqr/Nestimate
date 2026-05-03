# Tidy coefficients from a fitted mlvar model

Generic accessor for the tidy coefficient table stored on a
[`build_mlvar()`](https://mohsaqr.github.io/Nestimate/reference/build_mlvar.md)
result. Returns a `data.frame` with one row per `(outcome, predictor)`
pair and columns `outcome`, `predictor`, `beta`, `se`, `t`, `p`,
`ci_lower`, `ci_upper`, `significant`.

## Usage

``` r
coefs(x, ...)

# S3 method for class 'net_mlvar'
coefs(x, ...)

# Default S3 method
coefs(x, ...)
```

## Arguments

- x:

  A fitted model object — currently only `net_mlvar` is supported.

- ...:

  Unused.

## Value

A tidy `data.frame` of coefficient estimates.

## Details

Only the within-person (temporal) coefficients are tabulated — these are
the lagged fixed effects that populate `fit$temporal`. The
between-subjects effects that go into `fit$between` are handled via the
`D (I - Gamma)` transformation and are not exposed as a separate tidy
table.

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
