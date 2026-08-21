# Summary method for `net_entropy_trajectory`

Summary method for `net_entropy_trajectory`

## Usage

``` r
# S3 method for class 'net_entropy_trajectory'
summary(object, ...)
```

## Arguments

- object:

  A `net_entropy_trajectory` object.

- ...:

  Ignored.

## Value

Tidy per-group data.frame: windows, mean/sd/min/max entropy, entropy at
the first and last window, and their difference (negative =
routinization).
