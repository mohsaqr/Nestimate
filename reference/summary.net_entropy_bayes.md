# Summary method for `net_entropy_bayes`

Summary method for `net_entropy_bayes`

## Usage

``` r
# S3 method for class 'net_entropy_bayes'
summary(object, ...)
```

## Arguments

- object:

  A `net_entropy_bayes` object.

- ...:

  Ignored.

## Value

The tidy edge table (data.frame), one row per observed transition,
sorted by posterior mean contribution, returned invisibly. The
chain-level and per-edge tables are printed as a side effect.
