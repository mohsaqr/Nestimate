# Tidy per-state summary of a `chain_structure`.

Returns a single data.frame with one row per state, combining every
per-state metric
[`chain_structure()`](https://mohsaqr.github.io/Nestimate/reference/chain_structure.md)
computes. Always includes `state`, `classification`, `period`,
`return_probability` (the diagonal of the hitting matrix) and
`persistence` (the diagonal of the transition matrix); adds `sojourn`
whenever it is finite, the chain's `stationary_probability` when
irreducible, and absorption columns when the chain has any absorbing
states.

## Usage

``` r
# S3 method for class 'chain_structure'
summary(object, ...)
```

## Arguments

- object:

  A `chain_structure` object.

- ...:

  Ignored.

## Value

A `data.frame` with one row per state. Columns described above.

## Details

Columns are ordered for readability: identifiers first, classification
second, dynamic per-state metrics last.
