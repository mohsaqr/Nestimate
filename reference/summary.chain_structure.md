# Tidy per-state summary of a `chain_structure`

Returns a single data.frame with one row per state, combining every
per-state metric
[`chain_structure()`](https://saqr.me/Nestimate/reference/chain_structure.md)
computes. Always includes `state`, `classification`, `period`,
`persistence` (the diagonal of the transition matrix),
`return_probability` (the diagonal of the hitting matrix) and
`sojourn_steps` (`1 / (1 - persistence)`, which is `Inf` for an
absorbing state). Adds the chain's `stationary_probability` when the
chain is irreducible, and absorption columns when it has any absorbing
states: `absorption_probability` for a single absorbing state or one
`absorbed_in_<state>` column per state when there are several, plus
`mean_absorption_time`.

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

A `data.frame` with one row per state, of class
`c("summary_chain_structure", "data.frame")`, carrying the chain-level
flags (`is_regular`, `is_irreducible`, `is_aperiodic`, `is_reversible`,
`n_classes`, `absorbing_states`) as attributes, which its
[`print()`](https://rdrr.io/r/base/print.html) method shows as a header.
Columns as described above.

## Details

Columns are ordered for readability: identifiers first, classification
second, dynamic per-state metrics last.
