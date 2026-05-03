# Summary method for `net_transition_entropy`

Returns a tidy per-state contribution table sorted by share of the
chain-level entropy rate (largest first), so the dominant contributors
to \\h(P)\\ are visible at a glance. Each row contains the stationary
mass, the raw and normalised row entropy, the additive contribution
\\\pi_i H(P\_{i\cdot})\\, and that contribution as a percentage of
\\h(P)\\.

## Usage

``` r
# S3 method for class 'net_transition_entropy'
summary(object, ...)
```

## Arguments

- object:

  A `net_transition_entropy` object.

- ...:

  Ignored.

## Value

A `summary.net_transition_entropy` containing

- table:

  tidy per-state data.frame, sorted by `contribution_pct` descending

- chain:

  tidy chain-level data.frame with raw and normalised \\h(P)\\,
  \\H(\pi)\\, redundancy, and ceiling

- base:

  logarithm base used
