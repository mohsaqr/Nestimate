# Build a Co-occurrence Network (CNA)

Convenience wrapper for `build_network(method = "co_occurrence")`.
Computes co-occurrence counts from binary or sequence data.

## Usage

``` r
build_cna(data, start = FALSE, end = FALSE, ...)
```

## Arguments

- data:

  Data frame (sequences or per-observation frequencies) or a square
  symmetric matrix (correlation or covariance).

- start:

  Boundary marker prepended to every sequence as an explicit start state
  (a pure source: no incoming edges, every sequence's first transition
  is `start -> first_observed`). `FALSE` (default) adds nothing; `TRUE`
  uses the label `"Start"`; a single string uses that string as the
  label. Only valid for the transition methods (`relative`, `frequency`,
  `co_occurrence`, `attention`); errors otherwise.

- end:

  Boundary marker placed in the single cell after each sequence's last
  observed (non-`NA`) state, as an explicit terminal state (a pure sink:
  no outgoing edges, no self-loop – distinct from
  [`mark_terminal_state`](https://saqr.me/Nestimate/reference/mark_terminal_state.md),
  which fills all trailing NAs into an absorbing state). `FALSE`
  (default) adds nothing; `TRUE` uses the label `"End"`; a single string
  uses that string as the label. Same method restriction as `start`.

- ...:

  Additional arguments passed to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).

## Value

A `netobject` (see
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)).

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`cooccurrence`](https://saqr.me/Nestimate/reference/cooccurrence.md)
for delimited-field, bipartite, and other non-sequence co-occurrence
formats.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
net <- build_cna(seqs)
```
