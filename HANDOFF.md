# Session Handoff — 2026-05-02

## Completed

### `.extract_edges_from_matrix()` self-loop fix (R/estimate_network.R:273)

The helper used by `.wrap_netobject()` to populate every netobject's
`$edges` data.frame was silently filtering the diagonal:

```r
# before
if (directed) {
  idx <- which(mat != 0 & row(mat) != col(mat), arr.ind = TRUE)
} else {
  idx <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
}
# after
if (directed) {
  idx <- which(mat != 0, arr.ind = TRUE)
} else {
  idx <- which(mat != 0 & row(mat) <= col(mat), arr.ind = TRUE)
}
```

Result: every netobject built by `build_network()`, `build_mcml()`,
`build_mmm()`, `bootstrap_network()`, `wtna()`, or `as_tna()` now has
`$edges` containing every non-zero matrix entry (including the
diagonal), matching `$weights`. Previously `$weights` and `$edges`
silently disagreed on any matrix with a non-zero diagonal — the
canonical case being MCML macros, where diagonal entries represent
intra-cluster retention.

Symptom that triggered the investigation:
`cograph::centrality_degree(MCMLL_tna$macro)` returned 10 while
`cograph::centrality_degree(MCMLL_tna$macro$weights)` returned 12 on a
6-cluster fully-connected macro. The matrix path goes through
`igraph::graph_from_adjacency_matrix(weighted = TRUE)` which keeps the
diagonal; the netobject path went through cograph's
`network_to_igraph()` which builds from `$edges` and so saw a
loop-free graph.

## Current State

- All 7 internal callers of `.extract_edges_from_matrix()` re-checked;
  none assumed loop-free edges.
- Nestimate test suite: 781 tests, 0 failures, 0 errors, 18
  environment-gated equiv skips (set `NESTIMATE_EQUIV_TESTS=true` to
  run those).
- cograph regression test added on its side
  (`tests/testthat/test-validate-nestimate-bootstrap-permutation.R`)
  that asserts `centrality(netobj) == centrality(netobj$weights)` when
  `diag(weights) != 0`. Passes with this Nestimate.

## Open Issues

- Existing on-disk netobjects saved before this fix still have
  loop-free `$edges`. They will keep producing wrong centrality counts
  until rebuilt. No migration helper provided — users should rebuild
  from sequence/matrix data.

## Next Steps

- Optional: have cograph add a defensive fallback in
  `network_to_igraph()` that prefers `$weights` over `$edges` when
  both are present and the diagonals disagree. This was deliberately
  *not* added in this session because the source-of-truth fix here is
  preferred. Reconsider only if users hit the legacy-object case
  often.
- Bump `Nestimate` dev version when ready and note the bug fix in the
  next CRAN release line.

## Context

- File touched: `R/estimate_network.R` (one helper, ~5 lines).
- Verified end-to-end via synthetic 6-cluster MCML reproduction:
  `nrow($macro$edges)` went from 30 → 36 (the missing 6 self-loops),
  `centrality_degree($macro)` went from 10 → 12 matching the matrix
  path.
