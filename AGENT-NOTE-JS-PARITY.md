# Agent Note — Downstream Numerical Pin

Heads-up for anyone editing Nestimate's reliability / stability / Markov functions:

A sibling project (Dynalytics Desktop, via `tnaj`) maintains **cross-language parity tests** against these R functions at machine-ε tolerance. If you change numerical behaviour (formula, default argument, diagonal handling, factor-level ordering, summation order, rounding, scaling) on any of the functions below, the downstream parity suite will break.

**Pinned functions (scalar / `netobject` input):**
- `passage_time()`
- `markov_stability()`
- `markov_order_test()`
- `network_reliability()`
- `centrality_stability()`
- `casedrop_reliability()`

**Pending parity target (`netobject_group` input):**
The `.netobject_group` S3 methods that fan these functions out over a named list of networks are **not yet ported to JS**. The JS side plans thin wrappers taking `Record<string, TNA>` and returning `Record<string, Result>`, matching your `*_group` list naming / ordering / per-group error behaviour exactly. Changes to the group dispatch semantics (how groups are named, how a single-group failure propagates, aggregate fields on the `_group` S3 class) will be adopted on the JS side when Stage 6 lands — heads-up, not a blocker.

**Not a blocker for changes** — just notify the Dynalytics side so the JS port can follow, or expect a parity-test red and update tolerances deliberately.

**Safe changes** (no downstream impact):
- `print` / `summary` / `plot` S3 methods (display rounding does not affect raw computation; the parity suite recomputes raw values).
- Docs, examples, new arguments with backward-compatible defaults.
- New functions.

**Breaking changes** (parity suite will red):
- Formula changes in any of `.split_half_metrics`, `.scale_matrix`, `.calculate_cs`, `.mot_g2_*`, `.mogen_*`, `.mpt_*`.
- Default argument flips (`method`, `threshold`, `certainty`, `scale`, `loops`, `include_diag`).
- Changes to internal sampling / permutation loop structure.
- Output list-element renames or reorderings in `net_reliability`, `net_casedrop_reliability`, `net_stability`, `net_markov_stability`, `net_mpt`, `net_markov_order`.

**Where the tests live:**
`~/Documents/Github/Dynalytics_Desktop/verify/R/run_nestimate_*.R` — R wrappers that reproduce each function's internals and dump ground truth. Useful as reading material when refactoring any of the pinned functions, because the wrapper loops are a line-by-line shadow of the corresponding Nestimate internal.

No action required if you aren't touching these functions.
