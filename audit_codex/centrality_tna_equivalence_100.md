# Centrality tna Equivalence Audit

Date: 2026-06-29

Scope: 100 transition networks generated with `simulate_sequences()` from `tests/testthat/helper-simulate.R`, the local Saqrlab-style simulation stand-in used by Nestimate tests.

Tolerance: `1e-10`

Compared measures:
- `OutStrength`
- `InStrength`
- `ClosenessIn`
- `ClosenessOut`
- `Closeness`
- `Betweenness`
- `BetweennessRSP`
- `Diffusion`
- `Clustering`

Checks:
- Raw Nestimate centralities used `normalize_diffusion = FALSE` and were compared to `tna::centralities(normalize = FALSE)`.
- Nestimate's default `Diffusion` was compared separately to `tna::centralities(measures = "Diffusion", normalize = TRUE)`.

Results:
- Networks checked: `100`
- Passed: `100`
- Failed: `0`
- Maximum absolute delta, all checks: `2.7755575615628914e-17`
- Maximum absolute delta, raw centrality parity: `2.7755575615628914e-17`
- Maximum absolute delta, default normalized Diffusion: `0`

CSV detail: `audit_codex/centrality_tna_equivalence_100.csv`
