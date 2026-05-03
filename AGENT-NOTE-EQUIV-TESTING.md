# Agent Note: Equivalence Testing for Nestimate

**Date**: 2026-04-12 **From**: cograph session (audit-driven equivalence
testing) **For**: Next Nestimate agent session

## Context

The cograph package just completed a comprehensive cross-validation
audit. An external review (Gemini CLI) identified that high line
coverage (100%) masked integration-level gaps — 3 real bugs survived
19,736 tests because the tests validated formulas without testing
plumbing (directedness handling, R6/S3 dispatch, metadata flow). A 4th
bug (`.compute_modularity` undirected formula) was caught during the
equivalence testing itself.

The fix: 8 new equivalence test files that generate 100 random networks
each and compare every value element-by-element against igraph ground
truth. Results: **137 functions, 235,629 values checked, 0 failures, max
delta 9.09e-11**. Each test emits CSV reports with mean/median/p95/max
deltas and vitest-compatible JSON to the CVS validation platform.

## Task: Replicate for Nestimate

Nestimate has 32 test files and 32 R source modules but **zero
equivalence tests** (`test-*equiv*` files). The package wraps igraph,
qgraph, glasso, and mlVAR — each wrapping layer is a potential
divergence point.

### What needs equivalence testing

#### Priority 1: Network estimation (core pipeline)

- `build_network(data, method = "relative")` → verify `$weights` matrix
  matches manual row-normalized transition counts
- `build_network(data, method = "frequency")` → verify against manual
  [`table()`](https://rdrr.io/r/base/table.html) counts
- `build_network(data, method = "co_occurrence")` → verify against
  manual co-occurrence matrix
- `build_network(data, method = "cor")` → verify against
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)
- `build_network(data, method = "pcor")` → verify against
  `ppcor::pcor()` or manual inverse correlation
- `build_network(data, method = "glasso")` → verify against
  [`qgraph::EBICglasso()`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html)
  directly
- `build_network(data, method = "ising")` → verify against
  [`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html)
- `build_network(data, method = "attention")` → verify against manual
  decay-weighted formula
- Ground truth: manual matrix computation for simple methods,
  qgraph/IsingFit for regularized

#### Priority 2: Bootstrap & permutation

- [`bootstrap_network()`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)
  → verify CI bounds contain true edge weights at stated coverage
- `permutation_test()` → verify p-values against manual permutation
  distribution
- [`boot_glasso()`](https://mohsaqr.github.io/Nestimate/reference/boot_glasso.md)
  → verify against
  [`qgraph::EBICglasso()`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html)
  on bootstrap samples
- Ground truth: manual bootstrap resampling + independent qgraph calls

#### Priority 3: Higher-order methods

- [`build_hon()`](https://mohsaqr.github.io/Nestimate/reference/build_hon.md)
  → verify HON adjacency against manual state-space expansion
- [`build_hypa()`](https://mohsaqr.github.io/Nestimate/reference/build_hypa.md)
  → verify p-values against manual multi-hypergeometric formula
- [`build_honem()`](https://mohsaqr.github.io/Nestimate/reference/build_honem.md)
  → verify embedding dimensions against manual matrix factorization
- Ground truth: manual formula implementation with explicit loops

#### Priority 4: Specialized models

- [`build_mlvar()`](https://mohsaqr.github.io/Nestimate/reference/build_mlvar.md)
  → verify temporal/contemporaneous/between matrices against
  `mlVAR::mlVAR()` directly
- [`build_gimme()`](https://mohsaqr.github.io/Nestimate/reference/build_gimme.md)
  → verify against
  [`gimme::gimme()`](https://rdrr.io/pkg/gimme/man/gimmeSEM.html)
  directly
- [`association_rules()`](https://mohsaqr.github.io/Nestimate/reference/association_rules.md)
  → verify against
  [`arules::apriori()`](https://rdrr.io/pkg/arules/man/apriori.html)
  (test file exists but may not do numerical equiv)
- [`predict_links()`](https://mohsaqr.github.io/Nestimate/reference/predict_links.md)
  → verify link prediction scores against manual formula
- Ground truth: direct calls to wrapped packages

#### Priority 5: Utility functions

- [`centrality_stability()`](https://mohsaqr.github.io/Nestimate/reference/centrality_stability.md)
  CS coefficient → verify against manual case-dropping bootstrap
- [`cluster_data()`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md)
  → verify cluster assignments against `igraph::cluster_*()` directly
- `prepare_data()` → verify sequence parsing against manual state
  extraction

### Pattern to follow

Use the exact infrastructure from cograph’s
`test-equiv-network-summary.R`:

``` r

skip_on_cran()
skip_coverage_tests()  # Add this helper to Nestimate's helper-simulate.R
skip_if_not_installed("igraph")

# Report infrastructure
.equiv_log <- new.env(parent = emptyenv())
.equiv_log$rows <- list()

.log_result <- function(func, config, n_checked, n_passed, n_failed,
                        max_abs_err, mean_abs_err, median_abs_err,
                        p95_abs_err, reference_package, notes = "") { ... }

.write_report <- function() { ... }       # CSV to tmp/
.write_cvs_report <- function() { ... }   # vitest JSON to validation inbox
```

Key elements: - **100 random datasets per function** (vary n_actors,
n_states, sequence lengths, densities) - **Per-value comparison** —
compare every matrix cell, not just summaries - **Delta statistics**:
mean, median, p95, max absolute error per (function, config) pair -
**NA/NaN consistency check** before numeric comparison - **TOL = 1e-8**
for floating point, exact for integer values (membership, counts) -
**Console report** with per-function PASS/FAIL and delta stats during
test run - **CSV report** to `tmp/{module}_equivalence_report.csv` -
**CVS JSON** to
`../validation/data/inbox/nestimate-{module}-{timestamp}.json`

### CVS integration

1.  Register Nestimate in
    `/Users/mohammedsaqr/Documents/Github/validation/.validationrc.json`:

    ``` json
    { "name": "nestimate", "path": "../Nestimate", "testsDir": "tests/testthat" }
    ```

2.  Each test file emits vitest-compatible JSON with
    `_cvs.target = "nestimate"` and per-assertion delta/tolerance fields

3.  The validation dashboard at `localhost:3847` will then track
    Nestimate’s numerical precision alongside cograph, JStats, tna-js,
    and other subscribers

### Test data generation

Nestimate works with sequence data, not adjacency matrices. Use
`helper-simulate.R` which already provides: - `simulate_mtna()` —
multi-actor TNA sequences - `simulate_data("mlvar")` — panel data for
mlVAR

For 100-network configs:

``` r

set.seed(2026)
N <- 100L
configs <- lapply(seq_len(N), function(i) {
  list(
    n_actors = sample(c(10, 20, 30, 50), 1),
    n_states = sample(c(3, 5, 7, 9), 1),
    seq_length = sample(c(20, 50, 100), 1),
    seed = sample.int(100000, 1)
  )
})
```

### File naming

Follow cograph convention: `test-equiv-{module}.R` -
`test-equiv-estimators.R` — all 8 estimation methods -
`test-equiv-bootstrap.R` — bootstrap_network + boot_glasso -
`test-equiv-permutation.R` — permutation_test - `test-equiv-hon.R` —
HON/HYPA/HONEM - `test-equiv-mlvar.R` — mlVAR (if build_mlvar is
exported) - `test-equiv-cluster.R` — cluster_data -
`test-equiv-association-rules.R` — vs arules

### Skip helper

Add to `helper-simulate.R`:

``` r

skip_coverage_tests <- function() {
  run_coverage <- Sys.getenv("NESTIMATE_COVERAGE_TESTS", unset = "true")
  if (!identical(run_coverage, "true")) {
    skip("Coverage tests skipped (set NESTIMATE_COVERAGE_TESTS=true to run)")
  }
}
```

### What cograph found that Nestimate should watch for

1.  **Directedness plumbing**: Functions that accept both directed and
    undirected input may canonicalize edges incorrectly (the `pmin/pmax`
    bug)
2.  **Wrapper overhead**: `build_network` → `.estimator_*` →
    igraph/qgraph call — each layer could silently transform weights or
    drop metadata
3.  **Formula variants**: Different packages implement “the same” metric
    differently (e.g., igraph’s `assortativity_nominal` vs
    `assortativity` with integer types — both valid but different
    formulas)
4.  **Default rounding**: Functions with `digits` parameters that
    default to non-NULL values silently round before comparison — always
    test with `digits = NULL`
5.  **Stochastic reproducibility**: Boot/permutation methods need fixed
    seeds at BOTH the cograph/Nestimate level AND the internal RNG
    restoration level

### Expected outcome

- 7-10 new test files
- ~50,000+ individual value comparisons
- CSV + CVS JSON reports for every module
- Full integration with validation dashboard
- Any discovered bugs documented in LEARNINGS.md
