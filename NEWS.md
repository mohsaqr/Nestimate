# Nestimate (development)

## Bug fixes
* `.extract_edges_from_matrix()` no longer drops the diagonal. Netobjects
  built via `.wrap_netobject()` (and therefore everything from
  `build_network()`, `build_mcml()`, `bootstrap_network()`, `build_mmm()`,
  `wtna()`, `as_tna()`) now have `$edges` containing every non-zero
  matrix entry, including self-loops. Previously `$weights` and `$edges`
  were silently inconsistent on any matrix with a non-zero diagonal,
  causing downstream consumers (e.g. `cograph::centrality()` on an
  MCML macro) to under-count node degree by 2.

## New features
* `plot_state_frequencies()` — native S3 generic for state-frequency
  plots across `netobject`, `netobject_group`, `mcml`, and `htna`.
  Defaults to a marimekko (mosaic) layout where column widths reflect
  per-group totals and segment heights reflect within-group state
  proportions; also supports a colored-bars style and a per-group
  faceted marimekko. Uses the package Okabe-Ito palette throughout.
* `plot_mosaic()` — exported low-level marimekko primitive built on
  `geom_rect()` with cumulative-width / cumulative-height geometry.
  Reusable for any tidy `data.frame(group, state, weight)` input.

# Nestimate 0.4.4

## Bug fixes
* `passage_time()` and `markov_stability()` now raise an explicit error
  naming the dead state when a transition-matrix row sums to zero,
  instead of silently propagating `NaN` through `eigen`/`solve`. Zero
  rows mean the chain is not ergodic; mean first passage times are
  undefined. Shared helper `.mpt_normalize_rows()` factored out of both
  entry points.
* `.prepare_association_input()` no longer hard-rejects non-square
  numeric matrices. For association methods (glasso, pcor, cor) the
  netobject's `$data` slot is a numeric matrix (not a data.frame). Any
  downstream caller that row-subsetted `$data` and re-invoked the
  estimator (`centrality_stability()`, `bootstrap_network()`,
  `reliability()`) was silently producing NULL centralities caught by
  `tryCatch`, which surfaced as an "all centrality measures have zero
  variance" warning or all-`NaN` correlations. The matrix branch now
  recognises non-square input as raw observation data and recursively
  re-enters through the data-frame branch. Square symmetric matrices
  (pre-computed correlation / covariance) still go through the
  symmetric-matrix path with the symmetry check intact.

## New parameters
* `build_network()` gains `state_cols` and `metadata_cols` parameters
  (both default `NULL`). Explicit overrides for the state-vs-metadata
  column classifier, which previously used a "values-in-nodes"
  heuristic that silently misclassifies metadata columns whose values
  coincide with node labels (e.g. a `condition` column with levels
  `"A","B","C"` when nodes are `"A","B","C"`). Validation: error on
  overlap between the two vectors, error on column names not present
  in the input data. Forwarded through the `group = ...` recursive
  dispatch so per-group calls honour the override.

## Removed
* `plot.net_link_prediction()` and `plot.mcml()` removed. Nestimate is
  a computation engine — visualization is the user's concern.
  Previously both methods called `cograph::` directly, violating the
  stated dependency invariant (Nestimate -> cograph direction
  forbidden). Users call `cograph::splot(net)` or
  `cograph::plot_mcml(fit)` directly.

## Documentation
* `wtna()` `@param type` now flags that `type = "relative"` combined
  with `method = "cooccurrence"` produces an asymmetric matrix
  (conditional co-occurrence given row state), not a symmetric
  undirected weight matrix. Use `type = "frequency"` if symmetric
  counts are required.

## Testing infrastructure
* New numerical-equivalence tests (gated by `NESTIMATE_EQUIV_TESTS=true`):
  `test-equiv-permutation.R` (vs. `stats::p.adjust` + hand-coded base-R
  permutation loop), `test-equiv-mlvar.R` (vs. `mlVAR::mlVAR` at
  machine precision), `test-equiv-association-rules.R` (vs.
  `arules::apriori`), `test-equiv-link-prediction.R` (vs. clean-room
  matrix algebra + `igraph::similarity`), `test-equiv-centrality-stability.R`
  (vs. `bootnet::corStability`). Total ~162k per-value comparisons;
  all within machine precision except centrality-stability which uses
  a documented drop-grid tolerance because bootnet uses `igraph`
  path-based centrality and Nestimate uses Floyd-Warshall.
* New HON-family equivalence tests under `local_testing_and_equivalence/`
  validating HON, HONEM, HYPA, MOGen, and hypergraph against pathpy
  2.2.0 (via reticulate), `BiasedUrn`, `RSpectra`, and `HyperG`. Not
  shipped in the R-package `tests/` directory; added to `.Rbuildignore`.
* Branch-matrix coverage added for the four many-mode APIs (`wtna`,
  `bootstrap_network`, `build_clusters`, `sequence_plot`) — systematic
  cross-product tests over all combinations of mode parameters to
  catch regressions where one branch silently diverges.

# Nestimate 0.4.3

## CRAN resubmission (addresses incoming-check NOTEs on 0.4.2)
* Tarball now ships `build/vignette.rds` (the vignette index). Previous 0.4.2 build used `R CMD build --no-build-vignettes`, which preserved pre-built `inst/doc/*.html` but stripped the index — CRAN flagged "VignetteBuilder field but no prebuilt vignette index."
* `test-gimme.R` now `skip_on_cran()`. GIMME tests fit a lavaan SEM per subject and took ~50s locally (2-3× on Windows), pushing total check time to 11 min on win-devel. Full test suite still runs in CI and local dev.

# Nestimate 0.4.2

## CRAN resubmission
* Full `--as-cran --run-donttest` audit pass.
* Purged stale `.Rcheck/` and `Meta/` build artifacts from working tree; added explicit `^Nestimate\.Rcheck$` and `^\.\.Rcheck$` entries to `.Rbuildignore` as belt-and-suspenders against repeat-submission contamination.

# Nestimate 0.4.1

## CRAN resubmission
* Pre-built vignettes included in `inst/doc/` as required by CRAN.
* Fixed 301-redirect URLs in README.
* Added `skip_on_cran()` to slow test block to keep check time under 10 minutes.

# Nestimate 0.4.0

## New functions
* `build_mlvar()` — multilevel VAR networks from ESM/EMA panel data. Estimates temporal (directed), contemporaneous (undirected), and between-subjects (undirected) networks matching `mlVAR::mlVAR()` at machine precision.
* `build_mmm()` / `compare_mmm()` — mixture of Markov models via EM, with BIC/AIC/ICL model selection and optional covariate regression.
* `cooccurrence()` — standalone co-occurrence network builder supporting 6 input formats and 8 similarity methods.
* `sequence_compare()` — k-gram pattern comparison across groups with optional permutation testing.
* `sequence_plot()` / `distribution_plot()` — base-R sequence index and state distribution plots with clustering integration.
* `build_simplicial()`, `persistent_homology()`, `q_analysis()` — topological analysis of networks via simplicial complexes.
* `nct()` — Network Comparison Test matching `NetworkComparisonTest::NCT()` at machine precision.
* `build_gimme()` — group iterative mean estimation for idiographic networks via lavaan.
* `passage_time()`, `markov_stability()` — Markov chain passage times and stability analysis.
* `predict_links()` / `evaluate_links()` — link prediction with 6 structural similarity methods.
* `association_rules()` — Apriori association rule mining from sequences or binary matrices.
* `predictability()` — node predictability for glasso/pcor/cor networks.
* `build_hon()`, `build_honem()`, `build_hypa()`, `build_mogen()` — higher-order network methods (HON, HONEM, HYPA, MOGen) now `cograph_network`-compatible.

## New datasets
* `human_long`, `ai_long` — canonical long-format human–AI pair programming interaction sequences (10,796 turns, 429 sessions).
* `chatgpt_srl` — ChatGPT-generated SRL scale scores for psychological network analysis.
* `trajectories` — 138-student engagement trajectory matrix (15 timepoints, 3 states).

## API
* `build_clusters()`, `network_reliability()`, `permutation()`, and `prepare()` replace earlier internal names for consistency with the `build_*` naming convention.
* `mgm` estimator added (`method = "mgm"`) for mixed continuous + categorical data via nodewise lasso, matching `mgm::mgm()` at machine precision.

## Bug fixes
* `build_mmm()` no longer crashes on platforms where `parallel::detectCores()` returns `NA` (macOS ARM64 CRAN check failure).
* `gimme` convergence filter now correctly handles all typed `NA` variants (`NA_character_`, `NA_real_`, etc.).
* `NaN` values in numeric metadata aggregation (all-`NA` sessions) normalized to `NA_real_`.
* HYPA p-values corrected; `hypa_score` column renamed to `p_value`.

## CRAN compliance
* `.data` pronoun added to `globalVariables()`.
* `base::.rowSums()` / `base::.colSums()` replaced with `rowSums()` / `colSums()`.
* `dev.new()` guarded by `interactive()` — no side effects under knitr or CI.
* Equivalence test files excluded from the built tarball.

## Performance
* `do.call(rbind, ...)` replaced with `data.table::rbindlist()` in `mcml.R` and `sequence_compare.R`.

# Nestimate 0.3.4

* HYPA: Renamed `hypa_score` column to `p_value` for clarity. Added `$over`, `$under`, `$n_over`, `$n_under` fields to `net_hypa` objects. Scores are now pre-sorted with anomalous paths first.
* HYPA: `summary.net_hypa()` now shows over/under-represented paths separately with a configurable `n` parameter.
* `pathways.netobject()`: New S3 method to extract higher-order pathways directly from a netobject (builds HON or HYPA internally).
* `path_counts()`: Now handles NAs in trajectories by stripping them before k-gram counting.

# Nestimate 0.2.15

* Preparing for publication

# Nestimate 0.2.0

* Reduced hard dependencies from 6 to 4 Imports (ggplot2, glasso, data.table, cluster).
* Removed igraph from Imports — `centrality_stability()` and `boot_glasso()` now accept a `centrality_fn` parameter for external centrality computation.
* Removed tna from Imports — moved to Suggests (only used for input class detection).
* Implemented `graphical_var()` from scratch using coordinate descent lasso + graphical lasso with EBIC model selection, eliminating the graphicalVAR dependency.
* Dropped `ml_graphical_var()` — users should use `mlvar()` for multilevel VAR.
* Removed cograph plot wrappers — `plot.netobject()`, `plot.net_bootstrap()`, `plot.net_permutation()`, `plot.net_hon()`, `plot.net_hypa()` and `as_cograph()` removed. Users call cograph plotting functions directly on netobjects.
* Added `attention` estimator for decay-weighted transition networks.
* Increased test coverage from 84.5% to 96.1% (2780 tests).
* Passes R CMD check with 0 errors, 0 warnings, 0 notes.

# Nestimate 0.1.0

* Initial release. Split from Saqrlab v0.3.0.
* Core estimation via `build_network()` with 8 built-in estimators.
* Bootstrap inference (`bootstrap_network()`), permutation testing (`permutation()`), EBICglasso bootstrap (`boot_glasso()`).
* Higher-order networks: HON, HONEM, HYPA, MOGen.
* GIMME, MCML, multilevel VAR, graphical VAR.
* Temporal network analysis and velocity TNA.
* Dual-class `c("netobject", "cograph_network")` output for cograph compatibility.
