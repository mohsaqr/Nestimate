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
