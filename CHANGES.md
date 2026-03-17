# Changes

### 2026-03-17 — Add covariate-integrated MMM to build_mmm()
- R/mmm.R: Added `covariates` parameter to `build_mmm()`. Covariates influence cluster membership within EM: M-step replaces `colMeans(post)` with `nnet::multinom(post ~ covariates)` for covariate-dependent mixing proportions. New internals: `.mmm_softmax_mstep()`, `.mmm_predict_log_pi()`. Updated param count, return object ($covariates), print/summary/plot S3 methods. Strips covariate columns from raw_data before state extraction.
- tests/testthat/test-mmm.R: 9 new tests: covariate return structure, coefficient direction recovery, nested model LL, param count, S3 methods, plot, data.frame input, NA handling. 80/80 pass.
- Numerical equivalence: beta recovery within 0.4 of true (N=200), covariate LL > plain LL, BIC favors covariate model when effect is real, integrated captures covariate effects that post-hoc misses.

### 2026-03-17 — Add post-hoc covariate analysis to cluster_data()
- R/cluster_data.R: Added `covariates` parameter to `cluster_data()`. Accepts formula, character vector, string, or data.frame. Runs `nnet::multinom()` after clustering, produces cluster profiles (mean/SD/median for numeric, counts/% for categorical) and odds ratio table with z/p/CI. Stored in `$covariates` field. Internal functions: `.resolve_covariates()`, `.compute_cluster_profiles()`, `.run_covariate_analysis()`. Updated print/summary/plot S3 methods.
- DESCRIPTION: Added `nnet` to Suggests.
- tests/testthat/test-cluster_data.R: 105 new tests total (including earlier fix). Covariate tests: 5 input forms, profiles, k=2 glm cross-validation, k>2 manual verification, S3 methods, edge cases (NAs, constant, row mismatch), netobject metadata, build_network round-trip. 280/280 pass.

### 2026-03-17 — Fix cluster_data tests and add input extraction / dispatch coverage
- tests/testthat/test-cluster_data.R: Fixed 3 bugs (expect_s3_class/expect_null don't accept `info`, wrong encoding assertion, tna ward.D2 lowercasing). Added 10 new tests: input extraction (netobject, tna, cograph_network, association rejection) and build_network dispatch (class, default method, method override, sub-network sizes, clustering metadata). 175/175 pass.

### 2026-03-16 — Add cograph_network support across all downstream functions
- R/utils.R: Added `.as_netobject()` internal converter — coerces `cograph_network` to `netobject`. Decodes integer-encoded sequence data for transition methods, preserves numeric data for association methods. Infers method from `$meta$tna$method` or matrix symmetry.
- R/bootstrap_network.R: Added `cograph_network` coercion before netobject check.
- R/permutation_test.R: Added `cograph_network` coercion for both `x` and `y` before validation.
- R/reliability.R: Added `cograph_network` coercion in the variadic argument loop.
- R/centrality_stability.R: Added `cograph_network` coercion before netobject check.
- R/boot_glasso.R: Added `cograph_network` coercion before input dispatch.
- R/mcml.R: Added `cograph_network` case to `cluster_summary()` matrix extraction. Added coercion in `build_mcml()` before input type detection.
- R/mmm.R: Replaced inline `cograph_network` handling with `.as_netobject()` call.
- tests/testthat/test-cograph_network.R: 28 tests covering converter, bootstrap, permutation, reliability, centrality_stability, build_mmm, cluster_summary, boot_glasso.
- Tests: 1936 total (all pass, 0 failures).

### 2026-03-15 — Add reliability(), centrality_stability(); fix void/missing state handling
- R/reliability.R: New file — `reliability()` function for split-half network reliability assessment. Supports single and multi-model comparison, 4 metrics (mean/median/max abs dev, correlation), optional scaling (minmax/standardize/proportion) for cross-method comparability. Reuses `.precompute_per_sequence()` fast path for transitions. Includes `print.net_reliability` and `plot.net_reliability` S3 methods (density facets with dashed mean lines).
- R/estimators.R: Added `.void_markers` and `.clean_states()` helper. Applied to all 6 state-extraction paths (wide transitions, long transitions, wide co-occurrence, long co-occurrence, wide attention, long attention). TraMineR void (`%`), missing (`*`), empty strings, `"NA"`, `"NaN"` are now treated as `NA` before state detection.
- R/build_network.R: State/metadata column split now cleans void markers before checking `vals %in% nodes`. `$data` in netobject has void markers replaced with `NA` (character/factor columns only, numeric untouched).
- tests/testthat/test-reliability.R: 63 tests covering structure, reproducibility, metrics bounds, multi-model stacking, scaling options, warnings, edge cases, print/plot methods.
- R/centrality_stability.R: New file — `centrality_stability()` for CS-coefficient estimation. Computes strength centralities via matrix ops (no igraph). Betweenness/closeness use igraph. Pre-computed per-sequence counts for transitions. `loops` parameter (default FALSE) for diagonal exclusion in centrality only. `print`, `summary`, `plot` S3 methods. 3.7x faster than tna's `estimate_cs()`.
- tests/testthat/test-centrality_stability.R: 41 tests covering structure, reproducibility, all measures, loops parameter, custom parameters, tna equivalence, print/summary/plot methods.
- R/pathways.R: New file — `pathways()` S3 generic with methods for `net_hon`, `net_hypa`, `net_mogen`. Extracts higher-order pathway strings in `"A B -> C"` format compatible with `cograph::plot_simplicial()`. HON: higher-order edges (from_order > 1). HYPA: anomalous paths. MOGen: transitions at optimal/specified order. Supports `min_prob`, `order`, `type` filtering.
- R/hypa.R: Added `$edges` field to `net_hypa` output (alias of `$scores`) for consistency with HON/MOGen edge data structure.
- tests/testthat/test-pathways.R: 27 tests covering all three pathway methods, filtering, edge cases, HYPA $edges consistency, cograph parser compatibility.
- NAMESPACE: Exported `reliability`, `centrality_stability`, `pathways`, S3methods for print/plot/summary of `net_reliability`, `net_stability`, and pathways methods.

### 2026-03-14 — Add attention, wtna, import_onehot; remove tidyverse
- R/data_conversion.R: Rewrote `wide_to_long()` and `long_to_wide()` from dplyr/tidyr to base R; added `import_onehot()` for converting binary indicator data to sequences
- R/frequencies.R: Rewrote 6 internal helpers (`.standardize_long`, `.pivot_wide_to_long`, `.fmt_frequency`, `.fmt_onehot`, `.fmt_edgelist`, `.fmt_follows`) from dplyr/tidyr to base R
- R/utils.R: Deleted `safe_bind_rows()` (dead code using `dplyr::bind_rows`); cleaned globalVariables
- R/Nestimate-package.R: Removed `@import dplyr`; pruned globalVariables of all dplyr/tidyselect names
- R/estimators.R: Added `.count_attention_wide()`, `.count_attention_long()`, `.estimator_attention()` — decay-weighted attention transitions
- R/estimator_registry.R: Registered `"attention"` estimator
- R/estimate_network.R: Added `atna = "attention"` alias
- R/build_network.R: Added `attention` label in `print.netobject()`
- R/wtna.R: New file — `wtna()` for window-based TNA on one-hot binary matrices (transition, cooccurrence, both)
- DESCRIPTION: Removed dplyr and tidyr from Imports
- man/safe_bind_rows.Rd: Deleted
- Tests: Added test-estimator_attention.R (23 tests), test-import_onehot.R (22 tests), test-wtna.R (29 tests). All 2211 tests pass (1 pre-existing stochastic failure in mlvar).

### 2026-03-14 — Initial package creation (split from Saqrlab v0.3.0)
- Created Nestimate package with 22 R source files from Saqrlab computation modules
- R/build_network.R, estimators.R, estimator_registry.R: core estimation pipeline
- R/bootstrap_network.R, boot_glasso.R, permutation_test.R: inference
- R/gimme.R, hon.R, honem.R, hypa.R, mogen.R: advanced methods
- R/mcml.R, mlvar.R: multi-level methods
- R/temporal_network.R, velocity_tna.R: temporal analysis
- R/frequencies.R, data_conversion.R, extraction.R, network_comparison.R: utilities
- R/utils.R, estimate_network.R, fit_network_model.R: infrastructure/deprecated
- R/Nestimate-package.R: package init with .onLoad() for estimator registration
- tests/testthat/helper-simulate.R: standalone test data generators
- 18 test files, 2137/2138 passing
- docs/ARCHITECTURE.md, docs/STATUS.md: project documentation

### 2026-03-14 — Proximity timeline redesign
- R/temporal_network.R: Rewrote `.plot_proximity_timeline()` — smooth variable-width micro-segments via spline interpolation (~200 points), pre-scaled linewidth (0.5–3px range), Okabe-Ito colorblind-safe default palette, direct endpoint labels (no legend box). Added `highlight` parameter. Removed `render_edges`, `smooth` parameters and `.build_edge_segments()`.
- R/temporal_network.R: Rewrote `.compute_proximity_mds()` — strength-based interleaved slot positioning (hub at y=0, rank-spread to ±1).
