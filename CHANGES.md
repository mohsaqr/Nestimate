# Changes

### 2026-03-27 — CRAN extra check fixes: 0 errors / 0 warnings / 0 notes
- R/simplicial.R: Fixed `print.persistent_homology` and `plot.persistent_homology` examples — removed incorrect `sc <- build_simplicial(net, type = "clique")` step and changed `persistent_homology(sc)` to `persistent_homology(net)` (function takes a netobject/matrix, not a simplicial_complex); plot example converted to `\dontrun{}` to avoid Unicode rendering failure with β label in non-UTF8 locales
- R/centrality_stability.R: Removed non-existent `n_subsets = 3` parameter from all 3 examples (print/summary/plot.net_stability)
- R/mmm.R: Fixed `compare_mmm` examples — changed `k_range = 2:3` to `k = 2:3` (correct parameter name)

### 2026-03-26 — Informative bootstrap print methods
- R/bootstrap_network.R: Rewrote `print.net_bootstrap` to show method label, directionality, iterations, nodes, significant/total edge counts, CI level, inference method, CR range, and top 5 significant edges with mean weight, CI, and significance stars
- R/bootstrap_network.R: Rewrote `print.net_bootstrap_group` to show all groups, per-group sig/total edge counts, shared significant edge count, and top 5 shared edges with per-group means
- R/bootstrap_network.R: Fixed R quirk — `paste0(character(0), "→", character(0))` returns `"→"` with length 1; fixed by guarding with `which()` before `paste0`
- R/bootstrap_network.R: Fixed `[[` subscript-out-of-bounds for unknown method names; now uses `[` with `is.na()` guard
- tests/testthat/test-bootstrap_network.R: Updated print tests to match new output format
- tests/testthat/test-boot_glasso.R: Updated edge_diff fill-column test to match new `fill_val` continuous encoding
- Tests: 2739 pass, 0 failures

### 2026-03-20 — CLAUDE.md refresh, temporal spec, vibcoding datasets
- CLAUDE.md: Removed phantom modules, added simplicial complex, fixed reliability() name, added sidelined section
- docs/temporal/SPEC.md: Enterprise-grade temporal network build spec with cograph-compatible output
- docs/temporal/ref-*.R: Copied sidelined implementations as reference
- data-raw/vibcoding.R: Build script for 12 datasets from raw human-AI vibe coding CSV
- R/data.R: Documentation for 12 datasets (9 long, 2 wide, 1 edge list)
- data/*.rda: 12 datasets bundled (human_ai, human_ai_cat, human_ai_super, human_detailed, human_cat, human_super, ai_detailed, ai_cat, ai_super, human_wide, ai_wide, human_ai_edges)
- DESCRIPTION: Added LazyData: true, Depends: R (>= 3.5)
- Tests: 2779 pass, 1 pre-existing failure (test-prepare_data.R:382)

### 2026-03-18 — v0.2.0: Dependency reduction, coverage, CRAN prep
- DESCRIPTION: Reduced Imports from 6 to 4 (removed igraph, tna). Version 0.2.0.
- R/graphical_var.R: Reimplemented from scratch — own coordinate descent lasso + glasso + EBIC. Dropped graphicalVAR dependency and ml_graphical_var().
- R/centrality_stability.R: Added `centrality_fn` parameter, removed igraph dependency. Built-in strength via rowSums/colSums.
- R/boot_glasso.R: Added `centrality_fn` parameter, removed igraph dependency. Removed cograph plot branch.
- R/build_network.R: Removed plot.netobject, plot.netobject_ml (cograph wrappers).
- R/bootstrap_network.R: Removed plot.net_bootstrap (cograph wrapper).
- R/permutation_test.R: Removed plot.net_permutation (cograph wrapper).
- R/gimme.R: Removed cograph branches from plot.net_gimme.
- R/hon.R: Removed plot.net_hon (cograph wrapper).
- R/hypa.R: Removed plot.net_hypa (cograph wrapper).
- R/mmm.R: Removed cograph/patchwork from plot.net_mmm. Added _R_CHECK_LIMIT_CORES_ guard.
- R/mogen.R: Removed cograph branch from plot.net_mogen.
- R/utils.R: Removed as_cograph() and all methods.
- R/Nestimate-package.R: Removed @import tna, @importFrom igraph. Added @return to 39 exported functions.
- R/mcml.R: Fixed example, @noRd ordering, added print/summary roxygen.
- Tests: 2780 total (up from 2235), coverage 96.1%. Code review fixes applied.
- CRAN: NEWS.md, cran-comments.md, .Rbuildignore (.claude), DESCRIPTION (cph, acronyms, URL). R CMD check --as-cran: 0/0/0.
- docs/cograph-interface.md: Object format spec for cograph integration.

### 2026-03-17 — Add coverage tests for hon/honem/hypa/mogen/pathways
- tests/testthat/test-hon.R: Added Section 11 (16 tests) covering .hon_parse_input all-NA rows, invalid type stop, collapse_repeats edge cases, .hon_kld empty/Inf paths, .hon_get_extensions missing order, .hon_sequence_to_node empty input, .hon_graph_to_edgelist empty graph, .hon_assemble_output empty rules, build_hon no-valid-trajectories stop, plot.net_hon branches (no-cograph, no-HO-pathways, successful). Added Section 12 (6 tests) for pathways.net_hon: character return, empty for first-order HON, min_count/top/min_prob/order filters. 237 pass, 1 skip.
- tests/testthat/test-honem.R: Added Section 7 (1 test) for plot.net_honem dim<2 message. 31 pass.
- tests/testthat/test-hypa.R: Added Section 7 (7 tests) covering build_hypa no-valid-trajectories, no-edges stop, summary with anomalies, plot.net_hypa no-cograph/no-anomalies/successful. Added Section 8 (4 tests) for pathways.net_hypa: type=all/over/under, empty on no-anomalies. 41 pass, 1 skip.
- tests/testthat/test-mogen.R: Added Section 9 (15 tests) covering .mogen_log_likelihood empty/single/order0/missing-key paths, build_mogen no-valid-trajectories, mogen_transitions default order/empty return, path_counts data.frame/list/short/top/k-validation, state_frequencies data.frame/list, print.net_mogen LRT branch, plot.net_mogen pathways no-cograph/no-pathways. Added Section 10 (5 tests) for pathways.net_mogen: character return, order=0 empty, min_count/top/min_prob filters. 80 pass, 1 skip.

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
