# Learnings

### 2026-03-14
- [package split]: When splitting an R package, the `.onLoad()` function must be duplicated — it registers estimators in the global registry environment. Without it, `build_network()` fails with "Estimator not found."
- [test helpers]: Tests that depend on simulation functions from the parent package need standalone helpers. `helper-simulate.R` files in `tests/testthat/` are auto-loaded by testthat before test execution.
- [simulate_data("mlvar")]: The simplified test helper generates structurally correct VAR data but doesn't perfectly replicate the statistical properties of Saqrlab's version. One test (`test-mlvar.R:694`) fails marginally because of this. Fix: either improve the helper's data generation or adjust the test tolerance.
- [igraph::strength]: On unweighted temporal network snapshots, `strength()` equals `degree()`. For meaningful strength-based proximity positioning, the network needs weighted edges (multiple overlapping edges between the same pair).
- [proximity timeline]: Rank-based interleaved slot assignment (0, +s, -s, +2s, -2s...) guarantees one node at y=0 every bin. Rolling average smoothing destroys exact zeros — use spline interpolation only. Separate sign assignment creates two communities with a gap at center.
- [tidyverse removal]: `wide_to_long()` and `long_to_wide()` in `data_conversion.R` and 6 helpers in `frequencies.R` rewrote from dplyr/tidyr to base R. Key patterns: `do.call(rbind, lapply(...))` replaces `pivot_longer/pivot_wider`, `order()` replaces `arrange()`, `ave()` replaces `group_by+row_number()`, `table()+as.data.frame.matrix()` replaces `group_by+summarise+pivot_wider`. All 2137+ existing tests pass unchanged.
- [import_onehot rbind alignment]: When actors have different numbers of time points, `do.call(rbind, wide_list)` fails. Must align columns across groups by padding missing columns with `NA_character_` before binding.
- [wtna auto-detect]: `.wtna_auto_detect_codes()` must check `is.numeric(x)` before testing values against `c(0, 1)` — character columns cause `NAs introduced by coercion` warnings from `as.numeric()`.
- [attention estimator]: Returns raw decay-weighted counts (NOT row-normalized). Uses `tapply(d, pair_idx, sum)` for float accumulation — can't use `tabulate()` since weights are non-integer. Long format uses nested loops per group (position pairs), not vectorized like transitions.

### 2026-03-15
- [tna::group_regulation]: Has no group column — only T1-T26 sequence columns. For testing grouped operations, create a synthetic group column.
- [reliability split-half]: Named list assignment with duplicate keys silently overwrites — must collect objects and labels separately, then use `make.unique()` + `setNames()` to avoid losing models with the same method name.
- [data.frame row.names warning]: When building a data.frame inside `do.call(rbind, lapply(...))`, row names from the original data frame carry over and cause "row names were found from a short variable" warnings. Fix: set `row.names = NULL` in the `data.frame()` call.
- [TraMineR void markers]: `%` = void (padding), `*` = missing (non-response). tna silently filters these; Nestimate must do the same via `.clean_states()`. The `tna::engagement` dataset has 445 `%` values in T24/T25 — without cleaning, Nestimate detected a phantom 4th state.
- [build_network state/metadata split]: The `is_state_col` check must clean void markers before testing `vals %in% nodes`, otherwise columns with mixed valid+void values (like T24 with 85.5% valid + 14.5% `%`) get classified as metadata and dropped from `$data`. Also, cleaning `$data` values must skip numeric columns — one-hot (0/1) data should not be character-converted.
- [centrality stability]: tna's `estimate_cs` uses `loops = FALSE` (zeros diagonal before centrality). For relative networks, OutStrength = rowSums = 1 always WITH diagonal, but meaningful WITHOUT diagonal. The `loops` parameter controls this for centrality computation only — does not modify the stored matrix.
- [centrality fast path]: InStrength = colSums(mat), OutStrength = rowSums(mat) — no igraph needed. Betweenness/Closeness still require igraph. Strength-only runs at 3.8s vs full (with betweenness) at 10s for 9000 network builds.

### 2026-03-16
- [cograph_network conversion]: `.as_netobject()` in utils.R converts `cograph_network` → `netobject`. Key detail: only decode integer-encoded data for sequence methods (relative/frequency/co_occurrence/attention). Association methods (glasso/cor/pcor) store numeric observation data that must stay numeric. Detected via inferred method.
- [cograph_network structure]: `$weights` = matrix, `$nodes$label` = state names, `$data` = raw data (may be integer-encoded), `$directed` = logical, `$meta$tna$method` = estimation method (may be NULL). Class is `c("cograph_network", "list")`.
- [cluster_summary class]: `cluster_summary()` in mcml.R returns class `"mcml"`, not `"cluster_summary"`. The function name ≠ the class name.

### 2026-03-17
- [testthat info param]: `expect_s3_class()` and `expect_null()` do NOT accept `info` parameter. Use `expect_true(inherits(...), info=)` and `expect_true(is.null(...), info=)` instead.
- [tna ward.D2 bug]: tna v1.2.1 lowercases hclust method names, so `"ward.D2"` becomes `"ward.d2"` which errors. Skip ward.D2 in tna cross-validation tests. Our own code passes it correctly to `hclust()`.
- [test-mmm parallel fork]: `build_mmm()` tests fail during R CMD check with `"N simultaneous processes spawned"` — the check sandbox restricts `parallel::mclapply`. Pre-existing, not a regression.
- [build_network metadata]: `build_network()` classifies columns as state vs metadata by checking if all non-void values are in the network's node names. Character columns like Gender (M/F) can be misclassified as state columns if letters overlap. Numeric columns (Age, Score) always go to `$metadata` since they can't be state names.
- [nnet::multinom k=2]: For k=2 clusters, `summary(multinom)$coefficients` returns a named vector, not a matrix. Must check `is.null(dim(coefs))` and wrap into a 1-row matrix.
- [nnet vs glm tolerance]: `nnet::multinom` (BFGS) and `stats::glm(binomial)` (IRLS) use different optimizers. Coefficients match to ~1e-3, p-values to ~1e-2.
- [formula environment]: Building `stats::as.formula("cluster ~ X")` inside the fitting function ensures `cluster` resolves in `fit_df`. Building it in a different function causes `object 'fit_df' not found` errors due to formula environment scoping.
- [mmm covariate columns]: When `build_mmm()` receives a data.frame with covariate columns (numeric like Age), `.select_state_cols()` picks them up as state columns, inflating state count. Must strip covariate column names from `raw_data` before state extraction.
- [nnet::multinom matrix response]: `nnet::multinom()` accepts an N x k matrix as response — treats rows as multinomial proportions. This is the correct approach for soft-EM M-step. No weights needed; the posteriors already encode soft assignments.
- [mmm label switching]: EM may permute cluster labels relative to ground truth. Beta coefficients may have opposite sign. Tests must check both orderings (or use absolute values).
