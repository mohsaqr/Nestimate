# TODO

## Sidelined: velocity_tna (edge velocity/acceleration analysis)

**Status:** Moved to `sidelined/` — not exported, not tested, not documented in NAMESPACE.

**Files moved:**
- `sidelined/velocity_tna.R` — main implementation (1372 lines)
- `sidelined/test-velocity_tna.R` — test suite
- `sidelined/velocity_tna.Rd` — function documentation
- `sidelined/print.tna_velocity.Rd` — print method docs
- `sidelined/summary.tna_velocity.Rd` — summary method docs
- `sidelined/plot.tna_velocity.Rd` — plot method docs

**What it does:**
`velocity_tna()` takes a list of transition matrices ordered over time (e.g., weekly networks from weeks 1-15) and estimates how each edge weight changes — its velocity (rate of change) and acceleration (rate of change of the rate of change). This tells you which connections are strengthening, weakening, or stable, and whether those trends are accelerating or decelerating.

**Three estimation methods:**
1. **Regression** (`method = "regression"`): Fits OLS or beta regression per edge across time points. Returns slope (velocity), standard error, t-value, p-value, R-squared, standardized coefficient, percent change, and total change. Beta regression (`regression_type = "beta"`) handles [0,1]-bounded transition probabilities properly via Smithson-Verkuilen squeeze transform. Falls back to OLS if betareg fails.
2. **GLLA** (`method = "glla"`, Generalized Local Linear Approximation): Embeds the time series in a delay matrix and estimates derivatives via a weight matrix. No statistical inference (no p-values) but captures nonlinear local trends. Controlled by `n_embed` (embedding dimension) and `delta` (time step).
3. **Finite difference** (`method = "finite_diff"`): Simple numerical differentiation — forward, backward, or central differences. Fast, no assumptions, but noisy. No acceleration output.

**Key parameters:**
- `data`: list of square numeric matrices (transition matrices over time), OR raw long-format data with `time_col`/`id_col`/`action_col`
- `method`: `"regression"` (default), `"glla"`, `"finite_diff"`
- `regression_type`: `"ols"` (default) or `"beta"` (requires betareg in Suggests)
- `weights`: optional numeric vector of time-point weights for weighted regression
- `order`: polynomial order — 1 = linear (velocity only), 2 = quadratic (velocity + acceleration)
- `delta`: time step between matrices (default 1)
- `scaling`: optional post-hoc scaling of input matrices (`"minmax"`, `"max"`, `"rank"`, `"normalize"`)

**Returns `tna_velocity` S3 object with:**
- `velocity_matrix`: mean velocity per edge (p x p matrix, positive = strengthening)
- `velocity_series`: list of velocity matrices at each time point
- `acceleration_matrix`: mean acceleration per edge (NULL for finite_diff)
- `acceleration_series`: list of acceleration matrices at each time point
- `smoothed_matrices`: fitted/smoothed transition matrices (regression/GLLA only)
- `original_matrices`: the input matrices
- `edge_stats`: data.frame with per-edge regression statistics (regression only): from, to, slope, se, t_value, p_value, r_squared, standardized, pct_change, total_change
- `nodes`, `n_timepoints`, `n_nodes`, `n_edges`, `edges`

**S3 methods:**
- `print()`: shows method, top 5 strengthening/weakening edges, significance summary (regression)
- `summary()`: returns data.frame with velocity, direction, consistency per edge; includes regression stats if available
- `plot(type="network")`: velocity network via cograph::splot — green edges strengthening, red weakening
- `plot(type="series")`: matplot of observed edge weights over time with regression fit lines
- `plot(type="heatmap")`: from x to heatmap colored blue (negative velocity) to red (positive)

**Internal functions (18):**
- `.velocity_regression()`, `.velocity_glla()`, `.velocity_finite_diff()` — core estimators
- `.fit_ols_edge()`, `.fit_beta_edge()` — per-edge regression fitting
- `.glla_weight_matrix()` — GLLA weight matrix construction
- `.velocity_logistic()` — logistic regression variant (individual transitions, not matrix aggregation)
- `.detect_and_extract_matrices()` — auto-detect input type (list of matrices, netobject list, raw data)
- `.validate_matrix_list()` — validate matrix dimensions/names
- `.matrices_from_raw_data()` — build transition matrices from long-format data by time window
- `.plot_velocity_network()`, `.plot_velocity_series()`, `.plot_velocity_heatmap()` — plot dispatch functions

**Why sidelined:**
Velocity analysis operates on a fundamentally different input (list of matrices over time) than the rest of Nestimate (single dataset → single network). It needs its own cross-validation, its own integration tests, and possibly its own API design review before being included in a CRAN release. The regression internals (OLS, beta, logistic) are complex and the GLLA implementation needs validation against published reference implementations.

**To restore:**
1. Move files back: `mv sidelined/velocity_tna.R R/` and `mv sidelined/test-velocity_tna.R tests/testthat/`
2. Move man pages back: `mv sidelined/*.tna_velocity.Rd sidelined/velocity_tna.Rd man/`
3. Run `devtools::document()` to rebuild NAMESPACE
4. Run `devtools::test()` to verify tests pass
5. Remove references to velocity in CLAUDE.md class table and FEATURES.md if they were removed
