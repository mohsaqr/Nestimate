# Changelog

## Nestimate 0.7.3

### New features

- [`certainty()`](https://saqr.me/Nestimate/reference/certainty.md) —
  analytic Bayesian counterpart of
  [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
  for transition networks. Models each state’s outgoing transitions as a
  Dirichlet-Multinomial process (Jeffreys prior) and returns posterior
  mean, sd, credible interval and a stability decision per edge in
  closed form (no resampling). Returns the exact `net_bootstrap` object
  layout and carries class `c("net_certainty", "net_bootstrap")`, so it
  is a drop-in: every `net_bootstrap` method works on it. Completes the
  assessment trio certainty / stability (`bootstrap_network`) /
  reliability (`reliability`).

### Enhancements

- [`bayes_compare()`](https://saqr.me/Nestimate/reference/bayes_compare.md)
  results are now 100% compatible with the
  [`permutation()`](https://saqr.me/Nestimate/reference/permutation.md)
  format: the object carries class `c("net_bayes", "net_permutation")`
  with all `net_permutation` slots (`diff_sig`, `p_values`,
  `effect_size`, `iter`, `alpha`, `paired`, `adjust`), and its `summary`
  is a superset of `summary.net_permutation`
  (`from, to, weight_x, weight_y, diff, effect_size, p_value, sig` plus
  the Bayesian extras
  `count_x, count_y, ci_lower, ci_upper, ci_width, pd`). A
  [`bayes_compare()`](https://saqr.me/Nestimate/reference/bayes_compare.md)
  result is now a drop-in wherever a `net_permutation` is consumed.

## Nestimate 0.7.2

### New features

- [`bayes_compare()`](https://saqr.me/Nestimate/reference/bayes_compare.md)
  — Bayesian Dirichlet-Multinomial comparison of two transition
  networks, a complement to
  [`permutation()`](https://saqr.me/Nestimate/reference/permutation.md).
  Models each source state’s outgoing transitions as a
  Dirichlet-Multinomial process (Jeffreys prior) and returns, per edge,
  a posterior mean difference, a credible interval, the probability of
  direction (`pd`) and its two-sided p-equivalent. Adds
  `print`/`summary`/`plot` methods and `netobject_group` dispatch
  (all-pairwise or matched). Method source: Johnston & Jendoubi (2026),
  *How Delivery Mode Reshapes Resource Engagement: A Bayesian
  Differential Network Analysis*, TNA Workshop 2026.

## Nestimate 0.7.1

### New features

- [`as_networks()`](https://saqr.me/Nestimate/reference/as_networks.md)
  — promote a
  [`build_mcml_pc()`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
  result into a `netobject_group` (the psychometric-network counterpart
  of [`as_tna()`](https://saqr.me/Nestimate/reference/as_tna.md)).
  Singleton clusters with no within-network are dropped with a warning;
  an existing `netobject_group` passes through unchanged.

### Documentation

- Vignettes and articles now call package verbs directly instead of
  hand-assembled base-R subsetting rituals:
  [`markov_order_test()`](https://saqr.me/Nestimate/reference/markov_order_test.md)
  reads sequences straight from a fitted network
  (`markov_order_test(net)`); HYPA anomaly tables use
  `summary(hypa, order_by = "ratio")`; higher-order pathways use
  `pathways(hon, top = )`; grouped-clustering inspection uses
  [`cluster_diagnostics()`](https://saqr.me/Nestimate/reference/cluster_diagnostics.md).

## Nestimate 0.7.0

### New features (experimental)

- [`build_mcml_pc()`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
  — MCML aggregation for psychometric networks (cor / pcor /
  EBICglasso). Five aggregation methods with explicitly different
  statuses: `"average"` (descriptive block-mean; works without raw
  data), `"composite"` (cluster scores re-estimated with the chosen
  estimator — a genuine cluster-level network), `"loadings"` (composites
  weighted by mean within-cluster connection strength — Nestimate’s own
  weighting, not an EGA reimplementation), `"rv"` (Escoufier’s RV matrix
  correlation between blocks), and `"canonical"` (first canonical
  correlation — the upper bound for composite methods). Within-cluster
  networks re-estimated by default (`within = "reestimate"`) since a
  pcor submatrix is not the subsystem’s pcor network. Item diagnostics
  in `$loadings`: signed loadings (reverse-keyed items detected via the
  leading eigenvector of the within-block matrix and flipped in
  composites), cross-cluster strengths, and a `misfit` flag when an item
  is more connected to another cluster than its own (warned). Composites
  tolerate missing data (row-wise renormalized weighted means);
  Composite item weights are selectable via `weighting` — ten built-in
  schemes spanning three views of the cluster: the scale as scored
  (`"equal"`, `"item_total"`), the network’s view (`"strength"`,
  `"eigen"`, `"closeness"`, `"betweenness"`, `"expected_influence"`,
  `"specificity"` — the misfit margin as a weighting, zeroing items that
  belong as much to another cluster), and the latent-variable view
  (`"pca"`, `"factor"`); plus fully custom weighting via a named numeric
  vector or a `function(W_block, data_block, nodes)`.
  `aggregation = "loadings"` is the alias for composite + strength. The
  `"factor"` weighting exposes its extraction method via `fa_method`:
  `"ml"` (factanal), `"paf"` (iterated principal axis), `"minres"`
  (ULS), or `"cfa"` (one-factor lavaan model; with
  `cor_method = "polychoric"` the categorical DWLS factor model) — all
  operating on the `cor_method`-consistent correlation structure.
  Reverse-keyed handling works under every sign-carrying scheme
  (item-total correlations are computed on eigen-sign-pre-oriented
  columns so a reversed member cannot contaminate small clusters).
  `cor_method = "polychoric"` (via lavaan) supports ordinal items;
  `id_col` drops identifier columns so
  `convert_sequence_format(format = "frequency")` actor-profiles feed
  the function directly (the within-person co-occurrence view of event
  data). Returns class `mcml_pc` (macro + within netobjects, all
  undirected) with print/summary/plot; the composite/loadings macro is a
  full netobject, so
  [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
  [`vertex_bootstrap()`](https://saqr.me/Nestimate/reference/vertex_bootstrap.md),
  and
  [`vertex_compare()`](https://saqr.me/Nestimate/reference/vertex_compare.md)
  apply to it directly.
  [`cograph::plot_mcml()`](https://sonsoles.me/cograph/reference/plot_mcml.html)
  (\>= 2.3.8) renders the two-layer undirected MCML view. Experimental:
  API and formulas may change.
- [`loading_stability()`](https://saqr.me/Nestimate/reference/loading_stability.md)
  — case-bootstrap stability of the
  [`build_mcml_pc()`](https://saqr.me/Nestimate/reference/build_mcml_pc.md)
  composite weights (percentile CIs, sign-flip rates), with print and
  forest-style plot.

## Nestimate 0.6.5

### New features

- [`vertex_bootstrap()`](https://saqr.me/Nestimate/reference/vertex_bootstrap.md)
  — Snijders & Borgatti (1999) vertex bootstrap for network-level
  statistics (density, mean weight, strength centralization, weighted
  reciprocity, plus custom `statistic_fn`). Needs only the weight
  matrix, so it works on data-less netobjects
  ([`build_mlvar()`](https://saqr.me/Nestimate/reference/build_mlvar.md)
  constituents, `as_tna(mcml)` elements, plain matrices) where
  [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
  cannot run. Returns a tidy one-row-per-statistic
  `net_vertex_bootstrap` with print/summary/plot. Self-loops are
  preserved (diagonal carries the resampled vertex’s own self-weight);
  undirected replicates stay symmetric.
- [`vertex_compare()`](https://saqr.me/Nestimate/reference/vertex_compare.md)
  — the Snijders & Borgatti two-network test the vertex bootstrap was
  originally proposed for: z-tests and normal-approximation CIs for
  differences in network-level statistics between two networks
  (netobjects, matrices, or precomputed `net_vertex_bootstrap` objects).
  Tidy `net_vertex_comparison` result with print/summary/plot (forest
  plot of differences).
- [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
  and
  [`vertex_bootstrap()`](https://saqr.me/Nestimate/reference/vertex_bootstrap.md)
  gain `ci_method = c("percentile", "basic")`: basic intervals (Davison
  & Hinkley 1997, eq. 5.6) reflect the percentile bounds around the
  observed estimate, correcting first-order bootstrap bias. Default
  remains `"percentile"`.

## Nestimate 0.6.4

### Bug fixes

- [`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md)
  (sequence and edge-list paths) now records the *effective*
  directedness in `$meta$directed`: `FALSE` when
  `type = "cooccurrence"`, whose weights are symmetrized, instead of
  echoing the `directed` argument unchanged. Renderers that auto-detect
  directedness (e.g.,
  [`cograph::plot_mcml()`](https://sonsoles.me/cograph/reference/plot_mcml.html)
  with `directed = NULL`) now draw co-occurrence MCML objects as
  undirected networks automatically.

## Nestimate 0.6.0

CRAN release: 2026-05-31

### New features

- [`magnitude_difference()`](https://saqr.me/Nestimate/reference/magnitude_difference.md)
  compares the frequency (FTNA) and probability (TNA) views of a
  transition network and quantifies the per-edge discrepancy on a common
  scale, with five metrics, four scalings, and two polar
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) portraits
  (stacked and circular).
- Full persistent homology with a Vietoris-Rips filtration
  ([`persistent_homology()`](https://saqr.me/Nestimate/reference/persistent_homology.md),
  `build_simplicial(type = "vr")`) plus diagram tools
  [`bottleneck_distance()`](https://saqr.me/Nestimate/reference/bottleneck_distance.md)
  and
  [`persistence_landscape()`](https://saqr.me/Nestimate/reference/persistence_landscape.md).
- Network comparison:
  [`compare_model()`](https://saqr.me/Nestimate/reference/compare_model.md)
  (with `netobject_group` dispatch),
  [`summary.netobject()`](https://saqr.me/Nestimate/reference/summary.netobject.md),
  [`plot.net_comparison()`](https://saqr.me/Nestimate/reference/plot.net_comparison.md),
  and
  [`rename_models()`](https://saqr.me/Nestimate/reference/rename_models.md)
  for relabelling grouped network objects.

### Documentation

- The pkgdown reference index now lists every exported function;
  [`magnitude_difference()`](https://saqr.me/Nestimate/reference/magnitude_difference.md),
  [`casedrop_reliability()`](https://saqr.me/Nestimate/reference/casedrop_reliability.md),
  [`build_hypergraph()`](https://saqr.me/Nestimate/reference/build_hypergraph.md),
  [`hypergraph_measures()`](https://saqr.me/Nestimate/reference/hypergraph_measures.md),
  and
  [`cluster_data()`](https://saqr.me/Nestimate/reference/cluster_data.md)
  were previously absent.

### Packaging

- Removed the `Remotes:` field; `cograph` and `tna` are available from
  CRAN, so no non-CRAN source pin is needed.

## Nestimate 0.5.1

### Audit follow-up (clustering + MCML)

Followed `codex_docs/audit_clustering` and `codex_docs/audit_mcml`
recommendations across two modules. Eleven of thirteen findings
addressed; two deferred pending design decisions on numeric semantics
(`directed = FALSE` raw-data MCML, MMM first-non-NA initial state).

#### Bug fixes

- [`cluster_network()`](https://saqr.me/Nestimate/reference/cluster_network.md)
  now forwards distance-clustering arguments (`na_syms`, `weighted`,
  `lambda`, `seed`, `q`, `p`, `covariates`) to
  [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md)
  instead of silently passing them to
  [`build_network()`](https://saqr.me/Nestimate/reference/build_network.md).
  The split runs on caller `...` only — netobject `build_args` continue
  to flow only to the
  [`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
  step, protecting attention-method (`atna`) network history from being
  re-routed to weighted Hamming. (audit_clustering
  [\#1](https://github.com/mohsaqr/Nestimate/issues/1))
- [`.auto_detect_clusters()`](https://saqr.me/Nestimate/reference/dot-auto_detect_clusters.md)
  (used by
  [`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md)
  and
  [`cluster_summary()`](https://saqr.me/Nestimate/reference/cluster_summary.md))
  now requires `node_groups` to carry a node identifier column when
  shaped as a data.frame, or be a named atomic vector keyed by node
  label. Previously, a bare `cluster`-only data.frame was read
  positionally — silently mis-assigning nodes whenever `node_groups`
  rows were in a different order than `x$nodes`. (audit_mcml
  [\#1](https://github.com/mohsaqr/Nestimate/issues/1))
- [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md)
  now rejects all-missing input early with a clear message instead of
  failing indirectly downstream in pam/hclust. (audit_clustering
  [\#4](https://github.com/mohsaqr/Nestimate/issues/4))

#### New parameters

- `compare_mmm(return_fits = FALSE)` — when `TRUE`, the fitted `net_mmm`
  models are attached as `attr(result, "fits")` keyed by `k`, so users
  can pick the chosen model without re-running EM. Default behaviour
  unchanged. (audit_clustering
  [\#6](https://github.com/mohsaqr/Nestimate/issues/6))

#### Improvements

- [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md)
  validation messages now name the offending argument
  (`"'k' must be at least 2 (got k = 1)"`) rather than dumping the
  failing predicate. Top-level type checks switched to named-condition
  [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) for the same
  reason. (audit_clustering
  [\#2](https://github.com/mohsaqr/Nestimate/issues/2))

#### Documentation

- [`summary.mcml()`](https://saqr.me/Nestimate/reference/summary.mcml.md)
  roxygen corrected — was claiming a printing side effect that doesn’t
  exist. (audit_mcml
  [\#5](https://github.com/mohsaqr/Nestimate/issues/5))
- [`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md)
  `clusters = "<col>"` mode now documents its narrow contract: assigns
  each row’s group label to both endpoints, so it only makes sense for
  within-group edge lists. (audit_mcml
  [\#2](https://github.com/mohsaqr/Nestimate/issues/2))
- [`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md)
  `method` parameter doc now steers raw sequence / event-log inputs to
  `"sum"`, since the function counts observed transitions. Other methods
  are for weighted edge lists or pre-existing matrices. (audit_mcml
  [\#4](https://github.com/mohsaqr/Nestimate/issues/4))
- [`as_tna.mcml()`](https://saqr.me/Nestimate/reference/as_tna.md)
  “Excluded Clusters” section corrected — drop emits a
  [`warning()`](https://rdrr.io/r/base/warning.html) (was claimed
  silent) and only fires for `relative` method (was claimed
  unconditional). (audit_mcml
  [\#6](https://github.com/mohsaqr/Nestimate/issues/6))
- [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md)
  `na_syms` doc adds an explicit “Missing-value distance rule”
  subsection: NA becomes a comparable sentinel state, not pairwise
  deletion. (audit_clustering
  [\#3](https://github.com/mohsaqr/Nestimate/issues/3))
- [`build_mmm()`](https://saqr.me/Nestimate/reference/build_mmm.md) adds
  an “Initial states” section explaining first-column-verbatim init and
  that build_mmm does NOT honor build_clusters-style `na_syms` — only
  actual `NA` cells become NA inits. (audit_clustering
  [\#5](https://github.com/mohsaqr/Nestimate/issues/5), doc-only path)

#### Tests

- +12 new tests pinning the corrected contracts and the documented edge
  cases (misordered `node_groups` alignment, label propagation through
  [`state_distribution()`](https://saqr.me/Nestimate/reference/state_distribution.md),
  [`as_tna.mcml()`](https://saqr.me/Nestimate/reference/as_tna.md)
  drop-warning fixture, MMM first-column NA behaviour, and the four-way
  [`cluster_network()`](https://saqr.me/Nestimate/reference/cluster_network.md)
  arg-routing contract). Full sweep: 1628 / 1628 pass, 0 fail.

## Nestimate 0.5.0

### Bug fixes

- `.extract_edges_from_matrix()` no longer drops the diagonal.
  Netobjects built via `.wrap_netobject()` (and therefore everything
  from
  [`build_network()`](https://saqr.me/Nestimate/reference/build_network.md),
  [`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md),
  [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
  [`build_mmm()`](https://saqr.me/Nestimate/reference/build_mmm.md),
  [`wtna()`](https://saqr.me/Nestimate/reference/wtna.md),
  [`as_tna()`](https://saqr.me/Nestimate/reference/as_tna.md)) now have
  `$edges` containing every non-zero matrix entry, including self-loops.
  Previously `$weights` and `$edges` were silently inconsistent on any
  matrix with a non-zero diagonal, causing downstream consumers
  (e.g. [`cograph::centrality()`](https://sonsoles.me/cograph/reference/centrality.html)
  on an MCML macro) to under-count node degree by 2.

### New features

- [`plot_state_frequencies()`](https://saqr.me/Nestimate/reference/plot_state_frequencies.md)
  — native S3 generic for state-frequency plots across `netobject`,
  `netobject_group`, `mcml`, and `htna`. Defaults to a marimekko
  (mosaic) layout where column widths reflect per-group totals and
  segment heights reflect within-group state proportions; also supports
  a colored-bars style and a per-group faceted marimekko. Uses the
  package Okabe-Ito palette throughout.
- [`plot_mosaic()`](https://saqr.me/Nestimate/reference/plot_mosaic.md)
  — exported low-level marimekko primitive built on `geom_rect()` with
  cumulative-width / cumulative-height geometry. Reusable for any tidy
  `data.frame(group, state, weight)` input.

## Nestimate 0.4.4

### Bug fixes

- [`passage_time()`](https://saqr.me/Nestimate/reference/passage_time.md)
  and
  [`markov_stability()`](https://saqr.me/Nestimate/reference/markov_stability.md)
  now raise an explicit error naming the dead state when a
  transition-matrix row sums to zero, instead of silently propagating
  `NaN` through `eigen`/`solve`. Zero rows mean the chain is not
  ergodic; mean first passage times are undefined. Shared helper
  `.mpt_normalize_rows()` factored out of both entry points.
- `.prepare_association_input()` no longer hard-rejects non-square
  numeric matrices. For association methods (glasso, pcor, cor) the
  netobject’s `$data` slot is a numeric matrix (not a data.frame). Any
  downstream caller that row-subsetted `$data` and re-invoked the
  estimator
  ([`centrality_stability()`](https://saqr.me/Nestimate/reference/centrality_stability.md),
  [`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
  `reliability()`) was silently producing NULL centralities caught by
  `tryCatch`, which surfaced as an “all centrality measures have zero
  variance” warning or all-`NaN` correlations. The matrix branch now
  recognises non-square input as raw observation data and recursively
  re-enters through the data-frame branch. Square symmetric matrices
  (pre-computed correlation / covariance) still go through the
  symmetric-matrix path with the symmetry check intact.

### New parameters

- [`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
  gains `state_cols` and `metadata_cols` parameters (both default
  `NULL`). Explicit overrides for the state-vs-metadata column
  classifier, which previously used a “values-in-nodes” heuristic that
  silently misclassifies metadata columns whose values coincide with
  node labels (e.g. a `condition` column with levels `"A","B","C"` when
  nodes are `"A","B","C"`). Validation: error on overlap between the two
  vectors, error on column names not present in the input data.
  Forwarded through the `group = ...` recursive dispatch so per-group
  calls honour the override.

### Removed

- `plot.net_link_prediction()` and `plot.mcml()` removed. Nestimate is a
  computation engine — visualization is the user’s concern. Previously
  both methods called `cograph::` directly, violating the stated
  dependency invariant (Nestimate -\> cograph direction forbidden).
  Users call `cograph::splot(net)` or `cograph::plot_mcml(fit)`
  directly.

### Documentation

- [`wtna()`](https://saqr.me/Nestimate/reference/wtna.md) `@param type`
  now flags that `type = "relative"` combined with
  `method = "cooccurrence"` produces an asymmetric matrix (conditional
  co-occurrence given row state), not a symmetric undirected weight
  matrix. Use `type = "frequency"` if symmetric counts are required.

### Testing infrastructure

- New numerical-equivalence tests (gated by
  `NESTIMATE_EQUIV_TESTS=true`): `test-equiv-permutation.R`
  (vs. [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html) +
  hand-coded base-R permutation loop), `test-equiv-mlvar.R`
  (vs. [`mlVAR::mlVAR`](https://rdrr.io/pkg/mlVAR/man/mlVAR.html) at
  machine precision), `test-equiv-association-rules.R` (vs.
  [`arules::apriori`](https://rdrr.io/pkg/arules/man/apriori.html)),
  `test-equiv-link-prediction.R` (vs. clean-room matrix algebra +
  [`igraph::similarity`](https://r.igraph.org/reference/similarity.html)),
  `test-equiv-centrality-stability.R`
  (vs. [`bootnet::corStability`](https://rdrr.io/pkg/bootnet/man/corStability.html)).
  Total ~162k per-value comparisons; all within machine precision except
  centrality-stability which uses a documented drop-grid tolerance
  because bootnet uses `igraph` path-based centrality and Nestimate uses
  Floyd-Warshall.
- New HON-family equivalence tests under
  `local_testing_and_equivalence/` validating HON, HONEM, HYPA, MOGen,
  and hypergraph against pathpy 2.2.0 (via reticulate), `BiasedUrn`,
  `RSpectra`, and `HyperG`. Not shipped in the R-package `tests/`
  directory; added to `.Rbuildignore`.
- Branch-matrix coverage added for the four many-mode APIs (`wtna`,
  `bootstrap_network`, `build_clusters`, `sequence_plot`) — systematic
  cross-product tests over all combinations of mode parameters to catch
  regressions where one branch silently diverges.

## Nestimate 0.4.3

CRAN release: 2026-04-20

### CRAN resubmission (addresses incoming-check NOTEs on 0.4.2)

- Tarball now ships `build/vignette.rds` (the vignette index). Previous
  0.4.2 build used `R CMD build --no-build-vignettes`, which preserved
  pre-built `inst/doc/*.html` but stripped the index — CRAN flagged
  “VignetteBuilder field but no prebuilt vignette index.”
- `test-gimme.R` now `skip_on_cran()`. GIMME tests fit a lavaan SEM per
  subject and took ~50s locally (2-3× on Windows), pushing total check
  time to 11 min on win-devel. Full test suite still runs in CI and
  local dev.

## Nestimate 0.4.2

### CRAN resubmission

- Full `--as-cran --run-donttest` audit pass.
- Purged stale `.Rcheck/` and `Meta/` build artifacts from working tree;
  added explicit `^Nestimate\.Rcheck$` and `^\.\.Rcheck$` entries to
  `.Rbuildignore` as belt-and-suspenders against repeat-submission
  contamination.

## Nestimate 0.4.1

### CRAN resubmission

- Pre-built vignettes included in `inst/doc/` as required by CRAN.
- Fixed 301-redirect URLs in README.
- Added `skip_on_cran()` to slow test block to keep check time under 10
  minutes.

## Nestimate 0.4.0

### New functions

- [`build_mlvar()`](https://saqr.me/Nestimate/reference/build_mlvar.md)
  — multilevel VAR networks from ESM/EMA panel data. Estimates temporal
  (directed), contemporaneous (undirected), and between-subjects
  (undirected) networks matching
  [`mlVAR::mlVAR()`](https://rdrr.io/pkg/mlVAR/man/mlVAR.html) at
  machine precision.
- [`build_mmm()`](https://saqr.me/Nestimate/reference/build_mmm.md) /
  [`compare_mmm()`](https://saqr.me/Nestimate/reference/compare_mmm.md)
  — mixture of Markov models via EM, with BIC/AIC/ICL model selection
  and optional covariate regression.
- [`cooccurrence()`](https://saqr.me/Nestimate/reference/cooccurrence.md)
  — standalone co-occurrence network builder supporting 6 input formats
  and 8 similarity methods.
- [`sequence_compare()`](https://saqr.me/Nestimate/reference/sequence_compare.md)
  — k-gram pattern comparison across groups with optional permutation
  testing.
- [`sequence_plot()`](https://saqr.me/Nestimate/reference/sequence_plot.md)
  /
  [`distribution_plot()`](https://saqr.me/Nestimate/reference/distribution_plot.md)
  — base-R sequence index and state distribution plots with clustering
  integration.
- [`build_simplicial()`](https://saqr.me/Nestimate/reference/build_simplicial.md),
  [`persistent_homology()`](https://saqr.me/Nestimate/reference/persistent_homology.md),
  [`q_analysis()`](https://saqr.me/Nestimate/reference/q_analysis.md) —
  topological analysis of networks via simplicial complexes.
- [`nct()`](https://saqr.me/Nestimate/reference/nct.md) — Network
  Comparison Test matching
  [`NetworkComparisonTest::NCT()`](https://rdrr.io/pkg/NetworkComparisonTest/man/NCT.html)
  at machine precision.
- [`build_gimme()`](https://saqr.me/Nestimate/reference/build_gimme.md)
  — group iterative mean estimation for idiographic networks via lavaan.
- [`passage_time()`](https://saqr.me/Nestimate/reference/passage_time.md),
  [`markov_stability()`](https://saqr.me/Nestimate/reference/markov_stability.md)
  — Markov chain passage times and stability analysis.
- [`predict_links()`](https://saqr.me/Nestimate/reference/predict_links.md)
  /
  [`evaluate_links()`](https://saqr.me/Nestimate/reference/evaluate_links.md)
  — link prediction with 6 structural similarity methods.
- [`association_rules()`](https://saqr.me/Nestimate/reference/association_rules.md)
  — Apriori association rule mining from sequences or binary matrices.
- [`predictability()`](https://saqr.me/Nestimate/reference/predictability.md)
  — node predictability for glasso/pcor/cor networks.
- [`build_hon()`](https://saqr.me/Nestimate/reference/build_hon.md),
  [`build_honem()`](https://saqr.me/Nestimate/reference/build_honem.md),
  [`build_hypa()`](https://saqr.me/Nestimate/reference/build_hypa.md),
  [`build_mogen()`](https://saqr.me/Nestimate/reference/build_mogen.md)
  — higher-order network methods (HON, HONEM, HYPA, MOGen) now
  `cograph_network`-compatible.

### New datasets

- `human_long`, `ai_long` — canonical long-format human–AI pair
  programming interaction sequences (10,796 turns, 429 sessions).
- `chatgpt_srl` — ChatGPT-generated SRL scale scores for psychological
  network analysis.
- `trajectories` — 138-student engagement trajectory matrix (15
  timepoints, 3 states).

### API

- [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md),
  [`network_reliability()`](https://saqr.me/Nestimate/reference/network_reliability.md),
  [`permutation()`](https://saqr.me/Nestimate/reference/permutation.md),
  and [`prepare()`](https://saqr.me/Nestimate/reference/prepare.md)
  replace earlier internal names for consistency with the `build_*`
  naming convention.
- `mgm` estimator added (`method = "mgm"`) for mixed continuous +
  categorical data via nodewise lasso, matching
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) at machine
  precision.

### Bug fixes

- [`build_mmm()`](https://saqr.me/Nestimate/reference/build_mmm.md) no
  longer crashes on platforms where
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  returns `NA` (macOS ARM64 CRAN check failure).
- `gimme` convergence filter now correctly handles all typed `NA`
  variants (`NA_character_`, `NA_real_`, etc.).
- `NaN` values in numeric metadata aggregation (all-`NA` sessions)
  normalized to `NA_real_`.
- HYPA p-values corrected; `hypa_score` column renamed to `p_value`.

### CRAN compliance

- `.data` pronoun added to
  [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html).
- [`base::.rowSums()`](https://rdrr.io/r/base/colSums.html) /
  [`base::.colSums()`](https://rdrr.io/r/base/colSums.html) replaced
  with [`rowSums()`](https://rdrr.io/r/base/colSums.html) /
  [`colSums()`](https://rdrr.io/r/base/colSums.html).
- [`dev.new()`](https://rdrr.io/r/grDevices/dev.html) guarded by
  [`interactive()`](https://rdrr.io/r/base/interactive.html) — no side
  effects under knitr or CI.
- Equivalence test files excluded from the built tarball.

### Performance

- `do.call(rbind, ...)` replaced with
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
  in `mcml.R` and `sequence_compare.R`.

## Nestimate 0.3.4

- HYPA: Renamed `hypa_score` column to `p_value` for clarity. Added
  `$over`, `$under`, `$n_over`, `$n_under` fields to `net_hypa` objects.
  Scores are now pre-sorted with anomalous paths first.
- HYPA:
  [`summary.net_hypa()`](https://saqr.me/Nestimate/reference/summary.net_hypa.md)
  now shows over/under-represented paths separately with a configurable
  `n` parameter.
- [`pathways.netobject()`](https://saqr.me/Nestimate/reference/pathways.md):
  New S3 method to extract higher-order pathways directly from a
  netobject (builds HON or HYPA internally).
- [`path_counts()`](https://saqr.me/Nestimate/reference/path_counts.md):
  Now handles NAs in trajectories by stripping them before k-gram
  counting.

## Nestimate 0.2.15

- Preparing for publication

## Nestimate 0.2.0

- Reduced hard dependencies from 6 to 4 Imports (ggplot2, glasso,
  data.table, cluster).
- Removed igraph from Imports —
  [`centrality_stability()`](https://saqr.me/Nestimate/reference/centrality_stability.md)
  and
  [`boot_glasso()`](https://saqr.me/Nestimate/reference/boot_glasso.md)
  now accept a `centrality_fn` parameter for external centrality
  computation.
- Removed tna from Imports — moved to Suggests (only used for input
  class detection).
- Implemented `graphical_var()` from scratch using coordinate descent
  lasso + graphical lasso with EBIC model selection, eliminating the
  graphicalVAR dependency.
- Dropped `ml_graphical_var()` — users should use `mlvar()` for
  multilevel VAR.
- Removed cograph plot wrappers — `plot.netobject()`,
  `plot.net_bootstrap()`, `plot.net_permutation()`, `plot.net_hon()`,
  `plot.net_hypa()` and
  [`as_cograph()`](https://sonsoles.me/cograph/reference/as_cograph.html)
  removed. Users call cograph plotting functions directly on netobjects.
- Added `attention` estimator for decay-weighted transition networks.
- Increased test coverage from 84.5% to 96.1% (2780 tests).
- Passes R CMD check with 0 errors, 0 warnings, 0 notes.

## Nestimate 0.1.0

- Initial release. Split from Saqrlab v0.3.0.
- Core estimation via
  [`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
  with 8 built-in estimators.
- Bootstrap inference
  ([`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md)),
  permutation testing
  ([`permutation()`](https://saqr.me/Nestimate/reference/permutation.md)),
  EBICglasso bootstrap
  ([`boot_glasso()`](https://saqr.me/Nestimate/reference/boot_glasso.md)).
- Higher-order networks: HON, HONEM, HYPA, MOGen.
- GIMME, MCML, multilevel VAR, graphical VAR.
- Temporal network analysis and velocity TNA.
- Dual-class `c("netobject", "cograph_network")` output for cograph
  compatibility.
