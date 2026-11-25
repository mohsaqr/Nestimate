# Changelog

## Nestimate 1.1.2

Changes since the CRAN release 0.9.24.

### Breaking: results with a new shape

These results are now tidy data frames, one row per observation, with a
`group` column where there are groups. `format = "wide"` returns the old
one-column-per-group table;
[`summary()`](https://rdrr.io/r/base/summary.html) returns every row as
a plain data frame, and [`print()`](https://rdrr.io/r/base/print.html)
shows the first `n` rows.

- `sequence_metrics(group = )`: one row per index and group, with the
  difference, Hedges’ g and its interval, eta squared and p-values.
- `sequence_patterns(group = )`: one row per pattern and group, with
  `group_support`.
- [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  of a group of networks: one data frame with `group` and `state`
  columns; it returned a list of data frames.
- [`summary()`](https://rdrr.io/r/base/summary.html) of a grouped
  network: one row per metric and group.
- [`summary()`](https://rdrr.io/r/base/summary.html) of
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
  [`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  and
  [`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  returns one data frame; `table =` gives the second table (`"pairs"`,
  `"chain"`, `"absorption"`). They returned a list. Grouped results
  stack the table with a `group` column.
- [`summary()`](https://rdrr.io/r/base/summary.html) of
  [`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  returns the table and prints nothing.
- [`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  names its columns `observed_x` and `observed_y`.
- `summary(nct(), table = "edges")` returns the edge table;
  [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md) stores
  the M and S tests as `$structure` and `$strength`.
- [`print()`](https://rdrr.io/r/base/print.html) of a network lists its
  strongest edges instead of the weight matrix.

### Rerun an analysis made with 0.9.24 if it uses

- data whose column roles came from their names, or long-log times not
  in a whole format with a four-digit year. Declare `action`, `actor`,
  `time`, `state_cols`, `metadata_cols` or `cols`; one-hot data needs
  `format = "onehot"` or `codes`.
- wide data with metadata, `actor`, `time` or `order` columns, which
  were read as states.
- [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  or `compare_networks(test = "permutation")` on networks built with
  different estimator settings (say Pearson and Spearman): the
  difference was reported as a difference between the groups. Such pairs
  are now refused.
- [`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md),
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  or
  [`loading_stability()`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md):
  they now use the network’s own estimator settings.
- [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  or
  [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md):
  failed subsets count as below `threshold`, and
  [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  compares only the existing edges of a transition network
  (`target = "all"` gives the earlier values).
- [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  or
  [`markov_order_test()`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md):
  `dof = "paths"` is the default (`dof = "observed"` gives the earlier
  values).
- `build_hon(method = "hon")`, which now follows pyHON, or
  [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md),
  whose null is strength-preserving (the earlier null flagged anomalies
  in first-order data).
- [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md), whose
  `$structure` and `$strength` held each other’s values in 0.9.24, or
  strengths from [`summary()`](https://rdrr.io/r/base/summary.html),
  whose in and out spreads were swapped.
- [`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  temporal edges, which now run from predictor to outcome.
- [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  effect sizes, which depended on `min_freq`.
- [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  ties,
  [`evaluate_links()`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md),
  `scaling = "normalize"` on an undirected network (now directed when
  the result is asymmetric) or the kernel width of
  [`landscape_surface()`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  and
  [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md).
- corrected estimates in centrality (`BetweennessRSP`, betweenness ties,
  min-max normalisation),
  [`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md),
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md),
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md),
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
  [`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md),
  [`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md),
  [`entropy_bayes()`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md),
  [`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md),
  [`q_analysis()`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md),
  [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md),
  [`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md),
  [`bottleneck_distance()`](https://pak.dynasite.org/Nestimate/reference/bottleneck_distance.md),
  [`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md),
  [`mogen_transitions()`](https://pak.dynasite.org/Nestimate/reference/mogen_transitions.md),
  [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  and
  [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md).

### New features

- `build_network(method = "latency")`: a reply network weighted by reply
  speed (Templeton et al. 2022).
- `tna` and `group_tna` models are accepted wherever a network or
  sequences are.
- [`stability_landscape()`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md),
  [`landscape_surface()`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  and
  [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  draw a transition network as a ball-in-cup landscape.
- [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md),
  [`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md),
  [`sequence_metrics()`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  and
  [`pattern_outcome()`](https://pak.dynasite.org/Nestimate/reference/pattern_outcome.md):
  the sequences, patterns and indices behind any input.
  [`sequence_metrics()`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  compares every index between groups, within actors when an actor
  appears in several groups.
- [`mark_phases()`](https://pak.dynasite.org/Nestimate/reference/mark_phases.md)
  splits each actor’s events into `n` phases, for
  `build_network(group = "phase")`.
- [`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md),
  and `nested`, `actor` and `report` in
  [`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
  and
  [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md),
  for counts nested in actors.
- [`node_clusters()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md),
  [`label_predictions()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  and
  [`label_scores()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md).
- `actor` drops whole actors in
  [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  and
  [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md),
  and takes nested levels in
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md).
- [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  reads sequences as well as an event log, and its
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws a
  bootstrap confidence band around each trend.
- [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  gains `periods` (weeks or months), `period_unit` and `fill`; it and
  [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md)
  gain `cols` and a two-part `trim`.
- [`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
  and
  [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  read a group of networks;
  [`summary()`](https://rdrr.io/r/base/summary.html) of
  [`magnitude_difference()`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  returns its edges.
- [`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  stars the tested differences in `plot(type = "difference")`.
- [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md) gains
  `na_method` and accepts missing values and constant variables.
- `build_network(params = list(window_size = ))` counts sequences
  between blocks of positions, as
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
  does (it was ignored).
- The `build_*()` shortcuts take `group`, `scaling`, `threshold` and
  `params` as named arguments;
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  gains `from`, `to`, `weight` and `cols`;
  [`cooccurrence()`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md)
  and
  [`association_rules()`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  gain `format = "binary"`.
- New dataset `daily_actions`.

### Other changes

- [`print()`](https://rdrr.io/r/base/print.html) of
  [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  shows one row per pattern and group, rounded, with a `group` column;
  `summary(format = "long")` returns that table for every pattern. The
  result and [`summary()`](https://rdrr.io/r/base/summary.html) keep
  their shape.
- One colour rule for states: with no colours set, every figure draws
  the states in Okabe-Ito colours dealt in alphabetical order, and
  objects derived from a network keep its palette. The clusters of an
  mcml and the groups of an
  [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  plot take the colours the states leave.
  `plot_state_frequencies(colors = )` is now `state_colors =` (`colors`
  still works, with a warning).
- Plots that used ggplot’s default hues, viridis or green and red use
  Okabe-Ito and the package’s blue-red diverging colours, with shapes or
  linetypes where colour alone told groups apart.
- Sequence verbs refuse one-hot and association networks; resampling
  verbs refuse one-sequence data and networks built with `level`, and
  all but
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  refuse latency networks; Markov measures refuse correlation networks.
- [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  refuses a `time` column with missing values.
- `prepare_onehot()` and `plot_mosaic()` are no longer exported.
- [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  arguments after `group` must be named.
- `seed =` no longer moves the session’s random stream.

### Bug fixes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  dropped a group with few windows; it is now drawn as its windows and
  named in a warning.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  on a group of networks failed.
- [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  stored averaged metadata as a one-dimensional array, so
  `landscape_trend(along = )` failed.
- Colours set with
  [`set_state_colors()`](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
  were lost by
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md),
  [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  and other derived objects, and were matched against the wrong name
  column.
- [`subset()`](https://rdrr.io/r/base/subset.html) of
  [`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  tables printed at two decimals.
- One-state data failed in
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
  [`path_dependence()`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  and
  [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md).
- `labels =` did not relabel stored data.
- [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  and
  [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md)
  refused a `data.table` or tibble.
- [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  did not pass its arguments.
- `build_hon(method = "hon+")` missed some higher-order rules.
- The pyramid plot of
  [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  cut off its counts.

## Nestimate 0.9.24

CRAN release: 2026-10-07

### Documentation

- Fixes to the 0.9.23 help-page merge, found in review:
  - Pages are named after their function again: `build_mcml` and
    `plot_state_frequencies` had been named after a
    [`print()`](https://rdrr.io/r/base/print.html) method, and their
    usage now lists the function before its methods.
  - The descriptions of 51 methods (what a plot draws, what a print
    shows) are restored, as a “Methods” section on each function’s page.
  - Each page’s `...` entry states what each method does with extra
    arguments: methods that reject them say so,
    [`plot.net_bayes()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
    names
    [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
    and the rest say they are ignored.
  - `compare_networks`: the `...` and `labels` entries no longer repeat;
    `net_edge_betweenness` and `net_pruning_details` no longer list a
    non-existent class.
  - Every S3 method has a documented return value.

### Behaviour notes

- The `mgm` estimator rejects `scale = FALSE` with a classed error
  (`nestimate_mgm_unscaled`); this started in 0.9.23 and was not listed
  there. Missing values are handled pairwise, now passed explicitly.

### Tests

- Two more reference-implementation tests (a hand-written linear solve
  for
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
  a manual `multinom` re-computation) and unused reference helpers moved
  to the local suite; section headers that claimed equivalence tests are
  renamed. New tests: mgm on data with missing values, and the
  `scale = FALSE` error.

## Nestimate 0.3.0

CRAN release: 2026-04-08

- First release on CRAN.
