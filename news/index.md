# Changelog

## Nestimate 1.1.4

This version collects everything since 0.9.24, the last release on CRAN.
Some results are refined and some functions now ask you to say what your
columns are; see “Refined results” and “What works differently”.

### New

- **Latency networks.** `build_network(method = "latency")` builds a
  reply network in which faster replies carry more weight (Templeton et
  al. 2022).
- **Stability landscapes.**
  [`stability_landscape()`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  draws each state as a cup whose depth shows how long people stay in
  it.
  [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  follows the depths over time or along a covariate, and
  [`landscape_surface()`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  draws them along a measure of stability.
- **Sequence patterns and indices.**
  [`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
  finds the patterns in sequences and tests them between groups.
  [`sequence_metrics()`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  computes 23 indices for each sequence (entropy, turbulence,
  complexity, spells and others) and compares them between groups.
  [`pattern_outcome()`](https://pak.dynasite.org/Nestimate/reference/pattern_outcome.md)
  relates patterns to an outcome, and
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)
  returns the sequences behind any input.
- **Phases.**
  [`mark_phases()`](https://pak.dynasite.org/Nestimate/reference/mark_phases.md)
  splits each person’s events into phases, so
  `build_network(group = "phase")` gives one network per phase.
- **Nested data.**
  [`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md)
  tests counts nested in people or teams, and
  [`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
  and
  [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  take `nested` and `actor`.
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
  [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  and
  [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  take `actor` to resample whole people.
- **Entropy over time.**
  [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  reads sequences as well as event logs, and its plot shows a bootstrap
  band around each trend.
- **Sequence plots by period.**
  [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  can show actions per day, week or month (`periods`), and both
  [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  and
  [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md)
  can drop rare columns (`trim`).
- **tna models** can be passed wherever a network or sequences are
  accepted.
- **New dataset** `daily_actions`: the learning actions of 30 students
  over seven days (simulated).
- Smaller additions:
  [`node_clusters()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md),
  [`label_predictions()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  and
  [`label_scores()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md);
  [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md) handles
  missing values and constant variables;
  [`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
  and
  [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  read a group of networks; the `build_*()` shortcuts take `group`,
  `scaling`, `threshold` and `params`;
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  reads edge lists;
  [`cooccurrence()`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md)
  and
  [`association_rules()`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  read 0/1 tables (`format = "binary"`);
  [`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  marks significant differences in its difference plot; `window_size` in
  `build_network(params = )` now works as in tna.

### Refined results

These functions give refined results since 0.9.24.

- [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  reported the structure test as the strength test and the other way
  round, and [`summary()`](https://rdrr.io/r/base/summary.html) swapped
  the in- and out-strength spreads.
- [`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  drew temporal edges in the wrong direction (from outcome to
  predictor).
- [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  effect sizes changed with `min_freq`.
- [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  flagged anomalies in first-order data, which has none; it now uses a
  strength-preserving null. `build_hon(method = "hon")` now matches
  pyHON, and `method = "hon+"` missed some rules.
- [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  and `compare_networks(test = "permutation")` accepted two networks
  estimated in different ways (say Pearson and Spearman) and reported
  the difference as a difference between the groups. They now refuse.
- [`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md),
  [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  and
  [`loading_stability()`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  ignored the network’s own estimator settings.
- Numbers were corrected in centrality (`BetweennessRSP`, betweenness
  ties, min-max normalisation),
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
  [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md),
  [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md),
  [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  (ties),
  [`evaluate_links()`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md),
  and the smoothing of
  [`landscape_surface()`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  and
  [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md).

### What works differently

- **Say which column is which.** Nestimate no longer guesses a column’s
  role from its name. Give `action`, `actor` and `time` for an event
  log, or `state_cols` for wide data, and `format = "onehot"` for 0/1
  code columns. Extra columns in wide data (metadata, `actor`, `time`,
  `order`) are no longer read as states. Times must include a four-digit
  year.
- **Results are tables.** These used to return lists or wide tables and
  now return one data frame, one row per item, with a `group` column
  when groups are compared: `sequence_metrics(group = )`,
  `sequence_patterns(group = )`,
  [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  on several networks, and
  [`summary()`](https://rdrr.io/r/base/summary.html) of a grouped
  network,
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
  [`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md),
  [`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md),
  [`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  and [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md).
  `format = "wide"` gives the old layout; `table =` picks the second
  table where there is one. Printing shows the first rows,
  [`summary()`](https://rdrr.io/r/base/summary.html) all of them.
- **Printing a network** lists its strongest edges instead of the whole
  weight matrix.
- **New defaults.**
  [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  and
  [`markov_order_test()`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md)
  use `dof = "paths"` (`dof = "observed"` gives the old values).
  [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  compares only the edges a transition network has (`target = "all"`
  gives the old values), and it and
  [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  count a subset that fails as below the threshold.
- **Renamed.**
  [`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  columns are `observed_x` and `observed_y`;
  `plot_state_frequencies(colors = )` is `state_colors =` (the old name
  still works, with a warning).
- **Refused inputs.** Functions that read sequences refuse one-hot
  networks and association networks (correlation, partial correlation,
  glasso, Ising, mixed graphical models); resampling functions refuse a
  single sequence and networks built with `level`; latency networks can
  only be bootstrapped; Markov measures refuse correlation networks;
  [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  refuses missing times.
- `scaling = "normalize"` on an undirected network returns a directed
  network when the rows no longer match the columns.
- Arguments of
  [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  after `group` must be named.
- `prepare_onehot()` and `plot_mosaic()` are no longer exported.
- `seed =` no longer changes your session’s random numbers.

### Colours

Every figure colours states the same way: with no colours set, states
get colour-blind-safe Okabe-Ito colours in alphabetical order, and
anything derived from a network keeps its colours. Plots that used
rainbow, viridis or red-green now use these colours, with shapes or line
types where colour alone separated groups.

### Other fixes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  dropped groups with few windows.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  failed on a group of networks.
- `landscape_trend(along = )` failed on averaged metadata.
- Colours set with
  [`set_state_colors()`](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
  were lost after clustering.
- [`subset()`](https://rdrr.io/r/base/subset.html) of
  [`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  tables printed only two decimals.
- Data with a single state failed in
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
  [`path_dependence()`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  and
  [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md).
- `labels =` did not relabel the stored data.
- [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  and
  [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md)
  refused a data.table or tibble.
- [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  ignored its arguments.
- The pyramid plot of
  [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  cut off its counts, and its printout now shows one row per pattern and
  group.

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
