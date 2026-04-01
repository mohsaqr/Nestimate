# Package index

## Network Estimation

Core functions for building networks from data

- [`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
  : Build a Network
- [`estimate_network()`](https://mohsaqr.github.io/Nestimate/reference/estimate_network.md)
  : Estimate a Network (Deprecated)
- [`register_estimator()`](https://mohsaqr.github.io/Nestimate/reference/register_estimator.md)
  : Register a Network Estimator
- [`get_estimator()`](https://mohsaqr.github.io/Nestimate/reference/get_estimator.md)
  : Retrieve a Registered Estimator
- [`list_estimators()`](https://mohsaqr.github.io/Nestimate/reference/list_estimators.md)
  : List All Registered Estimators
- [`remove_estimator()`](https://mohsaqr.github.io/Nestimate/reference/remove_estimator.md)
  : Remove a Registered Estimator
- [`wtna()`](https://mohsaqr.github.io/Nestimate/reference/wtna.md) :
  Window-based Transition Network Analysis

## Higher-Order Networks

Methods for capturing higher-order dependencies

- [`build_hon()`](https://mohsaqr.github.io/Nestimate/reference/build_hon.md)
  : Build a Higher-Order Network (HON)
- [`build_honem()`](https://mohsaqr.github.io/Nestimate/reference/build_honem.md)
  : Build HONEM Embeddings for Higher-Order Networks
- [`build_hypa()`](https://mohsaqr.github.io/Nestimate/reference/build_hypa.md)
  : Detect Path Anomalies via HYPA
- [`build_mogen()`](https://mohsaqr.github.io/Nestimate/reference/build_mogen.md)
  : Build Multi-Order Generative Model (MOGen)
- [`pathways()`](https://mohsaqr.github.io/Nestimate/reference/pathways.md)
  : Extract Pathways from Higher-Order Network Objects
- [`mogen_transitions()`](https://mohsaqr.github.io/Nestimate/reference/mogen_transitions.md)
  : Extract Transition Table from a MOGen Model
- [`path_counts()`](https://mohsaqr.github.io/Nestimate/reference/path_counts.md)
  : Count Path Frequencies in Trajectory Data

## Bootstrap & Inference

Statistical inference for network estimation

- [`bootstrap_network()`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)
  : Bootstrap a Network Estimate
- [`boot_glasso()`](https://mohsaqr.github.io/Nestimate/reference/boot_glasso.md)
  : Bootstrap for Regularized Partial Correlation Networks
- [`permutation_test()`](https://mohsaqr.github.io/Nestimate/reference/permutation_test.md)
  : Permutation Test for Network Comparison

## Reliability & Stability

Assess reliability and stability of network estimates

- [`reliability()`](https://mohsaqr.github.io/Nestimate/reference/reliability.md)
  : Split-Half Reliability for Network Estimates
- [`centrality_stability()`](https://mohsaqr.github.io/Nestimate/reference/centrality_stability.md)
  : Centrality Stability Coefficient (CS-coefficient)

## Clustering & Grouping

Cluster-based and multilevel network analysis

- [`cluster_data()`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md)
  [`cluster_sequences()`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md)
  : Cluster Sequences by Dissimilarity
- [`cluster_summary()`](https://mohsaqr.github.io/Nestimate/reference/cluster_summary.md)
  [`csum()`](https://mohsaqr.github.io/Nestimate/reference/cluster_summary.md)
  : Cluster Summary Statistics
- [`cluster_mmm()`](https://mohsaqr.github.io/Nestimate/reference/cluster_mmm.md)
  : Cluster sequences using Mixed Markov Models
- [`cluster_network()`](https://mohsaqr.github.io/Nestimate/reference/cluster_network.md)
  : Cluster data and build per-cluster networks in one step
- [`build_mcml()`](https://mohsaqr.github.io/Nestimate/reference/build_mcml.md)
  : Build MCML from Raw Transition Data
- [`build_mmm()`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md)
  : Fit a Mixed Markov Model
- [`compare_mmm()`](https://mohsaqr.github.io/Nestimate/reference/compare_mmm.md)
  : Compare MMM fits across different k

## Simplicial Complex Analysis

Topological analysis of networks

- [`build_simplicial()`](https://mohsaqr.github.io/Nestimate/reference/build_simplicial.md)
  : Build a Simplicial Complex
- [`persistent_homology()`](https://mohsaqr.github.io/Nestimate/reference/persistent_homology.md)
  : Persistent Homology
- [`q_analysis()`](https://mohsaqr.github.io/Nestimate/reference/q_analysis.md)
  : Q-Analysis
- [`betti_numbers()`](https://mohsaqr.github.io/Nestimate/reference/betti_numbers.md)
  : Betti Numbers
- [`euler_characteristic()`](https://mohsaqr.github.io/Nestimate/reference/euler_characteristic.md)
  : Euler Characteristic
- [`simplicial_degree()`](https://mohsaqr.github.io/Nestimate/reference/simplicial_degree.md)
  : Simplicial Degree
- [`verify_simplicial()`](https://mohsaqr.github.io/Nestimate/reference/verify_simplicial.md)
  : Verify Simplicial Complex Against igraph

## Data Preparation

Convert and prepare data for network estimation

- [`prepare_data()`](https://mohsaqr.github.io/Nestimate/reference/prepare_data.md)
  : Prepare Event Log Data for Network Estimation
- [`prepare_for_tna()`](https://mohsaqr.github.io/Nestimate/reference/prepare_for_tna.md)
  : Prepare Data for TNA Analysis
- [`action_to_onehot()`](https://mohsaqr.github.io/Nestimate/reference/action_to_onehot.md)
  : Convert Action Column to One-Hot Encoding
- [`prepare_onehot()`](https://mohsaqr.github.io/Nestimate/reference/prepare_onehot.md)
  : Import One-Hot Encoded Data into Sequence Format
- [`wide_to_long()`](https://mohsaqr.github.io/Nestimate/reference/wide_to_long.md)
  : Convert Wide Sequences to Long Format
- [`long_to_wide()`](https://mohsaqr.github.io/Nestimate/reference/long_to_wide.md)
  : Convert Long Format to Wide Sequences
- [`convert_sequence_format()`](https://mohsaqr.github.io/Nestimate/reference/convert_sequence_format.md)
  : Convert Sequence Data to Different Formats

## Utilities

Helper functions and extractors

- [`predictability()`](https://mohsaqr.github.io/Nestimate/reference/predictability.md)
  : Compute Node Predictability
- [`frequencies()`](https://mohsaqr.github.io/Nestimate/reference/frequencies.md)
  : Sequence Data Conversion Functions
- [`state_frequencies()`](https://mohsaqr.github.io/Nestimate/reference/state_frequencies.md)
  : Compute State Frequencies from Trajectory Data
- [`aggregate_weights()`](https://mohsaqr.github.io/Nestimate/reference/aggregate_weights.md)
  [`wagg()`](https://mohsaqr.github.io/Nestimate/reference/aggregate_weights.md)
  : Aggregate Edge Weights
- [`centrality()`](https://mohsaqr.github.io/Nestimate/reference/centrality.md)
  : Compute Centrality Measures for a Network
- [`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md)
  : Convert cluster_summary to tna Objects
- [`extract_edges()`](https://mohsaqr.github.io/Nestimate/reference/extract_edges.md)
  : Extract Edge List with Weights
- [`extract_initial_probs()`](https://mohsaqr.github.io/Nestimate/reference/extract_initial_probs.md)
  : Extract Initial Probabilities from Model
- [`extract_transition_matrix()`](https://mohsaqr.github.io/Nestimate/reference/extract_transition_matrix.md)
  : Extract Transition Matrix from Model

## Other Methods

Additional network estimation approaches

- [`build_gimme()`](https://mohsaqr.github.io/Nestimate/reference/build_gimme.md)
  : GIMME: Group Iterative Multiple Model Estimation

## Data

Example datasets

- [`human_ai`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`human_ai_cat`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`human_ai_super`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`human_detailed`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`human_cat`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`human_super`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`ai_detailed`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`ai_cat`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`ai_super`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`human_wide`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  [`ai_wide`](https://mohsaqr.github.io/Nestimate/reference/vibcoding-data.md)
  : Human-AI Vibe Coding Interaction Data
- [`human_ai_edges`](https://mohsaqr.github.io/Nestimate/reference/human_ai_edges.md)
  : Human-AI Vibe Coding Edge List
- [`srl_strategies`](https://mohsaqr.github.io/Nestimate/reference/srl_strategies.md)
  : Self-Regulated Learning Strategy Frequencies
- [`learning_activities`](https://mohsaqr.github.io/Nestimate/reference/learning_activities.md)
  : Online Learning Activity Indicators
- [`group_regulation_long`](https://mohsaqr.github.io/Nestimate/reference/group_regulation_long.md)
  : Group Regulation in Collaborative Learning (Long Format)

## S3 Methods

Print, summary, and plot methods for package classes

- [`print(`*`<boot_glasso>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.boot_glasso.md)
  : Print Method for boot_glasso
- [`print(`*`<mcml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.mcml.md)
  : Print Method for mcml
- [`print(`*`<mmm_compare>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.mmm_compare.md)
  : Print Method for mmm_compare
- [`print(`*`<nestimate_data>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.nestimate_data.md)
  : Print Method for nestimate_data
- [`print(`*`<net_bootstrap>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_bootstrap.md)
  : Print Method for net_bootstrap
- [`print(`*`<net_bootstrap_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_bootstrap_group.md)
  : Print Method for net_bootstrap_group
- [`print(`*`<net_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_clustering.md)
  : Print Method for net_clustering
- [`print(`*`<net_gimme>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_gimme.md)
  : Print Method for net_gimme
- [`print(`*`<net_hon>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_hon.md)
  : Print Method for net_hon
- [`print(`*`<net_honem>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_honem.md)
  : Print Method for net_honem
- [`print(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_hypa.md)
  : Print Method for net_hypa
- [`print(`*`<net_mmm>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mmm.md)
  : Print Method for net_mmm
- [`print(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mogen.md)
  : Print Method for net_mogen
- [`print(`*`<net_permutation>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_permutation.md)
  : Print Method for net_permutation
- [`print(`*`<net_permutation_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_permutation_group.md)
  : Print Method for net_permutation_group
- [`print(`*`<net_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_reliability.md)
  : Print Method for net_reliability
- [`print(`*`<net_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_stability.md)
  : Print Method for net_stability
- [`print(`*`<netobject>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.netobject.md)
  : Print Method for Network Object
- [`print(`*`<netobject_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.netobject_group.md)
  : Print Method for Group Network Object
- [`print(`*`<netobject_ml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.netobject_ml.md)
  : Print Method for Multilevel Network Object
- [`print(`*`<persistent_homology>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.persistent_homology.md)
  : Print persistent homology results
- [`print(`*`<q_analysis>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.q_analysis.md)
  : Print Q-analysis results
- [`print(`*`<simplicial_complex>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.simplicial_complex.md)
  : Print a simplicial complex
- [`print(`*`<wtna_boot_mixed>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.wtna_boot_mixed.md)
  : Print Method for wtna_boot_mixed
- [`print(`*`<wtna_mixed>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.wtna_mixed.md)
  : Print Method for wtna_mixed
- [`summary(`*`<boot_glasso>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.boot_glasso.md)
  : Summary Method for boot_glasso
- [`summary(`*`<mcml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.mcml.md)
  : Summary Method for mcml
- [`summary(`*`<net_bootstrap>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_bootstrap.md)
  : Summary Method for net_bootstrap
- [`summary(`*`<net_bootstrap_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_bootstrap_group.md)
  : Summary Method for net_bootstrap_group
- [`summary(`*`<net_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_clustering.md)
  : Summary Method for net_clustering
- [`summary(`*`<net_gimme>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_gimme.md)
  : Summary Method for net_gimme
- [`summary(`*`<net_hon>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_hon.md)
  : Summary Method for net_hon
- [`summary(`*`<net_honem>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_honem.md)
  : Summary Method for net_honem
- [`summary(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_hypa.md)
  : Summary Method for net_hypa
- [`summary(`*`<net_mmm>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_mmm.md)
  : Summary Method for net_mmm
- [`summary(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_mogen.md)
  : Summary Method for net_mogen
- [`summary(`*`<net_permutation>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_permutation.md)
  : Summary Method for net_permutation
- [`summary(`*`<net_permutation_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_permutation_group.md)
  : Summary Method for net_permutation_group
- [`summary(`*`<net_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_stability.md)
  : Summary Method for net_stability
- [`summary(`*`<wtna_boot_mixed>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.wtna_boot_mixed.md)
  : Summary Method for wtna_boot_mixed
- [`plot(`*`<boot_glasso>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.boot_glasso.md)
  : Plot Method for boot_glasso
- [`plot(`*`<mcml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.mcml.md)
  : Plot Method for mcml
- [`plot(`*`<mmm_compare>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.mmm_compare.md)
  : Plot Method for mmm_compare
- [`plot(`*`<net_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_clustering.md)
  : Plot Sequence Clustering Results
- [`plot(`*`<net_gimme>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_gimme.md)
  : Plot Method for net_gimme
- [`plot(`*`<net_honem>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_honem.md)
  : Plot Method for net_honem
- [`plot(`*`<net_mmm>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mmm.md)
  : Plot Method for net_mmm
- [`plot(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mogen.md)
  : Plot Method for net_mogen
- [`plot(`*`<net_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_reliability.md)
  : Plot Method for net_reliability
- [`plot(`*`<net_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_stability.md)
  : Plot Method for net_stability
- [`plot(`*`<persistent_homology>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.persistent_homology.md)
  : Plot Persistent Homology
- [`plot(`*`<q_analysis>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.q_analysis.md)
  : Plot Q-Analysis
- [`plot(`*`<simplicial_complex>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.simplicial_complex.md)
  : Plot a Simplicial Complex
