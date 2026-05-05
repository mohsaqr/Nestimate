# Codex Documentation Rewrite

This folder contains rebuilt, knit-ready documentation for Nestimate. The
documents are written as tutorials: they load data, run the relevant functions,
show objects and plots, and explain how to read the output.

## Tutorials

- `01_network_estimation.Rmd` - transition, frequency, co-occurrence, WTNA,
  grouped networks, psychological networks, and validation.
- `02_sequence_clustering.Rmd` - distance clustering, model choice,
  diagnostics, per-cluster networks, MMM, and covariates.
- `03_sequence_visualization_and_comparison.Rmd` - sequence heatmaps, index
  plots, distribution plots, and k-gram comparison.
- `04_state_frequency_and_mosaic_plots.Rmd` - state composition plots,
  tidy frequency tables, and chi-square mosaics.
- `05_mcml.Rmd` - MCML from event logs, matrix aggregation, layer conversion,
  bootstrap, centrality, and interpretation.
- `06_higher_order_markov_simplicial.Rmd` - Markov passage times, stability,
  MOGen, HON, HYPA, and simplicial analysis.

## Audits

- `audit_clustering/report.Rmd`
- `audit_mcml/report.Rmd`

The audits are technical reviews of the code and tests. They are not tutorials.

