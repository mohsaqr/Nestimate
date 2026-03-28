# Changelog

## Nestimate 0.2.15

- Preparing for publication

## Nestimate 0.2.0

- Reduced hard dependencies from 6 to 4 Imports (ggplot2, glasso,
  data.table, cluster).
- Removed igraph from Imports —
  [`centrality_stability()`](https://mohsaqr.github.io/Nestimate/reference/centrality_stability.md)
  and
  [`boot_glasso()`](https://mohsaqr.github.io/Nestimate/reference/boot_glasso.md)
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
  `plot.net_hypa()` and `as_cograph()` removed. Users call cograph
  plotting functions directly on netobjects.
- Added `attention` estimator for decay-weighted transition networks.
- Increased test coverage from 84.5% to 96.1% (2780 tests).
- Passes R CMD check with 0 errors, 0 warnings, 0 notes.

## Nestimate 0.1.0

- Initial release. Split from Saqrlab v0.3.0.
- Core estimation via
  [`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
  with 8 built-in estimators.
- Bootstrap inference
  ([`bootstrap_network()`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)),
  permutation testing
  ([`permutation_test()`](https://mohsaqr.github.io/Nestimate/reference/permutation_test.md)),
  EBICglasso bootstrap
  ([`boot_glasso()`](https://mohsaqr.github.io/Nestimate/reference/boot_glasso.md)).
- Higher-order networks: HON, HONEM, HYPA, MOGen.
- GIMME, MCML, multilevel VAR, graphical VAR.
- Temporal network analysis and velocity TNA.
- Dual-class `c("netobject", "cograph_network")` output for cograph
  compatibility.
