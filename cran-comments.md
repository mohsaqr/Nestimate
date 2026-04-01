## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS (Darwin 25.3.0), R 4.5.2
* R CMD check --as-cran --run-donttest

## Resubmission

This is a resubmission. In this version I have:

* Added authors and year to all DOI references in DESCRIPTION per CRAN
  format requirements (e.g., `Saqr et al. (2025) <doi:...>`).

* Unwrapped `\donttest{}` for all examples that execute in under 5 seconds
  (data utilities, extraction, frequencies, link prediction, association
  rules, centrality, estimator registry). Added small fast toy examples
  (tiny data, iter=10) before `\donttest{}` blocks for slower functions
  (bootstrap, permutation, reliability, centrality stability). Most Rd
  files now have at least one automatically tested example.

## Downstream dependencies

This is a new CRAN submission. The package has no reverse dependencies.

The companion package `cograph` (visualization) uses Nestimate's `cograph_network` class objects but does not formally depend on Nestimate.
