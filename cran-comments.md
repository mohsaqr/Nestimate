## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS (Darwin 25.3.0), R 4.5.2
* R CMD check --as-cran --run-donttest

## Resubmission

This is a resubmission. In this version I have:

* Added authors and year to all DOI references in DESCRIPTION per CRAN
  format requirements (e.g., `Saqr et al. (2025) <doi:...>`).

* Unwrapped `\donttest{}` for all non-compute-intensive examples
  (simplicial complex functions, data utilities, extraction, frequencies,
  link prediction, association rules, centrality, estimator registry).
  Compute-intensive functions (bootstrap, permutation, clustering, MMM)
  retain `\donttest{}` with small fast toy examples before the wrapped block.

## Note

The system flags author names in DESCRIPTION references as misspelled.
We removed them previously to resolve the spelling notes. We received
feedback from CRAN that author names should be included, so we added
them again.

## Downstream dependencies

This is a new CRAN submission. The package has no reverse dependencies.

The companion package `cograph` (visualization) uses Nestimate's `cograph_network` class objects but does not formally depend on Nestimate.
