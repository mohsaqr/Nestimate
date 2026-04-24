## R CMD check results

0 errors | 0 warnings | 0 notes (local macOS `--as-cran --run-donttest`, win-builder devel/release, macOS builder, R-hub linux/windows/macos/macos-arm64).

## Test environments

* local macOS (Darwin 25.3.0), R 4.5.2 — `R CMD check --as-cran --run-donttest`
* win-builder R-devel
* win-builder R-release
* macOS builder (mac.r-project.org)
* R-hub: linux, windows, macos, macos-arm64

## Resubmission

This is a resubmission of 0.4.2 as 0.4.3, addressing the two NOTEs raised by
the CRAN incoming pre-tests on 0.4.2:

* **Vignette index.** The previous tarball was built with
  `R CMD build --no-build-vignettes`, which shipped pre-built `inst/doc/*.html`
  but omitted `build/vignette.rds`. 0.4.3 is built with plain `R CMD build`,
  so the vignette index is included.
* **Check time.** Windows check time was 11 min on 0.4.2. `tests/testthat/test-gimme.R`
  (which fits a `lavaan` SEM per subject via the `gimme` package) now calls
  `skip_on_cran()`. Local `--as-cran` check time dropped from 2m 41s to 2m 0.8s;
  Windows is projected to be ~8 min.

Earlier (0.4.2) cleanup, retained in 0.4.3:

* Full `--as-cran --run-donttest` audit, 0 errors / 0 warnings / 0 notes.
* `.Rbuildignore` strengthened with explicit `^Nestimate\.Rcheck$` and
  `^\.\.Rcheck$` entries.
* Stale build artifacts purged from the working tree.

## Notes on URL checks

The DOI `https://doi.org/10.1145/3706468.3706513` (ACM) returns HTTP 403
when accessed by automated URL checkers. The DOI is valid and resolves
correctly in a browser. ACM restricts automated access to their landing
pages. The same DOI format (`<doi:10.1145/3706468.3706513>`) is used in
the DESCRIPTION field where it is validated by the CRAN submission
infrastructure.

## Note on spelling

The DESCRIPTION references include author surnames (Saqr, López-Pernas)
which are flagged by the spell checker as misspelled words. These are
proper names and correct as written.

## Downstream dependencies

This is a new CRAN submission. The package has no reverse dependencies.

The companion package `cograph` (visualization) uses Nestimate's
`cograph_network` class objects but does not formally depend on Nestimate.
