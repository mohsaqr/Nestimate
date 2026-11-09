# Development-only CRAN review. Run through mini-ci, never on the laptop.
# Modes: test (full suite), check (build + --as-cran), revdeps (CRAN sources).
mode <- commandArgs(trailingOnly = TRUE)
stopifnot(length(mode) == 1L, mode %in% c("test", "check", "revdeps"))
out_dir <- file.path(".ci", "cran-review")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
release_lib <- file.path(out_dir, "released-library")
dir.create(release_lib, showWarnings = FALSE)
.libPaths(c(release_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 300)
Sys.setenv(OMP_THREAD_LIMIT = 2, `_R_CHECK_LIMIT_CORES_` = "true")

# Local sibling builds can hide incompatibilities with the published packages.
available <- available.packages()
stopifnot(all(c("idiographic", "psychnets", "cograph", "tna", "Nestimate") %in%
                rownames(available)))
release_packages <- c("idiographic", "psychnets", "cograph", "tna")
needs_install <- vapply(release_packages, function(pkg) {
  desc <- file.path(release_lib, pkg, "DESCRIPTION")
  !file.exists(desc) || read.dcf(desc)[1, "Version"] != available[pkg, "Version"]
}, logical(1))
if (any(needs_install)) {
  install.packages(release_packages[needs_install], lib = release_lib,
                   type = "source", dependencies = NA)
}
versions <- vapply(release_packages, function(pkg) {
  actual <- as.character(packageVersion(pkg, lib.loc = release_lib))
  stopifnot(actual == available[pkg, "Version"],
            normalizePath(find.package(pkg)) == normalizePath(file.path(release_lib, pkg)))
  actual
}, character(1))
print(versions)
saveRDS(available, file.path(out_dir, "available.rds"))
writeLines(c(R.version.string, paste(names(versions), versions)),
           file.path(out_dir, "environment.txt"))

if (mode == "test") {
  results <- testthat::test_local(stop_on_failure = FALSE)
  saveRDS(results, file.path(out_dir, "test-results.rds"))
  tab <- as.data.frame(results)
  totals <- c(failed = sum(tab$failed), errors = sum(tab$error),
              warnings = sum(tab$warning), skipped = sum(tab$skipped),
              passed = sum(tab$passed))
  print(totals)
  writeLines(capture.output(print(totals)), file.path(out_dir, "test-summary.txt"))
  if (totals["failed"] > 0 || totals["errors"] > 0) quit(status = 1L)
}

if (mode == "check") {
  tarball <- pkgbuild::build(".", dest_path = out_dir, binary = FALSE,
                            vignettes = TRUE, manual = FALSE, quiet = FALSE)
  result <- rcmdcheck::rcmdcheck(tarball, args = c("--as-cran", "--no-manual"),
                                 check_dir = file.path(out_dir, "check"),
                                 error_on = "never")
  print(result)
  saveRDS(result, file.path(out_dir, "check-result.rds"))
  writeLines(capture.output(print(result)), file.path(out_dir, "check-summary.txt"))
  if (result$status != 0L || length(result$errors) || length(result$warnings)) {
    quit(status = 1L)
  }
}

if (mode == "revdeps") {
  revdeps <- tools::package_dependencies("Nestimate", db = available,
                                          which = c("Depends", "Imports", "LinkingTo",
                                                    "Suggests"), reverse = TRUE)[[1]]
  print(revdeps)
  revdir <- file.path(out_dir, "revdeps")
  dir.create(revdir, showWarnings = FALSE)
  candidate_lib <- file.path(revdir, "candidate-library")
  baseline_lib <- file.path(revdir, "baseline-library")
  invisible(lapply(c(candidate_lib, baseline_lib), dir.create, showWarnings = FALSE))
  install.packages(".", repos = NULL, type = "source", lib = candidate_lib)
  install.packages("Nestimate", lib = baseline_lib, type = "source", dependencies = NA)
  stopifnot(as.character(packageVersion("Nestimate", lib.loc = candidate_lib)) ==
              read.dcf("DESCRIPTION")[1, "Version"],
            as.character(packageVersion("Nestimate", lib.loc = baseline_lib)) ==
              available["Nestimate", "Version"])
  summaries <- lapply(revdeps, function(pkg) {
    message("Checking reverse dependency: ", pkg)
    deps <- tools::package_dependencies(pkg, db = available,
                                        which = c("Depends", "Imports", "LinkingTo"),
                                        recursive = TRUE)[[1]]
    deps <- setdiff(deps, c("R", "Nestimate", rownames(installed.packages())))
    if (length(deps)) install.packages(deps, lib = release_lib, type = "source")
    src <- download.packages(pkg, destdir = revdir, type = "source")[1, 2]
    run_check <- function(lib, suffix) {
      result <- rcmdcheck::rcmdcheck(src, args = c("--no-manual", "--no-build-vignettes"),
                                     libpath = c(lib, .libPaths()),
                                     check_dir = file.path(revdir, paste0(pkg, "-", suffix)),
                                     env = c(`_R_CHECK_FORCE_SUGGESTS_` = "false",
                                             NOT_CRAN = "false"),
                                     error_on = "never")
      print(result)
      saveRDS(result, file.path(revdir, paste0(pkg, "-", suffix, ".rds")))
      result
    }
    candidate <- run_check(candidate_lib, "candidate")
    baseline <- if (candidate$status != 0L || length(candidate$errors) ||
                    length(candidate$warnings)) {
      run_check(baseline_lib, "baseline")
    } else NULL
    data.frame(package = pkg, version = available[pkg, "Version"],
               status = candidate$status,
               errors = length(candidate$errors), warnings = length(candidate$warnings),
               notes = length(candidate$notes),
               baseline_status = if (is.null(baseline)) NA_integer_ else baseline$status,
               baseline_errors = if (is.null(baseline)) NA_integer_ else length(baseline$errors),
               baseline_warnings = if (is.null(baseline)) NA_integer_ else length(baseline$warnings))
  })
  summary <- do.call(rbind, summaries)
  print(summary)
  write.csv(summary, file.path(revdir, "summary.csv"), row.names = FALSE)
  # Nonzero results require manual comparison with the saved baseline evidence.
  if (any(summary$status != 0L | summary$errors > 0 | summary$warnings > 0)) {
    quit(status = 1L)
  }
}
