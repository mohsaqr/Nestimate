#!/usr/bin/env Rscript

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
source("tests/testthat/helper-simulate.R")

if (!requireNamespace("tna", quietly = TRUE)) {
  stop("Package 'tna' is required for this equivalence audit.", call. = FALSE)
}

out_dir <- "audit_codex"
csv_path <- file.path(out_dir, "centrality_tna_equivalence_100.csv")
md_path <- file.path(out_dir, "centrality_tna_equivalence_100.md")

measures <- c(
  "OutStrength", "InStrength", "ClosenessIn", "ClosenessOut",
  "Closeness", "Betweenness", "BetweennessRSP", "Diffusion", "Clustering"
)
tolerance <- 1e-10
set.seed(20260629)

configs <- lapply(seq_len(100L), function(i) {
  list(
    id = i,
    seed = sample.int(1e6, 1L),
    n_actors = sample(c(15L, 20L, 30L, 45L, 60L), 1L),
    n_states = sample(3L:8L, 1L),
    seq_length = sample(c(8L, 10L, 12L, 16L, 20L, 25L, 30L), 1L)
  )
})

.max_abs_diff <- function(a, b) {
  same_missing <- (is.na(a) & is.na(b)) | (is.nan(a) & is.nan(b))
  keep <- !same_missing
  if (!any(keep)) return(0)
  if (any(is.na(a[keep]) != is.na(b[keep]))) return(Inf)
  max(abs(a[keep] - b[keep]), na.rm = TRUE)
}

.measure_diff <- function(ours, ref, measure) {
  .max_abs_diff(as.numeric(ours[[measure]]), as.numeric(ref[[measure]]))
}

rows <- do.call(rbind, lapply(configs, function(cfg) {
  dat <- simulate_sequences(
    n_actors = cfg$n_actors,
    n_states = cfg$n_states,
    seq_length = cfg$seq_length,
    seed = cfg$seed
  )

  nest_net <- build_network(dat, method = "relative")
  tna_net <- tna::tna(dat)

  nest_raw <- suppressMessages(net_centrality(
    nest_net, measures = measures, normalize_diffusion = FALSE
  ))
  tna_raw <- tna::centralities(tna_net, measures = measures)

  # The one intentional default difference: Nestimate defaults Diffusion to
  # range-normalized. Compare that default against tna's normalized Diffusion.
  nest_diff_default <- suppressMessages(net_centrality(
    nest_net, measures = "Diffusion"
  ))
  tna_diff_norm <- tna::centralities(
    tna_net, measures = "Diffusion", normalize = TRUE
  )

  diffs <- vapply(measures, function(m) {
    .measure_diff(nest_raw, tna_raw, m)
  }, numeric(1))
  default_diffusion_delta <- .measure_diff(
    nest_diff_default, tna_diff_norm, "Diffusion"
  )

  data.frame(
    id = cfg$id,
    seed = cfg$seed,
    n_actors = cfg$n_actors,
    n_states = cfg$n_states,
    seq_length = cfg$seq_length,
    max_abs_delta = max(c(diffs, default_diffusion_delta), na.rm = TRUE),
    raw_centrality_max_abs_delta = max(diffs, na.rm = TRUE),
    default_diffusion_max_abs_delta = default_diffusion_delta,
    status = if (max(c(diffs, default_diffusion_delta), na.rm = TRUE) <= tolerance) {
      "PASS"
    } else {
      "FAIL"
    },
    t(as.data.frame(diffs)),
    check.names = FALSE,
    row.names = NULL
  )
}))

utils::write.csv(rows, csv_path, row.names = FALSE)

.status_n <- function(x, label) {
  idx <- match(label, names(x))
  if (is.na(idx)) 0L else unname(as.integer(x[[idx]]))
}
status_counts <- table(rows$status)
pass_n <- .status_n(status_counts, "PASS")
fail_n <- .status_n(status_counts, "FAIL")
max_delta <- max(rows$max_abs_delta, na.rm = TRUE)
raw_max_delta <- max(rows$raw_centrality_max_abs_delta, na.rm = TRUE)
default_diffusion_max_delta <- max(
  rows$default_diffusion_max_abs_delta, na.rm = TRUE
)

lines <- c(
  "# Centrality tna Equivalence Audit",
  "",
  sprintf("Date: %s", format(Sys.Date())),
  "",
  "Scope: 100 transition networks generated with `simulate_sequences()` from `tests/testthat/helper-simulate.R`, the local Saqrlab-style simulation stand-in used by Nestimate tests.",
  "",
  sprintf("Tolerance: `%g`", tolerance),
  "",
  "Compared measures:",
  paste0("- `", measures, "`"),
  "",
  "Checks:",
  "- Raw Nestimate centralities used `normalize_diffusion = FALSE` and were compared to `tna::centralities(normalize = FALSE)`.",
  "- Nestimate's default `Diffusion` was compared separately to `tna::centralities(measures = \"Diffusion\", normalize = TRUE)`.",
  "",
  "Results:",
  sprintf("- Networks checked: `%d`", nrow(rows)),
  sprintf("- Passed: `%d`", pass_n),
  sprintf("- Failed: `%d`", fail_n),
  sprintf("- Maximum absolute delta, all checks: `%.17g`", max_delta),
  sprintf("- Maximum absolute delta, raw centrality parity: `%.17g`", raw_max_delta),
  sprintf("- Maximum absolute delta, default normalized Diffusion: `%.17g`", default_diffusion_max_delta),
  "",
  sprintf("CSV detail: `%s`", csv_path)
)
writeLines(lines, md_path)

print(rows[, c(
  "id", "n_actors", "n_states", "seq_length", "max_abs_delta", "status"
)])
cat("\n")
cat(sprintf("PASS=%d FAIL=%d max_abs_delta=%.17g\n", pass_n, fail_n, max_delta))

if (fail_n > 0L || !is.finite(max_delta) || max_delta > tolerance) {
  stop("Centrality equivalence audit failed.", call. = FALSE)
}
