# Regenerate nct-baseline.rds, the frozen complete-data nct() output that
# test-nct.R checks against: input-handling changes (missing data, column
# alignment, degenerate permutations) must not move a complete-data number.
#
# Regenerate ONLY when a change of numbers is intentional (record it in NEWS.md).
# Run from the package root with the package loaded:
#   Rscript -e 'devtools::load_all("."); source("tests/testthat/fixtures/make-nct.R")'

source(file.path("tests", "testthat", "helper-nct.R"))

baseline <- lapply(seq_len(nrow(nct_configs)), function(i) {
  cfg <- nct_configs[i, ]
  g <- nct_groups(cfg$seed)
  set.seed(cfg$seed)
  nct_numbers(nct(g$a, g$b, iter = nct_iter, paired = cfg$paired,
                  abs = cfg$abs, weighted = cfg$weighted))
})
saveRDS(baseline, file.path("tests", "testthat", "fixtures", "nct-baseline.rds"))
message("wrote tests/testthat/fixtures/nct-baseline.rds")
