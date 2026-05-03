# Equivalence test_that() blocks extracted from
# tests/testthat/test-pathways.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- Tests for pathways() ----

# Shared test data with clear higher-order dependencies
.make_ho_seqs <- function(n = 100, seed = 42) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    s <- character(20)
    s[1] <- sample(LETTERS[1:4], 1)
    for (j in 2:20) {
      if (j >= 3 && s[j - 2] == "A" && s[j - 1] == "B") {
        s[j] <- sample(c("C", "D"), 1, prob = c(0.95, 0.05))
      } else if (j >= 3 && s[j - 2] == "C" && s[j - 1] == "B") {
        s[j] <- sample(c("A", "D"), 1, prob = c(0.95, 0.05))
      } else {
        s[j] <- sample(LETTERS[1:4], 1)
      }
    }
    s
  })
}


# ---- pathways.net_hon ----

test_that("pathways output is parseable by plot_simplicial", {
  skip_if_not_installed("cograph")
  skip_if(!exists(".parse_pathways", envir = asNamespace("cograph")),
          "cograph:::.parse_pathways not available in this version")
  seqs <- .make_ho_seqs()
  hon <- build_hon(seqs, max_order = 3)
  pw <- pathways(hon, min_prob = 0.5)
  skip_if(length(pw) == 0, "No higher-order pathways found")

  parsed <- cograph:::.parse_pathways(pw, LETTERS[1:4])
  expect_length(parsed, length(pw))
  for (p in parsed) {
    expect_true(length(p$source) >= 1)
    expect_true(length(p$target) == 1)
  }
})
