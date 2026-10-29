# ---- pattern_outcome() ----

# A->B occurs in rows 1, 2, 3, 6 (y = p, p, p, f); absent rows 4, 5, 7, 8
# have y = f, p, f, f. B->A and B->B occur under one outcome value only.
.po_data <- function() {
  data.frame(T1 = c("A", "A", "A", "B", "B", "A", "B", "B"),
             T2 = c("B", "B", "B", "A", "B", "B", "A", "A"),
             y = c("p", "p", "p", "f", "p", "f", "f", "f"),
             stringsAsFactors = FALSE)
}

.po_fit <- function(...) {
  pattern_outcome(.po_data(), outcome = "y", cols = c("T1", "T2"),
                  min_freq = 1, min_support = 0, ...)
}

test_that("one pattern: estimate is the 2 x 2 log odds ratio", {
  res <- .po_fit(len = 2, n = 1)
  expect_identical(res$pattern, "A->B")
  # present: 3 p, 1 f; absent: 1 p, 3 f -> OR = (3 * 3) / (1 * 1)
  expect_equal(res$estimate, log(9), tolerance = 1e-8)
  expect_equal(res$odds_ratio, exp(res$estimate))
  expect_equal(res$support_p, 3 / 4)
  expect_equal(res$support_f, 1 / 4)
})

test_that("chisq is the equal-split statistic codyna ranks by", {
  res <- .po_fit(len = 1:2, n = 2)
  # A->B: 3 vs 1 containing sequences, expected 2 each
  expect_equal(res$chisq[res$pattern == "A->B"], 1)
  # A: 3 vs 4, expected 3.5 each
  expect_equal(res$chisq[res$pattern == "A"], 1 / 7)
})

test_that("the top n by priority are selected, listed by name", {
  res <- .po_fit(len = 1:2, n = 2)
  expect_identical(res$pattern, c("A", "A->B"))
  # start = "A" leaves A (support 7/8) and A->B (4/8)
  by_support <- .po_fit(len = 1:2, n = 1, priority = "support", start = "A")
  expect_identical(by_support$pattern, "A")
})

test_that("patterns absent under one outcome value are never selected", {
  res <- .po_fit(len = 2, n = 10)
  expect_identical(res$pattern, "A->B")
})

test_that("swapping the reference flips the log odds ratio", {
  ref_f <- .po_fit(len = 2, n = 1)
  ref_p <- .po_fit(len = 2, n = 1, reference = "p")
  expect_equal(ref_p$estimate, -ref_f$estimate, tolerance = 1e-8)
})

test_that("occurrences equal presence when no pattern repeats in a sequence", {
  presence <- .po_fit(len = 2, n = 1)
  counts <- .po_fit(len = 2, n = 1, occurrences = TRUE)
  expect_equal(counts$estimate, presence$estimate)
})

test_that("a pattern present in every sequence warns and is NA", {
  expect_warning(res <- .po_fit(len = 1:2, n = 3),
                 class = "nestimate_aliased_pattern")
  expect_true(is.na(res$estimate[res$pattern == "B"]))
})

test_that("missing outcomes are dropped with a classed warning", {
  d <- .po_data()
  d <- rbind(d, data.frame(T1 = "A", T2 = "B", y = NA))
  expect_warning(pattern_outcome(d, outcome = "y", cols = c("T1", "T2"),
                                 len = 2, n = 1, min_freq = 1,
                                 min_support = 0),
                 class = "nestimate_missing_outcome")
})

test_that("bad outcomes and empty selections are refused by class", {
  d <- .po_data()
  d$three <- rep(c("a", "b", "c"), length.out = nrow(d))
  expect_error(pattern_outcome(d, outcome = "three", cols = c("T1", "T2")),
               class = "nestimate_bad_input")
  expect_error(pattern_outcome(d, outcome = "nope", cols = c("T1", "T2")),
               class = "nestimate_bad_input")
  expect_error(pattern_outcome(d, outcome = "y", cols = c("T1", "T2"),
                               len = 2, min_freq = 99),
               class = "nestimate_no_patterns")
})

test_that("lift and support agree with sequence_patterns()", {
  res <- .po_fit(len = 2, n = 1)
  pat <- sequence_patterns(.po_data()[c("T1", "T2")], len = 2, min_freq = 1,
                           min_support = 0)
  ref <- pat[pat$pattern == "A->B", ]
  expect_equal(res$lift, ref$lift)
  expect_equal(res$support, ref$support)
  expect_identical(res$occurrences, ref$occurrences)
})

test_that("an outcome vector works for a network that stores no covariates", {
  d <- .po_data()
  from_column <- .po_fit(len = 2, n = 1)
  net <- build_tna(d[c("T1", "T2")])
  from_vector <- pattern_outcome(net, outcome = d$y, len = 2, n = 1,
                                 min_freq = 1, min_support = 0)
  expect_identical(from_vector, from_column)
  expect_error(pattern_outcome(net, outcome = "y"),
               class = "nestimate_bad_input")
  expect_error(pattern_outcome(net, outcome = c("p", "f")),
               class = "nestimate_bad_input")
})

test_that("a network keeps the outcome when the positions are named T1, T2", {
  net <- build_network(.po_data(), method = "relative")
  expect_identical(pattern_outcome(net, outcome = "y", len = 2, n = 1,
                                   min_freq = 1, min_support = 0)$pattern,
                   "A->B")
  # positions without numbered names: the outcome is read as a position
  unnamed <- .po_data()
  names(unnamed) <- c("first", "second", "y")
  absorbed <- build_network(unnamed, method = "relative")
  expect_error(pattern_outcome(absorbed, outcome = "y"),
               class = "nestimate_bad_input")
  kept <- build_network(.po_data(), method = "relative", metadata_cols = "y")
  res <- pattern_outcome(kept, outcome = "y", len = 2, n = 1, min_freq = 1,
                         min_support = 0)
  expect_identical(res$pattern, "A->B")
})
