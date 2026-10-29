# ---- sequence_patterns() ----

# S1 = A B A B, S2 = A B, S3 = B B B
.sp_tiny <- function() {
  matrix(c("A", "B", "A", "B",
           "A", "B", NA,  NA,
           "B", "B", "B", NA), nrow = 3, byrow = TRUE)
}

.sp_row <- function(tab, pattern) {
  tab[tab$pattern == pattern, , drop = FALSE]
}

test_that("n-gram counts, support and lift match hand values", {
  tab <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 1, min_support = 0)
  ab <- .sp_row(tab, "A->B")
  expect_identical(ab$occurrences, 3L)
  expect_identical(ab$n_sequences, 2L)
  expect_equal(ab$support, 2 / 3)
  expect_equal(ab$proportion, 3 / 6)
  # support(A) = 2/3, support(B) = 1
  expect_equal(ab$lift, (2 / 3) / (2 / 3 * 1))
  bb <- .sp_row(tab, "B->B")
  expect_identical(bb$occurrences, 2L)
  expect_identical(bb$n_sequences, 1L)
  expect_equal(bb$lift, 1 / 3)
  expect_identical(sort(tab$pattern), c("A->B", "B->A", "B->B"))
})

test_that("gapped patterns skip one position and never span a missing value", {
  tab <- sequence_patterns(.sp_tiny(), type = "gapped", gap = 1,
                           min_freq = 1, min_support = 0)
  expect_identical(sort(tab$pattern), c("A->*->A", "B->*->B"))
  bxb <- .sp_row(tab, "B->*->B")
  expect_identical(bxb$occurrences, 2L)
  expect_identical(bxb$n_sequences, 2L)
  expect_identical(bxb$length, 3L)
  expect_equal(bxb$lift, 2 / 3)
})

test_that("repeated patterns are runs of one state", {
  tab <- sequence_patterns(.sp_tiny(), type = "repeated", len = 2:3,
                           min_freq = 1, min_support = 0)
  expect_identical(sort(tab$pattern), c("B->B", "B->B->B"))
  expect_identical(.sp_row(tab, "B->B->B")$occurrences, 1L)
})

test_that("start, end and contain filter by whole states", {
  tab <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 1,
                           min_support = 0, start = "B")
  expect_identical(sort(tab$pattern), c("B->A", "B->B"))
  tab <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 1,
                           min_support = 0, end = "A")
  expect_identical(tab$pattern, "B->A")
  tab <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 1,
                           min_support = 0, contain = "A")
  expect_identical(sort(tab$pattern), c("A->B", "B->A"))
})

test_that("min_freq and min_support remove patterns but not proportions", {
  all_rows <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 1,
                                min_support = 0)
  kept <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 2,
                            min_support = 0.5)
  expect_identical(kept$pattern, "A->B")
  expect_equal(kept$proportion, .sp_row(all_rows, "A->B")$proportion)
})

test_that("support never rises when a pattern is extended", {
  tab <- sequence_patterns(trajectories, len = 2:3, min_freq = 1,
                           min_support = 0)
  tri <- tab[tab$length == 3L, , drop = FALSE]
  prefix <- sub("->[^>]*$", "", tri$pattern)
  prefix_support <- tab$support[match(prefix, tab$pattern)]
  expect_true(all(tri$support <= prefix_support))
})

test_that("row order of the sequences does not change the result", {
  data <- as.data.frame(trajectories, stringsAsFactors = FALSE)
  shuffled <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  expect_identical(sequence_patterns(data, len = 2:3),
                   sequence_patterns(shuffled, len = 2:3))
})

test_that("no pattern gives an empty table with the full column set", {
  tab <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 99)
  expect_identical(nrow(tab), 0L)
  expect_identical(names(tab), c("pattern", "length", "occurrences",
                                 "proportion", "n_sequences", "support",
                                 "lift"))
})

test_that("a state containing the separator is refused", {
  bad <- matrix(c("A->1", "B", "B", "A->1"), nrow = 2)
  expect_error(sequence_patterns(bad), class = "nestimate_bad_input")
})


# ---- group test ----

# A->B occurs in 3 of 4 "g2" sequences and 1 of 4 "g1" sequences
.sp_two <- function() {
  data.frame(T1 = c("A", "A", "A", "B", "A", "B", "B", "B"),
             T2 = c("B", "B", "B", "A", "B", "A", "A", "A"),
             g = rep(c("g2", "g1"), each = 4), stringsAsFactors = FALSE)
}

# A->B occurs in 3, 1 and 0 of 4 sequences in groups x, y, z
.sp_three <- function() {
  data.frame(T1 = c("A", "A", "A", "B", "A", "B", "B", "B", "B", "B", "B", "B"),
             T2 = rep("B", 12), g = rep(c("x", "y", "z"), each = 4),
             stringsAsFactors = FALSE)
}

.sp_grouped <- function(d, ...) {
  sequence_patterns(d, cols = c("T1", "T2"), len = 1:2, min_freq = 1,
                    min_support = 0, group = "g", ...)
}

test_that("two groups: Fisher p-value and per-group support", {
  res <- .sp_grouped(.sp_two())
  ab <- res[res$pattern == "A->B", ]
  # hypergeometric 4/4/4/4 margins: P(X <= 1 or X >= 3) = 34 / 70
  expect_equal(ab$p_value, 34 / 70)
  expect_equal(ab$support_g2, 3 / 4)
  expect_equal(ab$support_g1, 1 / 4)
  expect_true(ab$odds_ratio > 1)
  expect_equal(res$p_adj, stats::p.adjust(res$p_value, method = "BH"))
})

test_that("swapping the reference inverts the odds ratio and its interval", {
  ref_g1 <- .sp_grouped(.sp_two())
  ref_g2 <- .sp_grouped(.sp_two(), reference = "g2")
  ref_g2 <- ref_g2[match(ref_g1$pattern, ref_g2$pattern), ]
  finite <- is.finite(ref_g1$odds_ratio) & ref_g1$odds_ratio > 0
  # fisher.test() finds the conditional MLE by root-finding (uniroot), so
  # the estimate inverts only to the solver tolerance; the limits are exact
  expect_equal(ref_g2$odds_ratio[finite], 1 / ref_g1$odds_ratio[finite],
               tolerance = 1e-4)
  expect_equal(ref_g2$or_lower[finite], 1 / ref_g1$or_upper[finite],
               tolerance = 1e-6)
  expect_equal(ref_g2$p_value, ref_g1$p_value)
})

test_that("three groups: chi-square of presence by group and Cramer's V", {
  res <- .sp_grouped(.sp_three())
  ab <- res[res$pattern == "A->B", ]
  expect_equal(ab$chisq, 5.25)
  expect_equal(ab$cramers_v, sqrt(5.25 / 12))
  expect_equal(ab$p_value, stats::pchisq(5.25, 2, lower.tail = FALSE))
  expect_false("odds_ratio" %in% names(res))
})

test_that("a pattern in every sequence has no contrast", {
  three <- .sp_grouped(.sp_three())
  expect_equal(three$chisq[three$pattern == "B"], 0)
  two <- .sp_grouped(.sp_two())
  expect_true(is.na(two$odds_ratio[two$pattern == "B"]))
})

test_that("group as column, vector and grouped network agree", {
  d <- .sp_two()
  from_column <- .sp_grouped(d)
  from_vector <- sequence_patterns(d[c("T1", "T2")], len = 1:2, min_freq = 1,
                                   min_support = 0, group = d$g)
  expect_identical(from_vector, from_column)
  nets <- build_network(d, method = "relative", group = "g")
  from_network <- sequence_patterns(nets, len = 1:2, min_freq = 1,
                                    min_support = 0)
  expect_identical(sort(from_network$pattern), sort(from_column$pattern))
  expect_equal(from_network$p_value[order(from_network$pattern)],
               from_column$p_value[order(from_column$pattern)])
})

test_that("top returns the first rows in output order", {
  all_rows <- sequence_patterns(trajectories, len = 2:3)
  expect_identical(sequence_patterns(trajectories, len = 2:3, top = 3),
                   utils::head(all_rows, 3))
})

test_that("bad groups are refused by class, missing ones warn", {
  d <- .sp_two()
  expect_error(.sp_grouped(d[c("T1", "T2")]), class = "nestimate_bad_input")
  expect_error(sequence_patterns(d[c("T1", "T2")], group = c("a", "b")),
               class = "nestimate_bad_input")
  d$one <- "same"
  expect_error(sequence_patterns(d, cols = c("T1", "T2"), group = "one"),
               class = "nestimate_bad_input")
  d$g[1] <- NA
  expect_warning(.sp_grouped(d), class = "nestimate_missing_group")
})


# ---- pattern = (template search) ----

test_that("a template returns every concrete instance", {
  # S1 = A B A B, S2 = A B, S3 = B B B
  hits <- sequence_patterns(.sp_tiny(), pattern = "A->*->A", min_freq = 1,
                            min_support = 0)
  expect_identical(hits$pattern, "A->B->A")
  expect_identical(hits$occurrences, 1L)
  following_b <- sequence_patterns(.sp_tiny(), pattern = "B->*", min_freq = 1,
                                   min_support = 0)
  expect_identical(sort(following_b$pattern), c("B->A", "B->B"))
  expect_identical(following_b$occurrences[following_b$pattern == "B->B"], 2L)
})

test_that("a template without wildcards equals its n-gram row", {
  ngrams <- sequence_patterns(.sp_tiny(), len = 2, min_freq = 1,
                              min_support = 0)
  searched <- sequence_patterns(.sp_tiny(), pattern = "A->B", min_freq = 1,
                                min_support = 0)
  expect_identical(searched$occurrences,
                   ngrams$occurrences[ngrams$pattern == "A->B"])
  expect_equal(searched$lift, ngrams$lift[ngrams$pattern == "A->B"])
})

test_that("a malformed template is refused by class", {
  expect_error(sequence_patterns(.sp_tiny(), pattern = "A->->B"),
               class = "nestimate_bad_input")
})


# ---- group = "last_obs" ----

# Each sequence ends in its outcome: pass, fail, pass, fail
.sp_terminal <- function() {
  data.frame(T1 = c("A", "A", "B", "B"), T2 = c("B", "A", "A", "fail"),
             T3 = c("pass", "fail", "B", NA), T4 = c(NA, NA, "pass", NA),
             stringsAsFactors = FALSE)
}

test_that("last_obs groups by the final state and removes it", {
  res <- sequence_patterns(.sp_terminal(), group = "last_obs", len = 1:2,
                           min_freq = 1, min_support = 0)
  expect_false(any(grepl("pass|fail", res$pattern)))
  ab <- res[res$pattern == "A->B", ]
  # A->B occurs in sequences 1 and 3, both ending in pass
  expect_identical(ab$n_sequences, 2L)
  expect_equal(ab$support_pass, 1)
  expect_equal(ab$support_fail, 0)
})

test_that("pattern_outcome takes the final state as the outcome", {
  rows <- rep(c("A B A pass", "A B B pass", "B A A fail", "B B A fail",
                "A A B pass", "B A B fail"), 6)
  d <- as.data.frame(do.call(rbind, strsplit(rows, " ")),
                     stringsAsFactors = FALSE)
  d$V3[c(1, 8, 15)] <- "B"
  res <- pattern_outcome(d, outcome = "last_obs", len = 2, n = 3,
                         min_freq = 1)
  expect_identical(names(res)[7:8], c("support_fail", "support_pass"))
  expect_false(any(grepl("pass|fail", res$pattern)))
})

test_that("a gapped window with a missing middle position is skipped", {
  d <- data.frame(T1 = "A", T2 = NA_character_, T3 = "B")
  expect_identical(nrow(sequence_patterns(d, type = "gapped", gap = 1,
                                          min_freq = 1)), 0L)
  expect_identical(nrow(sequence_patterns(d, pattern = "A->*->B",
                                          min_freq = 1)), 0L)
})

test_that("a repeated length or gap counts each occurrence once", {
  d <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"), T3 = c("A", "B"))
  expect_identical(sequence_patterns(d, len = c(2, 2), min_freq = 1),
                   sequence_patterns(d, len = 2, min_freq = 1))
  expect_identical(sequence_patterns(d, type = "gapped", gap = c(1, 1),
                                     min_freq = 1),
                   sequence_patterns(d, type = "gapped", gap = 1,
                                     min_freq = 1))
})

test_that("no matching pattern keeps the grouped columns, silently", {
  d <- data.frame(T1 = c("A", "B", "A"), T2 = c("B", "A", "A"))
  two <- expect_silent(sequence_patterns(d[1:2, ], group = c("x", "y"),
                                         len = 3, min_freq = 1))
  expect_identical(nrow(two), 0L)
  expect_true(all(c("support_x", "support_y", "odds_ratio", "or_lower",
                    "or_upper", "p_value", "p_adj") %in% names(two)))
  three <- expect_silent(sequence_patterns(d, group = c("x", "y", "z"),
                                           len = 3, min_freq = 1))
  expect_true(all(c("chisq", "cramers_v", "p_value", "p_adj") %in%
                    names(three)))
})

test_that("a template with a trailing separator is refused", {
  d <- data.frame(T1 = "A", T2 = "B")
  expect_error(sequence_patterns(d, pattern = "A->", min_freq = 1),
               class = "nestimate_bad_input")
  expect_error(sequence_patterns(d, pattern = "->B", min_freq = 1),
               class = "nestimate_bad_input")
})
