# print() and summary(format =) of sequence_compare(). The object keeps its
# list shape (`$patterns` one row per pattern), which htna on CRAN reads.

.scp_fit <- function(test = "none") {
  seqs <- data.frame(
    T1 = c("A", "A", "B", "A", "B", "B", "A", "B"),
    T2 = c("B", "A", "B", "B", "A", "A", "B", "A"),
    T3 = c("A", "B", "A", "A", "B", "A", "A", "B"),
    T4 = c("B", "B", "A", "B", "A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )
  set.seed(1)
  sequence_compare(seqs, group = rep(c("x", "y"), each = 4L), sub = 2L,
                   min_freq = 1L, test = test, iter = 49L)
}

test_that("summary(format = 'long') has one row per pattern and group", {
  fit <- .scp_fit("permutation")
  wide <- summary(fit)
  long <- summary(fit, format = "long")
  expect_identical(wide, fit$patterns)
  expect_identical(names(long),
                   c("pattern", "length", "group", "count", "proportion",
                     "residual", "effect_size", "p_value"))
  expect_identical(nrow(long), nrow(wide) * length(fit$groups))
  # invariant: each group's long rows are that group's wide columns
  x_rows <- subset(long, group == "x")
  y_rows <- subset(long, group == "y")
  expect_identical(x_rows$pattern, wide$pattern)
  expect_identical(x_rows$count, wide$freq_x)
  expect_identical(y_rows$count, wide$freq_y)
  expect_identical(y_rows$residual, wide$resid_y)
  expect_identical(x_rows$p_value, wide$p_value)
})

test_that("summary(format = 'long') follows the test", {
  expect_true(all(c("statistic", "p_value") %in%
                    names(summary(.scp_fit("chisq"), format = "long"))))
  expect_identical(names(summary(.scp_fit("none"), format = "long")),
                   c("pattern", "length", "group", "count", "proportion",
                     "residual"))
  expect_error(summary(.scp_fit(), format = "tall"), "should be one of")
})

test_that("print() shows rounded rows with a group column", {
  fit <- .scp_fit("permutation")
  shown <- capture.output(print(fit, n = 1L, digits = 2L))
  expect_match(shown[1L], "Sequence Comparison  \\[3 patterns \\| 2 groups: x, y\\]")
  header <- shown[grepl("^ pattern", shown)]
  expect_length(header, 1L)
  expect_match(header, "group")
  expect_false(any(grepl("freq_x", shown)))
  expect_length(shown[grepl("^ [AB]->[AB] ", shown)], 2L)
  expect_true(any(grepl("2 more patterns \\(4 rows\\)", shown)))
  every <- capture.output(print(fit, n = Inf))
  expect_false(any(grepl("more patterns", every)))
  capture.output(shown_value <- withVisible(print(fit)))
  expect_false(shown_value$visible)
  expect_error(print(fit, n = 0), class = "nestimate_bad_input")
  expect_error(print(fit, digits = 1.5), class = "nestimate_bad_input")
})
