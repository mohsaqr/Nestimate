testthat::skip_on_cran()

test_that(".period_bin groups dates into Monday weeks, keeping empty weeks", {
  # 2026-02-04 is a Wednesday, 2026-02-08 a Sunday, 2026-02-23 a Monday;
  # the week of 2026-02-16 has no date and must stay on the axis
  d <- as.Date(c("2026-02-04", "2026-02-08", "2026-02-23", "2026-02-09"))
  wk <- .period_bin(d, "week", "day")
  expect_true(is.factor(wk))
  expect_identical(as.character(wk),
                   c("2026-02-02", "2026-02-02", "2026-02-23", "2026-02-09"))
  expect_identical(levels(wk),
                   c("2026-02-02", "2026-02-09", "2026-02-16", "2026-02-23"))
})

test_that(".period_bin reads text dates and POSIXct, and groups months", {
  txt <- c("2026-01-31", "2026-03-02")
  mo <- .period_bin(txt, "month", "day")
  expect_identical(as.character(mo), c("2026-01-01", "2026-03-01"))
  expect_identical(levels(mo), c("2026-01-01", "2026-02-01", "2026-03-01"))
  ct <- as.POSIXct(c("2026-02-04 10:00:00", "2026-02-10 09:00:00"), tz = "UTC")
  expect_identical(as.character(.period_bin(ct, "week", "day")),
                   c("2026-02-02", "2026-02-09"))
  # "day" leaves the column untouched
  expect_identical(.period_bin(txt, "day", "day"), txt)
})

test_that("period_unit refuses a column that is not a date", {
  expect_error(.period_bin(c("w1", "w2"), "week", "p"),
               class = "nestimate_bad_input")
  log <- data.frame(who = c("a", "a", "b"), act = c("x", "y", "x"),
                    p = c("w1", "w2", "w1"), stringsAsFactors = FALSE)
  expect_error(sequence_plot(log, actor = "who", action = "act",
                             periods = "p", period_unit = "week"),
               class = "nestimate_bad_input")
})
