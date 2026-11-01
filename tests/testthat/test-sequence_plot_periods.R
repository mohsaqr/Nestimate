# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

# Three students, three days, 1-4 actions per student-day.
.daily_long <- function() {
  data.frame(
    student = rep(c("s1", "s2", "s3"), times = c(5, 9, 6)),
    day     = c("d1", "d2", "d2", "d3", "d3",
                "d1", "d1", "d2", "d2", "d2", "d3", "d3", "d3", "d3",
                "d1", "d1", "d1", "d2", "d3", "d3"),
    action  = c("plan", "read", "write", "read", "plan",
                "read", "read", "plan", "write", "write", "read", "plan",
                "write", "read", "plan", "plan", "write", "read", "write",
                "read"))
}
.daily_net <- function(d = .daily_long()) {
  build_network(d, method = "relative", actor = "student",
                action = "action", session = "day")
}
.daily_layout <- function(fill, net = .daily_net()) {
  .period_layout(net$data, net$metadata$student, net$metadata$day, fill)
}

test_that("stretch: each actor-period fills its block in equal shares", {
  cells <- .daily_layout("stretch")$cells
  # s3 (actor 3) on d1 has three actions: thirds of [0, 1].
  s3_d1 <- subset(cells, actor == 3L & period == 1L)
  expect_equal(s3_d1$xleft,  c(0, 1, 2) / 3)
  expect_equal(s3_d1$xright, c(1, 2, 3) / 3)
  # s1 on d1 has one action: it fills the whole block.
  s1_d1 <- subset(cells, actor == 1L & period == 1L)
  expect_equal(c(s1_d1$xleft, s1_d1$xright), c(0, 1))
  # s2 on d3 has four: quarters of [2, 3].
  s2_d3 <- subset(cells, actor == 2L & period == 3L)
  expect_equal(s2_d3$xright - s2_d3$xleft, rep(0.25, 4))
  expect_equal(range(c(s2_d3$xleft, s2_d3$xright)), c(2, 3))
})

test_that("stretch: widths in every actor-period sum to one block", {
  cells <- .daily_layout("stretch")$cells
  total <- tapply(cells$xright - cells$xleft,
                  paste(cells$actor, cells$period), sum)
  expect_equal(as.vector(total), rep(1, length(total)))
  first <- tapply(cells$xleft, paste(cells$actor, cells$period), min)
  per   <- tapply(cells$period, paste(cells$actor, cells$period), `[`, 1L)
  expect_equal(as.vector(first), as.vector(per) - 1)
})

test_that("pad: slots are 1 / max steps wide and match the aligned grid", {
  lay <- .daily_layout("pad")
  expect_identical(lay$width, 4L)
  expect_equal(lay$cells$xright - lay$cells$xleft,
               rep(0.25, nrow(lay$cells)))
  expect_identical(dim(lay$grid), c(3L, 12L))
  expect_identical(rownames(lay$grid), c("s1", "s2", "s3"))
  levs <- lay$levels
  # s2's d3 actions: read, plan, write, read.
  expect_identical(levs[lay$grid["s2", c("d3.1", "d3.2", "d3.3", "d3.4")]],
                   c("read", "plan", "write", "read"))
  # s1 has one action on d1, the other slots are empty.
  expect_true(all(is.na(lay$grid["s1", c("d1.2", "d1.3", "d1.4")])))
  expect_identical(sum(!is.na(lay$grid)), nrow(.daily_long()))
})

test_that("the layout does not depend on the order actors arrive in", {
  d <- .daily_long()
  shuffled <- d[order(factor(d$student, levels = c("s3", "s1", "s2")),
                      seq_len(nrow(d))), ]
  a <- .daily_layout("pad")
  b <- .daily_layout("pad", .daily_net(shuffled))
  expect_identical(a$grid, b$grid)
})

test_that("factor levels set the period order and keep an empty period", {
  d <- .daily_long()
  d$day <- factor(d$day, levels = c("d3", "d0", "d1", "d2"))
  lay <- .daily_layout("pad", .daily_net(d))
  expect_identical(lay$period_levels, c("d3", "d0", "d1", "d2"))
  expect_true(all(is.na(lay$grid[, paste0("d0.", 1:4)])))
  # d3 is now the first block.
  expect_identical(lay$levels[lay$grid["s2", "d3.1"]], "read")
})

test_that("sequence_plot(periods =) draws both types and both fills", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  net <- .daily_net()
  heat <- sequence_plot(net, periods = "day")
  expect_named(heat, c("ord", "codes", "palette", "levels", "periods", "fill"))
  expect_setequal(heat$ord, 1:3)
  expect_identical(heat$fill, "stretch")
  expect_identical(heat$periods, c("d1", "d2", "d3"))
  idx <- sequence_plot(net, periods = "day", fill = "pad", type = "index",
                       row_gap = 0.15, sort = "start")
  expect_identical(idx$fill, "pad")
  expect_identical(idx$codes, .daily_layout("pad")$grid)
})

test_that("periods errors are classed", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  net <- .daily_net()
  expect_error(sequence_plot(net$data, periods = "day"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(net, periods = "week"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(net, fill = "pad"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(net, periods = "day", type = "distribution"),
               class = "nestimate_bad_input")
  no_actor <- build_network(.daily_long(), method = "relative",
                            action = "action", session = "day")
  expect_error(sequence_plot(no_actor, periods = "day"),
               class = "nestimate_bad_input")
  # Two weeks under one day label: two sequences per student-day.
  d <- rbind(transform(.daily_long(), week = "w1"),
             transform(.daily_long(), week = "w2"))
  two_weeks <- build_network(d, method = "relative", actor = "student",
                             action = "action", session = c("week", "day"))
  expect_error(sequence_plot(two_weeks, periods = "day"),
               class = "nestimate_bad_input")
})

test_that("a missing actor or period in long data names the column", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- .daily_long()
  d$day[2L] <- NA
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "day"),
               "`day`", class = "nestimate_bad_input")
  d <- .daily_long()
  d$student[2L] <- NA
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "day"),
               "`student`", class = "nestimate_bad_input")
})

test_that("a tibble without a stored palette plots without warnings", {
  skip_if_not_installed("tibble")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- tibble::as_tibble(.daily_long())
  expect_no_warning(
    sequence_plot(d, actor = "student", action = "action", periods = "day"))
  expect_no_warning(
    sequence_plot(d, actor = "student", action = "action", periods = "day",
                  type = "index"))
})

test_that("trim with periods cuts period blocks, not time points", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  net  <- .daily_net()
  full <- .daily_layout("pad")
  two  <- .period_layout(net$data, net$metadata$student, net$metadata$day,
                         "pad", trim = 2)
  expect_identical(two$period_levels, c("d1", "d2"))
  expect_identical(two$actors, full$actors)
  expect_true(all(two$cells$period <= 2L))
  # Widest kept actor-period (s2 on d2, 3 actions) sets the slot count.
  expect_identical(two$width, 3L)
  expect_identical(dim(two$grid), c(3L, 6L))
  # A count past the number of periods changes nothing.
  expect_identical(.period_layout(net$data, net$metadata$student,
                                  net$metadata$day, "pad", trim = 10)$grid,
                   full$grid)
  # Fraction: every student is active on d3, so 0.5 keeps all three days.
  expect_identical(.period_layout(net$data, net$metadata$student,
                                  net$metadata$day, "pad",
                                  trim = 0.5)$period_levels,
                   c("d1", "d2", "d3"))
  # Through the public verb, both types and both routes.
  idx <- sequence_plot(.daily_long(), actor = "student", action = "action",
                       periods = "day", type = "index", trim = 2)
  expect_identical(idx$periods, c("d1", "d2"))
  heat <- sequence_plot(net, periods = "day", fill = "pad", trim = 1)
  expect_identical(heat$periods, "d1")
  expect_identical(dim(heat$codes), c(3L, 3L))
  # An empty first factor level under trim = 1 leaves nothing to draw.
  d <- .daily_long()
  d$day <- factor(d$day, levels = c("d0", "d1", "d2", "d3"))
  expect_error(sequence_plot(.daily_net(d), periods = "day", trim = 1),
               class = "nestimate_bad_input")
})

test_that("group with periods facets by a per-actor column", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- .daily_long()
  d$klass <- ifelse(d$student == "s2", "B", "A")
  res <- sequence_plot(d, actor = "student", action = "action",
                       periods = "day", type = "index", group = "klass")
  expect_identical(res$groups, c("A", "B"))
  # Panel A holds s1 and s3 (actors 1, 3), panel B holds s2.
  expect_setequal(res$orders[[1L]], c(1L, 3L))
  expect_identical(res$orders[[2L]], 2L)
  expect_identical(sort(res$ord), 1:3)
  # The period axis is shared and the grid is the same as without groups.
  expect_identical(res$periods, c("d1", "d2", "d3"))
  expect_identical(res$codes, .daily_layout("stretch")$grid)
  # Through a netobject whose metadata carries the column; trim combines.
  # prepare() aggregates covariates to character, so a factor's level order
  # is honoured only when the metadata column itself is a factor.
  net <- .daily_net(d)
  res2 <- sequence_plot(net, periods = "day", type = "index",
                        group = "klass", trim = 2, combined = FALSE)
  expect_identical(res2$groups, c("A", "B"))
  expect_identical(res2$periods, c("d1", "d2"))
  net$metadata$klass <- factor(net$metadata$klass, levels = c("B", "A"))
  expect_identical(sequence_plot(net, periods = "day", type = "index",
                                 group = "klass")$groups, c("B", "A"))
  # Heatmap ignores group, as it does without periods.
  heat <- sequence_plot(net, periods = "day", group = "klass")
  expect_null(heat$groups)
})

test_that("group with periods rejects what it cannot place", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- .daily_long()
  d$klass <- ifelse(d$student == "s2", "B", "A")
  # A vector, not a column name.
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "day", type = "index",
                             group = c("A", "B", "A")),
               class = "nestimate_bad_input")
  # An unknown column.
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "day", type = "index", group = "nope"),
               class = "nestimate_bad_input")
  # A column that varies within actor.
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "day", type = "index", group = "day"),
               class = "nestimate_bad_input")
  # Missing group values.
  d$klass[1L] <- NA
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "day", type = "index", group = "klass"),
               class = "nestimate_bad_input")
  # Long data without periods still refuses group.
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             type = "index", group = "klass"),
               class = "nestimate_bad_input")
})

test_that("long data goes straight in and matches the network route", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  from_long <- sequence_plot(.daily_long(), actor = "student",
                             action = "action", periods = "day")
  from_net  <- sequence_plot(.daily_net(), periods = "day")
  expect_identical(from_long$codes, from_net$codes)
  expect_identical(from_long$ord, from_net$ord)
  expect_identical(from_long$codes, .daily_layout("stretch")$grid)
})

test_that("long data without periods plots one row per actor", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- .daily_long()
  res <- sequence_plot(d, actor = "student", action = "action",
                       type = "index", sort = "start")
  expect_identical(nrow(res$codes), 3L)
  # s2 has 9 actions in all, the longest sequence.
  expect_identical(ncol(res$codes), 9L)
  expect_identical(res$levels[res$codes[2, seq_len(3)]],
                   c("read", "read", "plan"))
  dist <- sequence_plot(d, actor = "student", action = "action",
                        type = "distribution")
  expect_type(dist, "list")
})

test_that("timestamps order long data and do not split periods", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- .daily_long()
  # Reverse the rows but give timestamps hours apart in the true order:
  # without periods these gaps would split sessions.
  d$time <- as.POSIXct("2026-03-02 09:00", tz = "UTC") +
    seq_len(nrow(d)) * 3 * 3600
  shuffled <- d[rev(seq_len(nrow(d))), ]
  res <- sequence_plot(shuffled, actor = "student", action = "action",
                       time = "time", periods = "day")
  expect_identical(res$codes, .daily_layout("stretch")$grid)
})

test_that("long-data argument errors are classed", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- .daily_long()
  expect_error(sequence_plot(d, actor = "student", action = "act"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(d, actor = "student", periods = "day"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(d, action = "action", periods = "day"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             group = rep("a", nrow(d))),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(as.matrix(d), actor = "student",
                             action = "action"),
               class = "nestimate_bad_input")
})
