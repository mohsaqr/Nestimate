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
  # Fraction: the busiest periods that hold that share of the actions.
  counts <- tabulate(full$cells$period, nbins = 3L)
  expect_gte(sum(sort(counts, decreasing = TRUE)[1:2]) / sum(counts), 0.5)
  expect_lt(max(counts) / sum(counts), 0.5)
  expect_identical(.period_layout(net$data, net$metadata$student,
                                  net$metadata$day, "pad",
                                  trim = 0.5)$period_levels,
                   full$period_levels[rank(-counts, ties.method = "first") <= 2L])
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

test_that("a Date or day counter keeps the periods nobody acted in", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- data.frame(student = rep(c("s1", "s2"), each = 3),
                  day = as.Date(c("2026-01-05", "2026-01-05", "2026-01-08",
                                  "2026-01-05", "2026-01-06", "2026-01-08")),
                  action = c("a", "b", "a", "b", "b", "c"))
  res <- sequence_plot(d, actor = "student", action = "action",
                       periods = "day", type = "index")
  expect_identical(res$periods, as.character(seq(as.Date("2026-01-05"),
                                                 as.Date("2026-01-08"),
                                                 by = "day")))
  # 7 January is on the axis with no cell; s1 has none on the 6th either.
  expect_true(all(is.na(res$codes[, startsWith(colnames(res$codes), "2026-01-07")])))
  expect_true(all(is.na(res$codes["s1", startsWith(colnames(res$codes), "2026-01-06")])))
  # Same through a network, whose metadata keeps the Date class.
  net <- build_network(d, method = "relative", actor = "student",
                       action = "action", session = "day")
  expect_identical(sequence_plot(net, periods = "day")$periods, res$periods)
  # A whole-number counter fills the gaps too; a label column does not.
  d$n <- as.integer(d$day - min(d$day)) + 1L
  expect_identical(sequence_plot(d, actor = "student", action = "action",
                                 periods = "n", type = "index")$periods,
                   as.character(1:4))
  d$lab <- paste0("d", d$n)
  expect_identical(sequence_plot(d, actor = "student", action = "action",
                                 periods = "lab", type = "index")$periods,
                   c("d1", "d2", "d4"))
  # A huge day-number range is not expanded past a count cut, and is
  # refused without one.
  d$far <- c(1, 1, 1, 1, 1, 1e12)
  expect_identical(sequence_plot(d, actor = "student", action = "action",
                                 periods = "far", type = "index",
                                 trim = 3)$periods, c("1", "2", "3"))
  expect_error(sequence_plot(d, actor = "student", action = "action",
                             periods = "far", type = "index"),
               class = "nestimate_bad_input")
  # trim counts the empty blocks as well.
  expect_identical(sequence_plot(d, actor = "student", action = "action",
                                 periods = "day", type = "index",
                                 trim = 3)$periods, res$periods[1:3])
})

test_that("the row half of trim ends each actor's row earlier and removes none", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  net <- .daily_net()
  layout_with <- function(trim, fill = "pad") {
    .period_layout(net$data, net$metadata$student, net$metadata$day, fill, trim)
  }
  full <- layout_with(NULL)
  per_actor <- function(lay) as.vector(tabulate(lay$cells$actor, nbins = 3L))
  # s1 has 5 actions, s2 9, s3 6.
  expect_identical(per_actor(full), c(5L, 9L, 6L))
  expect_identical(layout_with(c(NA, 4))$actors, full$actors)
  expect_identical(per_actor(layout_with(c(NA, 4))), c(4L, 4L, 4L))
  # ceiling(0.5 * 5) = 3, ceiling(0.5 * 9) = 5, ceiling(0.5 * 6) = 3.
  expect_identical(per_actor(layout_with(c(NA, 0.5))), c(3L, 5L, 3L))
  # Counted over the kept periods: on d1 alone s1 has 1, s2 2, s3 3.
  d1 <- layout_with(c(1, 2))
  expect_identical(d1$actors, full$actors)
  expect_identical(d1$period_levels, "d1")
  expect_identical(per_actor(d1), c(1L, 2L, 2L))
  # Kept cells sit exactly where they were: nothing spreads to fill the row.
  for_key <- function(lay) paste(lay$cells$actor, lay$cells$period,
                                 lay$cells$xleft, lay$cells$xright)
  stretched_full <- layout_with(NULL, "stretch")
  stretched_half <- layout_with(c(NA, 0.5), "stretch")
  expect_true(all(for_key(stretched_half) %in% for_key(stretched_full)))
  # They are the earliest actions of each actor.
  first_of <- function(lay, who) {
    cells <- lay$cells[lay$cells$actor == who, ]
    paste(cells$period, cells$step)[order(cells$period, cells$step)]
  }
  expect_identical(first_of(stretched_half, 2L),
                   utils::head(first_of(stretched_full, 2L), 5L))
  # Through the verb, with group: every actor keeps its row and panel.
  d <- .daily_long()
  d$klass <- ifelse(d$student == "s1", "B", "A")
  res <- sequence_plot(d, actor = "student", action = "action",
                       periods = "day", type = "index", group = "klass",
                       trim = c(NA, 2))
  expect_setequal(res$groups, c("A", "B"))
  expect_identical(rownames(res$codes), c("s1", "s2", "s3"))
  expect_error(sequence_plot(net, periods = "day", trim = c(NA, -1)))
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

test_that("a fractional trim removes the weakest periods wherever they sit", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  # Actions per day: d1 6, d2 1, d3 none, d4 6, d5 1, d6 4.
  ids <- c("a", "b", "c")
  d <- rbind(
    data.frame(id = rep(ids, each = 2L), day = 1L, act = "x"),
    data.frame(id = "a", day = 2L, act = "y"),
    data.frame(id = rep(ids, each = 2L), day = 4L, act = "x"),
    data.frame(id = "b", day = 5L, act = "y"),
    data.frame(id = rep(c("a", "c"), each = 2L), day = 6L, act = "x")
  )
  d$day <- as.Date("2024-01-01") + d$day - 1L
  days_with <- function(trim) {
    sequence_plot(d, actor = "id", action = "act", periods = "day",
                  type = "index", trim = trim)$periods
  }
  day <- function(i) as.character(as.Date("2024-01-01") + i - 1L)
  # No fraction: the empty day stays on the axis.
  expect_identical(days_with(NULL), day(1:6))
  # Actions per day: d1 6, d2 1, d4 6, d5 1, d6 4; 18 in all, d3 empty.
  # Busiest first they add up to 6, 12, 16, 17, 18.
  # 0.6 needs 10.8 actions: d1 and d4 hold 12.
  expect_identical(days_with(0.6), day(c(1, 4)))
  # 0.75 needs 13.5: d6 joins.
  expect_identical(days_with(0.75), day(c(1, 4, 6)))
  # 0.9 needs 16.2: d2 joins, the earlier of the two days with one action.
  expect_identical(days_with(0.9), day(c(1, 2, 4, 6)))
  # All of it: only the empty day goes.
  expect_identical(days_with(0.999), day(c(1, 2, 4, 5, 6)))
  # Equal activity: the earlier day stays.
  expect_identical(days_with(0.4), day(c(1, 4)))
  expect_identical(days_with(0.2), day(1))
  # A count is still the first n days, the empty one included.
  expect_identical(days_with(3), day(1:3))

  # The kept cells sit in the closed-up blocks; rows are counted over them.
  # Under 0.75, a has 2 + 2 + 2 kept actions, b 2 + 2, c 2 + 2 + 2.
  prep <- prepare(d, actor = "id", action = "act", session = "day",
                  time_threshold = FALSE)
  lay <- .period_layout(prep$sequence_data, prep$meta_data$id,
                        prep$meta_data$day, "pad", trim = c(0.75, 5))
  expect_identical(lay$actors, ids)
  expect_identical(as.vector(tabulate(lay$cells$actor, nbins = 3L)),
                   c(5L, 4L, 5L))
  expect_setequal(lay$cells$period, 1:3)
  # A row cut that empties a period clears it under a fractional column cut:
  # one action each leaves only d1.
  one <- .period_layout(prep$sequence_data, prep$meta_data$id,
                        prep$meta_data$day, "pad", trim = c(0.6, 1))
  expect_identical(one$period_levels, day(1))
})

test_that("the row half of trim keeps cell widths and row order in periods", {
  d <- rbind(c("x", "y", NA, NA, NA, NA), rep("x", 6),
             c("y", NA, NA, NA, NA, NA))
  actor <- c("a", "a", "b")
  period <- c(1, 2, 1)
  full <- .period_layout(d, actor, period, "pad", trim = NULL)
  cut <- .period_layout(d, actor, period, "pad", trim = c(NA, 1))
  # The row of six ends after one action; the width it set stays.
  expect_identical(cut$width, full$width)
  kept <- merge(cut$cells, full$cells, by = c("actor", "period", "step"))
  expect_identical(nrow(kept), nrow(cut$cells))
  expect_equal(kept$xleft.x, kept$xleft.y)
  expect_equal(kept$xright.x, kept$xright.y)
  # Rows are ordered from what they held before the cut.
  expect_identical(cut$trajectory, full$trajectory)
})
