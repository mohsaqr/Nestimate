# ---- sequence_metrics() ----

# Rows: A A B A | B B B | A B B B B
.sm_data <- function() {
  data.frame(T1 = c("A", "B", "A"), T2 = c("A", "B", "B"),
             T3 = c("B", "B", "B"), T4 = c("A", NA, "B"),
             T5 = c(NA, NA, "B"), stringsAsFactors = FALSE)
}

test_that("composition, spells and transitions match hand values", {
  r <- sequence_metrics(.sm_data())
  expect_identical(r$valid_n, c(4L, 3L, 5L))
  expect_equal(r$mean_spell_duration, c(4 / 3, 3, 2.5))
  # A A B A: proportions 3/4 and 1/4, two states
  expect_equal(r$longitudinal_entropy[1],
               -(0.75 * log(0.75) + 0.25 * log(0.25)) / log(2))
  expect_equal(r$simpson_diversity[1], 1 - (0.75^2 + 0.25^2))
  # transitions AA, AB, BA: one self-loop, two changes over three steps
  expect_equal(r$self_loop_tendency[1], 1 / 3)
  expect_equal(r$transition_rate[1], 2 / 3)
  expect_equal(r$transition_complexity, c(1, 0, 0.5))
})

test_that("initial-state, dominant and cyclic indices match hand values", {
  r <- sequence_metrics(.sm_data())
  expect_equal(r$initial_state_persistence, c(2 / 4, 1, 1 / 5))
  expect_equal(r$initial_state_influence_decay, c(1 - 1 / 2, 0, 1 / 2 - 0))
  expect_equal(r$cyclic_feedback_strength, c(1, 1, 2 / 3))
  expect_identical(r$dominant_state, c("A", "B", "B"))
  expect_equal(r$dominant_max_spell, c(2, 3, 4))
})

test_that("an emergent state needs a long spell that is not the opening one", {
  r <- sequence_metrics(.sm_data())
  expect_identical(r$emergent_state, c(NA, NA, "B"))
  expect_equal(r$emergent_state_persistence[3], 4)
  expect_equal(r$emergent_state_proportion[3], 4 / 5)
})

test_that("complexity and integrative potential match hand values", {
  r <- sequence_metrics(.sm_data(), favorable = "B")
  entropy <- -(0.75 * log(0.75) + 0.25 * log(0.25)) / log(2)
  spells <- c(2, 1, 1)
  expect_equal(r$complexity_index[1],
               0.4 * entropy + 0.4 * (2 / 3) + 0.2 * sd(spells) / mean(spells))
  # one spell: one state, no change, so the least complex sequence
  expect_identical(r$complexity_index[2], 0)
  expect_identical(sequence_metrics(data.frame(T1 = "A"))$complexity_index, 0)
  # position weights 1..4, B at position 3
  expect_equal(r$integrative_potential, c(3 / 10, 1, 14 / 15))
})

test_that("an empty sequence gives NA indices without warnings", {
  d <- rbind(.sm_data(), data.frame(T1 = NA, T2 = NA, T3 = NA, T4 = NA,
                                    T5 = NA))
  expect_no_warning(r <- sequence_metrics(d))
  expect_identical(r$valid_n[4], 0L)
  expect_true(is.na(r$transition_rate[4]))
  expect_true(is.na(r$dominant_state[4]))
})

test_that("rows follow the sequences, covariates on request", {
  d <- .sm_data()
  d$grp <- c("x", "y", "z")
  reversed <- d[3:1, ]
  expect_equal(sequence_metrics(reversed)$valid_n, c(5L, 3L, 4L))
  with_meta <- sequence_metrics(d, meta = TRUE)
  expect_identical(names(with_meta)[1], "grp")
  expect_identical(with_meta$grp, c("x", "y", "z"))
})

test_that("missing positions are not observations", {
  # - - A A B | A A - - B B B | A - A A
  d <- data.frame(T1 = c(NA, "A", "A"), T2 = c(NA, "A", NA),
                  T3 = c("A", NA, "A"), T4 = c("A", NA, "A"),
                  T5 = c("B", "B", NA), T6 = c(NA, "B", NA),
                  T7 = c(NA, "B", NA), stringsAsFactors = FALSE)
  r <- sequence_metrics(d, favorable = "A")
  # opening spell starts at the first observed state, ends at a gap
  expect_equal(r$initial_state_persistence, c(2 / 5, 2 / 7, 1 / 4))
  # transitions are adjacent observed pairs: A->B across the gap is none
  expect_equal(r$transition_rate, c(1 / 2, 0, 0))
  expect_equal(r$self_loop_tendency + r$transition_rate, c(1, 1, 1))
  expect_equal(r$mean_spell_duration, c(1.5, 2.5, 1.5))
  # weights over observed positions only: row 2 has A at 1, 2 of 1, 2, 5, 6, 7
  expect_equal(r$integrative_potential, c(7 / 12, 3 / 21, 1))
})

test_that("favourable states that never occur are refused by class", {
  expect_error(sequence_metrics(.sm_data(), favorable = "Z"),
               class = "nestimate_bad_input")
})

test_that("a return to the first state emerges above half the opening spell", {
  # opening spell 5, return spell 3: 2 * 3 > 5, so A emerges (codyna's rule)
  emerges <- as.data.frame(t(c(rep("A", 5), "B", rep("A", 3))))
  expect_identical(sequence_metrics(emerges)$emergent_state, "A")
  expect_identical(sequence_metrics(emerges)$emergent_state_persistence, 3)
  # opening spell 7, return spell 3: 2 * 3 < 7, no emergent state
  stays <- as.data.frame(t(c(rep("A", 7), "B", rep("A", 3))))
  expect_true(is.na(sequence_metrics(stays)$emergent_state))
})

# ---- group = : comparing the indices between groups ----

.sm_groups_data <- function() {
  data.frame(T1 = c("A", "A", "B", "A", "B", "B", "A", "B"),
             T2 = c("A", "B", "B", "A", "A", "B", "B", "A"),
             T3 = c("B", "B", "A", "A", "B", "A", "B", "A"),
             T4 = c("B", "A", "A", "B", "A", "A", "A", "B"),
             cond = rep(c("x", "y"), each = 4),
             stringsAsFactors = FALSE)
}

test_that("sequence_metrics(group =) reproduces hand-computed means and Hedges' g", {
  d <- .sm_groups_data()
  per_seq <- sequence_metrics(d, meta = TRUE)
  set.seed(1)
  cmp <- sequence_metrics(d, group = "cond", iter = 50, format = "wide")
  expect_true(is.data.frame(cmp))
  row <- cmp[cmp$index == "transition_rate", ]
  x <- per_seq$transition_rate[per_seq$cond == "x"]
  y <- per_seq$transition_rate[per_seq$cond == "y"]
  expect_equal(row$mean_x, mean(x))
  expect_equal(row$sd_y, sd(y))
  expect_equal(row$difference, mean(y) - mean(x))
  # F against stats::t.test(): tests/equivalence/test-sequence_metrics.R
  # Hedges' g by hand
  sp <- sqrt(((4 - 1) * var(x) + (4 - 1) * var(y)) / 6)
  g <- (1 - 3 / (4 * 8 - 9)) * (mean(y) - mean(x)) / sp
  expect_equal(row$hedges_g, g)
  expect_true(row$g_lower < g && g < row$g_upper)
  # p-values are permutation p-values, never below 1 / (iter + 1)
  ok <- !is.na(cmp$p_value)
  expect_true(all(cmp$p_value[ok] >= 1 / 51 & cmp$p_value[ok] <= 1))
  expect_equal(cmp$p_adj, stats::p.adjust(cmp$p_value, "BH"))
})

test_that("sequence_metrics(group =) leaves the per-sequence output unchanged without group", {
  d <- .sm_groups_data()
  expect_identical(names(sequence_metrics(d)),
                   names(sequence_metrics(d[, 1:4])))
  # invariant: eta squared lies in [0, 1]
  set.seed(2)
  cmp <- sequence_metrics(d, group = "cond", iter = 20)
  ok <- !is.na(cmp$eta_squared)
  expect_true(all(cmp$eta_squared[ok] >= 0 & cmp$eta_squared[ok] <= 1))
})

test_that("a grouped network is compared by its own groups, nested in its actor", {
  halves <- mark_phases(human_long, actor = "session_id",
                        time = "order_in_session", labels = c("early", "late"))
  net <- build_network(halves, actor = "session_id", action = "cluster",
                       group = "phase", method = "relative")
  set.seed(3)
  by_default <- sequence_metrics(net, iter = 20)
  set.seed(3)
  explicit <- sequence_metrics(net, group = "phase", actor = "session_id",
                               iter = 20)
  # the defaults are the network's group column and its actor (recorded as
  # the network's session label, which names the same sessions)
  expect_identical(summary(by_default), summary(explicit))
  expect_identical(by_default$index[1], "valid_n")
  expect_true(all(c("index", "group", "n", "mean", "sd", "hedges_g",
                    "p_adj") %in% names(by_default)))
  unique_states <- subset(by_default, index == "unique_states")
  expect_identical(unique_states$group, c("early", "late"))
  expect_identical(unique_states$n, c(425L, 429L))
  # meta = TRUE still lists the sequences, one row each
  per_seq <- sequence_metrics(net, meta = TRUE)
  expect_equal(nrow(per_seq), 854L)
  expect_true(all(c("phase", "session_id", "transition_rate") %in%
                    names(per_seq)))
})

test_that("sequence_metrics(group =) refuses an unknown or single group", {
  d <- .sm_groups_data()
  expect_error(sequence_metrics(d, group = "nope"), class = "nestimate_bad_input")
  d$cond <- "x"
  expect_error(sequence_metrics(d, group = "cond"), class = "nestimate_bad_input")
})

test_that("indices = returns the named indices only, per sequence and compared", {
  d <- .sm_groups_data()
  per_seq <- sequence_metrics(d, indices = c("transition_rate", "valid_n"))
  expect_identical(names(per_seq), c("transition_rate", "valid_n"))
  expect_equal(per_seq$valid_n, sequence_metrics(d)$valid_n)
  set.seed(4)
  cmp <- sequence_metrics(d, group = "cond", indices = "transition_rate",
                          iter = 9)
  expect_identical(cmp$index, c("transition_rate", "transition_rate"))
  expect_error(sequence_metrics(d, indices = "nope"), class = "nestimate_bad_input")
})

test_that("the group comparison is a data frame, one row per index and group", {
  d <- .sm_groups_data()
  set.seed(5)
  long <- sequence_metrics(d, group = "cond", iter = 20)
  set.seed(5)
  wide <- sequence_metrics(d, group = "cond", iter = 20, format = "wide")
  expect_s3_class(long, "net_sequence_metrics")
  expect_s3_class(long, "data.frame")
  expect_identical(names(long)[1:5], c("index", "group", "n", "mean", "sd"))
  expect_identical(nrow(long), 2L * nrow(wide))
  expect_identical(attr(long, "groups"), c("x", "y"))
  # invariant: each group's long rows are that group's wide columns
  x_rows <- subset(long, group == "x")
  y_rows <- subset(long, group == "y")
  expect_identical(x_rows$index, wide$index)
  expect_identical(x_rows$n, wide$n_x)
  expect_identical(y_rows$mean, wide$mean_y)
  expect_identical(x_rows$sd, wide$sd_x)
  expect_identical(y_rows$p_adj, wide$p_adj)
  expect_identical(x_rows$hedges_g, wide$hedges_g)
  # every row as a plain data frame; a subset is plain too
  expect_identical(class(summary(long)), "data.frame")
  expect_identical(dim(summary(long)), dim(long))
  expect_identical(class(head(long, 2)), "data.frame")
  expect_identical(class(x_rows), "data.frame")
})

test_that("print() of the group comparison is rounded and has a group column", {
  d <- .sm_groups_data()
  set.seed(6)
  cmp <- sequence_metrics(d, group = "cond", iter = 20)
  shown <- capture.output(print(cmp, n = 2L))
  expect_match(shown[1L], "Sequence metrics by cond  \\[\\d+ indices \\| 2 groups: x, y\\]")
  header <- shown[grepl("^ index", shown)]
  expect_length(header, 1L)
  expect_match(header, "group")
  rows <- shown[grepl("^ valid_n ", shown)]
  expect_length(rows, 2L)
  expect_true(any(grepl("more indices", shown)))
  capture.output(shown_value <- withVisible(print(cmp)))
  expect_false(shown_value$visible)
  expect_error(print(cmp, n = 0), class = "nestimate_bad_input")
  expect_error(sequence_metrics(d, group = "cond", format = "tall"),
               "should be one of")
})
