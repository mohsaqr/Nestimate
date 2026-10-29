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
