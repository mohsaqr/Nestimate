# Tests for R/data_conversion.R
# wide_to_long, long_to_wide, prepare_for_tna, action_to_onehot, prepare_onehot

# --- Fixtures -----------------------------------------------------------

make_wide <- function() {
  data.frame(
    V1 = c("A", "B", "C"),
    V2 = c("B", "C", "A"),
    V3 = c("C", "A", "B"),
    stringsAsFactors = FALSE
  )
}

make_long <- function() {
  data.frame(
    Actor = rep(1:3, each = 3),
    Time  = rep(1:3, 3),
    Action = c("A", "B", "C", "B", "C", "A", "C", "A", "B"),
    stringsAsFactors = FALSE
  )
}

# --- wide_to_long --------------------------------------------------------

test_that("wide_to_long basic conversion", {
  w <- make_wide()
  l <- wide_to_long(w)
  expect_s3_class(l, "data.frame")
  expect_equal(nrow(l), 9)
  expect_true(all(c("id", "Time", "Action") %in% names(l)))
  # Row 1 sequence: A, B, C
  sub1 <- l[l$id == 1, ]
  expect_equal(sub1$Action, c("A", "B", "C"))
  expect_equal(sub1$Time, 1:3)
})

test_that("wide_to_long auto-generates IDs when id_col is NULL", {
  w <- make_wide()
  l <- wide_to_long(w)
  expect_equal(sort(unique(l$id)), 1:3)
})

test_that("wide_to_long uses existing id_col", {
  w <- make_wide()
  w$person <- c("x", "y", "z")
  l <- wide_to_long(w, id_col = "person")
  expect_true("person" %in% names(l))
  expect_equal(sort(unique(l$person)), c("x", "y", "z"))
  expect_false("id" %in% names(l))
})

test_that("wide_to_long drops NA by default", {
  w <- data.frame(V1 = c("A", "B"), V2 = c("B", NA), stringsAsFactors = FALSE)
  l <- wide_to_long(w, drop_na = TRUE)
  expect_equal(nrow(l), 3)
  expect_false(any(is.na(l$Action)))
})

test_that("wide_to_long preserves NA when drop_na = FALSE", {
  w <- data.frame(V1 = c("A", "B"), V2 = c("B", NA), stringsAsFactors = FALSE)
  l <- wide_to_long(w, drop_na = FALSE)
  expect_equal(nrow(l), 4)
  expect_equal(sum(is.na(l$Action)), 1)
})

test_that("wide_to_long preserves extra columns", {
  w <- make_wide()
  w$group <- c("g1", "g1", "g2")
  l <- wide_to_long(w)
  expect_true("group" %in% names(l))
  expect_equal(l$group[l$id == 1], rep("g1", 3))
})

test_that("wide_to_long custom column names", {
  w <- data.frame(T1 = c("X", "Y"), T2 = c("Y", "X"), stringsAsFactors = FALSE)
  l <- wide_to_long(w, time_prefix = "T", action_col = "State", time_col = "Step")
  expect_true(all(c("State", "Step") %in% names(l)))
  expect_equal(nrow(l), 4)
})

test_that("wide_to_long errors on missing time prefix", {
  w <- data.frame(x = 1:3, y = 4:6)
  expect_error(wide_to_long(w), "No columns matching")
})

test_that("wide_to_long errors on invalid id_col", {
  w <- make_wide()
  expect_error(wide_to_long(w, id_col = "nonexistent"), "not found")
})

# --- long_to_wide --------------------------------------------------------

test_that("long_to_wide basic conversion", {
  l <- make_long()
  w <- long_to_wide(l)
  expect_s3_class(w, "data.frame")
  expect_equal(nrow(w), 3)
  expect_true(all(c("V1", "V2", "V3") %in% names(w)))
  # Check first actor's sequence
  row1 <- w[w$Actor == 1, ]
  expect_equal(as.character(row1$V1), "A")
  expect_equal(as.character(row1$V3), "C")
})

test_that("long_to_wide handles unequal sequence lengths", {
  l <- data.frame(
    Actor  = c(1, 1, 1, 2, 2),
    Time   = c(1, 2, 3, 1, 2),
    Action = c("A", "B", "C", "X", "Y"),
    stringsAsFactors = FALSE
  )
  w <- long_to_wide(l)
  expect_equal(nrow(w), 2)
  expect_true("V3" %in% names(w))
  expect_true(is.na(w$V3[w$Actor == 2]))
})

test_that("long_to_wide without Time column uses row order", {
  l <- data.frame(
    Actor  = c(1, 1, 2, 2),
    Action = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )
  w <- long_to_wide(l, time_col = "Time")
  expect_equal(nrow(w), 2)
  expect_equal(as.character(w$V1[w$Actor == 1]), "A")
  expect_equal(as.character(w$V2[w$Actor == 1]), "B")
})

test_that("long_to_wide custom prefix", {
  l <- make_long()
  w <- long_to_wide(l, time_prefix = "T")
  expect_true(all(c("T1", "T2", "T3") %in% names(w)))
})

test_that("long_to_wide errors on missing columns", {
  l <- data.frame(x = 1:3, y = letters[1:3])
  expect_error(long_to_wide(l), "Missing required columns")
})

# --- Round-trip tests ----------------------------------------------------

test_that("wide -> long -> wide round-trip preserves data", {
  w_orig <- make_wide()
  l <- wide_to_long(w_orig)
  w_back <- long_to_wide(l, id_col = "id")
  # Compare just the V columns
  vcols <- grep("^V[0-9]+$", names(w_back), value = TRUE)
  recovered <- w_back[order(w_back$id), vcols]
  rownames(recovered) <- NULL
  expect_equal(recovered, w_orig)
})

test_that("long -> wide -> long round-trip preserves data", {
  l_orig <- make_long()
  w <- long_to_wide(l_orig)
  l_back <- wide_to_long(w, id_col = "Actor")
  l_back <- l_back[order(l_back$Actor, l_back$Time), ]
  rownames(l_back) <- NULL
  expect_equal(l_back$Actor, l_orig$Actor)
  expect_equal(l_back$Action, l_orig$Action)
})

# --- prepare_for_tna -----------------------------------------------------

test_that("prepare_for_tna sequences type returns only V columns", {
  w <- make_wide()
  w$extra <- "meta"
  result <- prepare_for_tna(w, type = "sequences")
  expect_equal(names(result), c("V1", "V2", "V3"))
  expect_false("extra" %in% names(result))
})

test_that("prepare_for_tna long type converts to wide", {
  l <- make_long()
  result <- prepare_for_tna(l, type = "long")
  expect_true(all(grepl("^V[0-9]+$", names(result))))
  expect_equal(nrow(result), 3)
})

test_that("prepare_for_tna auto-detects wide", {
  w <- make_wide()
  result <- prepare_for_tna(w, type = "auto")
  expect_equal(ncol(result), 3)
})

test_that("prepare_for_tna auto-detects long", {
  l <- make_long()
  result <- prepare_for_tna(l, type = "auto")
  expect_true(all(grepl("^V[0-9]+$", names(result))))
})

test_that("prepare_for_tna validates state names", {
  w <- make_wide()
  expect_warning(
    prepare_for_tna(w, type = "sequences", state_names = c("A", "B"),
                    validate = TRUE),
    "Unknown states"
  )
})

test_that("prepare_for_tna skips validation when validate = FALSE", {
  w <- make_wide()
  expect_no_warning(
    prepare_for_tna(w, type = "sequences", state_names = c("A", "B"),
                    validate = FALSE)
  )
})

test_that("prepare_for_tna converts factors to characters", {
  w <- data.frame(V1 = factor(c("A", "B")), V2 = factor(c("B", "A")))
  result <- prepare_for_tna(w, type = "sequences")
  expect_true(all(vapply(result, is.character, logical(1))))
})

test_that("prepare_for_tna errors on no V columns", {
  d <- data.frame(x = 1:3, y = letters[1:3])
  expect_error(prepare_for_tna(d, type = "sequences"), "No sequence columns")
})

test_that("prepare_for_tna auto errors when ambiguous", {
  d <- data.frame(x = 1:3)
  expect_error(prepare_for_tna(d, type = "auto"), "Cannot auto-detect")
})

test_that("prepare_for_tna long errors on missing action col", {
  d <- data.frame(Actor = 1:3, Time = 1:3, State = c("A", "B", "C"))
  expect_error(prepare_for_tna(d, type = "long"), "Action column")
})

# --- action_to_onehot ----------------------------------------------------

test_that("action_to_onehot basic encoding", {
  d <- data.frame(id = 1:4, Action = c("A", "B", "A", "C"),
                  stringsAsFactors = FALSE)
  result <- action_to_onehot(d)
  expect_true(all(c("A", "B", "C") %in% names(result)))
  expect_false("Action" %in% names(result))
  expect_equal(result$A, c(1L, 0L, 1L, 0L))
  expect_equal(result$C, c(0L, 0L, 0L, 1L))
})

test_that("action_to_onehot keeps action when drop_action = FALSE", {
  d <- data.frame(Action = c("X", "Y"), stringsAsFactors = FALSE)
  result <- action_to_onehot(d, drop_action = FALSE)
  expect_true("Action" %in% names(result))
})

test_that("action_to_onehot sorts states", {
  d <- data.frame(Action = c("C", "A", "B"), stringsAsFactors = FALSE)
  result <- action_to_onehot(d, sort_states = TRUE)
  state_cols <- setdiff(names(result), "Action")
  expect_equal(state_cols, c("A", "B", "C"))
})

test_that("action_to_onehot custom states includes extras as all-zero", {
  d <- data.frame(Action = c("A", "B"), stringsAsFactors = FALSE)
  result <- action_to_onehot(d, states = c("A", "B", "Z"))
  expect_true("Z" %in% names(result))
  expect_equal(result$Z, c(0L, 0L))
})

test_that("action_to_onehot prefix works", {
  d <- data.frame(Action = c("A", "B"), stringsAsFactors = FALSE)
  result <- action_to_onehot(d, prefix = "s_")
  expect_true(all(c("s_A", "s_B") %in% names(result)))
})

test_that("action_to_onehot errors on missing column", {
  d <- data.frame(State = c("A", "B"))
  expect_error(action_to_onehot(d), "action_col")
})

# --- prepare_onehot ------------------------------------------------------

test_that("prepare_onehot basic conversion", {
  d <- data.frame(A = c(1, 0, 1), B = c(0, 1, 0))
  result <- prepare_onehot(d, cols = c("A", "B"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  # Row 1 time 1: A active
  expect_equal(result[[1]], "A")
  expect_true(attr(result, "windowed") == FALSE)
})

test_that("prepare_onehot with actor grouping", {
  d <- data.frame(
    actor = c(1, 1, 2, 2),
    X = c(1, 0, 0, 1),
    Y = c(0, 1, 1, 0)
  )
  result <- prepare_onehot(d, cols = c("X", "Y"), actor = "actor")
  expect_equal(nrow(result), 2)
})

test_that("prepare_onehot non-overlapping window", {
  d <- data.frame(A = c(1, 0, 1, 0), B = c(0, 1, 0, 1))
  result <- prepare_onehot(d, cols = c("A", "B"), window_size = 2,
                           window_type = "non-overlapping")
  expect_true(attr(result, "windowed"))
  expect_equal(attr(result, "window_size"), 2L)
})

test_that("prepare_onehot overlapping window", {
  d <- data.frame(A = c(1, 0, 1, 0), B = c(0, 1, 0, 1))
  result <- prepare_onehot(d, cols = c("A", "B"), window_size = 2,
                           window_type = "overlapping")
  expect_true(attr(result, "windowed"))
})

test_that("prepare_onehot aggregate within windows", {
  d <- data.frame(A = c(1, 0, 1, 0), B = c(0, 1, 0, 1))
  result <- prepare_onehot(d, cols = c("A", "B"), window_size = 2,
                           window_type = "non-overlapping", aggregate = TRUE)
  expect_true(attr(result, "windowed"))
})

test_that("prepare_onehot errors on missing cols", {
  d <- data.frame(A = c(1, 0))
  expect_error(prepare_onehot(d, cols = c("A", "Z")), "cols")
})

test_that("prepare_onehot codes attribute set", {
  d <- data.frame(X = c(1, 0), Y = c(0, 1))
  result <- prepare_onehot(d, cols = c("X", "Y"))
  expect_equal(attr(result, "codes"), c("X", "Y"))
})
