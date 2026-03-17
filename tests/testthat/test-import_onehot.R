# ---- prepare_onehot() tests ----

test_that("prepare_onehot basic conversion", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 0, 0, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B", "C"))
  expect_true(is.data.frame(result))
  # Should have columns like W0_T1, W0_T2, etc.
  expect_true(any(grepl("^W\\d+_T\\d+$", names(result))))
  # A and B should appear; C should be NA (all zeros)
  vals <- unlist(result[, grepl("^W\\d+_T\\d+$", names(result))])
  expect_true("A" %in% vals)
  expect_true("B" %in% vals)
})

test_that("prepare_onehot with actor grouping", {
  df <- data.frame(
    actor = c(1, 1, 1, 2, 2),
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"), actor = "actor")
  expect_true(is.data.frame(result))
  # Should have at least 2 rows (one per actor)
  expect_true(nrow(result) >= 2)
})

test_that("prepare_onehot non-overlapping window", {
  df <- data.frame(
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           window_size = 2, window_type = "non-overlapping")
  expect_true(is.data.frame(result))
  expect_true(attr(result, "windowed"))
  expect_equal(attr(result, "window_size"), 2L)
})

test_that("prepare_onehot overlapping window", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           window_size = 2, window_type = "overlapping")
  expect_true(is.data.frame(result))
  expect_true(attr(result, "windowed"))
})

test_that("prepare_onehot aggregate mode", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           window_size = 2, window_type = "non-overlapping",
                           aggregate = TRUE)
  expect_true(is.data.frame(result))
})

test_that("prepare_onehot default actor/session (single group)", {
  df <- data.frame(
    A = c(1, 0, 1),
    B = c(0, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"))
  expect_true(is.data.frame(result))
  # Should have 1 row (single group)
  expect_equal(nrow(result), 1)
})

test_that("prepare_onehot attributes set correctly", {
  df <- data.frame(
    A = c(1, 0, 1),
    B = c(0, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"))
  expect_false(attr(result, "windowed"))
  expect_equal(attr(result, "window_size"), 1L)
  expect_equal(attr(result, "codes"), c("A", "B"))

  result2 <- prepare_onehot(df, cols = c("A", "B"), window_size = 2)
  expect_true(attr(result2, "windowed"))
})

test_that("prepare_onehot with session", {
  df <- data.frame(
    actor = c(1, 1, 1, 1),
    session = c("s1", "s1", "s2", "s2"),
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           actor = "actor", session = "session")
  expect_true(is.data.frame(result))
})

test_that("prepare_onehot output feeds into build_network", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1, 0),
    B = c(0, 1, 0, 1, 0, 1),
    C = c(0, 0, 0, 0, 0, 0),
    actor = c(1, 1, 1, 2, 2, 2)
  )

  result <- prepare_onehot(df, cols = c("A", "B", "C"), actor = "actor")
  # The wide format should be usable for build_network
  seq_cols <- grep("^W\\d+_T\\d+$", names(result), value = TRUE)
  if (length(seq_cols) >= 2) {
    net <- build_network(result, method = "relative",
                         params = list(format = "wide", cols = seq_cols))
    expect_s3_class(net, "netobject")
  }
})

test_that("prepare_onehot validates inputs", {
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  expect_error(prepare_onehot(df, cols = c("X", "Y")))
  expect_error(prepare_onehot(df, cols = "A", actor = "missing_col"))
})
