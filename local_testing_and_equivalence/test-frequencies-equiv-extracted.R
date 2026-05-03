# Equivalence test_that() blocks extracted from
# tests/testthat/test-frequencies.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- frequencies() tests ----

test_that("frequencies works with tna package data", {
  skip_if_not_installed("tna")

  freq_long <- frequencies(tna::group_regulation_long,
                           action = "Action", id = "Actor")
  expect_true(is.matrix(freq_long))
  expect_true(is.integer(freq_long))
  expect_equal(nrow(freq_long), 9)

  freq_wide <- frequencies(tna::group_regulation, format = "wide")
  expect_equal(nrow(freq_wide), 9)
  expect_identical(freq_long, freq_wide)
})


# ---- convert_sequence_format() tests ----

test_that("convert_sequence_format works with tna package data (wide)", {
  skip_if_not_installed("tna")

  result <- convert_sequence_format(
    cbind(id = seq_len(nrow(tna::group_regulation)), tna::group_regulation),
    id_col = "id",
    format = "frequency"
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2000)
  # Should have 9 state columns + id + rid
  state_cols <- setdiff(names(result), c("id", "rid"))
  expect_equal(length(state_cols), 9)
})

test_that("convert_sequence_format works with tna package data (long)", {
  skip_if_not_installed("tna")

  result <- convert_sequence_format(
    tna::group_regulation_long,
    action = "Action", id_col = "Actor", time = "Time",
    format = "frequency"
  )

  expect_true(is.data.frame(result))
  expect_true("Actor" %in% names(result))
  state_cols <- setdiff(names(result), c("Actor", "rid"))
  expect_equal(length(state_cols), 9)
})

