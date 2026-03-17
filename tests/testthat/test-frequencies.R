# ---- frequencies() tests ----

test_that("frequencies works with long format data", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2, 2),
    Time = c(1, 2, 3, 1, 2, 3),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )

  freq <- frequencies(long_data, action = "Action", id = "Actor")

  expect_true(is.matrix(freq))
  expect_equal(nrow(freq), 2)
  expect_equal(ncol(freq), 2)
  expect_equal(rownames(freq), c("A", "B"))
  expect_equal(colnames(freq), c("A", "B"))
  expect_equal(freq["A", "B"], 2L)
  expect_equal(freq["B", "A"], 2L)
  expect_equal(freq["A", "A"], 0L)
  expect_equal(freq["B", "B"], 0L)
})

test_that("frequencies works with wide format data", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  freq <- frequencies(wide_data, format = "wide")

  expect_true(is.matrix(freq))
  expect_equal(freq["A", "B"], 2L)
  expect_equal(freq["B", "A"], 2L)
})

test_that("frequencies auto-detects format", {
  long_data <- data.frame(
    Actor = c(1, 1, 2, 2),
    Time = c(1, 2, 1, 2),
    Action = c("X", "Y", "Y", "X"),
    stringsAsFactors = FALSE
  )
  freq <- frequencies(long_data, action = "Action", id = "Actor")
  expect_equal(freq["X", "Y"], 1L)
  expect_equal(freq["Y", "X"], 1L)

  wide_data <- data.frame(T1 = c("X", "Y"), T2 = c("Y", "X"),
                           stringsAsFactors = FALSE)
  freq2 <- frequencies(wide_data)
  expect_equal(freq2["X", "Y"], 1L)
  expect_equal(freq2["Y", "X"], 1L)
})

test_that("frequencies handles multiple ID columns", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 1),
    Group = c("G1", "G1", "G2", "G2"),
    Time = c(1, 2, 1, 2),
    Action = c("A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )

  freq <- frequencies(long_data, action = "Action", id = c("Actor", "Group"))
  expect_equal(freq["A", "B"], 1L)
  expect_equal(freq["B", "A"], 1L)
})

test_that("frequencies handles NAs", {
  wide_data <- data.frame(
    T1 = c("A", "B"), T2 = c("B", NA), T3 = c(NA, NA),
    stringsAsFactors = FALSE
  )
  freq <- frequencies(wide_data, format = "wide")
  expect_equal(freq["A", "B"], 1L)
  expect_equal(sum(freq), 1L)

  long_data <- data.frame(
    Actor = c(1, 1, 1), Time = c(1, 2, 3),
    Action = c("A", NA, "B"), stringsAsFactors = FALSE
  )
  freq2 <- frequencies(long_data, action = "Action", id = "Actor")
  expect_equal(sum(freq2), 0L)
})

test_that("frequencies with cols parameter selects specific columns", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"), T3 = c("C", "C"),
    extra = c("ignore", "this"), stringsAsFactors = FALSE
  )
  freq <- frequencies(wide_data, cols = c("T1", "T2", "T3"), format = "wide")
  expect_true("C" %in% rownames(freq))
  expect_equal(freq["A", "B"], 1L)
  expect_equal(freq["B", "C"], 1L)
})

test_that("frequencies returns integer matrix", {
  wide_data <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"),
                           stringsAsFactors = FALSE)
  freq <- frequencies(wide_data, format = "wide")
  expect_true(is.integer(freq))
})

test_that("frequencies long format without id treats all as one sequence", {
  long_data <- data.frame(
    Time = 1:4, Action = c("A", "B", "A", "B"), stringsAsFactors = FALSE
  )
  freq <- frequencies(long_data, action = "Action")
  expect_equal(freq["A", "B"], 2L)
  expect_equal(freq["B", "A"], 1L)
})

test_that("frequencies errors on missing columns", {
  long_data <- data.frame(x = 1:3, y = c("A", "B", "C"))
  expect_error(frequencies(long_data, action = "Action", format = "long"),
               "not found")
  expect_error(frequencies(long_data, action = "y", id = "missing", format = "long"),
               "not found")
})

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

test_that("convert_sequence_format frequency format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"), T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "frequency")

  expect_true(is.data.frame(result))
  expect_true("id" %in% names(result))
  expect_true("A" %in% names(result))
  expect_true("B" %in% names(result))
  # Row 1 (id=1): A, B, A => A=2, B=1
  expect_equal(result[result$id == 1, "A"], 2L)
  expect_equal(result[result$id == 1, "B"], 1L)
  # Row 2 (id=2): B, A, B => B=2, A=1
  expect_equal(result[result$id == 2, "A"], 1L)
  expect_equal(result[result$id == 2, "B"], 2L)
})

test_that("convert_sequence_format onehot format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "B"), T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "onehot")

  expect_true(is.data.frame(result))
  # Row 1 (id=1): A and B present
  expect_equal(result[result$id == 1, "A"], 1L)
  expect_equal(result[result$id == 1, "B"], 1L)
  # Row 2 (id=2): only B present
  expect_equal(result[result$id == 2, "A"], 0L)
  expect_equal(result[result$id == 2, "B"], 1L)
})

test_that("convert_sequence_format edgelist format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "edgelist")

  expect_true(is.data.frame(result))
  expect_true(all(c("id", "from", "to") %in% names(result)))
  expect_equal(nrow(result), 2)
  # Row 1: A->B, Row 2: B->A
  expect_true(any(result$from == "A" & result$to == "B"))
  expect_true(any(result$from == "B" & result$to == "A"))
})

test_that("convert_sequence_format follows format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "follows")

  expect_true(is.data.frame(result))
  expect_true(all(c("id", "act", "follows") %in% names(result)))
  expect_equal(nrow(result), 2)
  # B follows A, A follows B
  expect_true(any(result$act == "B" & result$follows == "A"))
  expect_true(any(result$act == "A" & result$follows == "B"))
})

test_that("convert_sequence_format frequency format from long", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2),
    Time = c(1, 2, 3, 1, 2),
    Action = c("A", "B", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(long_data, action = "Action",
                                    id_col = "Actor", time = "Time",
                                    format = "frequency")

  expect_true(is.data.frame(result))
  expect_true("Actor" %in% names(result))
  # Actor 1: A=2, B=1
  expect_equal(result[result$Actor == 1, "A"], 2L)
  expect_equal(result[result$Actor == 1, "B"], 1L)
  # Actor 2: A=0, B=2
  expect_equal(result[result$Actor == 2, "A"], 0L)
  expect_equal(result[result$Actor == 2, "B"], 2L)
})

test_that("convert_sequence_format edgelist format from long", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2),
    Time = c(1, 2, 3, 1, 2),
    Action = c("A", "B", "A", "B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(long_data, action = "Action",
                                    id_col = "Actor", time = "Time",
                                    format = "edgelist")

  expect_true(all(c("Actor", "from", "to") %in% names(result)))
  # Actor 1: A->B, B->A  Actor 2: B->A
  expect_equal(nrow(result), 3)
})

test_that("convert_sequence_format handles NAs in wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", NA), T3 = c(NA, NA),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "frequency")

  # id=1: A=1, B=1   id=2: B=1
  expect_equal(result[result$id == 1, "A"], 1L)
  expect_equal(result[result$id == 1, "B"], 1L)
  expect_equal(result[result$id == 2, "B"], 1L)
})

test_that("convert_sequence_format with seq_cols selects specific columns", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    extra = c("X", "Y"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, seq_cols = c("T1", "T2"),
                                    id_col = "id", format = "frequency")

  # "extra" column values should NOT appear as state columns
  expect_false("X" %in% names(result))
  expect_false("Y" %in% names(result))
})

test_that("convert_sequence_format defaults id_col to first column", {
  wide_data <- data.frame(
    Student = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, format = "frequency")

  expect_true("Student" %in% names(result))
})

test_that("convert_sequence_format with multiple id_col from long", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 1),
    Group = c("G1", "G1", "G2", "G2"),
    Time = c(1, 2, 1, 2),
    Action = c("A", "B", "C", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(long_data, action = "Action",
                                    id_col = c("Actor", "Group"),
                                    time = "Time", format = "edgelist")

  expect_true(all(c("Actor", "Group", "from", "to") %in% names(result)))
  expect_equal(nrow(result), 2)
})

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

test_that("convert_sequence_format errors on missing columns", {
  data <- data.frame(x = 1:3, y = c("A", "B", "C"))

  expect_error(
    convert_sequence_format(data, seq_cols = c("missing1", "missing2")),
    "not found"
  )
})
