# ---- Synthetic fixtures ----

make_events <- function() {
  data.frame(
    student = c("s1", "s1", "s1", "s2", "s2"),
    code = c("A", "B", "C", "A", "B"),
    timestamp = as.POSIXct(c(
      "2024-01-01 10:00:00", "2024-01-01 10:05:00",
      "2024-01-01 10:10:00", "2024-01-01 10:00:00",
      "2024-01-01 10:05:00"
    )),
    stringsAsFactors = FALSE
  )
}

# ---- Basic: actor + action, no time ----

test_that("basic actor + action returns nestimate_data with correct structure", {
  ev <- make_events()
  res <- prepare_data(ev, actor = "student", action = "code")
  expect_s3_class(res, "nestimate_data")
  expect_named(res, c("sequence_data", "long_data", "meta_data",
                       "time_data", "statistics"))
  expect_null(res$time_data)
  expect_equal(res$statistics$total_sessions, 2L)
  expect_equal(res$statistics$total_actions, 5L)
  expect_equal(res$statistics$unique_actors, 2L)
  # Wide columns named T1, T2, ...
  expect_true(all(grepl("^T\\d+$", names(res$sequence_data))))
  expect_equal(nrow(res$sequence_data), 2L)
})

test_that("sequence order follows row order when no time given", {
  ev <- make_events()
  res <- prepare_data(ev, actor = "student", action = "code")
  s1_row <- which(res$meta_data[[grep("student|actor", names(res$meta_data),
                                       value = TRUE)[1]]] == "s1")
  expect_equal(as.character(res$sequence_data[s1_row, 1:3]), c("A", "B", "C"))
})

# ---- With time ----

test_that("time column is parsed and time_data is returned", {
  ev <- make_events()
  res <- prepare_data(ev, actor = "student", action = "code",
                      time = "timestamp")
  expect_false(is.null(res$time_data))
  expect_equal(ncol(res$time_data), res$statistics$max_sequence_length)
  expect_true(all(grepl("^time_T\\d+$", names(res$time_data))))
  expect_s3_class(res$time_data[[1]], "POSIXct")
})

# ---- Session splitting by time_threshold ----

test_that("time_threshold splits sessions within same actor", {
  ev <- data.frame(
    student = rep("s1", 4),
    code = c("A", "B", "C", "D"),
    timestamp = as.POSIXct(c(
      "2024-01-01 10:00:00", "2024-01-01 10:05:00",
      "2024-01-01 12:00:00", "2024-01-01 12:05:00"
    )),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code",
                      time = "timestamp", time_threshold = 900)
  expect_equal(res$statistics$total_sessions, 2L)
  expect_equal(nrow(res$sequence_data), 2L)
  # First session: A, B; Second session: C, D
  expect_equal(as.character(res$sequence_data[1, 1:2]), c("A", "B"))
  expect_equal(as.character(res$sequence_data[2, 1:2]), c("C", "D"))
})

test_that("large time_threshold keeps everything in one session", {
  ev <- make_events()
  res <- prepare_data(ev, actor = "student", action = "code",
                      time = "timestamp", time_threshold = 1e6)
  expect_equal(res$statistics$total_sessions, 2L)
})

# ---- Missing actor ----

test_that("missing actor treats all data as one actor", {
  ev <- data.frame(code = c("A", "B", "C"), stringsAsFactors = FALSE)
  res <- prepare_data(ev, action = "code")
  expect_equal(res$statistics$total_sessions, 1L)
  expect_null(res$statistics$unique_actors)
  expect_equal(as.character(res$sequence_data[1, ]), c("A", "B", "C"))
})

# ---- Multiple actor columns (interaction) ----

test_that("multiple actor columns create interaction grouping", {
  ev <- data.frame(
    student = c("s1", "s1", "s2", "s2"),
    group = c("g1", "g1", "g1", "g2"),
    code = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = c("student", "group"), action = "code")
  expect_equal(res$statistics$unique_actors, 3L)
  expect_equal(res$statistics$total_sessions, 3L)
})

# ---- Explicit session column ----

test_that("session column creates separate sessions per actor-session combo", {
  ev <- data.frame(
    student = c("s1", "s1", "s1", "s1"),
    course = c("math", "math", "bio", "bio"),
    code = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code",
                      session = "course")
  expect_equal(res$statistics$total_sessions, 2L)
})

test_that("multiple session columns create interaction sessions", {
  ev <- data.frame(
    student = rep("s1", 4),
    course = c("math", "math", "bio", "bio"),
    semester = c("fall", "spring", "fall", "fall"),
    code = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code",
                      session = c("course", "semester"))
  expect_equal(res$statistics$total_sessions, 3L)
})

# ---- Order column for tie-breaking ----

test_that("order column controls sequence within tied timestamps", {
  ev <- data.frame(
    student = rep("s1", 3),
    code = c("C", "A", "B"),
    timestamp = rep(as.POSIXct("2024-01-01 10:00:00"), 3),
    priority = c(3, 1, 2),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code",
                      time = "timestamp", order = "priority")
  expect_equal(as.character(res$sequence_data[1, 1:3]), c("A", "B", "C"))
})

# ---- Time parsing: ISO8601 string ----

test_that("ISO8601 T-separator timestamps are parsed", {
  ev <- data.frame(
    student = rep("s1", 2),
    code = c("A", "B"),
    ts = c("2024-01-01T10:00:00", "2024-01-01T10:05:00"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code", time = "ts")
  expect_false(is.null(res$time_data))
  expect_equal(nrow(res$sequence_data), 1L)
})

# ---- Time parsing: numeric auto-detected as unix ----

test_that("numeric time column auto-detected as unix timestamps", {
  ev <- data.frame(
    student = rep("s1", 3),
    code = c("A", "B", "C"),
    ts = c(1704100000, 1704100300, 1704100600),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code", time = "ts")
  expect_false(is.null(res$time_data))
  expect_s3_class(res$time_data[[1]], "POSIXct")
})

# ---- Time parsing: explicit unix milliseconds ----

test_that("unix milliseconds are correctly converted", {
  ev <- data.frame(
    student = rep("s1", 2),
    code = c("A", "B"),
    ts = c(1704100000000, 1704100300000),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code", time = "ts",
                      is_unix_time = TRUE, unix_time_unit = "milliseconds")
  t1 <- res$time_data[[1]][1]
  expect_equal(as.numeric(t1), 1704100000, tolerance = 1)
})

# ---- Custom time format ----

test_that("custom_format parses non-standard timestamps", {
  ev <- data.frame(
    student = rep("s1", 2),
    code = c("A", "B"),
    ts = c("01-Jan-2024 10:00", "01-Jan-2024 10:05"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code", time = "ts",
                      custom_format = "%d-%b-%Y %H:%M")
  expect_false(is.null(res$time_data))
})

# ---- Extra columns aggregated per session ----

test_that("numeric extra columns are aggregated as mean", {
  ev <- data.frame(
    student = c("s1", "s1", "s1"),
    code = c("A", "B", "C"),
    score = c(10, 20, 30),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code")
  expect_true("score" %in% names(res$meta_data))
  expect_equal(res$meta_data$score, 20)
})

test_that("character extra columns are aggregated as mode", {
  ev <- data.frame(
    student = c("s1", "s1", "s1"),
    code = c("A", "B", "C"),
    level = c("high", "high", "low"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code")
  expect_equal(res$meta_data$level, "high")
})

# ---- print method ----

test_that("print.nestimate_data produces expected output", {
  ev <- make_events()
  res <- prepare_data(ev, actor = "student", action = "code",
                      time = "timestamp")
  out <- capture.output(print(res))
  expect_true(any(grepl("Prepared Data", out)))
  expect_true(any(grepl("Sessions:", out)))
  expect_true(any(grepl("Actors:", out)))
  expect_true(any(grepl("Time data: available", out)))
})

test_that("print without actors omits Actors line", {
  ev <- data.frame(code = c("A", "B"), stringsAsFactors = FALSE)
  res <- prepare_data(ev, action = "code")
  out <- capture.output(print(res))
  expect_false(any(grepl("Actors:", out)))
})

# ---- Error conditions ----

test_that("non-data.frame input errors", {
  expect_error(prepare_data("not_a_df", actor = "a", action = "b"))
})

test_that("missing action column errors", {
  ev <- make_events()
  expect_error(prepare_data(ev, actor = "student", action = "nonexistent"))
})

test_that("missing actor column errors", {
  ev <- make_events()
  expect_error(prepare_data(ev, actor = "nonexistent", action = "code"))
})

test_that("missing time column errors", {
  ev <- make_events()
  expect_error(prepare_data(ev, actor = "student", action = "code",
                            time = "nonexistent"))
})

test_that("invalid time_threshold errors", {
  ev <- make_events()
  expect_error(prepare_data(ev, actor = "student", action = "code",
                            time_threshold = -1))
})

test_that("missing session column errors", {
  ev <- make_events()
  expect_error(prepare_data(ev, actor = "student", action = "code",
                            session = "nonexistent"))
})

test_that("missing order column errors", {
  ev <- make_events()
  expect_error(prepare_data(ev, actor = "student", action = "code",
                            order = "nonexistent"))
})

# ---- NA padding in wide format ----

test_that("shorter sessions are NA-padded in wide format", {
  ev <- data.frame(
    student = c("s1", "s1", "s1", "s2"),
    code = c("A", "B", "C", "X"),
    stringsAsFactors = FALSE
  )
  res <- prepare_data(ev, actor = "student", action = "code")
  expect_equal(res$statistics$max_sequence_length, 3L)
  s2_row <- which(res$meta_data$student == "s2")
  expect_true(is.na(res$sequence_data[s2_row, 2]))
  expect_true(is.na(res$sequence_data[s2_row, 3]))
})
