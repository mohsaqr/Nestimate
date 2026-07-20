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
  res <- prepare(ev, actor = "student", action = "code")
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
  res <- prepare(ev, actor = "student", action = "code")
  s1_row <- which(res$meta_data[[grep("student|actor", names(res$meta_data),
                                       value = TRUE)[1]]] == "s1")
  expect_equal(as.character(res$sequence_data[s1_row, 1:3]), c("A", "B", "C"))
})

# ---- With time ----

test_that("time column is parsed and time_data is returned", {
  ev <- make_events()
  res <- prepare(ev, actor = "student", action = "code",
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
  res <- prepare(ev, actor = "student", action = "code",
                      time = "timestamp", time_threshold = 900)
  expect_equal(res$statistics$total_sessions, 2L)
  expect_equal(nrow(res$sequence_data), 2L)
  # First session: A, B; Second session: C, D
  expect_equal(as.character(res$sequence_data[1, 1:2]), c("A", "B"))
  expect_equal(as.character(res$sequence_data[2, 1:2]), c("C", "D"))
})

test_that("large time_threshold keeps everything in one session", {
  ev <- make_events()
  res <- prepare(ev, actor = "student", action = "code",
                      time = "timestamp", time_threshold = 1e6)
  expect_equal(res$statistics$total_sessions, 2L)
})

# ---- Missing actor ----

test_that("missing actor treats all data as one actor", {
  ev <- data.frame(code = c("A", "B", "C"), stringsAsFactors = FALSE)
  res <- prepare(ev, action = "code")
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
  res <- prepare(ev, actor = c("student", "group"), action = "code")
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
  res <- prepare(ev, actor = "student", action = "code",
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
  res <- prepare(ev, actor = "student", action = "code",
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
  res <- prepare(ev, actor = "student", action = "code",
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
  res <- prepare(ev, actor = "student", action = "code", time = "ts")
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
  res <- prepare(ev, actor = "student", action = "code", time = "ts")
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
  res <- prepare(ev, actor = "student", action = "code", time = "ts",
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
  res <- prepare(ev, actor = "student", action = "code", time = "ts",
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
  res <- prepare(ev, actor = "student", action = "code")
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
  res <- prepare(ev, actor = "student", action = "code")
  expect_equal(res$meta_data$level, "high")
})

# ---- print method ----

test_that("print.nestimate_data produces expected output", {
  ev <- make_events()
  res <- prepare(ev, actor = "student", action = "code",
                      time = "timestamp")
  out <- capture.output(print(res))
  expect_true(any(grepl("Prepared Data", out)))
  expect_true(any(grepl("Sessions:", out)))
  expect_true(any(grepl("Actors:", out)))
  expect_true(any(grepl("Time data: available", out)))
})

test_that("print without actors omits Actors line", {
  ev <- data.frame(code = c("A", "B"), stringsAsFactors = FALSE)
  res <- prepare(ev, action = "code")
  out <- capture.output(print(res))
  expect_false(any(grepl("Actors:", out)))
})

# ---- Error conditions ----

test_that("non-data.frame input errors", {
  expect_error(prepare("not_a_df", actor = "a", action = "b"))
})

test_that("missing action column errors", {
  ev <- make_events()
  expect_error(prepare(ev, actor = "student", action = "nonexistent"))
})

test_that("missing actor column errors", {
  ev <- make_events()
  expect_error(prepare(ev, actor = "nonexistent", action = "code"))
})

test_that("missing time column errors", {
  ev <- make_events()
  expect_error(prepare(ev, actor = "student", action = "code",
                            time = "nonexistent"))
})

test_that("invalid time_threshold errors", {
  ev <- make_events()
  expect_error(prepare(ev, actor = "student", action = "code",
                            time_threshold = -1))
})

test_that("missing session column errors", {
  ev <- make_events()
  expect_error(prepare(ev, actor = "student", action = "code",
                            session = "nonexistent"))
})

test_that("missing order column errors", {
  ev <- make_events()
  expect_error(prepare(ev, actor = "student", action = "code",
                            order = "nonexistent"))
})

# ---- NA padding in wide format ----

test_that("shorter sessions are NA-padded in wide format", {
  ev <- data.frame(
    student = c("s1", "s1", "s1", "s2"),
    code = c("A", "B", "C", "X"),
    stringsAsFactors = FALSE
  )
  res <- prepare(ev, actor = "student", action = "code")
  expect_equal(res$statistics$max_sequence_length, 3L)
  s2_row <- which(res$meta_data$student == "s2")
  expect_true(is.na(res$sequence_data[s2_row, 2]))
  expect_true(is.na(res$sequence_data[s2_row, 3]))
})


# ---- Unix timestamp numeric path (L333-338) ----

test_that(".parse_time numeric unix path with milliseconds divisor", {
  # Call .parse_time directly via prepare with is_unix_time + milliseconds
  ev <- data.frame(
    student = rep("s1", 2),
    code = c("A", "B"),
    ts = c(1704100000000, 1704100300000),
    stringsAsFactors = FALSE
  )
  res <- prepare(ev, actor = "student", action = "code", time = "ts",
                      is_unix_time = TRUE, unix_time_unit = "milliseconds")
  expect_false(is.null(res$time_data))
  # Should be about 2024 timestamps
  expect_s3_class(res$time_data[[1]], "POSIXct")
})

test_that(".parse_time numeric unix path with microseconds divisor", {
  ev <- data.frame(
    student = rep("s1", 2),
    code = c("A", "B"),
    ts = c(1704100000000000, 1704100300000000),
    stringsAsFactors = FALSE
  )
  res <- prepare(ev, actor = "student", action = "code", time = "ts",
                      is_unix_time = TRUE, unix_time_unit = "microseconds")
  expect_false(is.null(res$time_data))
  expect_s3_class(res$time_data[[1]], "POSIXct")
})

test_that(".parse_time errors on unparseable string timestamps (L341-343)", {
  # Provide completely unparseable strings that cannot be coerced to numeric
  ev <- data.frame(
    student = rep("s1", 2),
    code = c("A", "B"),
    ts = c("not-a-date-xyz", "also-garbage"),
    stringsAsFactors = FALSE
  )
  expect_error(
    prepare(ev, actor = "student", action = "code", time = "ts"),
    "Could not parse"
  )
})

# ---- .aggregate_metadata: all-NA extra column (L362) ----

test_that("extra column all-NA returns NA in meta_data", {
  ev <- data.frame(
    student = c("s1", "s1"),
    code = c("A", "B"),
    level = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  res <- prepare(ev, actor = "student", action = "code")
  # level is extra column; all NA -> aggregated to NA
  expect_true("level" %in% names(res$meta_data))
  expect_true(is.na(res$meta_data$level))
})

# ---- .aggregate_metadata: tied mode emits message (L366-373) ----

test_that("tied-mode character extra column emits message and returns first value", {
  ev <- data.frame(
    student = c("s1", "s1", "s1", "s1"),
    code = c("A", "B", "A", "B"),
    level = c("high", "low", "medium", "low"),
    stringsAsFactors = FALSE
  )
  # "low" appears twice, "high" and "medium" once each -> tie between all equal
  # Actually: high=1, low=2, medium=1 -> mode is "low" (unique), no tie
  # To get a tie: use equal counts
  ev2 <- data.frame(
    student = c("s1", "s1"),
    code = c("A", "B"),
    level = c("high", "low"),
    stringsAsFactors = FALSE
  )
  expect_message(
    res <- prepare(ev2, actor = "student", action = "code"),
    "ties resolved by first occurrence"
  )
  expect_true(res$meta_data$level %in% c("high", "low"))
})

# ---------------------------------------------------------------------------
# Wide time matrix: built in one shot, not column-by-column
# ---------------------------------------------------------------------------
# prepare() used to convert time_data one column at a time with
#   for (j in seq_len(ncol(time_data))) time_data[[j]] <- as.POSIXct(...)
# Each `[[<-` copies the whole data.frame, so the cost was quadratic in the
# sequence length. A single long sequence (no `actor`) has one row and one
# column per event, so 27k events meant 27k frame copies (~26 s).

test_that("time_data columns are POSIXct regardless of sequence shape", {
  ev <- data.frame(
    student = c("a", "a", "a", "b", "b"),
    code    = c("x", "y", "z", "x", "z"),
    t       = as.POSIXct("2024-01-01 00:00:00", tz = "UTC") + c(0, 60, 120, 0, 60)
  )
  res <- prepare(ev, actor = "student", action = "code", time = "t")
  expect_true(all(vapply(res$time_data, function(x) inherits(x, "POSIXct"), TRUE)))
  expect_true(all(grepl("^time_T\\d+$", names(res$time_data))))
  expect_identical(nrow(res$time_data), 2L)

  # Single pooled sequence: one row, one column per event.
  res1 <- prepare(ev, action = "code", time = "t")
  expect_identical(nrow(res1$time_data), 1L)
  expect_identical(ncol(res1$time_data), 5L)
  expect_true(all(vapply(res1$time_data, function(x) inherits(x, "POSIXct"), TRUE)))
})

test_that("a long single sequence does not scale quadratically", {
  skip_on_cran()
  n <- 4000L
  ev <- data.frame(
    code = rep(c("x", "y", "z", "w"), length.out = n),
    t    = as.POSIXct("2024-01-01", tz = "UTC") + seq_len(n)
  )
  half <- system.time(prepare(ev[seq_len(n / 2L), ], action = "code", time = "t"))[["elapsed"]]
  full <- system.time(prepare(ev, action = "code", time = "t"))[["elapsed"]]
  # Linear would be ~2x. Quadratic would be ~4x. Allow generous slack for
  # timer noise on loaded CI machines, but a quadratic regression blows past it.
  expect_lt(full, max(0.5, half * 3))
})

test_that("time_threshold = FALSE switches session-interval splitting off", {
  # One actor, three events an hour apart.
  ev <- data.frame(
    s = "u1",
    a = c("A", "B", "C"),
    t = as.POSIXct("2024-01-01", tz = "UTC") + c(0, 3600, 7200)
  )

  # The default 900s interval splits every gap into its own session.
  expect_identical(
    nrow(prepare(ev, actor = "s", action = "a", time = "t")$sequence_data),
    3L
  )

  # FALSE keeps the actor as a single sequence, matching an infinite interval.
  off <- prepare(ev, actor = "s", action = "a", time = "t",
                 time_threshold = FALSE)$sequence_data
  expect_identical(nrow(off), 1L)
  expect_identical(unname(unlist(off[1, ])), c("A", "B", "C"))
  expect_equal(
    off,
    prepare(ev, actor = "s", action = "a", time = "t",
            time_threshold = Inf)$sequence_data
  )

  # build_network() forwards it: three lone events carry no transition,
  # one pooled sequence carries A -> B -> C.
  expect_identical(
    build_network(ev, actor = "s", action = "a", time = "t",
                  method = "tna")$n_edges,
    0L
  )
  expect_identical(
    build_network(ev, actor = "s", action = "a", time = "t",
                  time_threshold = FALSE, method = "tna")$n_edges,
    2L
  )

  # TRUE stays invalid: only FALSE is meaningful as a switch.
  expect_error(prepare(ev, actor = "s", action = "a", time = "t",
                       time_threshold = TRUE))
  expect_error(prepare(ev, actor = "s", action = "a", time = "t",
                       time_threshold = 0))
  expect_error(prepare(ev, actor = "s", action = "a", time = "t",
                       time_threshold = NA))
})

test_that("high-cardinality actor/session grouping does not overflow", {
  skip_on_cran()
  # 50,000 x 50,000 marginal levels: interaction() codes this as 2.5e9, past
  # .Machine$integer.max, and used to drop ~14% of sessions into a single "NA"
  # pseudo-session that manufactured transitions no input sequence contained.
  n <- 50000L
  ev <- data.frame(
    student = sprintf("stu-%05d", seq_len(n)),
    step    = sprintf("step-%05d", seq_len(n)),
    stringsAsFactors = FALSE
  )
  ev <- ev[rep(seq_len(n), each = 3L), ]
  ev$tm  <- 1700000000 + rep(c(0, 1, 2), n)
  ev$act <- rep(c("A", "B", "C"), n)

  expect_no_warning(
    prepared <- prepare(ev, actor = "student", session = "step", time = "tm",
                        action = "act", time_threshold = 3600)
  )
  expect_identical(nrow(prepared$sequence_data), n)

  # Every session is A -> B -> C, so C is terminal and its row must be empty.
  w <- build_network(ev, actor = "student", session = "step", time = "tm",
                     action = "act", time_threshold = 3600,
                     method = "tna")$weights
  expect_equal(unname(w["A", "B"]), 1)
  expect_equal(unname(w["B", "C"]), 1)
  expect_equal(unname(w["C", ]), c(0, 0, 0))
})

test_that("separator characters inside identifiers cannot collide", {
  # ("a | b", "c") and ("a", "b | c") both paste to "a | b | c".
  ev <- data.frame(
    student = c("a | b", "a | b", "a", "a"),
    session = c("c", "c", "b | c", "b | c"),
    tm      = 1700000000 + c(0, 1, 0, 1),
    act     = c("A", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  res <- prepare(ev, actor = "student", session = "session", time = "tm",
                 action = "act", time_threshold = 3600)

  expect_identical(nrow(res$sequence_data), 2L)
  # The merge key must stay distinct even though the display label does not.
  expect_identical(anyDuplicated(res$meta_data$.session_id), 0L)
  expect_identical(nrow(res$meta_data), 2L)
})

test_that("missing actor or session identifiers fail fast", {
  ev <- data.frame(
    student = c("u1", "u1", NA, NA),
    session = c("s1", "s1", "s2", "s2"),
    tm      = 1700000000 + c(0, 1, 0, 1),
    act     = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  # Previously this silently produced three sessions from two real ones.
  expect_error(
    prepare(ev, actor = "student", session = "session", time = "tm",
            action = "act", time_threshold = 3600),
    "Missing values in actor/session column"
  )
  ev2 <- ev
  ev2$student <- "u1"
  ev2$session[3:4] <- NA
  expect_error(
    prepare(ev2, actor = "student", session = "session", time = "tm",
            action = "act", time_threshold = 3600),
    "Missing values in actor/session column"
  )
})

test_that("grouping keeps interaction()'s row order so seeds stay reproducible", {
  set.seed(9)
  k <- 60L
  ev <- data.frame(
    student = sprintf("u%02d", rep(seq_len(k), each = 6)),
    session = sprintf("s%d", rep(1:2, each = 3, times = k)),
    act     = sample(c("A", "B", "C"), k * 6, TRUE),
    stringsAsFactors = FALSE
  )
  ev$tm <- 1700000000 + seq_len(nrow(ev))

  sorted  <- prepare(ev, actor = "student", session = "session",
                     time = "tm", action = "act")
  shuffled <- prepare(ev[sample(nrow(ev)), ], actor = "student",
                      session = "session", time = "tm", action = "act")
  # Input order must not leak into prepared row order: a finite seeded
  # bootstrap indexes rows, so a permutation would change its output.
  expect_identical(sorted$sequence_data, shuffled$sequence_data)
})
