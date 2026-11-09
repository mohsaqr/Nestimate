# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

# ---- prepare(): timezone-aware timestamp parsing ----

.tz_strings <- c("2024-01-01 10:00:00", "2024-01-01T10:00:00Z",
                 "2024-01-01 10:00:00+0200", "2024-01-01T10:00:00+02:00")
# Reference epochs computed independently (tna::prepare_data, timezone UTC):
# naive-as-UTC 10:00 = 1704103200; Z = 1704103200; +0200 = 08:00 UTC =
# 1704096000.
.tz_epochs <- c(1704103200, 1704103200, 1704096000, 1704096000)

# Evaluate `code` with the process TZ temporarily set (base R only).
.with_tz <- function(tz, code) {
  old <- Sys.getenv("TZ", unset = NA)
  Sys.setenv(TZ = tz)
  on.exit(if (is.na(old)) Sys.unsetenv("TZ") else Sys.setenv(TZ = old),
          add = TRUE)
  force(code)
}

test_that(".parse_time is independent of the machine time zone", {
  for_tz <- function(tz) {
    .with_tz(tz, as.numeric(.parse_time(.tz_strings, timezone = "UTC")))
  }
  expect_equal(for_tz("UTC"), .tz_epochs)
  expect_equal(for_tz("America/New_York"), .tz_epochs)
  expect_equal(for_tz("Asia/Tokyo"), .tz_epochs)
})

test_that("naive timestamps are interpreted in `timezone`", {
  naive <- "2024-01-01 10:00:00"
  utc <- as.numeric(.parse_time(naive, timezone = "UTC"))
  helsinki <- as.numeric(.parse_time(naive, timezone = "Europe/Helsinki"))
  expect_equal(utc - helsinki, 2 * 3600)
  # Offset strings are unaffected by `timezone`.
  z <- "2024-01-01T10:00:00Z"
  expect_equal(as.numeric(.parse_time(z, timezone = "UTC")),
               as.numeric(.parse_time(z, timezone = "Europe/Helsinki")))
})

test_that("mixed layouts in one column all parse", {
  mixed <- c("2024-01-01 10:00:00", "01/02/2024 10:00", "2024-01-03T10:00:00Z",
             "20240104100000", "1704448800")
  # 01/02/2024 alone cannot say which number is the day.
  expect_warning(parsed <- .parse_time(mixed, timezone = "UTC"),
                 class = "nestimate_ambiguous_dates")
  expect_false(anyNA(parsed))
  expect_s3_class(parsed, "POSIXct")
})

test_that("unix timestamps honour unit and timezone", {
  expect_equal(as.numeric(.parse_time(1704103200, is_unix_time = TRUE)),
               1704103200)
  expect_equal(as.numeric(.parse_time(1704103200000, is_unix_time = TRUE,
                                      unix_time_unit = "milliseconds")),
               1704103200)
  expect_equal(attr(.parse_time(1, is_unix_time = TRUE,
                                timezone = "Asia/Tokyo"), "tzone"),
               "Asia/Tokyo")
})

test_that("invalid timezone and unparsable values fail loudly", {
  expect_error(.parse_time("2024-01-01", timezone = "Mars/Olympus"),
               "timezone")
  expect_error(prepare(data.frame(a = "x", t = "2024-01-01"), action = "a",
                       time = "t", timezone = "Mars/Olympus"), "timezone")
  expect_error(.parse_time("not a date", timezone = "UTC"),
               "Could not parse")
})

test_that("prepare() session splitting is stable across machine time zones", {
  d <- data.frame(
    user = "A",
    t = c("2024-03-10 01:30:00", "2024-03-10 01:40:00",
          "2024-03-10 03:05:00", "2024-03-10 03:10:00"),
    action = c("read", "write", "read", "plan"),
    stringsAsFactors = FALSE
  )
  n_sessions <- function(tz) {
    .with_tz(tz, {
      prepare(d, actor = "user", action = "action", time = "t",
              time_threshold = 30 * 60, timezone = "UTC")$statistics$total_sessions
    })
  }
  # 01:40 -> 03:05 is 85 min > 30 min: two sessions, wherever the machine is
  # (2024-03-10 is the US DST switch; local-time parsing would move the gap).
  expect_equal(n_sessions("UTC"), 2L)
  expect_equal(n_sessions("America/New_York"), 2L)
})

test_that("build_network() forwards timezone to prepare()", {
  d <- data.frame(
    user = rep("A", 4),
    t = c("2024-01-01 10:00:00", "2024-01-01 10:01:00",
          "2024-01-01 12:00:00", "2024-01-01 12:01:00"),
    action = c("read", "write", "read", "plan"),
    stringsAsFactors = FALSE
  )
  net <- build_network(d, method = "frequency", actor = "user",
                       action = "action", time = "t", timezone = "UTC")
  # Two sessions of two events: exactly two transitions.
  expect_equal(sum(net$weights), 2)
})

test_that("time_data is expressed in `timezone`, not the machine's zone", {
  # 10:00 UTC must show as 10:00 UTC in time_data whatever TZ is set to
  old_tz <- Sys.getenv("TZ", unset = NA)
  Sys.setenv(TZ = "America/New_York")
  on.exit(if (is.na(old_tz)) Sys.unsetenv("TZ") else Sys.setenv(TZ = old_tz),
          add = TRUE)
  df <- data.frame(a = c("x", "y"),
                   t = c("2024-01-01 10:00:00", "2024-01-01 10:01:00"))
  p <- prepare(df, action = "a", time = "t", timezone = "UTC")
  expect_identical(attr(p$time_data$time_T1, "tzone"), "UTC")
  expect_identical(format(p$time_data$time_T1), "2024-01-01 10:00:00")
  p2 <- prepare(df, action = "a", time = "t", timezone = "Asia/Tokyo")
  expect_identical(attr(p2$time_data$time_T1, "tzone"), "Asia/Tokyo")
  expect_identical(format(p2$time_data$time_T1), "2024-01-01 10:00:00")
})

test_that("is_unix_time = TRUE reads text as Unix seconds, never as a layout", {
  # "202401011030" matches the compact YYYYMMDDHHMM layout; declared Unix
  # time reads it as 202401011030 seconds instead.
  d <- data.frame(a = c("x", "y"), t = c("1700000000", "1700000060"))
  p <- prepare(d, action = "a", time = "t", is_unix_time = TRUE)
  expect_equal(as.numeric(p$time_data$time_T1), 1700000000)
  d2 <- data.frame(a = c("x", "y"), t = c("202401011030", "202401011031"))
  as_date <- prepare(d2, action = "a", time = "t")
  expect_equal(format(as_date$time_data$time_T1, tz = "UTC"),
               "2024-01-01 10:30:00")
  as_unix <- suppressWarnings(
    prepare(d2, action = "a", time = "t", is_unix_time = TRUE))
  expect_equal(as.numeric(as_unix$time_data$time_T1), 202401011030)
  expect_error(prepare(data.frame(a = "x", t = "2024-01-01"), action = "a",
                       time = "t", is_unix_time = TRUE),
               class = "nestimate_bad_input")
})
