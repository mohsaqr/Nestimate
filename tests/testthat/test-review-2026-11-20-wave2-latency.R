testthat::skip_on_cran()

# Review 2026-11-20, wave 2: method = "latency".
# C12: a grouped build resolved its caps on reply pairs the groups' networks
#      did not read (a state outside `alphabet`).
# C13: an event without a time was sorted last into a session of its own,
#      silently dropping its reply.

test_that("a grouped build calibrates the cap on the replies in the alphabet (C12)", {
  # B -> A takes 1 s; C -> B (100 s) involves C, which is outside the alphabet
  d <- data.frame(sender = c("A", "B", "C"), stamp = c(0, 1, 101), grp = "g")
  build <- function(...) {
    suppressMessages(
      build_network(d, "latency", action = "sender", time = "stamp",
                    latency_cap = 1, params = list(alphabet = c("A", "B")),
                    ...),
      classes = "nestimate_single_sequence")
  }
  ungrouped <- build()
  grouped <- build(group = "grp")
  # cap = the 100th percentile of the one in-alphabet reply = 1 s, so B -> A
  # sits at the cap: weight = floor = 0.5
  expect_equal(ungrouped$latency_cap_seconds, 1)
  expect_equal(ungrouped$weights["B", "A"], 0.5)
  expect_equal(grouped$g$latency_cap_seconds, 1)
  expect_equal(grouped$g$weights, ungrouped$weights, tolerance = 1e-12)
})

test_that("grouped and pooled latency networks agree under an alphabet, every scope (C12)", {
  .local_seed(3)
  d <- data.frame(team = rep(c("r", "b"), each = 60),
                  chat = rep(paste0("c", 1:6), each = 20),
                  sender = sample(c("A", "B", "C", "D"), 120, replace = TRUE),
                  at = cumsum(stats::rexp(120, 1 / 30)))
  agree <- vapply(c("global", "actor", "session", "responder"), function(s) {
    build <- function(...) {
      suppressMessages(
        build_network(d, "latency", actor = "chat", action = "sender",
                      time = "at", latency_scope = s, latency_aggregate = "sum",
                      alphabet = c("A", "B", "C"), ...),
        classes = "nestimate_latency_fallback")
    }
    full <- build()
    grouped <- build(group = "team")
    total <- Reduce(`+`, lapply(grouped, function(n) n$weights))
    isTRUE(all.equal(grouped$r$latency_cap_seconds,
                     full$latency_cap_seconds)) &&
      isTRUE(all.equal(total, full$weights, tolerance = 1e-12))
  }, logical(1))
  expect_true(all(agree))
})

test_that("an event without a time is refused, naming its row (C13)", {
  d <- data.frame(sender = c("A", "B", "C"), stamp = c(0, 1, NA_real_))
  expect_error(
    build_network(d, "latency", action = "sender", time = "stamp",
                  latency_cap = "10s"),
    class = "nestimate_bad_input")
  expect_error(
    build_network(d, "latency", action = "sender", time = "stamp",
                  latency_cap = "10s"),
    "row\\(s\\) 3")
  # the same in a grouped build
  expect_error(
    build_network(transform(d, grp = "g"), "latency", action = "sender",
                  time = "stamp", group = "grp", latency_cap = "10s"),
    class = "nestimate_bad_input")
  # an empty time string is a missing time
  d_text <- data.frame(sender = c("A", "B", "C", "A"),
                       stamp = c("2024-01-01 10:00:00", "",
                                 "2024-01-01 10:00:05",
                                 "2024-01-01 10:00:09"))
  expect_error(
    build_network(d_text, "latency", action = "sender", time = "stamp",
                  latency_cap = "10s"),
    class = "nestimate_bad_input")
  # other methods still take a log with a missing time
  expect_no_error(suppressMessages(
    build_network(d, "relative", action = "sender", time = "stamp"),
    classes = "nestimate_single_sequence"))
})

test_that("a log with every time present is read as before (C13)", {
  # A 0, B 10, A 70, B 100 | gap | C 2100, A 2400, A 2410
  d <- data.frame(chat = "c1",
                  sender = c("A", "B", "A", "B", "C", "A", "A"),
                  at = c(0, 10, 70, 100, 2100, 2400, 2410))
  net <- build_network(d, "latency", actor = "chat", action = "sender",
                       time = "at", latency_cap = "60s", latency_aggregate = "sum")
  w <- function(t) 1 - 0.5 * min(t / 60, 1)
  expected <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"),
                                              c("A", "B", "C")))
  expected["B", "A"] <- w(10) + w(30)
  expected["A", "B"] <- w(60)
  expected["A", "C"] <- w(300)
  expected["A", "A"] <- w(10)
  expect_equal(net$weights, expected, tolerance = 1e-12)
  # an event with neither an action nor a time forms no reply: not refused
  d_na <- data.frame(sender = c("A", "B", NA, "A"), stamp = c(0, 1, NA, 3))
  na_net <- suppressMessages(
    build_network(d_na, "latency", action = "sender", time = "stamp",
                  latency_cap = "10s"),
    classes = "nestimate_single_sequence")
  expect_equal(na_net$weights["B", "A"], 1 - 0.5 * 1 / 10)
  expect_equal(na_net$weights["A", "B"], 1 - 0.5 * 2 / 10)
})
