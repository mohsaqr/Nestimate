test_that("mark_phases splits each actor's events into consecutive phases", {
  log <- data.frame(who = c("a", "a", "a", "a", "a", "b", "b", "b"),
                    t = c(5, 1, 3, 2, 4, 2, 1, 3),
                    act = c("e", "a", "c", "b", "d", "y", "x", "z"),
                    stringsAsFactors = FALSE)
  out <- mark_phases(log, actor = "who", time = "t", n = 2,
                     labels = c("early", "late"))
  # rows and columns unchanged, one column added
  expect_identical(out[, names(log)], log)
  expect_identical(names(out), c(names(log), "phase"))
  # by hand: actor a has 5 events (times 1..5): ceiling(k * 2 / 5) gives
  # 1, 1, 2, 2, 2; actor b has 3 (times 1..3): 1, 2, 2
  expect_identical(out$phase[order(out$who, out$t)],
                   c("early", "early", "late", "late", "late",
                     "early", "late", "late"))
})

test_that("mark_phases gives each actor phases that differ by at most one event", {
  log <- data.frame(who = rep(c("a", "b", "c"), times = c(7, 10, 4)),
                    t = c(7:1, 1:10, c(2, 4, 1, 3)),
                    stringsAsFactors = FALSE)
  out <- mark_phases(log, actor = "who", time = "t", n = 3)
  sizes <- table(out$who, out$phase)
  expect_true(all(apply(sizes, 1, \(r) max(r) - min(r)) <= 1))
  # invariant: within an actor, phases never go back in time order
  ordered <- out[order(out$who, out$t), ]
  steps <- tapply(as.integer(ordered$phase), ordered$who, \(p) all(diff(p) >= 0))
  expect_true(all(steps))
})

test_that("mark_phases without time keeps the row order and labels by number", {
  log <- data.frame(act = letters[1:4], stringsAsFactors = FALSE)
  out <- mark_phases(log)
  expect_identical(out$phase, c("1", "1", "2", "2"))
})

test_that("mark_phases refuses bad input with classed errors", {
  log <- data.frame(who = c("a", "a"), t = c(1, 2), phase = c("x", "y"),
                    stringsAsFactors = FALSE)
  expect_error(mark_phases(log, actor = "who"), class = "nestimate_bad_input")
  expect_error(mark_phases(log, actor = "nope", name = "p"),
               class = "nestimate_bad_input")
  expect_error(mark_phases(log, n = 1, name = "p"), class = "nestimate_bad_input")
  expect_error(mark_phases(log, n = 2, labels = c("x", "x"), name = "p"),
               class = "nestimate_bad_input")
  bad_time <- data.frame(t = c("1", "2"), stringsAsFactors = FALSE)
  expect_error(mark_phases(bad_time, time = "t"), class = "nestimate_bad_input")
  na_actor <- data.frame(who = c("a", NA), stringsAsFactors = FALSE)
  expect_error(mark_phases(na_actor, actor = "who"), class = "nestimate_bad_input")
})

test_that("mark_phases reproduces the first-half / second-half split", {
  # the split the sequence-comparison vignette used to build by hand
  dat <- human_long[order(human_long$session_id, human_long$order_in_session), ]
  n_per <- ave(dat$order_in_session, dat$session_id, FUN = length)
  pos <- ave(dat$order_in_session, dat$session_id, FUN = seq_along)
  by_hand <- ifelse(pos <= n_per %/% 2, "early", "late")
  out <- mark_phases(dat, actor = "session_id", time = "order_in_session",
                     labels = c("early", "late"))
  expect_identical(out$phase, by_hand)
})
