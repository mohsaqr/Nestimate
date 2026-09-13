# extract_pathways() -- cutting an event log into pathways.

make_log <- function() {
  data.frame(
    id  = c(rep("a", 10), rep("b", 4)),
    grp = c(rep("s1", 6), rep("s2", 4), rep("s3", 4)),
    act = c("Try", "Wrong", "Hint", "Retry", "Right", "Praise",
            "Try", "Wrong", "Retry", "Wrong",
            "Try", "Right", "Praise", "Quit"),
    stringsAsFactors = FALSE
  )
}
OUT <- c("Right", "Wrong")

test_that("unit gives one pathway per group and keeps every event", {
  p <- extract_pathways(make_log(), action = "act", group = c("id", "grp"))
  expect_s3_class(p, "data.frame")
  expect_equal(nrow(p), 3L)
  expect_equal(sum(p$length), nrow(make_log()))
  expect_equal(p$opens, c("Try", "Try", "Try"))
})

test_that("unit + terminal truncates at the group's LAST terminal state", {
  p <- extract_pathways(make_log(), action = "act", group = c("id", "grp"),
                        terminal = OUT)
  expect_true(all(p$closes %in% OUT))
  # group s1 runs Try Wrong Hint Retry Right Praise -> cut after Right
  expect_equal(p$path[p$grp == "s1"], "Try -> Wrong -> Hint -> Retry -> Right")
})

test_that("segments tile the group: every event in exactly one pathway", {
  p <- extract_pathways(make_log(), action = "act", group = c("id", "grp"),
                        type = "segments", terminal = OUT)
  expect_true(all(p$closes %in% OUT))
  # s1 = (Try Wrong)(Hint Retry Right); s2 = (Try Wrong)(Retry Wrong); s3 = (Try Right)
  expect_equal(nrow(p), 5L)
  expect_equal(p$path[1:2], c("Try -> Wrong", "Hint -> Retry -> Right"))
})

test_that("anchored opens one pathway per anchor and may overlap", {
  p <- extract_pathways(make_log(), action = "act", group = c("id", "grp"),
                        type = "anchored", anchor = "Wrong", terminal = OUT)
  expect_true(all(p$opens == "Wrong"))
  expect_equal(p$path[1], "Wrong -> Hint -> Retry -> Right")
  # the trailing Wrong of s2 has no following outcome, so it opens nothing
  expect_equal(nrow(p), 2L)
})

test_that("pathway ids restart within each group", {
  p <- extract_pathways(make_log(), action = "act", group = c("id", "grp"),
                        type = "segments", terminal = OUT)
  expect_equal(p$pathway[p$grp == "s1"], c(1L, 2L))
  expect_equal(p$pathway[p$grp == "s3"], 1L)
})

test_that("resolve appends a label, first rule winning, and keeps the raw end", {
  p <- extract_pathways(make_log(), action = "act", group = c("id", "grp"),
                        resolve = list(Exited = c("Praise", "Quit"),
                                       Solved = "Right"))
  expect_true("ends" %in% names(p))
  expect_equal(p$closes[p$grp == "s3"], "Exited")   # has Praise AND Quit
  expect_equal(p$ends[p$grp == "s3"], "Quit")       # raw final state kept
  expect_equal(p$closes[p$grp == "s1"], "Solved")   # Praise only, no Quit
  expect_true(grepl("Solved$", p$path[p$grp == "s1"]))
  expect_equal(p$closes[p$grp == "s2"], "Unresolved")
})

test_that("order sorts before cutting", {
  d <- make_log()[c(6:1, 10:7, 14:11), ]
  d$seq <- c(6:1, 10:7, 14:11)
  p <- extract_pathways(d, action = "act", group = c("id", "grp"),
                        order = "seq", terminal = OUT)
  expect_equal(p$path[p$grp == "s1"], "Try -> Wrong -> Hint -> Retry -> Right")
})

test_that("bad input is refused with a message naming the problem", {
  d <- make_log()
  expect_error(extract_pathways(d, action = "nope", group = "id"), "not found")
  expect_error(extract_pathways(d, action = "act", group = "id",
                                type = "anchored", anchor = "Wrong"),
               "needs both")
  expect_error(extract_pathways(d, action = "act", group = "id",
                                type = "segments"), "needs `terminal`")
  expect_error(
    extract_pathways(d, action = "act", group = "id", type = "anchored",
                     anchor = "Nonexistent", terminal = OUT),
    class = "nestimate_no_pathway")
})

test_that("a group whose rows are not contiguous keeps only its own events", {
  # Regression: the cuts are spans between a group's first and last row, so an
  # interleaved log (a learner leaving a unit and returning to it) made a
  # pathway swallow the events of whatever groups sat in between. Silent: the
  # path was simply longer and wrong.
  log <- data.frame(
    id  = c("a", "a", "b", "b", "b", "a", "a"),
    act = c("Try", "Wrong", "Try", "Hint", "Right", "Retry", "Right"),
    stringsAsFactors = FALSE
  )
  p <- extract_pathways(log, action = "act", group = "id")

  expect_equal(nrow(p), 2L)
  a <- p$path[p$id == "a"]
  b <- p$path[p$id == "b"]
  expect_equal(a, "Try -> Wrong -> Retry -> Right")
  expect_equal(b, "Try -> Hint -> Right")
  expect_equal(p$length[p$id == "a"], 4L)
  expect_equal(p$length[p$id == "b"], 3L)
  expect_false(grepl("Hint", a, fixed = TRUE))
})

test_that("interleaving does not change pathway order or per-group results", {
  # Same events, one contiguous and one interleaved: identical pathways, and
  # the groups still come back in order of first appearance.
  tidy <- data.frame(
    id  = c("x", "x", "x", "y", "y"),
    act = c("Try", "Wrong", "Right", "Try", "Right"),
    stringsAsFactors = FALSE
  )
  mixed <- data.frame(
    id  = c("x", "y", "x", "x", "y"),
    act = c("Try", "Try", "Wrong", "Right", "Right"),
    stringsAsFactors = FALSE
  )
  expect_equal(extract_pathways(tidy,  action = "act", group = "id"),
               extract_pathways(mixed, action = "act", group = "id"))
  expect_equal(extract_pathways(mixed, action = "act", group = "id")$id,
               c("x", "y"))
})

test_that("segments and anchored cuts also ignore interleaved groups", {
  log <- data.frame(
    id  = c("a", "a", "b", "a", "a", "b"),
    act = c("Try", "Wrong", "Try", "Retry", "Right", "Right"),
    stringsAsFactors = FALSE
  )
  seg <- extract_pathways(log, action = "act", group = "id",
                          type = "segments", terminal = c("Right", "Wrong"))
  expect_equal(seg$path[seg$id == "a"], c("Try -> Wrong", "Retry -> Right"))
  expect_equal(seg$path[seg$id == "b"], "Try -> Right")

  anc <- extract_pathways(log, action = "act", group = "id",
                          type = "anchored", anchor = "Wrong",
                          terminal = c("Right", "Wrong"))
  expect_equal(anc$path, "Wrong -> Retry -> Right")
  expect_equal(anc$id, "a")
})
