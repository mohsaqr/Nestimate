# Grouping safety for prepare_onehot()
#
# prepare_onehot() used to key its actor/session, aggregation, and reshape
# groups with interaction(), which fails in two independent ways:
#
#   1. Separator collision. interaction() pastes with ".", so ("a.b", "c") and
#      ("a", "b.c") both become "a.b.c" and two unrelated sequences merge.
#   2. Integer overflow. The combination code is built over the *marginal*
#      level space, so two columns of ~46,341 distinct values exceed 32-bit
#      integers, the codes become NA, and every affected row lands in one group.
#
# Both are silent -- the wrong answer looks like a valid one. The grouping now
# goes through .observed_group_id(), which numbers observed combinations only.
#
# These tests also pin the row ORDER, because a same-seed bootstrap built on
# prepared data depends on it.

test_that("separator-colliding actor/session ids stay separate sequences", {
  data <- data.frame(
    A = c(1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L),
    actor   = c("a.b", "a.b", "a",   "a"),
    session = c("c",   "c",   "b.c", "b.c"),
    stringsAsFactors = FALSE
  )
  keys <- subset(data, select = c(actor, session))

  # the old key really did collapse these -- guards the premise of the test
  expect_equal(length(unique(as.character(interaction(keys, drop = TRUE)))), 1L)
  expect_equal(length(unique(.observed_group_id(keys))), 2L)

  res <- prepare_onehot(data, cols = c("A", "B"), actor = "actor",
                        session = "session", window_size = 1L)
  expect_equal(nrow(res), 2L)
})


test_that("distinct actor/session pairs each yield one sequence", {
  n <- 40L
  data <- data.frame(
    A = rep(c(1L, 0L), length.out = n),
    B = rep(c(0L, 1L), length.out = n),
    actor   = sprintf("a%02d", seq_len(n)),
    session = sprintf("s%02d", seq_len(n)),
    stringsAsFactors = FALSE
  )
  res <- prepare_onehot(data, cols = c("A", "B"), actor = "actor",
                        session = "session", window_size = 1L)
  expect_equal(nrow(res), n)
})


test_that("sequences come back in first-appearance order of actor/session", {
  # ground truth: b/y appears first, so it must be the first output row.
  data <- data.frame(
    A = c(1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L),
    actor   = c("b", "b", "a", "a"),
    session = c("y", "y", "x", "x"),
    stringsAsFactors = FALSE
  )
  res <- prepare_onehot(data, cols = c("A", "B"), actor = "actor",
                        session = "session", window_size = 1L)
  expect_equal(nrow(res), 2L)
  # row 1 is actor b (A then B), row 2 is actor a (A then B)
  expect_equal(unname(unlist(head(res, 1), use.names = FALSE))[1L], "A")
})


test_that("missing actor/session identifiers error instead of merging", {
  data <- data.frame(
    A = c(1L, 0L, 1L),
    B = c(0L, 1L, 0L),
    actor   = c("a", NA, "b"),
    session = c("s", "s", "s"),
    stringsAsFactors = FALSE
  )
  expect_error(
    prepare_onehot(data, cols = c("A", "B"), actor = "actor",
                   session = "session", window_size = 1L),
    "Missing values in actor/session"
  )
})


test_that(".split_first_seen partitions completely, in first-appearance order", {
  ids <- c(3L, 1L, 3L, 2L, 1L)
  idx <- .split_first_seen(ids)

  # every row assigned exactly once
  expect_equal(sort(unlist(idx, use.names = FALSE)), seq_along(ids))
  # groups ordered by where they first appear: 3, then 1, then 2
  expect_equal(vapply(idx, function(i) ids[i[1L]], integer(1)),
               c(3L, 1L, 2L), ignore_attr = TRUE)
  # matches the unique()-scan idiom it replaced
  expect_equal(lapply(unique(ids), function(k) which(ids == k)),
               unname(idx))

  expect_equal(length(.split_first_seen(integer(0))), 0L)
  expect_equal(.split_first_seen(1L), list(`1` = 1L))
})


test_that("prepare_onehot stays linear in row count", {
  skip_on_cran()
  make <- function(n) {
    data.frame(
      A = rep(c(1L, 0L), length.out = n),
      B = rep(c(0L, 1L), length.out = n),
      actor   = sprintf("a%05d", seq_len(n)),
      session = sprintf("s%05d", seq_len(n)),
      stringsAsFactors = FALSE
    )
  }
  run <- function(n) {
    system.time(
      prepare_onehot(make(n), cols = c("A", "B"), actor = "actor",
                     session = "session", window_size = 1L)
    )[["elapsed"]]
  }
  t_small <- run(1000L)
  t_large <- run(8000L)
  # 8x the rows: linear would be ~8x, the old quadratic scan was ~64x.
  # Generous bound so a loaded CI machine cannot flake it.
  expect_lt(t_large, max(8 * t_small, 0.5) * 4)
})
