# ---- Categorical covariates in wide sequence data ----
# Positions are named as a run (T1, T2, ...); any other categorical column
# of such a frame is metadata, not a sequence position.

.cov_wide <- function() {
  data.frame(T1 = c("A", "B", "A", "C", "B", "A"),
             T2 = c("B", "C", "C", "A", "A", "B"),
             T3 = c("C", "A", "B", "B", "C", NA),
             outcome = c("pass", "fail", "pass", "fail", "pass", "fail"),
             course = c("x", "x", "y", "y", "z", "z"),
             stringsAsFactors = FALSE)
}

test_that("character covariates are metadata, not states", {
  net <- build_network(.cov_wide(), method = "relative")
  expect_identical(sort(net$nodes$label), c("A", "B", "C"))
  expect_identical(sort(names(net$metadata)), c("course", "outcome"))
  expect_identical(names(net$data), c("T1", "T2", "T3"))
})

test_that("inferred roles give the same network as declared ones", {
  inferred <- build_network(.cov_wide(), method = "relative")
  declared <- build_network(.cov_wide(), method = "relative",
                            metadata_cols = c("outcome", "course"))
  expect_identical(inferred$weights, declared$weights)
})

test_that("covariates sharing values with each other are still covariates", {
  d <- .cov_wide()
  d$pre <- c("low", "high", "low", "high", "low", "high")
  d$post <- c("high", "high", "low", "high", "low", "low")
  net <- build_network(d, method = "relative")
  expect_identical(sort(net$nodes$label), c("A", "B", "C"))
  expect_true(all(c("pre", "post") %in% names(net$metadata)))
})

test_that("roles are left alone without a single run of position names", {
  # Two runs of equal size: ambiguous
  tied <- data.frame(A1 = c("a", "b"), A2 = c("b", "a"),
                     B1 = c("x", "y"), B2 = c("y", "x"))
  net <- build_network(tied, method = "relative")
  expect_identical(sort(net$nodes$label), c("a", "b", "x", "y"))
  # No numbered names at all: every character column is a position
  plain <- data.frame(first = c("a", "b"), second = c("b", "a"),
                      label = c("x", "y"))
  net <- build_network(plain, method = "relative")
  expect_identical(sort(net$nodes$label), c("a", "b", "x", "y"))
})

test_that("a position whose states occur nowhere else stays a position", {
  # constant first state
  starts <- data.frame(V1 = c("A", "A", "A"), V2 = c("B", "C", "B"),
                       V3 = c("C", "B", "C"))
  net <- build_network(starts, method = "relative")
  expect_identical(sort(net$nodes$label), c("A", "B", "C"))
  # terminal state found only in the last column
  ends <- data.frame(T1 = c("A", "A", "B", "B"), T2 = c("A", "B", "A", "B"),
                     T3 = c("B", NA, "B", "A"), T4 = c("C", NA, NA, "C"))
  net <- build_network(ends, method = "relative")
  expect_true("C" %in% net$nodes$label)
  expect_identical(names(net$data), c("T1", "T2", "T3", "T4"))
})

test_that("an all-missing position column stays a position", {
  d <- .cov_wide()
  d$T4 <- NA_character_
  net <- build_network(d, method = "relative")
  expect_true("T4" %in% names(net$data))
  expect_false("T4" %in% names(net$metadata))
})

test_that("explicit state_cols still override the inference", {
  net <- build_network(.cov_wide(), method = "relative",
                       state_cols = c("T1", "T2", "T3", "course"))
  expect_true(all(c("x", "y", "z") %in% net$nodes$label))
})

test_that("grouped wide data keeps covariates out of every child", {
  nets <- build_network(.cov_wide(), method = "relative", group = "course")
  labels <- unique(unlist(lapply(nets, \(n) n$nodes$label)))
  expect_identical(sort(labels), c("A", "B", "C"))
})

test_that("the outcome column survives into pattern_outcome()", {
  net <- build_network(.cov_wide(), method = "relative")
  expect_identical(extract_sequences(net, meta = TRUE)$outcome,
                   .cov_wide()$outcome)
})

test_that("one long sequence pivoted wide keeps every position", {
  # Each column holds one value; with thousands of columns any per-pair
  # comparison of columns would not finish
  one_row <- as.data.frame(matrix(rep(c("A", "A", "B", "C", "A", "D"), 500),
                                  nrow = 1), stringsAsFactors = FALSE)
  elapsed <- system.time(
    net <- build_network(one_row, method = "relative")
  )[["elapsed"]]
  expect_identical(ncol(net$data), 3000L)
  expect_identical(sort(net$nodes$label), c("A", "B", "C", "D"))
  expect_lt(elapsed, 10)
})

test_that("an all-missing logical position stays in the run of positions", {
  d <- data.frame(T1 = c("A", "B"), T2 = NA, outcome = c("pass", "fail"))
  expect_true(is.logical(d$T2))
  expect_identical(names(extract_sequences(d)), c("T1", "T2"))
  net <- build_network(d, method = "relative")
  expect_false(any(c("pass", "fail") %in% net$nodes$label))
  expect_true("outcome" %in% names(net$metadata))
})
