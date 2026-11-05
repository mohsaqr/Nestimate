# Windows on sequence data come from params = list(window_size =).

.two_rows <- function() {
  data.frame(T1 = c("A", "A"), T2 = "B", T3 = "C", T4 = "D")
}

test_that("window_size counts between consecutive blocks of positions", {
  # Blocks {T1, T2} and {T3, T4}: each row gives A-C, A-D, B-C, B-D once.
  net <- build_network(.two_rows(), method = "frequency",
                       params = list(window_size = 2L))
  w <- net$weights[c("A", "B", "C", "D"), c("A", "B", "C", "D")]
  expect_equal(unname(w), rbind(c(0, 0, 2, 2), c(0, 0, 2, 2),
                                c(0, 0, 0, 0), c(0, 0, 0, 0)))
  rel <- build_network(.two_rows(), method = "relative",
                       params = list(window_size = 2L))
  expect_equal(rel$weights[["A", "C"]], 0.5)
  expect_equal(rel$weights[["B", "D"]], 0.5)
})

test_that("window_size counts co-occurrence inside each block", {
  # Blocks {T1, T2} and {T3, T4}: A-B and C-D only, never across blocks.
  w <- build_network(.two_rows(), method = "co_occurrence",
                     params = list(window_size = 2L))$weights
  expect_gt(w[["A", "B"]], 0)
  expect_gt(w[["C", "D"]], 0)
  expect_equal(w[["A", "C"]], 0)
  expect_equal(w[["B", "D"]], 0)
})

test_that("window_size = 1 and no window give the same transitions", {
  a <- build_network(.two_rows(), method = "frequency",
                     params = list(window_size = 1L))
  b <- build_network(.two_rows(), method = "frequency")
  expect_equal(unname(a$weights), unname(b$weights))
})

test_that("windowed = TRUE means a window of 2", {
  a <- build_network(.two_rows(), method = "frequency",
                     params = list(windowed = TRUE))
  b <- build_network(.two_rows(), method = "frequency",
                     params = list(window_size = 2L))
  expect_equal(a$weights, b$weights)
})

test_that("long, wide and grouped data window alike", {
  wide <- data.frame(T1 = c("A", "B", "A", "B"), T2 = c("B", "C", "C", "A"),
                     T3 = c("C", "A", "B", "C"), T4 = c("A", "B", "A", "A"))
  params <- list(window_size = 2L)
  ref <- build_network(wide, method = "frequency", params = params)$weights

  long <- data.frame(
    id = rep(seq_len(4L), each = 4L),
    step = rep(seq_len(4L), 4L),
    state = as.vector(t(as.matrix(wide)))
  )
  from_long <- build_network(long, method = "frequency", action = "state",
                             actor = "id", order = "step", params = params)
  expect_equal(from_long$weights, ref)

  grouped <- build_network(cbind(g = c("x", "x", "y", "y"), wide),
                           method = "frequency", group = "g", params = params)
  x <- build_network(wide[1:2, ], method = "frequency", params = params)
  expect_equal(grouped$x$weights, x$weights)
  expect_equal(grouped$x$weights + grouped$y$weights[rownames(grouped$x$weights),
                                                     colnames(grouped$x$weights)],
               ref)
})

test_that("an actor or covariate column does not shift the blocks", {
  params <- list(window_size = 2L)
  ref <- build_network(.two_rows(), method = "frequency", params = params)$weights
  d <- cbind(person = c("p1", "p2"), .two_rows(), score = c(1, 2))
  net <- build_network(d, method = "frequency", actor = "person",
                       metadata_cols = "score", params = params)
  expect_equal(net$weights, ref)
})

test_that("the window is recorded on the data, and resampling respects it", {
  set.seed(11)
  d <- as.data.frame(matrix(sample(c("A", "B", "C"), 240L, TRUE), 30L, 8L,
                            dimnames = list(NULL, paste0("T", seq_len(8L)))))
  net <- build_network(d, method = "frequency", params = list(window_size = 2L))
  plain <- build_network(d, method = "frequency")
  expect_true(isTRUE(attr(net$data, "windowed")))
  expect_equal(attr(net$data, "window_size"), 2L)
  # 3 block pairs of 2 x 2 per row, against 7 adjacent pairs per row
  expect_equal(sum(net$weights), 30 * 3 * 4)
  expect_equal(sum(plain$weights), 30 * 7)

  boot <- bootstrap_network(net, iter = 200L, seed = 1)
  expect_lt(abs(sum(boot$mean) - sum(net$weights)) / sum(net$weights), 0.05)
})

test_that("windows the data records are not replaced", {
  d <- .two_rows()
  attr(d, "windowed") <- TRUE
  attr(d, "window_size") <- 2L
  attr(d, "window_span") <- 1L
  a <- build_network(d, method = "frequency")
  b <- build_network(d, method = "frequency", params = list(window_size = 3L))
  expect_equal(a$weights, b$weights)
})

test_that("methods without a windowed count raise an error", {
  invisible(lapply(c("attention", "ngram", "gap", "reverse"), function(m) {
    expect_error(
      build_network(.two_rows(), method = m, params = list(window_size = 2L)),
      class = "nestimate_bad_input"
    )
  }))
})

test_that("a window that is not a positive whole number raises an error", {
  invisible(lapply(list(0, -1, 1.5, "2", c(2, 3)), function(size) {
    expect_error(
      build_network(.two_rows(), method = "frequency",
                    params = list(window_size = size)),
      class = "nestimate_bad_input"
    )
  }))
})

test_that("mode = 'overlapping' warns and counts blocks", {
  expect_warning(
    net <- build_network(.two_rows(), method = "frequency",
                         params = list(window_size = 2L, mode = "overlapping")),
    class = "nestimate_ignored_mode"
  )
  expect_null(net$params$mode)
  ref <- build_network(.two_rows(), method = "frequency",
                       params = list(window_size = 2L))
  expect_equal(net$weights, ref$weights)
})
