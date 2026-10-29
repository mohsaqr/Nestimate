testthat::skip_on_cran()

# plot() for grouped passage times and grouped Markov stability

.mk_groups <- function() {
  d <- data.frame(T1 = c("A", "B", "C", "A", "B", "C", "A", "C"),
                  T2 = c("B", "C", "A", "C", "A", "B", "B", "A"),
                  T3 = c("C", "A", "B", "B", "C", "A", "C", "B"),
                  g = rep(c("x", "y"), each = 4), stringsAsFactors = FALSE)
  build_network(d, method = "relative", group = "g")
}

test_that("grouped passage times plot one panel per group", {
  p <- plot(passage_time(.mk_groups()))
  expect_s3_class(p, "ggplot")
  built <- expect_silent(ggplot2::ggplot_build(p))
  expect_identical(nrow(built$layout$layout), 2L)
  # one fill scale shared by both panels, so equal colours mean equal steps
  fills <- vapply(p$scales$scales, \(sc) "fill" %in% sc$aesthetics,
                  logical(1))
  expect_identical(sum(fills), 1L)
})

test_that("grouped Markov stability plots every state for every group", {
  p <- plot(markov_stability(.mk_groups()))
  expect_s3_class(p, "ggplot")
  built <- expect_silent(ggplot2::ggplot_build(p))
  points <- built$data[[2]]
  # 3 states x 2 groups x 6 metrics
  expect_identical(nrow(points), 36L)
  expect_identical(length(unique(points$shape)), 2L)
})

test_that("the single-network stability plot gives every state a fill", {
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time")
  built <- ggplot2::ggplot_build(plot(markov_stability(net)))
  expect_false(anyNA(built$data[[1]]$fill))
})

test_that("passage-time tiles keep row = from, column = to", {
  P <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE,
              dimnames = list(c("A", "B"), c("A", "B")))
  # A -> B takes 1 / 0.1 = 10 steps, B -> A 1 / 0.2 = 5
  tiles <- plot(passage_time(P), log_scale = FALSE)$data
  expect_identical(trimws(tiles$label[tiles$from == "A" & tiles$to == "B"]),
                   "10.0")
  expect_identical(trimws(tiles$label[tiles$from == "B" & tiles$to == "A"]),
                   "5.0")
  grouped <- plot(passage_time(.mk_groups()), log_scale = FALSE)$data
  mt <- passage_time(.mk_groups())[["x"]]$matrix
  x_tiles <- grouped[grouped$group == "x", ]
  expect_equal(as.numeric(x_tiles$label[x_tiles$from == "A" &
                                          x_tiles$to == "B"]),
               round(mt["A", "B"], 1))
})

test_that("ten groups get ten distinct colour-shape pairs", {
  d <- data.frame(T1 = rep(c("A", "B"), 10), T2 = "A", T3 = "B",
                  g = rep(letters[1:10], each = 2), stringsAsFactors = FALSE)
  p <- plot(markov_stability(build_network(d, method = "relative",
                                            group = "g")))
  expect_silent(ggplot2::ggplot_build(p))
  marks <- Nestimate:::.group_marks(10L)
  expect_identical(anyDuplicated(paste(marks$colour, marks$shape)), 0L)
  # within the first nine groups colour and shape both differ
  expect_identical(anyDuplicated(marks$shape[1:9]), 0L)
  expect_error(Nestimate:::.group_marks(82L), class = "nestimate_bad_input")
})
