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
