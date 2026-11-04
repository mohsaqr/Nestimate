# summary() of a network: the 13 descriptive metrics on a tiny hand-checked
# weight matrix. Rows are sources, so out-strength is the row sum.

.sum_net <- function() {
  seqs <- data.frame(
    V1 = c("A", "B", "A", "C", "B", "A"), V2 = c("B", "C", "B", "A", "C", "C"),
    V3 = c("C", "A", "C", "B", "A", "B"), V4 = c("D", "D", "D", "D", "D", "D"),
    stringsAsFactors = FALSE
  )
  build_network(seqs, method = "frequency")
}

test_that("summary() reports the SD of out-strength from row sums and in-strength from column sums", {
  net <- .sum_net()
  s <- summary(net)
  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), 13L)
  value <- function(m) s$value[s$metric == m]
  # D is a sink: out-strengths 6, 6, 6, 0; in-strengths 3, 4, 5, 6
  expect_equal(value("Mean Out-Strength"), 4.5)
  expect_equal(value("Mean In-Strength"), 4.5)
  expect_equal(value("SD Out-Strength"), stats::sd(c(6, 6, 6, 0)))
  expect_equal(value("SD In-Strength"), stats::sd(c(3, 4, 5, 6)))
  expect_equal(value("SD Out-Strength"), stats::sd(rowSums(net$weights)))
  expect_equal(value("SD In-Strength"), stats::sd(colSums(net$weights)))
})

test_that("transposing the weights swaps the in- and out-strength SDs and keeps the means", {
  W <- matrix(c(0, 3, 1, 0,  2, 0, 4, 0,  1, 1, 0, 5,  0, 0, 0, 0), 4,
              byrow = TRUE, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  s  <- .summary_metrics_from_weights(W, directed = TRUE)
  st <- .summary_metrics_from_weights(t(W), directed = TRUE)
  value <- function(x, m) x$value[x$metric == m]
  expect_equal(value(s, "SD Out-Strength"), value(st, "SD In-Strength"))
  expect_equal(value(s, "SD In-Strength"),  value(st, "SD Out-Strength"))
  expect_equal(value(s, "Mean Out-Strength"), value(s, "Mean In-Strength"))
  expect_equal(value(s, "SD Out-Strength"), stats::sd(c(4, 6, 7, 0)))
})
