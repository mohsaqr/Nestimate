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

# ---- Undirected and signed networks ----

.metric <- function(W, directed, metric) {
  s <- .summary_metrics_from_weights(W, directed = directed)
  s$value[s$metric == metric]
}

test_that("an undirected edge counts once in Edge Count and Density", {
  star <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  star[1, 2:4] <- star[2:4, 1] <- c(1, 2, 3)
  expect_equal(.metric(star, FALSE, "Edge Count"), 3)
  expect_equal(.metric(star, FALSE, "Network Density"), 0.5)
  # The same matrix read as directed has six arcs out of twelve.
  expect_equal(.metric(star, TRUE, "Edge Count"), 6)
  expect_equal(.metric(star, TRUE, "Network Density"), 0.5)
})

test_that("negative edges count, and their magnitude is their distance", {
  # Path A - B - C with A-B negative: distances 1, 2 and 3, in both directions.
  W <- matrix(c(0, -1, 0,
                -1, 0, 2,
                0, 2, 0), 3, byrow = TRUE,
              dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  expect_equal(.metric(W, FALSE, "Edge Count"), 2)
  expect_equal(.metric(W, FALSE, "Network Density"), 2 / 3)
  expect_equal(.metric(W, FALSE, "Mean Distance"), 2)
  expect_equal(.metric(W, TRUE, "Reciprocity"), 1)
})

test_that("summary() Edge Count agrees with $n_edges on a glasso network", {
  skip_if_not_installed("glmnet")
  net <- build_network(srl_strategies, method = "glasso")
  s <- summary(net)
  expect_equal(s$value[s$metric == "Edge Count"], net$n_edges)
})

test_that("flipping every sign leaves the counts and distances unchanged", {
  W <- matrix(c(0, 2, 1,
                2, 0, 4,
                1, 4, 0), 3, byrow = TRUE,
              dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  metrics <- c("Edge Count", "Network Density", "Mean Distance")
  expect_equal(vapply(metrics, \(m) .metric(-W, FALSE, m), 0),
               vapply(metrics, \(m) .metric(W, FALSE, m), 0))
})
