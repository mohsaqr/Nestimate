test_that("bayes_compare returns a well-formed net_bayes object", {
  set.seed(1)
  d1 <- data.frame(V1 = sample(LETTERS[1:4], 60, TRUE),
                   V2 = sample(LETTERS[1:4], 60, TRUE),
                   V3 = sample(LETTERS[1:4], 60, TRUE))
  d2 <- data.frame(V1 = sample(LETTERS[1:4], 60, TRUE),
                   V2 = sample(LETTERS[1:4], 60, TRUE),
                   V3 = sample(LETTERS[1:4], 60, TRUE))
  n1 <- build_network(d1, method = "relative")
  n2 <- build_network(d2, method = "relative")
  b <- bayes_compare(n1, n2, draws = 2000, seed = 42)

  expect_true(inherits(b, "net_bayes"))
  expect_true(all(c("diff", "ci_lower", "ci_upper", "prob_x", "prob_y",
                    "sig", "summary") %in% names(b)))
  expect_true(is.data.frame(b$summary))
  expect_true(all(c("from", "to", "prob_x", "prob_y", "count_x", "count_y",
                    "mean_diff", "ci_lower", "ci_upper", "ci_width", "sig")
                  %in% names(b$summary)))
  # CI brackets the posterior mean difference
  expect_true(all(b$summary$ci_lower <= b$summary$mean_diff + 1e-8))
  expect_true(all(b$summary$ci_upper >= b$summary$mean_diff - 1e-8))
  expect_true(all(b$summary$ci_width >= 0))
})

test_that("posterior mean difference matches the closed-form Dirichlet mean", {
  d1 <- data.frame(V1 = c("A","A","B","B","C","C"),
                   V2 = c("B","C","A","C","A","B"))
  d2 <- data.frame(V1 = c("A","A","B","C","C","C"),
                   V2 = c("B","B","A","A","A","B"))
  n1 <- build_network(d1, method = "relative")
  n2 <- build_network(d2, method = "relative")
  b <- bayes_compare(n1, n2, prior = 0.5, draws = 200, seed = 1)

  nodes <- n1$nodes$label
  c1 <- n1$frequency_matrix[nodes, nodes]
  c2 <- n2$frequency_matrix[nodes, nodes]
  a1 <- c1 + 0.5; a2 <- c2 + 0.5
  mean1 <- a1 / rowSums(a1)
  mean2 <- a2 / rowSums(a2)
  expect_equal(b$diff, mean1 - mean2, tolerance = 1e-12)
  expect_equal(b$prob_x, mean1, tolerance = 1e-12)
})

test_that("seed makes credible intervals reproducible", {
  d1 <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","C"))
  d2 <- data.frame(V1 = c("A","C","B","A"), V2 = c("C","B","A","B"))
  n1 <- build_network(d1, method = "relative")
  n2 <- build_network(d2, method = "relative")
  b1 <- bayes_compare(n1, n2, draws = 1000, seed = 7)
  b2 <- bayes_compare(n1, n2, draws = 1000, seed = 7)
  expect_identical(b1$ci_lower, b2$ci_lower)
  expect_identical(b1$ci_upper, b2$ci_upper)
})

test_that("significance rule honours both magnitude and bound thresholds", {
  d1 <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","C"))
  d2 <- data.frame(V1 = c("A","C","B","A"), V2 = c("C","B","A","B"))
  n1 <- build_network(d1, method = "relative")
  n2 <- build_network(d2, method = "relative")
  b <- bayes_compare(n1, n2, draws = 4000, seed = 3)
  s <- b$summary
  ci_excl <- (s$ci_lower > 0) | (s$ci_upper < 0)
  nearest <- pmin(abs(s$ci_lower), abs(s$ci_upper))
  expected <- ci_excl & abs(s$mean_diff) > b$mean_threshold &
    nearest > b$bound_threshold
  expect_identical(s$sig, expected)
})

test_that("non-transition methods are rejected", {
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(120), ncol = 4))
  n1 <- build_network(d, method = "cor")
  n2 <- build_network(d + rnorm(120, 0, 0.1), method = "cor")
  expect_error(bayes_compare(n1, n2), "transition-count method")
})

test_that("mismatched methods and node sets error", {
  d1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
  n_rel <- build_network(d1, method = "relative")
  n_freq <- build_network(d1, method = "frequency")
  expect_error(bayes_compare(n_rel, n_freq), "Methods must match")

  d2 <- data.frame(V1 = c("A","B","D"), V2 = c("B","D","A"))
  n_other <- build_network(d2, method = "relative")
  expect_error(bayes_compare(n_rel, n_other), "Nodes must be the same")
})

test_that("grouped dispatch returns all pairwise comparisons", {
  s <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","A","C"),
    grp = c("X","X","Y","Y","Z","Z")
  )
  nets <- build_network(s, method = "relative", group = "grp")
  bg <- bayes_compare(nets, draws = 500, seed = 1)
  expect_true(inherits(bg, "net_bayes_group"))
  expect_length(bg, 3L)  # choose(3, 2)
  expect_true(all(c("comparison", "from", "to", "mean_diff") %in%
                    names(summary(bg))))
})

test_that("frequency-method networks are accepted", {
  d1 <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","C"))
  d2 <- data.frame(V1 = c("A","C","B","A"), V2 = c("C","B","A","B"))
  n1 <- build_network(d1, method = "frequency")
  n2 <- build_network(d2, method = "frequency")
  b <- bayes_compare(n1, n2, draws = 500, seed = 1)
  expect_true(inherits(b, "net_bayes"))
  # transition probabilities are row-normalised, hence in [0, 1]
  expect_true(all(b$prob_x >= 0 & b$prob_x <= 1))
})
