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
  expect_true(inherits(b, "netdifference"))
  expect_true(all(c("diff", "ci_lower", "ci_upper", "prob_x", "prob_y",
                    "sig", "summary", "difference_matrix") %in% names(b)))
  expect_equal(b$difference_matrix, b$diff)
  expect_true(is.data.frame(b$summary))
  expect_true(all(c("from", "to", "weight_x", "weight_y", "count_x", "count_y",
                    "diff", "effect_size", "ci_lower", "ci_upper", "ci_width",
                    "p_difference", "p_value", "sig")
                  %in% names(b$summary)))
  # CI brackets the posterior mean difference
  expect_true(all(b$summary$ci_lower <= b$summary$diff + 1e-8))
  expect_true(all(b$summary$ci_upper >= b$summary$diff - 1e-8))
  expect_true(all(b$summary$ci_width >= 0))
})

test_that("bayes_compare result is a net_permutation drop-in", {
  set.seed(2)
  d1 <- data.frame(V1 = sample(LETTERS[1:4], 60, TRUE),
                   V2 = sample(LETTERS[1:4], 60, TRUE),
                   V3 = sample(LETTERS[1:4], 60, TRUE))
  d2 <- data.frame(V1 = sample(LETTERS[1:4], 60, TRUE),
                   V2 = sample(LETTERS[1:4], 60, TRUE),
                   V3 = sample(LETTERS[1:4], 60, TRUE))
  b <- bayes_compare(build_network(d1, method = "relative"),
                     build_network(d2, method = "relative"),
                     draws = 1000, seed = 1)
  expect_true(inherits(b, "net_permutation"))
  expect_true(inherits(b, "netdifference"))
  # all net_permutation slots present
  expect_true(all(c("x", "y", "diff", "diff_sig", "p_values", "effect_size",
                    "summary", "method", "iter", "alpha", "paired", "adjust",
                    "difference_matrix")
                  %in% names(b)))
  # permutation summary columns present
  perm_cols <- c("from", "to", "weight_x", "weight_y", "diff",
                 "effect_size", "p_value", "sig")
  expect_true(all(perm_cols %in% names(summary(b))))
  # diff_sig is diff masked by significance
  expect_equal(b$diff_sig, b$diff * b$sig)
})

test_that("bayes_compare can be represented as a netdifference with CI metadata", {
  d1 <- data.frame(V1 = c("A","A","B","B","C","C"),
                   V2 = c("B","C","A","C","A","B"))
  d2 <- data.frame(V1 = c("A","A","B","C","C","C"),
                   V2 = c("B","B","A","A","A","B"))
  b <- bayes_compare(build_network(d1, method = "relative"),
                     build_network(d2, method = "relative"),
                     draws = 500, seed = 1)
  d <- as_netdifference(b)

  expect_s3_class(d, "netdifference")
  expect_s3_class(b, "netdifference")
  # the coercion is a PURE netdifference: net_permutation/net_bayes classes
  # dropped so cograph's splot routes it to the difference renderer, not the
  # permutation renderer
  expect_false(inherits(d, "net_permutation"))
  expect_false(inherits(d, "net_bayes"))
  expect_equal(unname(d$weights), unname(b$diff_sig))
  expect_equal(d$difference_matrix, b$diff)
  expect_equal(d$ci_lower, b$ci_lower)
  expect_equal(d$ci_upper, b$ci_upper)
  expect_equal(d$p_values, b$p_values)
  expect_equal(d$sig, b$sig)
  # source posterior means carried as plain matrices; print() works
  expect_equal(d$x, b$prob_x)
  expect_equal(d$y, b$prob_y)
  expect_equal(d$p_difference, b$p_difference)
  expect_silent(invisible(capture.output(print(d))))

  d_full <- as_netdifference(b, significant_only = FALSE)
  expect_equal(unname(d_full$weights), unname(b$diff))
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
  expected <- ci_excl & abs(s$diff) > b$mean_threshold &
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
  expect_true(all(c("comparison", "from", "to", "diff") %in%
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

test_that("bayes_compare works on edge-betweenness networks", {
  w1 <- data.frame(
    T1 = c("A", "A", "B", "C", "A", "B", "C", "A"),
    T2 = c("B", "B", "C", "A", "C", "A", "B", "B"),
    T3 = c("C", "C", "A", "B", "B", "C", "A", "C"),
    stringsAsFactors = FALSE
  )
  w2 <- data.frame(
    T1 = c("A", "C", "B", "C", "B", "A", "C", "B"),
    T2 = c("C", "A", "A", "B", "C", "B", "A", "A"),
    T3 = c("B", "B", "C", "A", "A", "C", "B", "C"),
    stringsAsFactors = FALSE
  )
  eb1 <- net_edge_betweenness(build_network(w1, method = "relative"))
  eb2 <- net_edge_betweenness(build_network(w2, method = "relative"))

  b <- bayes_compare(eb1, eb2, draws = 300, seed = 1)
  b_again <- bayes_compare(eb1, eb2, draws = 300, seed = 1)

  expect_true(inherits(b, "net_bayes"))
  expect_identical(b$method, "edge_betweenness")
  expect_identical(b$source_method, "relative")
  expect_equal(b$p_values, b_again$p_values)               # seed-reproducible
  expect_equal(b$observed_diff, eb1$weights - eb2$weights)
  expect_true(all(b$ci_lower <= b$ci_upper))
  expect_true(all(b$p_values >= 0 & b$p_values <= 1))
  expect_true(all(b$diff >= b$ci_lower - 1e-9 & b$diff <= b$ci_upper + 1e-9))
  expect_true(is.data.frame(summary(b)))
  expect_true(any(grepl("Edge-Betweenness", capture.output(print(b)))))
})

test_that("bayes_compare edge-betweenness self-comparison finds nothing credible", {
  seqs <- data.frame(
    T1 = c("A", "B", "C", "A", "B"), T2 = c("B", "C", "A", "B", "C"),
    T3 = c("C", "A", "B", "C", "A"), stringsAsFactors = FALSE
  )
  eb <- net_edge_betweenness(build_network(seqs, method = "relative"))
  b <- bayes_compare(eb, eb, draws = 300, seed = 2)
  expect_equal(sum(b$sig), 0)
  expect_true(all(abs(b$diff) < 1))    # posterior mean diff hovers near zero
})

test_that("bayes_compare edge-betweenness validates its inputs", {
  seqs <- data.frame(
    T1 = c("A", "B", "C", "A"), T2 = c("B", "C", "A", "B"),
    T3 = c("C", "A", "B", "C"), stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  eb_rel <- net_edge_betweenness(net)
  eb_freq <- net_edge_betweenness(build_network(seqs, method = "frequency"))
  eb_noinv <- net_edge_betweenness(net, invert = FALSE)

  expect_error(bayes_compare(eb_rel, net, draws = 50),
               "Both x and y")
  expect_error(bayes_compare(eb_rel, eb_freq, draws = 50),
               "relative")
  expect_error(bayes_compare(eb_rel, eb_noinv, draws = 50),
               "same edge-betweenness `invert` setting")
})
