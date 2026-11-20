# permutation(): the permutation test of two transition networks. The observed differences of edge weights and
# of centralities are deterministic and are compared with tna::permutation_test(); the multiplicity adjustment is
# compared with stats::p.adjust(). The permutation p-values themselves are a known difference with tna: tna
# permutes whole sequences, Nestimate permutes the per-sequence transition counts, so the two null distributions
# are drawn differently under one seed (FINDINGS.md).
# Re-expressed from local_testing_and_equivalence/test-equiv-permutation.R (p.adjust part),
# test-permutation-equiv-extracted.R, test-equiv-permutation-centralities.R and test-equiv-bootstrap.R (part 6).
# The hand-coded permutation loop of test-equiv-permutation.R replays the package's own random stream, so it is
# not independent evidence and stays a standing study.

# tna's edge statistics restricted to the edges present in either network (the rows summary() reports).
rs_tna_edge_stats <- function(test, model_x, model_y) {
  stats <- rs_split_edge_names(test$edges$stats)
  present <- rs_edge_table(weight_x = model_x$weights, weight_y = model_y$weights)
  present <- present[present$weight_x != 0 | present$weight_y != 0, ]
  merge(present, stats, by = c("from", "to"))
}

halves <- list(x = tna::group_regulation[seq_len(1000L), ], y = tna::group_regulation[1000L + seq_len(1000L), ])

test_that("observed edge differences equal tna::permutation_test(); p-values differ by the permutation unit", {
  skip_without_reference("tna")
  model_x <- tna::tna(halves$x)
  model_y <- tna::tna(halves$y)
  set.seed(1)
  reference <- rs_tna_edge_stats(tna::permutation_test(model_x, model_y, iter = 200L), model_x, model_y)
  found <- permutation(build_network(halves$x, method = "relative"), build_network(halves$y, method = "relative"),
                       iter = 200L, seed = 1)
  expected <- data.frame(from = reference$from, to = reference$to, weight_x = reference$weight_x,
                         weight_y = reference$weight_y, diff = reference$diff_true)
  expect_reference(summary(found), expected, reference = "tna",
                   covers = c("permutation()", "summary.net_permutation()"),
                   quantity = "group_regulation halves, observed differences", tolerance = 1e-12, key = c("from", "to"))
  expect_known_difference(summary(found), data.frame(from = reference$from, to = reference$to, p_value = reference$p_value),
                          reference = "tna", covers = "permutation()",
                          quantity = "group_regulation halves, permutation p-values", tolerance = 1e-10,
                          key = c("from", "to"),
                          reason = "tna permutes whole sequences between the groups; Nestimate permutes the per-sequence transition counts. The null distributions are drawn differently under one seed, so the p-values agree only up to Monte Carlo error.")
})

test_that("observed edge differences equal tna on three pairs of random sequence sets", {
  skip_without_reference("tna")
  pairs <- lapply(seq_len(3L), \(i) list(x = rs_random_sequences(20L, 4L, 15L, 200L + i),
                                         y = rs_random_sequences(20L, 4L, 15L, 300L + i)))
  invisible(lapply(seq_along(pairs), \(i) {
    model_x <- tna::tna(pairs[[i]]$x)
    model_y <- tna::tna(pairs[[i]]$y)
    reference <- rs_tna_edge_stats(tna::permutation_test(model_x, model_y, iter = 50L), model_x, model_y)
    found <- permutation(build_network(pairs[[i]]$x, method = "relative"),
                         build_network(pairs[[i]]$y, method = "relative"), iter = 50L, seed = 1)
    expect_reference(summary(found), data.frame(from = reference$from, to = reference$to, diff = reference$diff_true),
                     reference = "tna", covers = c("permutation()", "summary.net_permutation()"),
                     quantity = sprintf("random pair %d, observed differences", i), tolerance = 1e-12,
                     key = c("from", "to"))
  }))
})

test_that("observed centrality differences of permutation(measures =) equal tna::permutation_test(measures =)", {
  skip_without_reference("tna")
  measures <- c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut", "Closeness", "Betweenness",
                "BetweennessRSP", "Diffusion", "Clustering")
  x <- rs_random_sequences(60L, 5L, 6L, 401L)
  y <- rs_random_sequences(60L, 5L, 6L, 402L)
  reference <- tna::permutation_test(tna::tna(x), tna::tna(y), iter = 30L, measures = measures)
  found <- permutation(build_network(x, method = "relative"), build_network(y, method = "relative"),
                       iter = 30L, measures = measures, seed = 1)
  # The centrality table is reachable only as a field of the result (API gap).
  stats <- found$centralities$stats
  reference_stats <- data.frame(state = as.character(reference$centralities$stats$state),
                                centrality = as.character(reference$centralities$stats$centrality),
                                diff_true = reference$centralities$stats$diff_true)
  expect_reference(subset(stats, centrality != "Diffusion"), subset(reference_stats, centrality != "Diffusion"),
                   reference = "tna", covers = "permutation()", quantity = "centrality differences, eight measures",
                   tolerance = 1e-9, key = c("state", "centrality"))
  expect_known_difference(subset(stats, centrality == "Diffusion"), subset(reference_stats, centrality == "Diffusion"),
                          reference = "tna", covers = "permutation()", quantity = "centrality differences, Diffusion",
                          tolerance = 1e-9, key = c("state", "centrality"),
                          reason = "Nestimate scales Diffusion centrality to 0-1 before differencing; tna 1.3.1 returns it raw (tna >= 1.5 scales it as Nestimate does).")
})

test_that("adjusted p-values equal stats::p.adjust() of the unadjusted ones, seven methods", {
  x <- rs_random_sequences(8L, 4L, 15L, 7L)
  y <- rs_random_sequences(8L, 4L, 15L, 10L)
  net_x <- build_network(x, method = "relative")
  net_y <- build_network(y, method = "relative")
  raw <- permutation(net_x, net_y, iter = 50L, adjust = "none", seed = 124L)
  invisible(lapply(c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), \(method) {
    adjusted <- permutation(net_x, net_y, iter = 50L, adjust = method, seed = 124L)
    # The p-value matrices are fields of the result (API gap); the adjustment runs over every cell.
    expected <- matrix(stats::p.adjust(as.vector(raw$p_values), method = method), nrow(raw$p_values),
                       dimnames = dimnames(raw$p_values))
    expect_reference(adjusted$p_values, expected, reference = "stats", covers = "permutation()",
                     quantity = sprintf("p-values adjusted by %s", method), tolerance = 1e-8)
  }))
})

test_that("observed edge differences of a grouped network equal tna::permutation_test() of its two groups", {
  skip_without_reference("tna")
  sequences <- rs_random_sequences(40L, 4L, 12L, 501L)
  group <- rep(c("x", "y"), length.out = nrow(sequences))
  model_x <- tna::tna(sequences[group == "x", ])
  model_y <- tna::tna(sequences[group == "y", ])
  reference <- rs_tna_edge_stats(tna::permutation_test(model_x, model_y, iter = 20L), model_x, model_y)
  found <- permutation(build_network(cbind(sequences, grp = group), method = "relative", group = "grp"),
                       iter = 20L, seed = 1)
  expect_reference(summary(found),
                   data.frame(group = "x vs y", from = reference$from, to = reference$to, weight_x = reference$weight_x,
                              weight_y = reference$weight_y, diff = reference$diff_true),
                   reference = "tna", covers = c("permutation()", "summary.net_permutation_group()"),
                   quantity = "two groups, observed differences", tolerance = 1e-12, key = c("group", "from", "to"))
})
