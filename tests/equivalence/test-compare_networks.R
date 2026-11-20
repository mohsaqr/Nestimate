# compare_networks(), edge_differences(), node_differences(), global_differences(): the descriptive comparison
# of two transition networks, against tna (the two models' weights, tna::centralities() of each, and the summary
# metrics of tna::compare()). The tables are asked for unrounded (digits = 15).
# Re-expressed from local_testing_and_equivalence/test-equiv-moved-compare_networks.R (weights of tna inputs).

halves <- list(early = tna::group_regulation[seq_len(200L), ], late = tna::group_regulation[1000L + seq_len(200L), ])

test_that("edge weights and differences equal the weights of tna::tna() of each network", {
  skip_without_reference("tna")
  comparison <- compare_networks(early = build_network(halves$early, method = "relative"),
                                 late = build_network(halves$late, method = "relative"))
  weights <- rs_edge_table(weight_a = tna::tna(halves$early)$weights, weight_b = tna::tna(halves$late)$weights)
  expected <- data.frame(pair = "early vs late", from = weights$from, to = weights$to, weight_a = weights$weight_a,
                         weight_b = weights$weight_b, diff = weights$weight_a - weights$weight_b)
  expect_reference(edge_differences(comparison, digits = 15L), expected, reference = "tna",
                   covers = c("compare_networks()", "edge_differences()"), quantity = "edges, early vs late",
                   tolerance = 1e-12, key = c("pair", "from", "to"))
})

test_that("tna models as input give tna's weights", {
  skip_without_reference("tna")
  model_early <- tna::tna(halves$early)
  model_late <- tna::tna(halves$late)
  comparison <- compare_networks(early = model_early, late = model_late, measures = NULL)
  weights <- rs_edge_table(weight_a = model_early$weights, weight_b = model_late$weights)
  expect_reference(edge_differences(comparison, digits = 15L),
                   data.frame(pair = "early vs late", weights), reference = "tna",
                   covers = c("compare_networks()", "edge_differences()"), quantity = "edges of tna models",
                   tolerance = 1e-12, key = c("pair", "from", "to"))
})

test_that("node values equal tna::centralities() of each network", {
  skip_without_reference("tna")
  measures <- c("InStrength", "OutStrength", "Betweenness")
  comparison <- compare_networks(early = build_network(halves$early, method = "relative"),
                                 late = build_network(halves$late, method = "relative"), measures = measures)
  centrality <- lapply(halves, \(h) as.data.frame(tna::centralities(tna::tna(h), measures = measures)))
  expected <- do.call(rbind, lapply(measures, \(m) data.frame(
    pair = "early vs late", node = as.character(centrality$early$state), measure = m,
    value_a = centrality$early[[m]], value_b = centrality$late[[m]],
    diff = centrality$early[[m]] - centrality$late[[m]])))
  expect_reference(node_differences(comparison, digits = 15L), expected, reference = "tna",
                   covers = c("compare_networks()", "node_differences()"), quantity = "node centralities, early vs late",
                   tolerance = 1e-10, key = c("pair", "node", "measure"))
})

test_that("global metrics equal the summary metrics of tna::compare()", {
  skip_without_reference("tna")
  comparison <- compare_networks(early = build_network(halves$early, method = "relative"),
                                 late = build_network(halves$late, method = "relative"))
  metrics <- as.data.frame(tna::compare(tna::tna(halves$early), tna::tna(halves$late))$summary_metrics)
  conventions <- c("Distance", "Canberra")
  global <- global_differences(comparison, digits = 15L)
  expected <- data.frame(pair = "early vs late", metrics)
  expect_reference(subset(global, !metric %in% conventions), subset(expected, !metric %in% conventions),
                   reference = "tna", covers = c("compare_networks()", "global_differences()"),
                   quantity = "global metrics, early vs late", tolerance = 1e-10, key = c("pair", "category", "metric"))
  # The two conventions of test-compare_model.R: dCor against tna's dCor squared; Canberra leaving out only the
  # cells zero in both networks against tna leaving out the cells zero in either.
  expect_known_difference(subset(global, metric %in% conventions), subset(expected, metric %in% conventions),
                          reference = "tna", covers = "global_differences()",
                          quantity = "distance correlation and Canberra, early vs late", tolerance = 1e-10,
                          key = c("pair", "category", "metric"),
                          reason = "Nestimate returns dCor (tna 1.3.1 its square) and counts a cell present in one network only in the Canberra distance (tna 1.3.1 leaves it out).")
  pick <- \(name) metrics$value[metrics$metric == name]
  expect_reference(summary(comparison, digits = 15L),
                   data.frame(pair = "early vs late", mean_abs_diff = pick("Mean Abs. Diff."),
                              max_abs_diff = pick("Max Abs. Diff."), pearson = pick("Pearson"),
                              spearman = pick("Spearman"), cosine = pick("Cosine"), jaccard = pick("Jaccard")),
                   reference = "tna", covers = "summary.net_network_comparison()",
                   quantity = "summary, early vs late", tolerance = 1e-10, key = "pair")
})

test_that("network-level metrics equal the network metrics of tna::compare()", {
  skip_without_reference("tna")
  comparison <- compare_networks(early = build_network(halves$early, method = "relative"),
                                 late = build_network(halves$late, method = "relative"))
  metrics <- as.data.frame(tna::compare(tna::tna(halves$early), tna::tna(halves$late))$network_metrics)
  expected <- rbind(data.frame(network = "early", metric = metrics$metric, value = metrics$x),
                    data.frame(network = "late", metric = metrics$metric, value = metrics$y))
  sd_strength <- c("SD Out-Strength", "SD In-Strength")
  found <- network_metrics(comparison, digits = 15L)
  expect_reference(subset(found, !metric %in% sd_strength), subset(expected, !metric %in% sd_strength), reference = "tna",
                   covers = "network_metrics()", quantity = "network metrics, early and late", tolerance = 1e-10,
                   key = c("network", "metric"))
  expect_known_difference(subset(found, metric %in% sd_strength), subset(expected, metric %in% sd_strength),
                          reference = "tna", covers = "network_metrics()", quantity = "SD of strengths, early and late",
                          tolerance = 1e-10, key = c("network", "metric"),
                          reason = "tna 1.3.1 reports the SD of the column sums (in-strength) as 'SD Out-Strength' and of the row sums as 'SD In-Strength'; Nestimate reports them under their own names.")
  weights <- lapply(halves, \(h) tna::tna(h)$weights)
  expect_reference(subset(found, metric %in% sd_strength),
                   data.frame(network = rep(c("early", "late"), each = 2L), metric = rep(sd_strength, 2L),
                              value = c(stats::sd(rowSums(weights$early)), stats::sd(colSums(weights$early)),
                                        stats::sd(rowSums(weights$late)), stats::sd(colSums(weights$late)))),
                   reference = "independent: SD of row and column sums of tna's weights", covers = "network_metrics()",
                   quantity = "SD of strengths from the definition", tolerance = 1e-12, key = c("network", "metric"))
})
