# compare_model(): the comparison of two weight matrices (edge-level metrics, 22 summary metrics of deviation,
# correlation, dissimilarity, similarity and pattern, and network-level metrics) against tna::compare() 1.3.1, for
# matrices, networks, tna models (the default method converts them), cograph networks and two members of a group.
# Re-expressed from local_testing_and_equivalence/test-equiv-summary-compare.R (30 random draws of
# group_regulation, maxima of differences) and test-equiv-network-similarity.R (random matrices; its comparison of
# an internal worker with its pre-refactor formulas compares the package with itself and stays there).
# compare_model() returns a list of tables with no accessor, so they are read as fields (API gap).
#
# Three conventions differ from tna 1.3.1, each shown as a known difference beside the definition it follows:
# - "Distance" correlation: Nestimate returns dCor (Szekely, Rizzo & Bakirov 2007); tna returns dCor squared.
# - Canberra: Nestimate sums |x - y| / (|x| + |y|) over every cell but the 0 / 0 ones (Lance & Williams), so a cell
#   present in one network only adds 1; tna leaves out every cell that is zero in either network.
# - "SD Out-Strength" / "SD In-Strength": tna 1.3.1 reports the SD of the column sums as out-strength and of the
#   row sums as in-strength; Nestimate reports the SD of the row sums (out) and column sums (in).
# The edge metric difference_index, (x - y) / y, is not defined where y = 0 (tna Inf, Nestimate NA).

halves <- list(x = tna::group_regulation[seq_len(200L), ], y = tna::group_regulation[1000L + seq_len(200L), ])
conventions <- c("Distance", "Canberra")
sd_strength <- c("SD Out-Strength", "SD In-Strength")

rs_random_matrix <- function(n, seed) {
  set.seed(seed)
  states <- paste0("S", seq_len(n))
  matrix(stats::runif(n * n), n, n, dimnames = list(states, states))
}

# Distance correlation of two vectors from its definition: double-centred distance matrices.
rs_distance_correlation <- function(x, y) {
  centre <- \(d) d - outer(rowMeans(d), colMeans(d), `+`) + mean(d)
  a <- centre(as.matrix(stats::dist(x)))
  b <- centre(as.matrix(stats::dist(y)))
  sqrt(mean(a * b) / sqrt(mean(a * a) * mean(b * b)))
}

# Canberra distance with 0 / 0 terms counted as 0 (Lance & Williams 1966).
rs_canberra <- function(x, y) {
  keep <- abs(x) + abs(y) > 0
  sum(abs(x - y)[keep] / (abs(x) + abs(y))[keep])
}

rs_summary_check <- function(found, reference, covers, label, tolerance) {
  metrics <- as.data.frame(reference$summary_metrics)
  expect_reference(subset(found$summary_metrics, !metric %in% conventions), subset(metrics, !metric %in% conventions),
                   reference = "tna", covers = covers, quantity = sprintf("%s, summary metrics", label),
                   tolerance = tolerance, key = c("category", "metric"))
}

test_that("summary metrics of two random matrices equal tna::compare(); dCor and Canberra their definitions", {
  skip_without_reference("tna")
  invisible(lapply(c(4L, 7L, 12L), \(n) {
    x <- rs_random_matrix(n, 100L + n)
    y <- rs_random_matrix(n, 200L + n)
    found <- compare_model(x, y)
    reference <- tna::compare(x, y, scaling = "none")
    rs_summary_check(found, reference, c("compare_model()", "compare_model.matrix()"),
                     sprintf("%d x %d random matrices", n, n), 1e-12)
    expect_reference(subset(found$summary_metrics, metric %in% conventions),
                     data.frame(category = c("Correlations", "Dissimilarities"), metric = conventions,
                                value = c(rs_distance_correlation(as.vector(x), as.vector(y)), rs_canberra(as.vector(x), as.vector(y)))),
                     reference = "independent: distance correlation and Canberra distance from their definitions",
                     covers = c("compare_model()", "compare_model.matrix()"),
                     quantity = sprintf("%d x %d random matrices, dCor and Canberra", n, n), tolerance = 1e-12,
                     key = c("category", "metric"))
    expect_known_difference(subset(found$summary_metrics, metric == "Distance"),
                            subset(as.data.frame(reference$summary_metrics), metric == "Distance"),
                            reference = "tna", covers = "compare_model()",
                            quantity = sprintf("%d x %d random matrices, distance correlation", n, n),
                            tolerance = 1e-12, key = c("category", "metric"),
                            reason = "Nestimate returns the distance correlation dCor; tna 1.3.1 returns its square (Nestimate did so until 1.0.5).")
  }))
})

test_that("edge, summary and network metrics of two transition networks equal tna::compare() per scaling", {
  skip_without_reference("tna")
  model_x <- build_network(halves$x, method = "tna")
  model_y <- build_network(halves$y, method = "tna")
  tna_x <- tna::tna(halves$x)
  tna_y <- tna::tna(halves$y)
  scalings <- c("none", "minmax", "max", "rank", "log1p", "softmax", "quantile")
  covers <- c("compare_model()", "compare_model.netobject()")
  # One stacked table per part, keyed by the scaling: the two sides of every scaling in one comparison.
  runs <- lapply(scalings, \(scaling) list(found = compare_model(model_x, model_y, scaling = scaling),
                                           tna = tna::compare(tna_x, tna_y, scaling = scaling)))
  stack <- \(part, side) do.call(rbind, lapply(seq_along(scalings), \(i) {
    cbind(scaling = scalings[i], as.data.frame(runs[[i]][[side]][[part]]))
  }))
  summary_found <- stack("summary_metrics", "found")
  summary_tna <- stack("summary_metrics", "tna")
  summary_key <- c("scaling", "category", "metric")
  expect_reference(subset(summary_found, !metric %in% conventions), subset(summary_tna, !metric %in% conventions),
                   reference = "tna", covers = covers, quantity = "summary metrics, seven scalings", tolerance = 1e-10,
                   key = summary_key)
  expect_known_difference(subset(summary_found, metric == "Distance"), subset(summary_tna, metric == "Distance"),
                          reference = "tna", covers = "compare_model()", quantity = "distance correlation, seven scalings",
                          tolerance = 1e-10, key = summary_key,
                          reason = "Nestimate returns the distance correlation dCor; tna 1.3.1 returns its square.")
  # Canberra differs only where a scaling keeps zero weights; softmax and quantile leave none.
  with_zeros <- c("none", "minmax", "max", "rank", "log1p")
  expect_known_difference(subset(summary_found, metric == "Canberra" & scaling %in% with_zeros),
                          subset(summary_tna, metric == "Canberra" & scaling %in% with_zeros),
                          reference = "tna", covers = "compare_model()", quantity = "Canberra, scalings that keep zeros",
                          tolerance = 1e-10, key = summary_key,
                          reason = "Nestimate leaves out only the cells that are zero in both networks (a cell present in one network adds 1); tna 1.3.1 leaves out every cell that is zero in either.")
  expect_reference(subset(summary_found, metric == "Canberra" & !scaling %in% with_zeros),
                   subset(summary_tna, metric == "Canberra" & !scaling %in% with_zeros), reference = "tna",
                   covers = covers, quantity = "Canberra, scalings without zeros", tolerance = 1e-10, key = summary_key)
  # Three edge ratios are not defined where their denominator is 0 (tna Inf or NaN, Nestimate NA).
  ratios <- c("difference_index", "similarity_strength_index", "standardized_score_inflation")
  edges_found <- stack("edge_metrics", "found")
  edges_tna <- stack("edge_metrics", "tna")
  edge_key <- c("scaling", "source", "target")
  expect_reference(edges_found, edges_tna[setdiff(names(edges_tna), ratios)], reference = "tna", covers = covers,
                   quantity = "edge metrics, seven scalings", tolerance = 1e-10, key = edge_key)
  long <- \(edges) do.call(rbind, lapply(ratios, \(ratio) {
    data.frame(edges[c("scaling", "source", "target")], ratio = ratio, value = edges[[ratio]])
  }))
  ratios_found <- long(edges_found)
  ratios_tna <- long(edges_tna)
  defined <- is.finite(ratios_tna$value)
  ratio_key <- c("scaling", "source", "target", "ratio")
  found_defined <- do.call(paste, ratios_found[ratio_key]) %in% do.call(paste, ratios_tna[defined, ratio_key])
  expect_reference(ratios_found[found_defined, ], ratios_tna[defined, ], reference = "tna", covers = covers,
                   quantity = "edge ratios where defined", tolerance = 1e-10, key = ratio_key)
  expect_not_comparable(ratios_found[!found_defined, ], ratios_tna[!defined, ], reference = "tna",
                        covers = "compare_model()", quantity = "edge ratios where the denominator is 0", key = ratio_key,
                        reason = "A ratio with a zero denominator is not defined: tna returns Inf or NaN, Nestimate NA.")
  network_found <- stack("network_metrics", "found")
  network_tna <- stack("network_metrics", "tna")
  expect_reference(subset(network_found, !metric %in% sd_strength), subset(network_tna, !metric %in% sd_strength),
                   reference = "tna", covers = covers, quantity = "network metrics, seven scalings", tolerance = 1e-10,
                   key = c("scaling", "metric"))
  expect_known_difference(subset(network_found, metric %in% sd_strength), subset(network_tna, metric %in% sd_strength),
                          reference = "tna", covers = "compare_model()", quantity = "SD of strengths, seven scalings",
                          tolerance = 1e-10, key = c("scaling", "metric"),
                          reason = "tna 1.3.1 reports the SD of the column sums (in-strength) as 'SD Out-Strength' and of the row sums as 'SD In-Strength'; Nestimate reports them under their own names.")
})

test_that("the SDs of out- and in-strength equal the SDs of the row and column sums", {
  skip_without_reference("tna")
  x <- tna::tna(halves$x)$weights
  y <- tna::tna(halves$y)$weights
  found <- compare_model(build_network(halves$x, method = "tna"), build_network(halves$y, method = "tna"))
  expect_reference(subset(found$network_metrics, metric %in% sd_strength),
                   data.frame(metric = sd_strength, x = c(stats::sd(rowSums(x)), stats::sd(colSums(x))),
                              y = c(stats::sd(rowSums(y)), stats::sd(colSums(y)))),
                   reference = "independent: SD of row and column sums of tna's weights",
                   covers = c("compare_model()", "compare_model.netobject()"), quantity = "SD of strengths",
                   tolerance = 1e-12, key = "metric")
})

test_that("tna models, cograph networks and two group members give tna::compare()'s summary metrics", {
  skip_without_reference("tna")
  skip_without_reference("cograph")
  tna_x <- tna::tna(halves$x)
  tna_y <- tna::tna(halves$y)
  reference <- tna::compare(tna_x, tna_y)
  rs_summary_check(compare_model(tna_x, tna_y), reference, "compare_model.default()", "tna models", 1e-10)
  rs_summary_check(compare_model(cograph::as_cograph(tna_x$weights), cograph::as_cograph(tna_y$weights)), reference,
                   "compare_model.cograph_network()", "cograph networks", 1e-10)
  grouped <- build_network(rbind(cbind(halves$x, grp = "x"), cbind(halves$y, grp = "y")), method = "tna", group = "grp")
  rs_summary_check(compare_model(grouped, i = "x", j = "y"), reference, "compare_model.netobject_group()",
                   "two group members", 1e-10)
})
