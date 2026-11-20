# predict_links(), evaluate_links(), summary.net_link_prediction(): link prediction scores on an undirected binary
# graph (Liben-Nowell and Kleinberg 2007; Zhou, Lu and Zhang 2009 for resource allocation; Katz 1953) and their
# evaluation. References: the six scores written from their definitions in base R, igraph::similarity() for Jaccard
# and Adamic-Adar, and the area under the ROC curve as the Mann-Whitney statistic of stats::wilcox.test().
# Re-expressed from local_testing_and_equivalence/test-equiv-link-prediction.R (108 random graphs by their largest
# difference). The prediction table is read from the result's field `predictions` (API gap: summary() returns only
# per-method aggregates).

ho_graph <- function(n, density, seed) {
  set.seed(seed)
  adjacency <- matrix(0L, n, n, dimnames = list(paste0("n", seq_len(n)), paste0("n", seq_len(n))))
  upper <- which(upper.tri(adjacency))
  adjacency[sample(upper, max(1L, round(length(upper) * density)))] <- 1L
  adjacency + t(adjacency)
}

# The six scores of every unordered pair, one row per method and pair.
ho_link_scores <- function(adjacency, katz) {
  degree <- rowSums(adjacency)
  common <- adjacency %*% adjacency
  union <- outer(degree, degree, "+") - common
  scores <- list(
    common_neighbors = common,
    resource_allocation = adjacency %*% diag(ifelse(degree > 0, 1 / degree, 0)) %*% adjacency,
    adamic_adar = adjacency %*% diag(ifelse(degree > 1, 1 / log(degree), 0)) %*% adjacency,
    jaccard = ifelse(union > 0, common / union, 0),
    preferential_attachment = outer(degree, degree),
    katz = solve(diag(nrow(adjacency)) - katz * adjacency) - diag(nrow(adjacency)))
  cells <- which(upper.tri(adjacency), arr.ind = TRUE)
  nodes <- rownames(adjacency)
  do.call(rbind, lapply(names(scores), \(method) {
    data.frame(method = method, from = nodes[cells[, 1L]], to = nodes[cells[, 2L]], score = scores[[method]][cells])
  }))
}

inputs <- list(list(n = 8L, density = 0.3, seed = 1001L), list(n = 12L, density = 0.2, seed = 1002L),
               list(n = 16L, density = 0.4, seed = 1003L))

test_that("the six link-prediction scores equal their definitions", {
  invisible(lapply(inputs, \(input) {
    adjacency <- ho_graph(input$n, input$density, input$seed)
    found <- predict_links(adjacency, weighted = FALSE, exclude_existing = FALSE, include_self = FALSE, katz_damping = 0.05)
    expect_reference(subset(found$predictions, select = c(method, from, to, score)), ho_link_scores(adjacency, 0.05),
                     reference = "independent: link-prediction scores from their definitions, base R",
                     covers = "predict_links()", quantity = sprintf("scores, %d nodes, density %.1f", input$n, input$density),
                     tolerance = 1e-8, key = c("method", "from", "to"))
  }))
})

test_that("Jaccard and Adamic-Adar scores equal igraph::similarity()", {
  skip_without_reference("igraph")
  invisible(lapply(inputs, \(input) {
    adjacency <- ho_graph(input$n, input$density, input$seed)
    graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", diag = FALSE)
    cells <- which(upper.tri(adjacency), arr.ind = TRUE)
    nodes <- rownames(adjacency)
    expected <- rbind(
      data.frame(method = "jaccard", from = nodes[cells[, 1L]], to = nodes[cells[, 2L]],
                 score = igraph::similarity(graph, method = "jaccard")[cells]),
      data.frame(method = "adamic_adar", from = nodes[cells[, 1L]], to = nodes[cells[, 2L]],
                 score = igraph::similarity(graph, method = "invlogweighted")[cells]))
    found <- predict_links(adjacency, methods = c("jaccard", "adamic_adar"), weighted = FALSE, exclude_existing = FALSE)
    expect_reference(subset(found$predictions, select = c(method, from, to, score)), expected, reference = "igraph",
                     covers = "predict_links()", quantity = sprintf("Jaccard and Adamic-Adar, %d nodes", input$n),
                     tolerance = 1e-8, key = c("method", "from", "to"))
  }))
})

test_that("AUC of evaluate_links() equals the Mann-Whitney statistic of stats::wilcox.test()", {
  adjacency <- ho_graph(16L, 0.3, 1004L)
  held_out <- data.frame(from = c("n1", "n2", "n3", "n5", "n7"), to = c("n9", "n11", "n12", "n14", "n16"))
  found <- predict_links(adjacency, weighted = FALSE, exclude_existing = TRUE, katz_damping = 0.05)
  table <- found$predictions
  is_true <- paste(table$from, table$to) %in% paste(held_out$from, held_out$to) |
    paste(table$to, table$from) %in% paste(held_out$from, held_out$to)
  auc <- vapply(split(data.frame(score = table$score, is_true = is_true), table$method), \(rows) {
    test <- suppressWarnings(stats::wilcox.test(rows$score[rows$is_true], rows$score[!rows$is_true], exact = FALSE))
    unname(test$statistic) / (sum(rows$is_true) * sum(!rows$is_true))
  }, numeric(1))
  expect_reference(subset(evaluate_links(found, held_out), select = c(method, auc)),
                   data.frame(method = names(auc), auc = unname(auc)), reference = "stats", covers = "evaluate_links()",
                   quantity = "area under the ROC curve", tolerance = 1e-12, key = "method",
                   note = "Mann-Whitney U / (n_true n_false), ties counted one half")
})

test_that("the per-method summary equals the mean, sd, max and min of the scores", {
  adjacency <- ho_graph(12L, 0.3, 1005L)
  found <- predict_links(adjacency, weighted = FALSE, exclude_existing = TRUE, katz_damping = 0.05)
  scores <- ho_link_scores(adjacency, 0.05)
  absent <- adjacency[cbind(match(scores$from, rownames(adjacency)), match(scores$to, rownames(adjacency)))] == 0
  candidates <- subset(scores, absent)
  aggregate_by <- \(f) vapply(split(candidates$score, candidates$method), f, numeric(1))
  expected <- data.frame(method = names(aggregate_by(mean)), n_predictions = aggregate_by(length),
                         score_mean = aggregate_by(mean), score_sd = aggregate_by(stats::sd),
                         score_max = aggregate_by(max), score_min = aggregate_by(min))
  expect_reference(summary(found), expected, reference = "independent: link-prediction scores from their definitions, base R",
                   covers = c("predict_links()", "summary.net_link_prediction()"), quantity = "per-method summary",
                   tolerance = 1e-8, key = "method")
})
