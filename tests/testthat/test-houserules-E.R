# House-rule fixes in the comparison / pruning / permutation family:
# classed input errors with named messages, corrected printed text.

he_groups <- function() {
  build_network(group_regulation_long, method = "relative",
                actor = "Actor", action = "Action", time = "Time",
                group = "Achiever")
}

he_matrix <- function() {
  matrix(c(0.2, 0.5, 0.3, 0.1, 0.6, 0.3, 0.4, 0.4, 0.2), 3, 3, byrow = TRUE,
         dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
}

test_that("compare_model() group members are checked with a classed, named error", {
  nets <- he_groups()
  expect_error(compare_model(nets, i = 1, j = 1), "two different networks",
               class = "nestimate_bad_input")
  expect_error(compare_model(nets, i = "nope"), "name or index",
               class = "nestimate_bad_input")
  expect_error(compare_model(nets, i = 1, j = 9), class = "nestimate_bad_input")
})

test_that("compare_model() refuses a scaling or plot type typo by class", {
  a <- he_matrix()
  expect_error(compare_model(a, a * 2, scaling = "zzz"), "`scaling`",
               class = "nestimate_bad_input")
  cmp <- compare_model(a, a * 2, network = FALSE)
  expect_error(plot(cmp, type = "zzz"), "`type`",
               class = "nestimate_bad_input")
})

test_that("compare_model(scaling = 'log1p') refuses weights at or below -1", {
  a <- he_matrix()
  below <- a
  below["A", "B"] <- -2
  at <- a
  at["A", "B"] <- -1
  expect_error(compare_model(a, below, scaling = "log1p"), "above -1",
               class = "nestimate_bad_input")
  expect_error(compare_model(at, a, scaling = "log1p", network = FALSE),
               class = "nestimate_bad_input")
  # the domain it is defined on is unaffected
  ok <- compare_model(a, a * 2, scaling = "log1p", network = FALSE)
  expect_equal(ok$matrices$x, log1p(a))
})

test_that("compare_model() size and missing-weight errors are named", {
  a <- he_matrix()
  expect_error(compare_model(a, matrix(1, 2, 2)), "same number of nodes",
               class = "nestimate_bad_input")
  with_na <- a
  with_na[1, 2] <- NA
  expect_error(compare_model(a, with_na), "missing weights",
               class = "nestimate_bad_input")
  expect_error(compare_model(matrix("a", 2, 2), matrix("a", 2, 2)),
               "must be numeric", class = "nestimate_bad_input")
})

test_that("compare_model.default() raises a classed error", {
  # called directly: dispatch needs S3method(compare_model, default) in
  # NAMESPACE, which the next devtools::document() adds
  expect_error(compare_model.default(list(1)), "Got: list",
               class = "nestimate_bad_input")
})

test_that("print.net_comparison() has no empty centrality block", {
  a <- he_matrix()
  none <- suppressWarnings(compare_model(a, a * 2, measures = "foo",
                                         network = FALSE))
  expect_false(any(grepl("Centrality correlations",
                         capture.output(print(none)))))
  some <- suppressMessages(compare_model(a, a * 2, measures = "InStrength",
                                         network = FALSE))
  expect_true(any(grepl("Centrality correlations",
                        capture.output(print(some)))))
})

test_that(".network_similarity() length mismatch is classed", {
  expect_error(.network_similarity(1:3, 1:4), class = "nestimate_bad_input")
  expect_error(.network_similarity(diag(2), diag(3)),
               class = "nestimate_bad_input")
})

test_that("subtract_networks() names itself for an unsupported object", {
  err <- expect_error(subtract_networks(list(a = 1), diag(2)),
                      class = "nestimate_bad_input")
  expect_match(conditionMessage(err), "subtract_networks()", fixed = TRUE)
  expect_no_match(conditionMessage(err), "compare_model", fixed = TRUE)
  expect_error(subtract_networks(matrix(1:6, 2), matrix(1:6, 2)),
               "must be square", class = "nestimate_bad_input")
})

test_that("summary() network methods validate with classed errors", {
  expect_error(summary.netobject(list(weights = diag(2))),
               class = "nestimate_bad_input")
  net <- he_groups()[[1]]
  net$weights <- 1
  expect_error(summary(net), "square matrix", class = "nestimate_bad_input")
  expect_error(summary(he_groups(), combined = NA),
               class = "nestimate_bad_input")
})

test_that("print of an uncombined group summary lists every group", {
  nets <- he_groups()
  out <- capture.output(print(summary(nets, combined = FALSE)))
  expect_identical(grep("^Group:", out, value = TRUE),
                   paste0("Group: ", names(nets), " "))
})

test_that("nct() refuses a p_adjust or na_method typo by class", {
  d1 <- srl_strategies[1:100, ]
  d2 <- srl_strategies[101:200, ]
  expect_error(nct(d1, d2, iter = 2, p_adjust = "zzz"), "`p_adjust`",
               class = "nestimate_bad_input")
  expect_error(nct(d1, d2, iter = 2, na_method = "zzz"), "`na_method`",
               class = "nestimate_bad_input")
})

test_that("compare_networks() argument errors carry nestimate_compare_error", {
  nets <- he_groups()
  expect_error(compare_networks(nets, iter = 1),
               class = "nestimate_compare_error")
  expect_error(compare_networks(nets, alpha = 2),
               class = "nestimate_bad_input")
  cmp <- compare_networks(nets)
  expect_error(plot(cmp, type = "edges", top_n = 0),
               class = "nestimate_compare_error")
  expect_error(edge_differences(cmp, digits = -1),
               class = "nestimate_compare_error")
})

test_that("permutation() and bayes_compare() argument errors are named", {
  nets <- he_groups()
  expect_error(permutation(nets[[1]], nets[[2]], iter = 1), "`iter`",
               class = "nestimate_bad_input")
  expect_error(permutation(nets[[1]], nets[[2]], alpha = 2), "`alpha`",
               class = "nestimate_bad_input")
  expect_error(permutation(nets[[1]], nets[[2]], iter = 5, seed = "a"),
               "`seed`", class = "nestimate_bad_input")
  expect_error(permutation(nets[[1]], 1), "netobjects",
               class = "nestimate_bad_input")
  expect_error(bayes_compare(nets[[1]], nets[[2]], draws = 1), "`draws`",
               class = "nestimate_bad_input")
  expect_error(bayes_compare(nets[[1]], nets[[2]], ci = 2), "`ci`",
               class = "nestimate_bad_input")
})

test_that("the permutation print names the floor that binds", {
  nets <- he_groups()
  res <- permutation(nets[[1]], nets[[2]], iter = 10, seed = 1,
                     actor = "Actor")
  out <- capture.output(print(res))
  expect_true(any(grepl("with 10 iterations no p-value can fall below",
                        out, fixed = TRUE)))
  expect_false(any(grepl("actors no p-value", out, fixed = TRUE)))
})

test_that("permutation_diagnostics() says y must be NULL for a group", {
  nets <- he_groups()
  expect_error(
    permutation_diagnostics(nets, nets[[1]], actor = "Actor", iter = 10),
    "`y` must be NULL", class = "nestimate_bad_input")
  expect_error(
    permutation_diagnostics(nets[[1]], actor = "Actor", iter = 10),
    "must be netobjects", class = "nestimate_bad_input")
})

test_that("rename_models() errors are named", {
  nets <- he_groups()
  expect_error(rename_models(nets, c("a", "a")), "unique",
               class = "nestimate_bad_input")
  expect_error(rename_models(nets, "a"), "one name per network",
               class = "nestimate_bad_input")
  expect_error(rename_models(nets, c("a", NA)), "missing or empty",
               class = "nestimate_bad_input")
})

test_that("net_prune() refuses a missing or non-scalar setting by class", {
  net <- he_groups()[[1]]
  expect_error(net_prune(net, level = NA), "'level'",
               class = "nestimate_bad_input")
  expect_error(net_prune(net, lowest = c(0.1, 0.2)), "'lowest'",
               class = "nestimate_bad_input")
  expect_error(net_prune(net, threshold = "a"), "'threshold'",
               class = "nestimate_bad_input")
})
