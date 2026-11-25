# Regression tests for the 2026-11-20 review findings (area "misc"):
# A1-A10, B1-B3, B6, B10-B19 and P9. Every expected value is computed by hand
# from the definition (see the comment above each expectation).
testthat::skip_on_cran()

# ---- A1: shortest-path ties are relative to the path length ----

test_that("A1: a strictly shorter path is not a tie at any weight scale", {
  # A->B 1e-12, B->C 1e-12, A->C 3e-12 as lengths (invert = FALSE): the only
  # shortest A->C path is A-B-C (2e-12 < 3e-12), so B lies on 1 of 1 paths.
  W <- matrix(0, 3, 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  W["A", "B"] <- 1e-12
  W["B", "C"] <- 1e-12
  W["A", "C"] <- 3e-12
  btw <- .betweenness(W, directed = TRUE, invert = FALSE)
  expect_equal(unname(btw), c(0, 1, 0))
  sp <- .floyd_warshall_sp(W, invert = FALSE)
  expect_equal(sp$D[1, 3], 2e-12)
  expect_equal(sp$sigma[1, 3], 1)
  # scale equivariance: the same network in units of 1e-12
  expect_equal(.betweenness(W * 1e12, directed = TRUE, invert = FALSE), btw)
})

test_that("A1: a tie broken only by floating-point addition stays a tie", {
  # distances 1/10 + 1/15 and 1/6 are equal: two shortest A->C paths, one
  # through B, so B's betweenness is 1/2.
  W <- matrix(0, 3, 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  W["A", "B"] <- 10
  W["B", "C"] <- 15
  W["A", "C"] <- 6
  expect_equal(.floyd_warshall_sp(W)$sigma[1, 3], 2)
  expect_equal(unname(.betweenness(W, directed = TRUE)), c(0, 0.5, 0))
})

# ---- A2, A7: hypergraph membership is a positive incidence cell ----

test_that("A2: clique centrality ignores non-member (non-positive) cells", {
  # g = {A, B} (C's summed weight -1 is not membership), h = {A, D}: the
  # clique expansion is the path B - A - D plus isolated C, whose Perron
  # vector is (A, B, C, D) = (1/sqrt(2), 1/2, 0, 1/2).
  d <- data.frame(player = c("A", "B", "C", "A", "D"),
                  group = c("g", "g", "g", "h", "h"), w = c(1, 1, -1, 1, 1))
  hg <- bipartite_groups(d, "player", "group", "w")
  cc <- hypergraph_centrality(hg, type = "clique")$clique
  expect_equal(unname(cc), c(1 / sqrt(2), 0.5, 0, 0.5), tolerance = 1e-8)
})

test_that("A7: random-walk Laplacian uses members only", {
  # C's cell in g is -2 (not a member): the walk is that of g = {A, B},
  # h = {A, C} with unit cells, whose stationary distribution is the degree
  # share (2, 1, 1) / 4, and which equals the binary hypergraph's Laplacian.
  signed <- bipartite_groups(
    data.frame(p = c("A", "B", "C", "A", "C"), g = c("g", "g", "g", "h", "h"),
               w = c(1, 1, -2, 1, 1)), "p", "g", "w")
  binary <- bipartite_groups(
    data.frame(p = c("A", "B", "A", "C"), g = c("g", "g", "h", "h")), "p", "g")
  L_signed <- hypergraph_laplacian(signed, type = "random_walk")
  L_binary <- hypergraph_laplacian(binary, type = "random_walk")
  expect_equal(unname(attr(L_signed, "pi")), c(0.5, 0.25, 0.25))
  expect_equal(unclass(L_signed)[, ], unclass(L_binary)[, ])
  expect_equal(unname(L_signed["A", "B"]), -sqrt(2) / 4)
})

# ---- A3, A4, B2: direction is decided after conversion and alignment ----

test_that("A3: a symmetric matrix with reordered columns stays undirected", {
  nm <- c("A", "B", "C")
  W <- matrix(0, 3, 3, dimnames = list(nm, nm))
  W["A", "B"] <- W["B", "A"] <- 1
  res <- compare_model(W, W[, c("A", "C", "B")])
  edges <- subset(res$network_metrics, metric == "Edge Count")
  # one undirected edge A-B on both sides
  expect_identical(c(edges$x, edges$y), c(1, 1))
})

test_that("A4: a tna model keeps its direction as the second argument", {
  skip_if_not_installed("tna")
  x <- tna::tna(data.frame(V1 = c("A", "B"), V2 = c("B", "A")))
  res <- compare_model(x, x)
  edges <- subset(res$network_metrics, metric == "Edge Count")
  # two arcs, A->B and B->A, on both sides
  expect_identical(c(edges$x, edges$y), c(2, 2))
})

test_that("B2: a change in the lower triangle is never hidden", {
  x <- matrix(0, 2, 2)
  y <- x
  y[2, 1] <- 1
  d <- subtract_networks(x, y)
  expect_true(d$directed)
  expect_identical(d$n_edges, 1L)
  expect_equal(d$weights[2, 1], -1)
  # either order gives the same direction
  expect_true(subtract_networks(y, x)$directed)
  # a directed network on either side makes the difference directed
  net <- build_network(data.frame(T1 = c("A", "B"), T2 = c("B", "A")),
                       method = "relative")
  sym <- matrix(0.5, 2, 2, dimnames = dimnames(net$weights))
  expect_true(subtract_networks(sym, net)$directed)
})

# ---- B3, B6: pruning ----

test_that("B3: the disparity filter reads the positive edges only", {
  # node 1: positive weights 0.4 (to 2) and 0.2 (to 3), strength 0.6, degree
  # 2; alpha(1-2) = 1 - 0.4/0.6 = 1/3, alpha(1-3) = 1 - 0.2/0.6 = 2/3. Nodes
  # 2 and 3 have degree 1 (alpha 1). The -0.3 edge to 4 is never kept.
  w <- matrix(0, 4, 4)
  w[1, 2:4] <- w[2:4, 1] <- c(0.4, 0.2, -0.3)
  net <- subtract_networks(w, 0 * w)
  expect_false(net$directed)
  strict <- net_prune(net, "disparity", level = 0.01)
  expect_true(all(strict$weights == 0))
  loose <- net_prune(net, "disparity", level = 0.5)
  expect_equal(unname(loose$weights[1, ]), c(0, 0.4, 0, 0))
  expect_equal(unname(loose$weights[, 1]), c(0, 0.4, 0, 0))
})

test_that("B6: a bootstrap of other data with equal weights is refused", {
  d <- data.frame(T1 = c("A", "A"), T2 = c("A", "B"))
  x <- build_network(d, method = "relative")
  y <- build_network(d[rep(1:2, 200), ], method = "relative")
  expect_equal(unname(x$weights), unname(y$weights))
  by <- bootstrap_network(y, iter = 19, seed = 4)
  expect_error(net_prune(x, "bootstrap", boot = by),
               class = "nestimate_bad_input")
  bx <- bootstrap_network(x, iter = 19, seed = 4)
  expect_s3_class(net_prune(x, "bootstrap", boot = bx), "netobject")
})

# ---- A5: absorption on the thresholded chain ----

test_that("A5: absorption agrees with hitting probabilities under tol", {
  # A -> B (0.001) is below tol = 0.01: from A the chain reaches only C+D
  P <- rbind(c(0, .001, .999, 0), c(0, 1, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0))
  dimnames(P) <- list(LETTERS[1:4], LETTERS[1:4])
  cs <- chain_structure(P, tol = 0.01)
  expect_equal(unname(cs$hitting_probabilities["A", c("B", "C", "D")]),
               c(0, 1, 1))
  expect_equal(unname(cs$absorption_probabilities["A", ]), c(0, 1))
  expect_equal(colnames(cs$absorption_probabilities), c("B", "C+D"))
  # with the default tol the 0.001 edge counts
  expect_equal(unname(chain_structure(P)$absorption_probabilities["A", ]),
               c(0.001, 0.999))
})

# ---- A6, A8, A10: association rules ----

test_that("A6: two binary columns with one name are refused", {
  d <- data.frame(A = c(1, 0), A = c(0, 1), check.names = FALSE)
  expect_error(association_rules(d, format = "binary"),
               class = "nestimate_bad_input")
})

test_that("A8: rules from a wtna network built on an unnamed matrix", {
  # windows {V1, V2}, {V1}, {V2}: V1 -> V2 has support 1/3, confidence
  # (1/3) / (2/3) = 1/2 and lift (1/2) / (2/3) = 3/4; V2 -> V1 the same.
  X <- rbind(c(1, 1), c(1, 0), c(0, 1))
  net <- wtna(X, method = "cooccurrence", window_size = 1)
  r <- association_rules(net, min_support = 0.1, min_confidence = 0.1,
                         min_lift = 0)
  rules <- summary(r)
  expect_identical(nrow(rules), 2L)
  expect_setequal(paste(rules$antecedent, rules$consequent),
                  c("V1 V2", "V2 V1"))
  expect_equal(rules$support, rep(1 / 3, 2))
  expect_equal(rules$confidence, rep(0.5, 2))
  expect_equal(rules$lift, rep(0.75, 2))
})

test_that("A10: no complete window gives an empty result", {
  net <- wtna(data.frame(A = 1, B = 1), method = "cooccurrence",
              window_size = 2, mode = "overlapping")
  r <- association_rules(net)
  expect_identical(r$n_rules, 0L)
  expect_identical(r$n_transactions, 0L)
  expect_identical(nrow(summary(r)), 0L)
})

# ---- A9, P9: plots ----

test_that("A9: `cols` selects columns of a data.table", {
  skip_if_not_installed("data.table")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  df <- data.frame(id = 1:2, T1 = c("A", "B"), T2 = c("B", "A"))
  dt <- data.table::as.data.table(df)
  expect_equal(distribution_plot(dt, cols = c("T1", "T2")),
               distribution_plot(df, cols = c("T1", "T2")))
  s_dt <- sequence_plot(dt, cols = c("T1", "T2"))
  s_df <- sequence_plot(df, cols = c("T1", "T2"))
  expect_equal(s_dt$levels, c("A", "B"))
  expect_equal(s_dt$codes, s_df$codes)
})

test_that("P9: a group value equal to a column name does not exempt it", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"),
                  action = c("plan", "review"))
  expect_error(sequence_plot(d, type = "index", group = rep("action", 2),
                             legend = "none"),
               class = "nestimate_bad_input")
})

# ---- B1: split by group ----

test_that("B1: select = 'split' with a group puts each group in one half", {
  skip_if_not_installed("lme4")
  g <- rep(seq_len(11), each = 4)
  set.seed(4)
  d <- data.frame(g, x = rnorm(44), x2 = rnorm(44), y = rnorm(11)[g] + rnorm(44))
  fit <- outcome_model(d, "y", c("x", "x2"), group = "g", select = "split",
                       seed = 1, n_select = 1)
  used <- as.integer(rownames(stats::model.frame(fit$model)))
  fit_groups <- unique(d$g[used])
  sel_groups <- setdiff(unique(d$g), fit_groups)
  # half of the 11 groups (5) select, the other 6 are fitted, whole
  expect_length(sel_groups, 5L)
  expect_length(intersect(d$g[used], d$g[-used]), 0L)
  expect_identical(fit$n, 24L)
  one_group <- transform(d, g = 1L)
  expect_error(outcome_model(one_group, "y", c("x", "x2"), group = "g",
                             select = "split", seed = 1),
               class = "nestimate_bad_input")
})

# ---- B10, B19: link prediction ----

test_that("B10: every index runs on a one-node network", {
  methods <- c("common_neighbors", "resource_allocation", "adamic_adar",
               "jaccard", "preferential_attachment", "katz")
  invisible(lapply(methods, function(m) {
    r <- predict_links(matrix(0, 1, 1), methods = m)
    expect_identical(nrow(r$predictions), 0L)
    expect_identical(dim(r$scores[[m]]), c(1L, 1L))
  }))
})

test_that("B19: evaluate_links() without pairs gives a zero-row frame", {
  w <- matrix(c(0, 1, 1, 0), 2)
  ev <- evaluate_links(predict_links(w, methods = "jaccard"), w)
  expect_s3_class(ev, "data.frame")
  expect_identical(nrow(ev), 0L)
  expect_named(ev, c("method", "auc", "average_precision", "precision_at_5",
                     "precision_at_10", "precision_at_20"))
})

# ---- B11: path dependence on one state ----

test_that("B11: path_dependence() on a single state", {
  # one context A -> A, always followed by A: no entropy, no divergence
  r <- path_dependence(data.frame(T1 = "A", T2 = "A", T3 = "A"), min_count = 1)
  expect_identical(r$contexts$context, "A -> A")
  expect_equal(r$contexts$KL, 0)
  expect_equal(r$contexts$H_order1, 0)
})

# ---- B12, B15: magnitude difference ----

test_that("B12: equal weights on both views give zero difference", {
  # all-A sequences: counts 4, probabilities 1; both views are constant, so
  # they agree, at every scale that puts them on one range
  d <- data.frame(T1 = rep("A", 4), T2 = rep("A", 4))
  vals <- vapply(c("tna_range", "minmax", "rank_minmax"), function(s) {
    magnitude_difference(d, scale = s)$edges$value
  }, numeric(1))
  expect_equal(unname(vals), c(0, 0, 0))
  # every count 1, every probability 0.5: the counts map to 0.5
  d2 <- data.frame(T1 = c("A", "B", "A", "B"), T2 = c("A", "B", "B", "A"))
  md <- magnitude_difference(d2)
  expect_equal(md$edges$ftna, rep(0.5, 4))
  expect_equal(md$edges$value, rep(0, 4))
  # `scale = "none"` still compares the raw weights: |1 - 4| = 3
  expect_equal(magnitude_difference(d, scale = "none")$edges$value, 3)
})

test_that("B15: format = 'long' without action is refused", {
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  expect_error(magnitude_difference(d, format = "long"),
               class = "nestimate_bad_input")
})

# ---- B13: vertex bootstrap statistic names ----

test_that("B13: a custom statistic cannot reuse a statistic's name", {
  w <- matrix(1, 3, 3)
  diag(w) <- 0
  expect_error(vertex_bootstrap(w, iter = 5, seed = 1, statistics = "density",
                                statistic_fn = list(density = function(w) 5)),
               class = "nestimate_bad_input")
  expect_error(vertex_compare(w, w, iter = 5, seed = 1,
                              statistic_fn = list(mean_weight = function(w) 5)),
               class = "nestimate_bad_input")
  ok <- vertex_bootstrap(w, iter = 5, seed = 1, statistics = "density",
                         statistic_fn = list(five = function(w) 5))
  expect_s3_class(ok, "net_vertex_bootstrap")
})

# ---- B16: pathway evidence ranks by weight ----

test_that("B16: the strongest bridge is kept as evidence", {
  # A -> D is bridged by B (weights 1 + 1) and by C (5 + 5): C is kept
  w <- matrix(0, 4, 4, dimnames = rep(list(LETTERS[1:4]), 2))
  w["A", "B"] <- w["B", "D"] <- 1
  w["A", "C"] <- w["C", "D"] <- 5
  p <- predict_links(w, methods = "katz", weighted = TRUE)
  expect_identical(pathways(p, top = 1, max_evidence = 1), "A C -> D")
  # swapping the weights swaps the bridge
  w2 <- w
  w2["A", "B"] <- w2["B", "D"] <- 5
  w2["A", "C"] <- w2["C", "D"] <- 1
  p2 <- predict_links(w2, methods = "katz", weighted = TRUE)
  expect_identical(pathways(p2, top = 1, max_evidence = 1), "A B -> D")
})

# ---- B17: state colours follow labels ----

test_that("B17: a relabelled network is coloured by its labels", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  net <- build_network(d, method = "relative", labels = c(A = "Alpha"))
  expect_identical(state_colors(net)$state, c("Alpha", "B"))
  net <- set_state_colors(net, c(Alpha = "#FF0000", B = "#0000FF"))
  sp <- sequence_plot(net)
  expect_identical(sp$levels, c("Alpha", "B"))
  expect_identical(sp$palette, c("#FF0000", "#0000FF"))
  # a palette set before relabelling moves with its states
  plain <- set_state_colors(build_network(d, method = "relative"),
                            c(A = "#FF0000", B = "#0000FF"))
  relabelled <- .apply_node_labels(plain, c(A = "Alpha"))
  expect_identical(state_colors(relabelled)$color, c("#FF0000", "#0000FF"))
  expect_identical(state_colors(relabelled)$state, c("Alpha", "B"))
})

# ---- B18: integrative potential at extreme omega ----

test_that("B18: integrative potential stays finite for large omega", {
  # positions 1 (A) and 2 (B, favourable): weight share 2^w / (1 + 2^w)
  d <- data.frame(T1 = "A", T2 = "B")
  ip <- function(om) {
    sequence_metrics(d, favorable = "B", omega = om)$integrative_potential
  }
  expect_equal(ip(1), 2 / 3)
  expect_equal(ip(2), 4 / 5)
  expect_equal(ip(1100), 1)
  expect_equal(ip(-1100), 0)
})
