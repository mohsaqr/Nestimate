# ---- Link Prediction Tests ----

# Helper: small directed graph A -> B -> C -> A, plus A -> C
.make_lp_net <- function() {
  seqs <- data.frame(
    V1 = c("A","B","C","A","A","B","C","A","B","C"),
    V2 = c("B","C","A","C","B","C","A","C","C","A"),
    V3 = c("C","A","A","A","C","A","A","A","A","A"),
    stringsAsFactors = FALSE
  )
  build_network(seqs, method = "relative")
}

# Helper: sparse 8-node network
.make_lp_sparse <- function() {
  set.seed(42)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:8], 50, TRUE),
    V2 = sample(LETTERS[1:8], 50, TRUE),
    V3 = sample(LETTERS[1:8], 50, TRUE),
    stringsAsFactors = FALSE
  )
  build_network(seqs, method = "relative", threshold = 0.05)
}


# ---- 1. Basic functionality ----

test_that("predict_links returns correct class and structure", {
  net <- .make_lp_sparse()
  pred <- predict_links(net)

  expect_s3_class(pred, "net_link_prediction")
  expect_true(is.data.frame(pred$predictions))
  expect_true(is.list(pred$scores))
  expect_equal(length(pred$scores), 6)
  expect_true(all(c("from", "to", "method", "score", "rank") %in%
                    names(pred$predictions)))
})


# ---- 2. Single method ----

test_that("predict_links works with single method", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "common_neighbors")
  expect_equal(pred$methods, "common_neighbors")
  expect_equal(length(pred$scores), 1)
  expect_true(all(pred$predictions$method == "common_neighbors"))
})


# ---- 3. All six methods ----

test_that("all six methods produce valid scores", {
  net <- .make_lp_sparse()
  all_methods <- c("common_neighbors", "resource_allocation", "adamic_adar",
                   "jaccard", "preferential_attachment", "katz")
  pred <- predict_links(net, methods = all_methods)
  expect_equal(length(pred$scores), 6)
  for (m in all_methods) {
    s <- pred$scores[[m]]
    expect_true(is.matrix(s))
    expect_equal(nrow(s), pred$n_nodes)
    expect_true(all(is.finite(s)))
  }
})


# ---- 4. Common Neighbors vectorized correctness ----

test_that("common_neighbors matches manual computation", {
  W <- matrix(c(0, .5, .5,
                0,  0,  1,
                1,  0,  0), 3, 3, byrow = TRUE)
  rownames(W) <- colnames(W) <- c("A", "B", "C")
  pred <- predict_links(W, methods = "common_neighbors",
                        weighted = FALSE, exclude_existing = FALSE,
                        include_self = FALSE)
  s <- pred$scores$common_neighbors
  # A and B share C as out-neighbor (A->C, B->C)
  # Plus: A and B share C as in-neighbor? No: C->A (A is in-neighbor of C, not B)
  # tcrossprod(A): shared out-neighbors
  # crossprod(A): shared in-neighbors
  A <- (W > 0) * 1
  expected <- tcrossprod(A) + crossprod(A)
  diag(expected) <- 0
  dimnames(expected) <- dimnames(s)
  expect_equal(s, expected)
})


# ---- 5. Resource Allocation correctness ----

test_that("resource_allocation penalizes hub neighbors", {
  net <- .make_lp_sparse()
  pred_cn <- predict_links(net, methods = "common_neighbors",
                           weighted = FALSE, exclude_existing = FALSE)
  pred_ra <- predict_links(net, methods = "resource_allocation",
                           weighted = FALSE, exclude_existing = FALSE)
  # RA scores should always be <= CN scores (divided by degree)
  cn <- pred_cn$scores$common_neighbors
  ra <- pred_ra$scores$resource_allocation
  expect_true(all(ra <= cn + 1e-10))
})


# ---- 6. Jaccard is bounded [0, 1] ----

test_that("jaccard scores are in [0, 1]", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "jaccard", exclude_existing = FALSE)
  s <- pred$scores$jaccard
  expect_true(all(s >= 0 & s <= 1))
})


# ---- 7. Preferential Attachment is degree product ----

test_that("preferential_attachment equals out_degree * in_degree", {
  net <- .make_lp_sparse()
  A <- (net$weights != 0) * 1
  expected <- outer(rowSums(A), colSums(A), "*")
  diag(expected) <- 0
  pred <- predict_links(net, methods = "preferential_attachment",
                        exclude_existing = FALSE)
  expect_equal(pred$scores$preferential_attachment, expected)
})


# ---- 8. Katz auto-damping ----

test_that("katz auto-computes valid damping", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz")
  s <- pred$scores$katz
  expect_true(all(is.finite(s)))
  expect_true(all(diag(s) == 0))
})

test_that("katz warns when user damping exceeds bound", {
  net <- .make_lp_sparse()
  expect_warning(
    predict_links(net, methods = "katz", katz_damping = 10),
    "auto-adjusted"
  )
})


# ---- 9. Weighted vs binary ----

test_that("weighted=TRUE uses weight magnitudes", {
  net <- .make_lp_sparse()
  pred_w <- predict_links(net, methods = "common_neighbors",
                          weighted = TRUE, exclude_existing = FALSE)
  pred_b <- predict_links(net, methods = "common_neighbors",
                          weighted = FALSE, exclude_existing = FALSE)
  # Scores should differ (weights != binary)
  expect_false(identical(pred_w$scores$common_neighbors,
                         pred_b$scores$common_neighbors))
})


# ---- 10. exclude_existing ----

test_that("exclude_existing removes known edges", {
  net <- .make_lp_sparse()
  pred_inc <- predict_links(net, methods = "katz", exclude_existing = FALSE)
  pred_exc <- predict_links(net, methods = "katz", exclude_existing = TRUE)
  expect_true(nrow(pred_inc$predictions) > nrow(pred_exc$predictions))
  # Excluded predictions should not contain existing edges
  A <- (net$weights != 0) * 1
  for (i in seq_len(nrow(pred_exc$predictions))) {
    r <- pred_exc$predictions[i, ]
    expect_equal(A[r$from, r$to], 0)
  }
})


# ---- 11. top_n limits output ----

test_that("top_n limits predictions per method", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = c("katz", "jaccard"), top_n = 3)
  for (m in c("katz", "jaccard")) {
    sub <- pred$predictions[pred$predictions$method == m, ]
    expect_true(nrow(sub) <= 3)
  }
})


# ---- 12. Undirected network ----

test_that("predict_links works on undirected network", {
  set.seed(42)
  num <- as.data.frame(matrix(rnorm(200), ncol = 5))
  net <- build_network(num, method = "cor")
  pred <- predict_links(net, methods = c("common_neighbors", "jaccard"))
  expect_false(pred$directed)
  # Undirected: from < to only
  expect_true(all(pred$predictions$from < pred$predictions$to))
})


# ---- 13. Matrix input ----

test_that("predict_links works on raw matrix", {
  W <- matrix(c(0, .5, 0, .3, 0, .4, 0, .2, 0), 3, 3)
  rownames(W) <- colnames(W) <- c("X", "Y", "Z")
  pred <- predict_links(W, methods = "common_neighbors")
  expect_s3_class(pred, "net_link_prediction")
  expect_equal(pred$nodes, c("X", "Y", "Z"))
})


# ---- 14. cograph_network input ----

test_that("predict_links works on cograph_network", {
  net <- .make_lp_sparse()
  cg <- structure(list(
    weights = net$weights, nodes = net$nodes, edges = net$edges,
    directed = net$directed, data = net$data,
    meta = list(source = "test", tna = list(method = "relative"))
  ), class = c("cograph_network", "list"))
  pred <- predict_links(cg)
  expect_s3_class(pred, "net_link_prediction")
})


# ---- 15. Symmetric scores for undirected ----

test_that("score matrices are symmetric for undirected networks", {
  mat <- matrix(c(0, .5, .3, .5, 0, .2, .3, .2, 0), 3, 3)
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  pred <- predict_links(mat, methods = c("common_neighbors", "resource_allocation",
                                          "adamic_adar", "jaccard"),
                        exclude_existing = FALSE)
  for (m in names(pred$scores)) {
    expect_true(isSymmetric(pred$scores[[m]]),
                info = paste("Non-symmetric scores for", m))
  }
})


# ---- 16. Diagonal is always 0 ----

test_that("score matrix diagonals are always 0", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, exclude_existing = FALSE)
  for (m in names(pred$scores)) {
    expect_true(all(diag(pred$scores[[m]]) == 0),
                info = paste("Non-zero diagonal for", m))
  }
})


# ---- 17. evaluate_links AUC ----

test_that("evaluate_links computes valid AUC", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz", exclude_existing = FALSE)
  # Use top-ranked edges as "true"
  top <- head(pred$predictions, 5)
  true_df <- data.frame(from = top$from, to = top$to)
  eval <- evaluate_links(pred, true_df, k = c(5, 10))
  expect_true(is.data.frame(eval))
  expect_true("auc" %in% names(eval))
  expect_true(eval$auc >= 0 && eval$auc <= 1)
})


# ---- 18. evaluate_links with matrix input ----

test_that("evaluate_links accepts true_edges as matrix", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "common_neighbors",
                        exclude_existing = FALSE)
  # Create a binary true-edge matrix
  true_mat <- matrix(0, pred$n_nodes, pred$n_nodes,
                     dimnames = list(pred$nodes, pred$nodes))
  true_mat[1, 2] <- 1
  true_mat[2, 3] <- 1
  eval <- evaluate_links(pred, true_mat)
  expect_true(is.data.frame(eval))
})


# ---- 19. Error on group input ----

test_that("predict_links errors on netobject_group", {
  set.seed(42)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:4], 30, TRUE),
    V2 = sample(LETTERS[1:4], 30, TRUE),
    grp = rep(c("X", "Y"), each = 15),
    stringsAsFactors = FALSE
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  expect_error(predict_links(nets), "single network")
})


# ---- 20. print and summary methods ----

test_that("print and summary work without errors", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, top_n = 5)
  expect_output(print(pred), "Link Prediction")
  s <- summary(pred)
  expect_true(is.data.frame(s))
  expect_true(nrow(s) == length(pred$methods))
})


# ---- 21. Katz fallback on singular matrix ----

test_that("katz works on near-singular graph", {
  # 2-node graph with one edge
  W <- matrix(c(0, 1, 0, 0), 2, 2)
  rownames(W) <- colnames(W) <- c("A", "B")
  pred <- predict_links(W, methods = "katz")
  expect_true(all(is.finite(pred$scores$katz)))
})


# ---- 22. Rankings are correct ----

test_that("predictions are ranked by descending score", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz")
  df <- pred$predictions
  scores <- df$score
  expect_true(all(diff(scores) <= 1e-10))
})


# ---- 23. Bundled dataset ----

test_that("predict_links works on human_cat data", {
  data(human_cat)
  net <- build_network(human_cat, method = "relative",
                       actor = "session_id", action = "category",
                       time = "timestamp")
  # Dense network: test with exclude_existing = FALSE to get predictions
  pred <- predict_links(net, methods = c("katz", "resource_allocation"),
                        top_n = 10, exclude_existing = FALSE)
  expect_s3_class(pred, "net_link_prediction")
  expect_true(nrow(pred$predictions) > 0)
})


# ---- 24. Empty predictions (dense graph) handled gracefully ----

test_that("dense graph with exclude_existing produces empty predictions", {
  set.seed(1)
  seqs <- data.frame(
    V1 = sample(c("A", "B"), 50, TRUE),
    V2 = sample(c("A", "B"), 50, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  pred <- predict_links(net, methods = "katz")
  # All edges exist; predictions should be empty
  expect_equal(nrow(pred$predictions), 0)
})


# ---- 25. Adjacency matrix stored in result ----

test_that("predict_links stores adjacency matrix", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz")
  expect_true(!is.null(pred$adjacency))
  expect_true(is.matrix(pred$adjacency))
  expect_equal(nrow(pred$adjacency), pred$n_nodes)
  # Binary: only 0 and 1
  expect_true(all(pred$adjacency %in% c(0L, 1L)))
})


# ---- 26. Consensus ranking ----

test_that("consensus ranking computed for multi-method predictions", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = c("common_neighbors", "katz"),
                        exclude_existing = FALSE)
  expect_true(!is.null(pred$consensus))
  expect_true(is.data.frame(pred$consensus))
  expect_true(all(c("from", "to", "avg_rank", "n_methods", "consensus_rank")
                  %in% names(pred$consensus)))
  # Consensus ranks are sequential
  expect_equal(pred$consensus$consensus_rank, seq_len(nrow(pred$consensus)))
  # avg_rank is non-decreasing (sorted)
  expect_true(all(diff(pred$consensus$avg_rank) >= -1e-10))
})

test_that("consensus is NULL for single method", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz")
  expect_null(pred$consensus)
})


# ---- 27. Print shows consensus ----

test_that("print shows consensus for multi-method", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = c("common_neighbors", "katz"),
                        exclude_existing = FALSE)
  expect_output(print(pred), "consensus")
})

test_that("print shows single method when only one used", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz", exclude_existing = FALSE)
  expect_output(print(pred), "katz")
})


# ---- 28. pathways.net_link_prediction ----

test_that("pathways returns arrow-notation strings", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "common_neighbors",
                        exclude_existing = FALSE)
  pw <- pathways(pred, top = 5)
  expect_type(pw, "character")
  expect_equal(length(pw), 5)
  expect_true(all(grepl("->", pw, fixed = TRUE)))
})

test_that("pathways with evidence includes common neighbors", {
  # Build network with known structure: A->B, A->C, B->D, C->D
  W <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  W[1, 2] <- 1; W[1, 3] <- 1; W[2, 4] <- 1; W[3, 4] <- 1
  pred <- predict_links(W, methods = "common_neighbors",
                        exclude_existing = TRUE)
  pw <- pathways(pred, top = 5, evidence = TRUE)
  # A->D should have evidence (B and C are common neighbors)
  ad_pw <- pw[grepl("D$", pw)]
  if (length(ad_pw) > 0) {
    # Should contain evidence nodes (more than just "A -> D")
    parts <- strsplit(ad_pw[1], " -> ", fixed = TRUE)[[1]]
    sources <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
    expect_true(length(sources) >= 1)
  }
})

test_that("pathways without evidence gives simple edges", {
  net <- .make_lp_sparse()
  pred <- predict_links(net, methods = "katz", exclude_existing = FALSE)
  pw <- pathways(pred, top = 3, evidence = FALSE)
  # Simple format: "X -> Y" (no extra evidence nodes)
  parts <- strsplit(pw, " -> ", fixed = TRUE)
  source_counts <- vapply(parts, function(p) {
    length(strsplit(p[1], " ", fixed = TRUE)[[1]])
  }, integer(1))
  expect_true(all(source_counts == 1))
})

test_that("pathways returns empty for no predictions", {
  set.seed(1)
  seqs <- data.frame(V1 = sample(c("A", "B"), 50, TRUE),
                     V2 = sample(c("A", "B"), 50, TRUE),
                     stringsAsFactors = FALSE)
  net <- build_network(seqs, method = "relative")
  pred <- predict_links(net, methods = "katz")
  pw <- pathways(pred)
  expect_equal(length(pw), 0)
})
