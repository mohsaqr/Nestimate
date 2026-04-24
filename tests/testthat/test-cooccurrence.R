# ---- cooccurrence() tests ----

# === Shared test data ===
# Three transactions:  {A, B, C}, {B, C}, {A, C}
# Co-occurrence matrix (raw counts):
#      A  B  C
#  A   2  1  2
#  B   1  2  2
#  C   2  2  3
# Frequencies: A=2, B=2, C=3, N=3

.co_test_list <- list(c("A", "B", "C"), c("B", "C"), c("A", "C"))
.co_expected_raw <- matrix(
  c(2, 1, 2,
    1, 2, 2,
    2, 2, 3),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
)

# ========================================
# 1. Input format tests
# ========================================

test_that("list input produces correct co-occurrence matrix", {
  net <- cooccurrence(.co_test_list)
  expect_s3_class(net, "netobject")
  expect_equal(net$weights, .co_expected_raw)
  expect_false(net$directed)
  expect_equal(net$method, "co_occurrence_fn")
})

test_that("delimited input produces correct co-occurrence matrix", {
  df <- data.frame(
    id = 1:3,
    items = c("A; B; C", "B; C", "A; C"),
    stringsAsFactors = FALSE
  )
  net <- cooccurrence(df, field = "items", sep = ";")
  expect_equal(net$weights, .co_expected_raw)
})

test_that("multi-column delimited input pools across columns", {
  df <- data.frame(
    col1 = c("A; B", "B", "A"),
    col2 = c("C", "C", "C"),
    stringsAsFactors = FALSE
  )
  net <- cooccurrence(df, field = c("col1", "col2"), sep = ";")
  expect_equal(net$weights, .co_expected_raw)
})

test_that("long/bipartite input produces correct co-occurrence matrix", {
  df <- data.frame(
    doc = c(1, 1, 1, 2, 2, 3, 3),
    item = c("A", "B", "C", "B", "C", "A", "C"),
    stringsAsFactors = FALSE
  )
  net <- cooccurrence(df, field = "item", by = "doc")
  expect_equal(net$weights, .co_expected_raw)
})

test_that("binary matrix input produces correct co-occurrence matrix", {
  bin <- matrix(
    c(1, 1, 1,
      0, 1, 1,
      1, 0, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(NULL, c("A", "B", "C"))
  )
  net <- cooccurrence(bin)
  expect_equal(net$weights, .co_expected_raw)
})

test_that("binary data.frame input produces correct co-occurrence matrix", {
  df <- data.frame(A = c(1, 0, 1), B = c(1, 1, 0), C = c(1, 1, 1))
  net <- cooccurrence(df)
  expect_equal(net$weights, .co_expected_raw)
})

test_that("wide sequence input produces correct co-occurrence matrix", {
  df <- data.frame(
    V1 = c("A", "B", "A"),
    V2 = c("B", "C", "C"),
    V3 = c("C", NA, NA),
    stringsAsFactors = FALSE
  )
  net <- cooccurrence(df)
  expect_equal(net$weights, .co_expected_raw)
})

# ========================================
# 2. All formats produce identical results
# ========================================

test_that("all 6 formats produce the same co-occurrence matrix", {
  net_list <- cooccurrence(.co_test_list)

  df_del <- data.frame(items = c("A;B;C", "B;C", "A;C"),
                        stringsAsFactors = FALSE)
  net_del <- cooccurrence(df_del, field = "items", sep = ";")

  df_multi <- data.frame(c1 = c("A;B", "B", "A"),
                          c2 = c("C", "C", "C"),
                          stringsAsFactors = FALSE)
  net_multi <- cooccurrence(df_multi, field = c("c1", "c2"), sep = ";")

  df_long <- data.frame(
    doc = c(1, 1, 1, 2, 2, 3, 3),
    item = c("A", "B", "C", "B", "C", "A", "C"),
    stringsAsFactors = FALSE
  )
  net_long <- cooccurrence(df_long, field = "item", by = "doc")

  bin <- matrix(c(1, 1, 1, 0, 1, 1, 1, 0, 1), nrow = 3, byrow = TRUE,
                dimnames = list(NULL, c("A", "B", "C")))
  net_bin <- cooccurrence(bin)

  df_wide <- data.frame(V1 = c("A", "B", "A"), V2 = c("B", "C", "C"),
                         V3 = c("C", NA, NA), stringsAsFactors = FALSE)
  net_wide <- cooccurrence(df_wide)

  expect_equal(net_del$weights, net_list$weights)
  expect_equal(net_multi$weights, net_list$weights)
  expect_equal(net_long$weights, net_list$weights)
  expect_equal(net_bin$weights, net_list$weights)
  expect_equal(net_wide$weights, net_list$weights)
})

# ========================================
# 3. Similarity tests
# ========================================

test_that("similarity = 'none' returns raw counts", {
  net <- cooccurrence(.co_test_list, similarity = "none")
  expect_equal(net$weights, .co_expected_raw)
})

test_that("similarity = 'jaccard' is correct", {
  net <- cooccurrence(.co_test_list, similarity = "jaccard")
  W <- net$weights
  # Jaccard(A,B) = C[A,B] / (freq_A + freq_B - C[A,B]) = 1 / (2+2-1) = 1/3
  expect_equal(W["A", "B"], 1 / 3)
  # Jaccard(A,C) = 2 / (2+3-2) = 2/3
  expect_equal(W["A", "C"], 2 / 3)
  # Jaccard(B,C) = 2 / (2+3-2) = 2/3
  expect_equal(W["B", "C"], 2 / 3)
  # Diagonal: Jaccard(A,A) = 2 / (2+2-2) = 1
  expect_equal(W["A", "A"], 1)
  # Symmetric
  expect_equal(W["A", "B"], W["B", "A"])
})

test_that("similarity = 'cosine' is correct", {
  net <- cooccurrence(.co_test_list, similarity = "cosine")
  W <- net$weights
  # Cosine(A,B) = 1 / sqrt(2*2) = 0.5
  expect_equal(W["A", "B"], 0.5)
  # Cosine(A,C) = 2 / sqrt(2*3) = 2/sqrt(6)
  expect_equal(W["A", "C"], 2 / sqrt(6))
  # Diagonal: Cosine(A,A) = 2 / sqrt(2*2) = 1
  expect_equal(W["A", "A"], 1)
})

test_that("similarity = 'inclusion' is correct", {
  net <- cooccurrence(.co_test_list, similarity = "inclusion")
  W <- net$weights
  # Inclusion(A,B) = 1 / min(2,2) = 0.5
  expect_equal(W["A", "B"], 0.5)
  # Inclusion(A,C) = 2 / min(2,3) = 1
  expect_equal(W["A", "C"], 1)
  # Inclusion(B,C) = 2 / min(2,3) = 1
  expect_equal(W["B", "C"], 1)
})

test_that("similarity = 'association' is correct", {
  net <- cooccurrence(.co_test_list, similarity = "association")
  W <- net$weights
  # Association(A,B) = 1 / (2*2) = 0.25
  expect_equal(W["A", "B"], 0.25)
  # Association(A,C) = 2 / (2*3) = 1/3
  expect_equal(W["A", "C"], 1 / 3)
  # Association(B,C) = 2 / (2*3) = 1/3
  expect_equal(W["B", "C"], 1 / 3)
})

test_that("similarity = 'dice' is correct", {
  net <- cooccurrence(.co_test_list, similarity = "dice")
  W <- net$weights
  # Dice(A,B) = 2*1 / (2+2) = 0.5
  expect_equal(W["A", "B"], 0.5)
  # Dice(A,C) = 2*2 / (2+3) = 0.8
  expect_equal(W["A", "C"], 0.8)
  # Dice(A,A) = 2*2 / (2+2) = 1
  expect_equal(W["A", "A"], 1)
})

test_that("similarity = 'equivalence' is correct", {
  net <- cooccurrence(.co_test_list, similarity = "equivalence")
  W <- net$weights
  # Equivalence(A,B) = 1^2 / (2*2) = 0.25
  expect_equal(W["A", "B"], 0.25)
  # Equivalence(A,C) = 2^2 / (2*3) = 4/6 = 2/3
  expect_equal(W["A", "C"], 2 / 3)
  # Equivalence(A,A) = 2^2 / (2*2) = 1
  expect_equal(W["A", "A"], 1)
})

test_that("similarity = 'relative' row-normalizes correctly", {
  net <- cooccurrence(.co_test_list, similarity = "relative")
  W <- net$weights
  # Each row sums to 1
  expect_equal(rowSums(W), c(A = 1, B = 1, C = 1))
  # A row: 2/5, 1/5, 2/5
  expect_equal(W["A", ], c(A = 2 / 5, B = 1 / 5, C = 2 / 5))
})

# ========================================
# 4. Parameters: threshold, min_occur, diagonal, top_n
# ========================================

test_that("threshold filters low-weight edges after normalization", {
  # Raw counts: A-B = 1, A-C = 2, B-C = 2
  # With similarity = "none", threshold = 2 drops A-B
  net <- cooccurrence(.co_test_list, threshold = 2)
  W <- net$weights
  expect_equal(W["A", "B"], 0)
  expect_equal(W["B", "A"], 0)
  expect_equal(W["A", "C"], 2)
  expect_equal(W["B", "C"], 2)
})

test_that("threshold applies after normalization", {
  # Jaccard: A-B = 1/3, A-C = 2/3, B-C = 2/3
  net <- cooccurrence(.co_test_list, similarity = "jaccard", threshold = 0.5)
  W <- net$weights
  # A-B = 1/3 < 0.5 → zeroed
  expect_equal(W["A", "B"], 0)
  # A-C = 2/3 >= 0.5 → kept
  expect_equal(W["A", "C"], 2 / 3)
})

test_that("min_occur drops infrequent entities", {
  # A appears in 2 transactions, B in 2, C in 3
  # min_occur = 3 keeps only C → but single-item transactions → error
  # min_occur = 2 keeps A, B, C (all >= 2)
  trans <- list(c("A", "B", "C"), c("B", "C"), c("A", "C"), c("D"))
  net <- cooccurrence(trans, min_occur = 2L)
  expect_false("D" %in% colnames(net$weights))
  expect_true(all(c("A", "B", "C") %in% colnames(net$weights)))
})

test_that("min_occur = 1 keeps all entities (default)", {
  trans <- list(c("A", "B"), c("B", "C"), c("D"))
  net <- cooccurrence(trans, min_occur = 1L)
  expect_true("D" %in% colnames(net$weights))
})

test_that("diagonal = FALSE zeros the diagonal", {
  net <- cooccurrence(.co_test_list, diagonal = FALSE)
  expect_equal(diag(net$weights), c(A = 0, B = 0, C = 0))
  expect_equal(net$weights["A", "B"], 1)
})

test_that("diagonal = FALSE + similarity works", {
  net <- cooccurrence(.co_test_list, diagonal = FALSE, similarity = "cosine")
  expect_equal(diag(net$weights), c(A = 0, B = 0, C = 0))
  # freq still computed from column sums of B (not diagonal)
  expect_equal(net$weights["A", "B"], 0.5)
})

test_that("top_n keeps only top N edges", {
  # 3 unique off-diagonal pairs: A-B=1, A-C=2, B-C=2
  net <- cooccurrence(.co_test_list, diagonal = FALSE, top_n = 2L)
  W <- net$weights
  # A-B (weight=1) should be zeroed; A-C and B-C (weight=2) kept
  expect_equal(W["A", "B"], 0)
  expect_equal(W["A", "C"], 2)
  expect_equal(W["B", "C"], 2)
})

# ========================================
# 5. netobject integration
# ========================================

test_that("co_occurrence returns valid netobject structure", {
  net <- cooccurrence(.co_test_list)
  expect_s3_class(net, "netobject")
  expect_s3_class(net, "cograph_network")
  expect_false(net$directed)
  expect_equal(net$n_nodes, 3L)
  expect_true(net$n_edges > 0L)
  expect_equal(nrow(net$nodes), 3L)
  expect_true(all(c("id", "label", "name") %in% names(net$nodes)))
  expect_true(all(c("from", "to", "weight") %in% names(net$edges)))
})

test_that("co_occurrence params are stored correctly", {
  net <- cooccurrence(.co_test_list, similarity = "jaccard", threshold = 0.1)
  expect_equal(net$params$similarity, "jaccard")
  expect_equal(net$params$threshold, 0.1)
  expect_equal(net$params$n_transactions, 3L)
  expect_equal(net$params$n_items, 3L)
})


# ========================================
# 6. Edge cases
# ========================================

test_that("single-item transactions produce zero off-diagonal", {
  trans <- list(c("A"), c("B"), c("C"))
  net <- cooccurrence(trans)
  W <- net$weights
  expect_equal(W["A", "B"], 0)
  expect_equal(W["B", "C"], 0)
  expect_equal(diag(W), c(A = 1, B = 1, C = 1))
})

test_that("all-identical transactions give expected counts", {
  trans <- list(c("X", "Y"), c("X", "Y"), c("X", "Y"))
  net <- cooccurrence(trans)
  expect_equal(net$weights["X", "Y"], 3)
  expect_equal(net$weights["X", "X"], 3)
})

test_that("duplicate items within a transaction are de-duplicated", {
  trans <- list(c("A", "A", "B"), c("B", "B", "C"))
  net <- cooccurrence(trans)
  expect_equal(net$weights["A", "A"], 1)
  expect_equal(net$weights["A", "B"], 1)
})

test_that("empty strings and NAs are removed from transactions", {
  trans <- list(c("A", "", NA, "B"), c("B", "C"))
  net <- cooccurrence(trans)
  expect_true("A" %in% colnames(net$weights))
  expect_true("B" %in% colnames(net$weights))
})

test_that("delimited with extra whitespace is handled", {
  df <- data.frame(kw = c("  A ; B  ;  C  ", " B;C "), stringsAsFactors = FALSE)
  net <- cooccurrence(df, field = "kw", sep = ";")
  expect_equal(sort(colnames(net$weights)), c("A", "B", "C"))
  expect_equal(net$weights["A", "B"], 1)
})

test_that("error on empty input", {
  expect_error(cooccurrence(list()), "No non-empty transactions")
  expect_error(cooccurrence(list(character(0), character(0))),
               "No non-empty transactions")
})

test_that("error on invalid format", {
  expect_error(cooccurrence("not_a_valid_input"), "Cannot detect input format")
})

test_that("large number of items works", {
  items <- paste0("item_", seq_len(50))
  trans <- lapply(seq_len(20), function(i) sample(items, 10))
  net <- cooccurrence(trans)
  expect_equal(ncol(net$weights), length(unique(unlist(trans))))
  expect_true(isSymmetric(net$weights))
})

test_that("binary matrix without colnames gets auto-named", {
  bin <- matrix(c(1, 0, 1, 1), nrow = 2)
  net <- cooccurrence(bin)
  expect_equal(colnames(net$weights), c("V1", "V2"))
})

test_that("wide sequence with void markers removes them", {
  df <- data.frame(V1 = c("A", "%"), V2 = c("B", "*"),
                   stringsAsFactors = FALSE)
  net <- cooccurrence(df)
  expect_equal(net$params$n_transactions, 1L)
})

test_that("min_occur filters before co-occurrence computation", {
  # D appears only once → dropped with min_occur = 2
  trans <- list(c("A", "B"), c("A", "B"), c("A", "D"))
  net <- cooccurrence(trans, min_occur = 2L)
  expect_false("D" %in% colnames(net$weights))
  expect_equal(net$weights["A", "B"], 2)
})
