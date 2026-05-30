# Equivalence tests locking the correct behavior for three regressions caught
# by the 2026-05-10 Codex review of the working tree:
#
#   #1  prepare_onehot(window_size = 1L, window_type = "overlapping") drops
#       the first row unconditionally even though there is no overlap to remove.
#   #2  .count_transitions_wide() ignores `weighted = TRUE` whenever the input
#       is windowed (attr "windowed" is TRUE).
#   #3  .count_cooccurrence_wide() ignores `weighted = TRUE` whenever the
#       input is windowed.
#
# These tests are deterministic worked examples (every cell is hand-checkable),
# not random. They double as regression tests so future refactors cannot
# silently revert.

# ---------------------------------------------------------------------------
# #1: prepare_onehot — size-1 overlapping must preserve all rows
# ---------------------------------------------------------------------------

test_that("prepare_onehot(window_size = 1L, overlapping) keeps the first row per actor", {
  # Single-row actor: with the original bug, the leading row of each group
  # was unconditionally removed even when window_size = 1L (no overlap to
  # remove), leaving zero rows. After the fix the row survives and the
  # output mirrors size-1 non-overlapping (which is the semantically
  # equivalent configuration).
  df <- data.frame(
    A = c(1L, 0L),
    B = c(0L, 1L),
    actor = c("a1", "a2")  # one row per actor
  )
  ov <- prepare_onehot(df, cols = c("A", "B"), actor = "actor",
                       window_size = 1L, window_type = "overlapping")
  no <- prepare_onehot(df, cols = c("A", "B"), actor = "actor",
                       window_size = 1L, window_type = "non-overlapping")
  expect_equal(nrow(ov), nrow(no))
  expect_equal(nrow(ov), nrow(df))
})

test_that("prepare_onehot size-1 overlapping == size-1 non-overlapping", {
  # With window_size = 1L the two window types are semantically identical:
  # each row is its own window. They must produce the same state cells.
  df <- data.frame(
    A = c(1L, 0L, 1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L, 0L, 1L),
    actor = c(1L, 1L, 1L, 2L, 2L, 2L)
  )
  ov <- prepare_onehot(df, cols = c("A", "B"), actor = "actor",
                       window_size = 1L, window_type = "overlapping")
  no <- prepare_onehot(df, cols = c("A", "B"), actor = "actor",
                       window_size = 1L, window_type = "non-overlapping")
  expect_equal(nrow(ov), nrow(no))
  expect_equal(ov$A, no$A)
  expect_equal(ov$B, no$B)
})

# ---------------------------------------------------------------------------
# #2: .count_transitions_wide — windowed branch must honor `weighted`
# ---------------------------------------------------------------------------
#
# Hand calculation:
#   data has 2 rows; states A, B; window_size = 2, span = 1, nc = 4.
#   Effective window = 2; q = 4/2 - 1 = 1; loop runs once with
#   j_idx = c(1, 2), k_idx = c(3, 4).
#
#   Row S1 = (A, B, A, B), all valid → 4 windowed pairs:
#     (1,3) A→A, (1,4) A→B, (2,3) B→A, (2,4) B→B  [one of each]
#   Row S2 = (A, B, NA, NA) → all (j,k) pairs hit at least one NA → 0 pairs.
#
#   Unweighted: each cell count = 1 (4 transitions total).
#   Weighted (weight = 1/row_seq_length, matching the non-windowed branch):
#     S1 length = 4, contributes 1/4 to each of the 4 cells; S2 contributes 0.
#     Each cell = 0.25; total weighted mass = 1.
#
# So weighted MUST differ from unweighted, and must equal unweighted / 4 in
# this configuration.

test_that(".count_transitions_wide windowed honors weighted = TRUE", {
  data <- data.frame(
    id = c("S1", "S2"),
    T1 = c("A", "A"),
    T2 = c("B", "B"),
    T3 = c("A", NA_character_),
    T4 = c("B", NA_character_),
    stringsAsFactors = FALSE
  )
  attr(data, "windowed") <- TRUE
  attr(data, "window_size") <- 2L
  attr(data, "window_span") <- 1L
  cols <- c("T1", "T2", "T3", "T4")

  unweighted <- Nestimate:::.count_transitions_wide(
    data, id = "id", cols = cols, alphabet = c("A", "B"), weighted = FALSE)
  weighted <- Nestimate:::.count_transitions_wide(
    data, id = "id", cols = cols, alphabet = c("A", "B"), weighted = TRUE)

  expect_false(identical(unweighted, weighted),
               info = "weighted must change the windowed transition counts")
  # Exact hand-computed values:
  expect_equal(unname(unweighted),
               matrix(c(1, 1, 1, 1), 2, 2,
                      dimnames = list(c("A","B"), c("A","B"))) |> unname())
  expect_equal(unname(weighted),
               matrix(c(0.25, 0.25, 0.25, 0.25), 2, 2,
                      dimnames = list(c("A","B"), c("A","B"))) |> unname())
})

# ---------------------------------------------------------------------------
# #3: .count_cooccurrence_wide — windowed branch must honor `weighted`
# ---------------------------------------------------------------------------
#
# Hand calculation:
#   data has 2 rows; states A, B; window_size = 2, span = 1, nc = 2.
#   Effective window = 2; q = 2/2 - 1 = 0; loop runs once (q + 1 = 1) with
#   j_idx = c(1, 2). Inside, both j and k iterate over j_idx.
#
#   Row S1 = (A, B):
#     pairs (j,k) ∈ {(1,1), (1,2), (2,1), (2,2)} = {A→A, A→B, B→A, B→B}
#   Row S2 = (A, A):
#     pairs (j,k) = {A→A, A→A, A→A, A→A} = 4×A→A
#
#   Unweighted: A→A = 1 + 4 = 5; A→B = 1; B→A = 1; B→B = 1.
#   Weighted (weight 1/row_seq_length): S1 length = 2, S2 length = 2;
#     each pair contributes 1/2.
#     A→A = 1/2 + 4/2 = 2.5; A→B = 0.5; B→A = 0.5; B→B = 0.5.

test_that(".count_cooccurrence_wide windowed honors weighted = TRUE", {
  data <- data.frame(
    id = c("S1", "S2"),
    T1 = c("A", "A"),
    T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )
  attr(data, "windowed") <- TRUE
  attr(data, "window_size") <- 2L
  attr(data, "window_span") <- 1L
  cols <- c("T1", "T2")

  unweighted <- Nestimate:::.count_cooccurrence_wide(
    data, id = "id", cols = cols, alphabet = c("A", "B"), weighted = FALSE)
  weighted <- Nestimate:::.count_cooccurrence_wide(
    data, id = "id", cols = cols, alphabet = c("A", "B"), weighted = TRUE)

  expect_false(identical(unweighted, weighted),
               info = "weighted must change the windowed co-occurrence counts")
  expect_equal(unname(unweighted),
               matrix(c(5, 1, 1, 1), 2, 2) |> unname())
  expect_equal(unname(weighted),
               matrix(c(2.5, 0.5, 0.5, 0.5), 2, 2) |> unname())
})
