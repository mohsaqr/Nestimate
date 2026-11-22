# wtna(): window-based transition and co-occurrence networks of one-hot data, against tna::import_onehot()
# followed by ftna(), tna() and tna's co-occurrence model (tna: the owner's co-authored package, a separate
# implementation; open-source, to be confirmed by the owner). Re-expressed from
# local_testing_and_equivalence/test-equiv-wtna-overlapping.R (25 seeds x 4 window sizes; here 2 seeds x 3 sizes
# x 2 modes). The weights are read with extract_transition_matrix().

# One-hot data: four 0/1 code columns and 2-5 actors in blocks, as the old test drew it.
tn_wtna_data <- function(seed) {
  set.seed(seed)
  codes <- c("A", "B", "C", "D")
  n <- sample(c(30L, 45L, 60L), 1L)
  data <- as.data.frame(matrix(stats::rbinom(4L * n, 1L, stats::runif(1L, .2, .6)), ncol = 4L,
                               dimnames = list(NULL, codes)))
  data$actor <- sort(rep_len(seq_len(sample(2:5, 1L)), n))
  data
}

# wtna()'s three networks and tna's, for one input, as two named lists of matrices.
tn_wtna_pair <- function(seed, window, mode) {
  codes <- c("A", "B", "C", "D")
  data <- tn_wtna_data(seed)
  imported <- tna::import_onehot(data, cols = codes, actor = "actor", window_size = window,
                                 window_type = if (identical(mode, "overlapping")) "sliding" else "tumbling")
  found <- \(method, type) extract_transition_matrix(wtna(data, method = method, type = type, codes = codes,
                                                         actor = "actor", window_size = window, mode = mode))
  list(found = list(counts = found("transition", "frequency"), probabilities = found("transition", "relative"),
                    cooccurrence = found("cooccurrence", "frequency")),
       reference = list(counts = unclass(tna::ftna(imported)$weights),
                        probabilities = unclass(tna::tna(imported)$weights),
                        cooccurrence = unclass(tna::build_model(imported, type = "co-occurrence",
                                                                scaling = character(0))$weights)))
}

parts <- c(counts = "transition counts", probabilities = "transition probabilities",
           cooccurrence = "co-occurrence counts")

test_that("tumbling windows, and sliding windows of width 2, equal tna::import_onehot() networks", {
  skip_without_reference("tna")
  grid <- rbind(expand.grid(seed = c(3L, 17L), window = c(2L, 3L, 5L), mode = "non-overlapping", stringsAsFactors = FALSE),
                expand.grid(seed = c(3L, 17L), window = 2L, mode = "overlapping", stringsAsFactors = FALSE))
  invisible(Map(\(seed, window, mode) {
    pair <- tn_wtna_pair(seed, window, mode)
    lapply(names(parts), \(part) {
      expect_reference(pair$found[[part]], pair$reference[[part]], reference = "tna", covers = "wtna()",
                       quantity = sprintf("%s, seed %d, window %d, %s", parts[[part]], seed, window, mode),
                       tolerance = 1e-10)
    })
  }, grid$seed, grid$window, grid$mode))
})

test_that("sliding windows wider than 2 differ from tna 1.3.1, whose sliding rule tna 1.5.5 corrected", {
  skip_without_reference("tna")
  # tna 1.3.1 is the version on the ci machine. LEARNINGS.md: "installed tna here is 1.3.1 with the old sliding
  # window. tna 1.5.5 source has the corrected rule that Nestimate now matches." Width 2 is the same under both
  # rules (compared above).
  grid <- expand.grid(seed = c(3L, 17L), window = c(3L, 5L))
  invisible(Map(\(seed, window) {
    pair <- tn_wtna_pair(seed, window, "overlapping")
    lapply(names(parts), \(part) {
      expect_known_difference(pair$found[[part]], pair$reference[[part]], reference = "tna", covers = "wtna()",
                              quantity = sprintf("%s, seed %d, window %d, overlapping", parts[[part]], seed, window),
                              tolerance = 1e-10, agrees_with = "mode = \"non-overlapping\"",
                              reason = "tna 1.3.1 builds sliding windows by the rule tna 1.5.5 replaced; Nestimate follows the corrected rule.")
    })
  }, grid$seed, grid$window))
})

test_that("build_wtna() counts and probabilities equal tna::import_onehot() networks", {
  skip_without_reference("tna")
  codes <- c("A", "B", "C", "D")
  invisible(Map(\(seed, window) {
    data <- tn_wtna_data(seed)
    imported <- tna::import_onehot(data, cols = codes, actor = "actor", window_size = window, window_type = "tumbling")
    label <- sprintf("seed %d, window %d, non-overlapping", seed, window)
    expect_reference(extract_transition_matrix(build_wtna(data, codes = codes, actor = "actor", window_size = window)),
                     unclass(tna::ftna(imported)$weights), reference = "tna", covers = "build_wtna()",
                     quantity = paste("transition counts,", label), tolerance = 1e-10)
    expect_reference(extract_transition_matrix(build_wtna(data, codes = codes, actor = "actor", window_size = window,
                                                          type = "relative")),
                     unclass(tna::tna(imported)$weights), reference = "tna", covers = "build_wtna()",
                     quantity = paste("transition probabilities,", label), tolerance = 1e-10)
  }, c(3L, 17L, 3L, 17L), c(2L, 2L, 3L, 3L)))
})
